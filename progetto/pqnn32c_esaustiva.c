#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <limits.h>

#define	MATRIX		float*
#define	VECTOR		float*

#define DATASET		0
#define QUERYSET	1


typedef struct {
	char* filename; //
	MATRIX ds; // data set 
	MATRIX qs; // query set
	int n; // numero di punti del data set
	int d; // numero di dimensioni del data/query set
	int nq; // numero di punti del query set
	int knn; // numero di ANN approssimati da restituire per ogni query
	int m; // numero di gruppi del quantizzatore prodotto
	int k; // numero di centroidi di ogni sotto-quantizzatore
	int kc; // numero di centroidi del quantizzatore coarse
	int w; // numero di centroidi del quantizzatore coarse da selezionare per la ricerca non esaustiva
	int nr; // dimensione del campione dei residui nel caso di ricerca non esaustiva
	float eps; // 
	int tmin; //
	int tmax; //
	int exaustive; // tipo di ricerca: (0=)non esaustiva o (1=)esaustiva
	int symmetric; // tipo di distanza: (0=)asimmetrica ADC o (1=)simmetrica SDC
	int silent;
	int display;
	// nns: matrice row major order di interi a 32 bit utilizzata per memorizzare gli ANN
	// sulla riga i-esima si trovano gli ID (a partire da 0) degli ANN della query i-esima
	//
	int* ANN;
	VECTOR vq;
	int* pq;
	int* query_pq;

	MATRIX codebook; // per E. contiene quantizzatori prodotto. Per N.E. contiene quantizzatori grossolani
	
	MATRIX distanze_simmetriche;
	int nDist;
	MATRIX distanze_asimmetriche;
	// ---------------------------------------
	// Strutture ad-hoc ricerca non esaustiva

	// Vettore contenente alla posizione i l'indice di qc(Y_i) in codebook
	unsigned short * qc_indexes;

	// Residual product quantizators in non exhaustive search
	// matrix type: Row-major-order
	MATRIX residual_codebook;

	// Learning set nella ricerca non esastiva. Contiene i residui dei primi nr vettori
	MATRIX residual_set;

	// Lista di liste (secondo livello dell'inverted index)
	struct entry* v; 
} params;

int calcolaIndice(int i, int j){
	//funzione che calcola l'indice per la matrice delle distanze_simmetriche
	return i*(i-1)/2+j;
}

int dist_eI(params* input,int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	int i;
	int ret=0;
	for(i=start; i<end; i++){
		ret += pow(input->ds[punto1*input->d+i]-input->ds[punto2*input->d+i], 2.0);
	}
	return ret;
}

int dist_e(params* input,int punto1, int punto2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_eI(input, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

int calcolaPQ(params* input, int x, int start, int end){
	// estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=1.79E+308;
    int imin=-1;
    float temp;
    for(i=0; i<input->k; i++){
        temp=dist_eI(input, x, i, start, end);
        if(temp<min){ 
            min=temp;
            imin=i;
        }
    }
    return imin;
}

float dist_simmetricaI(params* input, int centroide1, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	for(i=start; i<end; i++){
		ret += pow( input->codebook[centroide1*input->d+i] - input->codebook[centroide2*input->d+i] , 2);
	}
	return ret;
}

float dist_simmetrica(params* input, int centroide1, int centroide2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_simmetricaI(input, centroide1, centroide2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

float dist_asimmetricaI(params* input, MATRIX set, int punto1, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int i, c;
	float ret=0;
	for(i=start; i<end; i++){
		ret += pow( set[punto1*input->d+i] - input->codebook[centroide2*input->d+i] , 2);
	}
	return ret;
}

float dist_asimmetrica(params* input, MATRIX set, int punto1, int punto2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_asimmetricaI(input, set, punto1, input->pq[punto2], i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

float distI(params* input, int* quantizer, int punto1, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int c1=quantizer[punto1*input->m+(start/input->m)];
	//printf("breakpoint Dist %d %d\n", c1, centroide2);
	if(c1==centroide2){
		return 0;
	}else{
		//row major order-------------------------------------------
//		if(c1<centroide2){
//			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(centroide2, c1)];
//		}else{
//			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(c1, centroide2)];
//		}
		//column major order-------------------------------------------
		if(c1<centroide2){
			return input->distanze_simmetriche[start/input->m+calcolaIndice(centroide2, c1)*input->m];
		}else{
			return input->distanze_simmetriche[start/input->m+calcolaIndice(c1, centroide2)*input->m];
		}
		//-------------------------------------------
	}
}

float dist(params* input, int* quantizer, int punto1, int punto2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(distI(input, quantizer, punto1, input->pq[punto2], i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

int calcolaQueryPQ(params* input, int x, int start, int end){
	// estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=1.79E+308;
    int imin=-1;
    float temp;
	//printf("breakpoint PQ\n");
	// valore 2 messo temporaneamente, inpossibile calcolare la distanza simmetrica per calcolare il centroide
	if(input->symmetric==2){
		for(i=0; i<input->k; i++){
			temp=distI(input, input->query_pq, x, i, start, end);
			if(temp<min){ 
				min=temp;
				imin=i;
			}	
			//printf("breakpoint PQ %d\n", i);
		}
	}else{
		for(i=0; i<input->k; i++){
			temp=dist_asimmetricaI(input, input->qs, x, i, start, end);
			if(temp<min){ 
				min=temp;
				imin=i;
			}
		}
	}
	return imin;
}

void kmeans(params* input, int start, int end, int n_centroidi){
	// estremi start incluso ed end escluso
	int i, j, k, t;
	int count;
	float fob1, fob2;
	MATRIX codebook;
	VECTOR min;

	codebook = alloc_matrix(n_centroidi, input->n); // row-major-order?
    if(codebook==NULL) exit(-1);
	
	//
	// Inizializzazione del codebook
	//		-Scelta dei k vettori casuali
	//
	
	for(i=0; i<n_centroidi; i++){
		k=rand()%input->n;
		for(j=start; j<end; j++){
			codebook[i*input->d+j]=input->ds[k*input->d+j];
		}
	}
    
//	for(i=0; i<input->n; i++){
//		input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
//	}
	//--------------------------------------------------------
	min=(VECTOR) alloc_matrix(input->n, 1);
	for(i=0; i<input->n; i++){
		min[i]=1.79E+308;
	}
	float temp;
	for(i=0; i<input->n; i++){
		for(j=0; j<input->k; j++){
			temp=dist_eI(input, i, j, start, end);
			if(temp<min[i]){ 
				min[i]=temp;
				input->pq[i*input->m+(start/input->m)]=j;
			}
		}
	}
	//--------------------------------------------------------
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && (fob2-fob1) > input->eps); t++){
		for(i=0; i<n_centroidi; i++){
			count=0;
			for(j=start; j<end; j++){
				codebook[i*input->d+j]=0; // con calloc forse è più veloce. 
			}
			
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			
			for(j=0; j<input->n; j++){
				if(input->pq[j*input->m+(start/input->m)]==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					for(k=start; k<end; k++){
						codebook[i*input->d+k]+=input->ds[j*input->d+k];
					}
				}
			}
			
			for(j=start; j<end; j++){
				if(count!=0){ 
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					codebook[i*input->d+j]=codebook[i*input->d+j]/count;
				}
			}
			
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
		}
		
//		for(i=0; i<input->n; i++){
//			input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
//		}

		for(i=0; i<input->n; i++){
			min[i]=1.79E+308;
		}
		float temp;
		for(i=0; i<input->n; i++){
			for(j=0; j<input->k; j++){
				temp=dist_eI(input, i, j, start, end);
				if(temp<min[i]){ 
					min[i]=temp;
					input->pq[i*input->m+(start/input->m)]=j;
				}
			}
		}

//-----------------------------------
		
		fob1=fob2;
		fob2=0;
		
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(i=0; i<input->n; i++){
			fob2+=pow(dist_eI(input, i, input->pq[i*input->m+(start/input->m)], start, end), 2.0);
		}
	}
	
	_mm_free(min);

	input->codebook=codebook;
}

void creaMatricedistanze(params* input){
	int i, j, k;
	MATRIX distanze_simmetriche;
	input->nDist=input->k*(input->k+1)/2;
	distanze_simmetriche = alloc_matrix(input->m, input->nDist);
	if(distanze_simmetriche==NULL) exit(-1);
	//row major order---------------------------------------------------------
//	for(k=0; k<input->m; k++){
//		for(i=1; i<input->k; i++){
//			for(j=0; j<i; j++){
//				distanze_simmetriche[k*input->nDist+calcolaIndice(i, j)] = dist_simmetricaI(input, i, j, k*input->m, (k+1)*input->m);
//				// verificare se qui va usata la distanza simmetrica o no
//			}
//		}
//	}
	//column major order---------------------------------------------------------
	for(i=1; i<input->k; i++){
		for(j=0; j<i; j++){
			for(k=0; k<input->m; k++){
				distanze_simmetriche[k+calcolaIndice(i, j)*input->m] = dist_simmetricaI(input, i, j, k*input->m, (k+1)*input->m);
				// verificare se qui va usata la distanza simmetrica o no
			}
		}
	}
	//---------------------------------------------------------
	input->distanze_simmetriche=distanze_simmetriche;
}

void bubbleSort(VECTOR arr, int* arr2, int n, int nit){ 
   int i, j, t1;
   float t2; 
   for (i = 0; i < nit; i++)
       for (j = n-2; j > i-1; j--)  
           if (arr[j] > arr[j+1]){
			   t2=arr[j];
			   arr[j]=arr[j+1];
			   arr[j+1]=t2;
			   t1=arr2[j];
			   arr2[j]=arr2[j+1];
			   arr2[j+1]=t1;
		   }
} 

void calcolaNN(params* input, int query){
	int i, j, k;
	VECTOR distanze=alloc_matrix(input->n, 1);
	VECTOR m;
	int* ind;

	if(input->knn<1050){
		if(input->symmetric==0){
			for(i=0; i<input->n; i++){
				distanze[i]=dist_asimmetrica(input, input->qs, query, i);
			}
		}else{
			for(i=0; i<input->n; i++){
				distanze[i]=dist(input, input->query_pq, query, i);
			}
		}

		m=(VECTOR) _mm_malloc(input->knn*sizeof(float),16);
		for(i=0; i<input->knn; i++){
			m[i]=1.79E+308;
			input->ANN[query*input->knn+i]=-1;
		}
		for(i=0; i<input->n; i++){
			for(j=0; j<input->knn; j++){
				if(distanze[i]<m[j]){
					for(k=input->knn-1; k>j; k--){
						if(m[k-1]!=-1){
							input->ANN[query*input->knn+k]=input->ANN[query*input->knn+k-1];
							m[k]=m[k-1];
						}
					}
					input->ANN[query*input->knn+j]=i;
					m[j]=distanze[i];
					printf("%d %d", i, j);
					break;
				}
			}
		}
		_mm_free(m);
	}else{
		ind=(int*) _mm_malloc(input->n*sizeof(int), 16);
		//printf("breakpoint 1\n");
		if(input->symmetric==0){
			for(i=0; i<input->n; i++){
				distanze[i]=dist_asimmetrica(input, input->qs, query, i);
				ind[i]=i;
			}
		}else{
			for(i=0; i<input->n; i++){
				distanze[i]=dist(input, input->query_pq, query, i);
				ind[i]=i;
			}
		}
		bubbleSort(distanze, ind, input->n, input->knn);
		for(i=0; i<input->knn; i++){
			input->ANN[query*input->knn+i]=ind[i];
		}
		_mm_free(ind);
	}
	dealloc_matrix(distanze);
}

void pqnn_index_esaustiva(params* input){
	int i, dStar;
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 16); 
	dStar=input->d/input->m;
	for(i=0; i<input->m; i++){
		kmeans(input, i*dStar, (i+1)*dStar, input->k);
	}
	if(input->symmetric==1){
		creaMatricedistanze(input);
	}
}

void pqnn_search_esaustiva(params* input){
	int i, j, c;
	if(input->symmetric==1){
		input->query_pq=(int*)_mm_malloc(input->nq*input->m*sizeof(int), 16);
		if(input->query_pq==NULL) exit(-1);
		c=input->d/input->m;
		for(i=0; i<input->nq; i++){
			for(j=0; j<input->m; j++){
				input->query_pq[i*input->m+j]=calcolaQueryPQ(input, i, j*c, (j+1)*c);
			}
		}
	}
	for(i=0; i<input->nq; i++){
		calcolaNN(input, i);
	}
	_mm_free(input->codebook);
	_mm_free(input->pq);
	if(input->symmetric==1){
		_mm_free(input->distanze_simmetriche);
	}else{
		_mm_free(input->query_pq);
		_mm_free(input->distanze_asimmetriche);
	}
}