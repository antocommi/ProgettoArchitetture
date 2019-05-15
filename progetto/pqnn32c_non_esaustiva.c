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

//Entry della s.d. multilivello
struct entry{
	int index;
	VECTOR q;
	//temporaneo
	//Serve per gestire liste a dimensione sconosciuta. 
	struct entry * next;
};

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}


void free_block(void* p) { 
	_mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(float),rows*cols);
}


void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

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
	// 	estremi start incluso ed end escluso
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
	if(input->exaustive==0)
		input->qc_indexes[x]= (unsigned char) imin;
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

float dist_asimmetricaI(params* input, MATRIX set, int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int i, c;
	float ret=0;
	c=input->pq[punto2];
	for(i=start; i<end; i++){
		ret += pow( set[punto1*input->d+i] - input->codebook[c*input->d+i] , 2);
	}
	return ret;
}

float dist_asimmetrica(params* input, MATRIX set, int punto1, int punto2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_asimmetricaI(input, set, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

float distI(params* input, int* quantizer, int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int c1, c2;
	if(punto1==punto2){
		return 0;
	}else{
		c1=quantizer[punto1*input->m+(start/input->m)];
		c2=input->pq[punto2*input->m+(start/input->m)];
		if(c1<c2){
			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(c2, c1)];
		}else{
			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(c1, c2)];
		}
	}
}
float dist(params* input, int* quantizer, int punto1, int punto2){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(distI(input, quantizer, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

int calcolaQueryPQ(params* input, int x, int start, int end){
	//  Estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=1.79E+308;
    int imin=-1;
    float temp;
	if(input->symmetric==1){
		for(i=0; i<input->k; i++){
			temp=distI(input, input->query_pq, x, i, start, end);
			if(temp<min){ 
				min=temp;
				imin=i;
			}
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

// ritorna indice del centroide c più vicino ad x.
int PQ_non_esaustiva(params* input, int x, int start, int end, int n_centroidi){
	// 	Estremi start incluso ed end escluso
    //	
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=1.79E+308;
    int imin=-1;
    float temp;
    for(i=0; i<n_centroidi; i++){
        temp=dist_eI(input, x, i, start, end);
        if(temp<min){ 
            min=temp;
            imin=i;
        }
    }
	
    return imin;
}

// source : Rappresenta la sorgente da cui calcolare il codebook (prodotto o vettoriale)
//			lo spazio deve essere già allocato a priori. 
// dest   : Rappresenta la destinazione dove dovranno essere inseriti i centroidi calcolati
void kmeans_from(params* input, int start, int end, int n_centroidi, float* source, float* dest, int dest_columns){
	
	// estremi start incluso ed end escluso

	int i, j, k, t;
	int count;
	float fob1, fob2;

	// residual_codebook = alloc_matrix(n_centroidi, input->nr); // row-major-order?
    // if(residual_codebook==NULL) exit(-1);
	// memset(residual_codebook,0,n_centroidi*input->nr); //azzera tutto il codebook

	
	//
	// Inizializzazione del codebook
	//		-Scelta dei k vettori casuali
	// 
	
    for(i=0; i<n_centroidi; i++){
		k=rand()%input->nr;
		for(j=start; j<end; j++){
			dest[i*dest_columns+j]=source[k*input->d+j];
		}
    }
	printf("--4--\n");

	// Assegnazione dei vettori ai centroidi casuali individuati
    for(i=0; i<input->nr; i++){
       dest[i*+(start/input->m)]=PQ_non_esaustiva(input, i, start, end, n_centroidi);
    }
	printf("--5--\n");
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && (fob2-fob1) > input->eps); t++){
		for(i=0; i<n_centroidi; i++){
			count=0; 
			memset(&dest[i*dest_columns+start], 0, end-start);
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			printf("--7--\n");
			for(j=0; j<input->nr; j++){
				printf("--8--\n");
				// pq[i][j]: centroide i-esimo del j-esimo sottogruppo
				if(input->pq[i*input->d+(start/input->m)]==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					printf("--9--\n");
					count++;
					for(k=start; k<end; k++){
						dest[i*input->d+k]+=sorg[j*input->d+k];
						printf("--11--\n");
					}
				}
			}
			
			for(j=start; j<end; j++){
				if(count!=0){ 
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					printf("--12--\n");
					dest[i*input->d+j]=dest[i*input->d+j]/count;
				}
			}
			
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
		}
		
		for(i=0; i<input->nr; i++){
			dest[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
		}
		
		fob1=fob2;
		fob2=0;
		
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(i=0; i<input->nr; i++){
			fob2+=pow(dist_eI(input, i, input->pq[i*input->m+(start/input->m)], start, end), 2.0);
		}
	}
}

void kmeans(params* input, int start, int end, int n_centroidi){
	// estremi start incluso ed end escluso
	int i, j, k, t;
	int count;
	float fob1, fob2;
	float* codebook;

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
    
    for(i=0; i<input->n; i++){
        input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
    }
    
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
		
		for(i=0; i<input->n; i++){
			input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
		}
		
		fob1=fob2;
		fob2=0;
		
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(i=0; i<input->n; i++){
			fob2+=pow(dist_eI(input, i, input->pq[i*input->m+(start/input->m)], start, end), 2.0);
		}
	}

	input->codebook=codebook;
}

void creaMatricedistanze(params* input){
	int i, j, k;
	MATRIX distanze_simmetriche;
	input->nDist=input->k*(input->k+1)/2;
	distanze_simmetriche = alloc_matrix(input->m, input->nDist);
	if(distanze_simmetriche==NULL) exit(-1);
	for(k=0; k<input->m; k++){
		for(i=1; i<input->k; i++){
			for(j=0; j<i; j++){
				distanze_simmetriche[k*input->nDist+calcolaIndice(i, j)] = dist_simmetricaI(input, i, j, k*input->m, (k+1)*input->m);
				// verificare se qui va usata la distanza simmetrica o no
			}
		}
	}
	input->distanze_simmetriche=distanze_simmetriche;
}

void bubbleSort(float* arr, int* arr2, int n){ 
   int i, j, temp; 
   for (i = 0; i < n-1; i++)    
       for (j = 0; j < n-i-1; j++)  
           if (arr[j] > arr[j+1]){
			   temp=arr[j];
			   arr[j]=arr[j+1];
			   arr[j+1]=temp;
			   temp=arr2[j];
			   arr2[j]=arr2[j+1];
			   arr2[j+1]=temp;
		   }
} 

void calcolaNN(params* input, int query){
	int i;
	VECTOR distanze=alloc_matrix(input->n, 1);
	int* d2=(int*) _mm_malloc(input->n*sizeof(int),16);
	if(input->symmetric==0){
		for(i=0; i<input->n; i++){
			distanze[i]=dist_asimmetrica(input, input->qs, query, i);
			d2[i]=i;
		}
	}else{
		for(i=0; i<input->n; i++){
			distanze[i]=dist(input, input->query_pq, query, i);
			d2[i]=i;
		}
	}
	bubbleSort(distanze, d2, input->n);
	dealloc_matrix(distanze);
	for(i=0; i<input->knn; i++){
		input->ANN[query*input->knn+i]=d2[i];
	}
	_mm_free(d2);
}

void inizializza_learning_set(params* input){
	//TODO: 
	//AL momento sceglie i primi nr come elementi del learning set. 
	input->residual_set = _mm_malloc(sizeof(float)*input->nr*input->d, 16);
	if(input->residual_set==NULL) exit(-1);
	
	//inizializza vettore
	input->qc_indexes = _mm_malloc(sizeof(unsigned char)*input->nr,16);
	if(input->qc_indexes==NULL) exit(-1);

}

// Ritorna il quantizzatore prodotto completo (con d dimensioni) del residuo r
VECTOR qp_of_r(params* input, int r){
	int qp_index, dStar;
	float* res;
	dStar = input->d/input->m;
	res = _mm_malloc(sizeof(float)*input->d, 16);
	for(int i=0;i<input->m;i++){
		qp_index = input->pq[r*input->d+i];
		for(int j=0;j<dStar;j++){
			res[i*input->m+j] = input->residual_codebook[qp_index*input->d+i*dStar+j];
		}
	}
	return res;
}

// Aggiunge a input.v la entry new alla posizione i-esima
void add (struct entry * new, int i, params* input){
	struct entry* vett;
	vett=input->v;
	if(vett[i].next== NULL){
		vett[i].next= new;
		new->next=NULL;
	}
	else{
		new->next = vett[i].next;
		vett[i].next = new;
	}
}

// Inizializza il vettore di entry v in modo tale da avere una lista di liste
// 
void inizializzaSecLiv(params* input){
	int qc_i;
	struct entry* new;
	input->v = _mm_malloc(sizeof(struct entry)*input->kc,16);
	if(input->v==NULL) return;
	for(int y= 0;y<input->nr;y++){
		qc_i = input->qc_indexes[y];
		new = _mm_malloc(sizeof(struct entry),16);
		if(new==NULL) exit(-1);
		new->index=y;
		new->q = qp_of_r(input, y);
		add(new,qc_i,input);
	}
}

float dist_coarse_and_residual(params* input, int qc, int y){
	// qc 		: indice del quantizzatore grossolano nel codebook in input
	// y	: puntatore al vettore residuo pari a r(y)=y-qc(y)
	//	
	//	<-------------------------------------------------------------->
	//	
	//	return -> distanza euclidea tra qc e residual, entrambi vettori a d coordinate
	int i; 
	float sum=0; //somma parziale
	for(i=0; i<input->m; i++){
		sum+=pow(input->codebook[qc*input->d+i]-input->ds[y*input->d+i], 2);
	}
	return sum;


}

// Calcola il centroide grossolano associato ad y.
int qc_index(params* input, int y){ 
	return input->qc_indexes[y];
}

void compute_residual(params* input, float* res, int qc_i, int y){
	// qc_i : corrisponde all' indice del quantizzatore grossolano nel codebook in input
	// y 	: indice del punto y appartenente al dataset ds in input
	//
	// -----------------------------------------
	// ritorna un puntatore al residuo r(y)
	for(int i=0; i<input->d;i++)
		res[i]=input->ds[y*input->d+i] - input->codebook[qc_i*input->d+i]; // r(y) = y - qc(y)
}

// Calcola tutti i residui dei vettori appartenenti al learning set
void calcola_residui(params* input){
	int qc_i; 
	float* ry; // puntatore al residuo corrente nel residual_codebook
	//ry = _mm_malloc(input->d*sizeof(float),16);
	for(int y=0;y<input->nr;y++){ // Per ogni y in Nr (learning-set):
		qc_i = qc_index(input,y); // Calcola il suo quantizzatore grossolano qc(y)
		ry = &input->residual_set[y*input->nr];
		compute_residual(input,ry,qc_i,y); // calcolo del residuo r(y) = y - qc(y)
	}

	
}


void pqnn_index_non_esaustiva(params* input){
	int i;
	float* tmp;
	printf("--1--\n");
	inizializza_learning_set(input);//selezionati i primi nr del dataset
	input->pq = (int*) _mm_malloc(input->nr*input->m*sizeof(int), 16);
	printf("--2--\n");
	tmp = input->residual_set;
	input->residual_set=input->qs;
	kmeans_from(input, 0, input->d, input->kc);//calcolo dei q. grossolani memorizzati messi in codebook
	input->residual_set=tmp; //scambio di puntatori per calcolare i centroidi grossolani dal learning set
	
	calcola_residui(input);
	//calcolo dei quantizzatori prodotto
	for(i=0;i<input->m;i++){
		kmeans_from(input, i*dStar, (i+1)*dStar, input->k);
	}
	inizializzaSecLiv(input);
}

void pqnn_search_non_esaustiva(params* input){

}