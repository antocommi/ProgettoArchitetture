/**************************************************************************************
 * 
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2018/19
 * 
 * Progetto dell'algoritmo di Product Quantization for Nearest Neighbor Search
 * in linguaggio assembly x86-32 + SSE
 * 
 * Fabrizio Angiulli, aprile 2019
 * 
 **************************************************************************************/

/*
 
 Software necessario per l'esecuzione:

     NASM (www.nasm.us)
     GCC (gcc.gnu.org)

 entrambi sono disponibili come pacchetti software 
 installabili mediante il packaging tool del sistema 
 operativo; per esempio, su Ubuntu, mediante i comandi:

     sudo apt-get install nasm
     sudo apt-get install gcc

 potrebbe essere necessario installare le seguenti librerie:

     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
     sudo apt-get install libc6-dev-i386

 Per generare il file eseguibile:

 nasm -f elf32 pqnn32.nasm && gcc -O0 -m32 -msse pqnn32.o pqnn32c.c -o pqnn32c && ./pqnn32c
 
 oppure
 
 ./runpqnn32

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <limits.h>
#include <float.h>

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

	MATRIX codebook; // per Ex contiene quantizzatori prodotto. Per NotEx contiene quantizzatori grossolani
	
	MATRIX distanze_simmetriche;
	int nDist;
	MATRIX distanze_asimmetriche;
	// ---------------------------------------
	// Strutture ad-hoc ricerca non esaustiva

	// Vettore contenente alla posizione i l'indice di qc(Y_i) in codebook
	int * qc_indexes;

	// Coarse q. dim: input.d x input.kc
	MATRIX qc;

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

struct kmeans_data{

	// Sorgente da cui si impara il codebook
	float* source;
	
	// Per ogni vettore contiene il centroide di appartenenza
	int* index; 
	
	// Per ogni riga contiene il centroide per intero
	float* dest;

	// Dimensione degli index
	int index_rows,index_colums;

	// Numero di centroidi da calcolare
	int n_centroidi;
};

/*
 * 
 *	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere 
 * 	memorizzate mediante array di array (float**).
 * 
 * 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
 * 	matrici per righe (row-major order) o per colonne (column major-order).
 *
 * 	L'assunzione corrente è che le matrici siano in row-major order.
 * 
 */

void stampa_matrice_flt(float* M, int rows, int col){
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<col;j++){
			printf(" %.2f",M[i*col+j]);
		}
		printf("\n---%d---\n",i);
	}
}

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


/*
 * 
 * 	load_data
 * 	=========
 * 
 *	Legge da file una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 * 
 * 	Codifica del file:
 * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
 * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
 * 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
 * 
 *****************************************************************************
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	della matrice. 
 *****************************************************************************
 * 
 */
MATRIX load_data(char* filename, int *n, int *d) {	
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
		
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;
	
	return data;
}


void save_ANN(char* filename, int* ANN, int nq, int knn) {	
	FILE* fp;
	int i, j;
	char fpath[256];
	
	sprintf(fpath, "%s.ann", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < nq; i++) {
		for (j = 0; j < knn; j++)
			fprintf(fp, "%d ", ANN[i*knn+j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}


extern void pqnn32_index(params* input);
extern int* pqnn32_search(params* input);

//funzioni fatte da noi

int calcolaIndice(int i, int j){
	//funzione che calcola l'indice per la matrice delle distanze_simmetriche
	return i*(i-1)/2+j;
}

float dist_eI(params* input, struct kmeans_data* data, int punto, int centroide, int start, int end){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	for(i=start; i<end; i++){
		ret += pow( data->source[punto*input->d+i]-data->dest[centroide*input->d+i], 2.0);
	}
	return ret;
}

/*




int dist_e(params* input, MATRIX set, int punto, int centroide){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_eI(input, set, punto, centroide, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}*/

int calcolaPQ(params* input, int x, int start, int end){
	/*
	// estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=1.79E+308;
    int imin=-1;
    float dist;
    for(i=0; i<input->k; i++){
        dist=dist_eI(input, input->ds, x, i, start, end);
        if(dist<min){ 
            min=dist;
            imin=i;
        }
    }
    return imin;
    */
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

float dist_asimmetricaI(params* input, MATRIX set, int punto, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	// centroide è un punto del dataset
	//
	// punto può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int i, c;
	float ret=0;
	for(i=start; i<end; i++){
		ret += pow( set[punto*input->d+i] - input->codebook[centroide2*input->d+i] , 2);
	}
	return ret;
}

float dist_asimmetrica(params* input, MATRIX set, int punto, int centroide){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_asimmetricaI(input, set, punto, input->pq[centroide], i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

float distI(params* input, int* quantizer, int punto, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	// centroide è un punto del dataset
	//
	// punto può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int c1=quantizer[punto*input->m+(start/input->m)];
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

float dist(params* input, int* quantizer, int punto, int centroide){
	int i;
	float sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(distI(input, quantizer, punto, input->pq[centroide], i*input->m, (i+1)*input->m), 2);
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
    float dist;
	//printf("breakpoint PQ\n");
	// valore 2 messo temporaneamente, inpossibile calcolare la distanza simmetrica per calcolare il centroide
	if(input->symmetric==2){
		for(i=0; i<input->k; i++){
			dist=distI(input, input->query_pq, x, i, start, end);
			if(dist<min){ 
				min=dist;
				imin=i;
			}	
			//printf("breakpoint PQ %d\n", i);
		}
	}else{
		for(i=0; i<input->k; i++){
			dist=dist_asimmetricaI(input, input->qs, x, i, start, end);
			if(dist<min){ 
				min=dist;
				imin=i;
			}
		}
	}
	return imin;
}


int PQ_non_esaustiva(params* input, int x, int start, int end, struct kmeans_data* data){
	// 	Estremi -> [start,end)
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=FLT_MAX;
    int imin=-1;
    float dist;
    for(i=0; i<data->n_centroidi; i++){
        dist=dist_eI(input, data, x, i, start, end);
        if(dist<min){ 
            min=dist;
            imin=i;
        }
    }
	// printf("trovato c=%d a distanza=%.2f di x=%d \n", imin, min, x);
    return imin;
}

// source : Rappresenta la sorgente da cui calcolare il codebook (prodotto o vettoriale)
//			lo spazio deve essere già allocato a priori. 
// dest   : Rappresenta la destinazione dove dovranno essere inseriti i centroidi calcolati
void kmeans_from(params* input, struct kmeans_data* data, int start, int end ){
	// estremi start incluso ed end escluso
	int i, j, k, t, c, imin;
	int count;
	float fob1, fob2;
	int* index;
	//
	// Inizializzazione del codebook
	//		-Scelta dei k vettori casuali
	//
	printf("\t --y--start=[%d]&end=[%d]\n",start,end);
    for(i=0; i<data->n_centroidi; i++){
		k = rand()%input->nr;
		for(j=start; j<end; j++){
			data->dest[i*input->d+j]=data->source[k*input->d+j];
		}
    }
	// Assegnazione dei vettori ai centroidi casuali individuati
	// stampa_matrice_flt(input->qc_indexes, input->nr, 1);
    for(i=0; i<data->index_rows; i++){
		if(start==16)
			printf("%d\n",PQ_non_esaustiva(input, i, start, end, data));
		data->index[i*data->index_colums+start/(data->index_colums)] = PQ_non_esaustiva(input, i, start, end, data);
    }
	printf("----------------\n");
	// stampa_matrice_flt(input->qc_indexes, input->nr, 1);
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && (fob2-fob1) > input->eps); t++){
		for(i=0; i<data->n_centroidi; i++){
			count=0; 
			memset(&data->dest[i*input->d+start], 0, (end-start)*sizeof(float));
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			for(j=0; j<input->nr; j++){
				if(data->index[j*data->index_colums+start]==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					for(k=start; k<end; k++){
						data->dest[i*input->d+k] += data->source[j*input->d+k];
					}
				}
			}
			for(j=start; j<end; j++){
				if(count!=0){ 
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					data->dest[i*input->d+j]=data->dest[i*input->d+j]/count;
				}
			}
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
		}
		for(i=0; i<input->nr; i++){
			// printf("prima: %d \n", data->index[i*data->index_colums+start]);
			if(start==16)
			printf("%d\n",PQ_non_esaustiva(input, i, start, end, data));
			data->index[i*data->index_colums+start/(data->index_colums)]=PQ_non_esaustiva(input, i, start, end, data);
			// printf("dopo: %d \n", data->index[i*data->index_colums+start]);
		}
		printf("----------------\n");
		fob1=fob2;
		fob2=0;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(i=0; i<input->nr; i++){
			fob2+=pow(dist_eI(input, data, i, data->index[i*data->index_colums+start/(data->index_colums)], start, end), 2.0);
		}
		// printf("delta=%.2f - %.2f - %.2f \n",fob2-fob1, fob2, fob1 );
	}
	printf("\t --x--\n");
}

void kmeans(params* input, int start, int end, int n_centroidi){

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
	int scambi=1; 
	for (i = 0; i < nit && scambi==1; i++){
		scambi=0;
    	for (j = n-2; j > i-1; j--)  
        	if (arr[j] > arr[j+1]){
				t2=arr[j];
				arr[j]=arr[j+1];
				arr[j+1]=t2;
				t1=arr2[j];
				arr2[j]=arr2[j+1];
				arr2[j+1]=t1;
				scambi=1;
			}
	}
}

void merge(VECTOR arr, int* arr2, int i, int j, int k){
	VECTOR a=(VECTOR) _mm_malloc(k-i+1 ,16);
	int* b=(int*) _mm_malloc(k-i+1 ,16);
	int c=0;
	int i2=i;
	int j1=j;
	while(i<j && j1<k){
		if(arr[i]<arr[j]){
			a[c]=arr[i];
			b[c]=arr2[i];
			i++;
		}else{
			a[c]=arr[j1];
			b[c]=arr2[j1];
			j1++;
		}
		c++;
	}
	while(i<j){
		a[c]=arr[i];
		b[c]=arr2[i];
		i++;
		c++;
	}
	while(j1<k){
		a[c]=arr[j1];
		b[c]=arr2[j1];
		j1++;
		c++;
	}
	c=i2;
	while(i2<k){
		arr[i2]=a[i2-c];
		arr2[i2]=b[i2-c];
		i2++;
	}
	_mm_free(a);
	_mm_free(b);
}

void mergesort(VECTOR arr, int* arr2, int i1, int i2){
	double t1;
	int t2;
	if(i2-i1<2) return;
	if(i2-i2==2){
		if(arr[i1]<=arr[i2]) return;
		t1=arr[i1];
		arr[i1]=arr[i2];
		arr[i2]=t1;
		t2=arr2[i1];
		arr2[i1]=arr2[i2];
		arr2[i2]=t2;
		return;
	}
	t2=(i2+i1)/2;
	mergesort(arr, arr2, i1, t2);
	mergesort(arr, arr2, t2, i2);
	merge(arr, arr2, i1, t2, i2);
}

void calcolaNN(params* input, int query){
	int i, j, k;
	VECTOR distanze=alloc_matrix(input->n, 1);
	VECTOR m;
	int* ind;

	if(input->knn<4500){
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
					//printf("%d %d", i, j);
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
		//mergesort(distanze, ind, 0, input->n);
		for(i=0; i<input->knn; i++){
			input->ANN[query*input->knn+i]=ind[i];
		}
		_mm_free(ind);
	}
	dealloc_matrix(distanze);
}



// Ritorna il quantizzatore prodotto completo (con d dimensioni) del residuo r
VECTOR qp_of_r(params* input, int r){
	int qp_index, dStar;
	float* res;
	dStar = input->d/input->m;
	res = _mm_malloc(sizeof(float)*input->d, 16);
	for(int i=0;i<input->m;i++){
		printf("\t1\n");
		qp_index = input->pq[r*input->d+i];
		printf("\t2\n");
		for(int j=0;j<dStar;j++){
			res[i*input->m+j] = input->residual_codebook[qp_index*input->d+i*dStar+j];
		}
		printf("\t3\n");
	}
	return res;
}

// Aggiunge a input.v la entry new alla posizione i-esima
void add (struct entry * new, int i, params* input){
	struct entry* vett;
	vett=input->v;
	printf("\t1\n");
	if(vett[i].next== NULL){
		vett[i].next= new;
		new->next=NULL;
		printf("\t2\n");
	}
	else{
		new->next = vett[i].next;
		vett[i].next = new;
		printf("\t2\n");
	}
}

// Inizializza il vettore di entry v in modo tale da avere una lista di liste
// 
void inizializzaSecLiv(params* input){
	int qc_i, y;
	struct entry* res;
	printf("\nxxxxx\n");
	input->v = malloc(sizeof(struct entry)*input->kc);
	if(input->v==NULL) exit(-1);
	printf("\nbbbbb\n");
	res = malloc(sizeof(struct entry)*input->nr);
	if(res==NULL) exit(-1);
	printf("\naaaaa\n");
	for(y=0;y<input->nr;y++){
		qc_i = input->qc_indexes[y];
		res[y].index=y;
		res[y].q = qp_of_r(input, y);
		add(&res[y],qc_i,input);
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
	int i;
	for(i=0; i<input->d;i++){
		res[i]=input->ds[y*input->d+i] - input->qc[qc_i*input->d+i]; // r(y) = y - qc(y)
	}
}

// Calcola tutti i residui dei vettori appartenenti al learning set
void calcola_residui(params* input){
	int qc_i; 
	float* ry; // puntatore al residuo corrente nel residual_codebook
	for(int y=0;y<input->nr;y++){ // Per ogni y in Nr (learning-set):
		qc_i = qc_index(input,y); // Calcola il suo quantizzatore grossolano qc(y)
		ry = &(input->residual_set[y*input->d]);
		compute_residual(input,ry,qc_i,y); // calcolo del residuo r(y) = y - qc(y)
	}
}

void pqnn_index_non_esaustiva(params* input){
	int i, dStar;
	float* tmp;
	struct kmeans_data* data;

	data = _mm_malloc(sizeof(struct kmeans_data),16);
	dStar = input->d/input->m;

	//TODO: 
	//AL momento sceglie i primi nr come elementi del learning set. 
	input->residual_set = (float*) _mm_malloc(sizeof(float)*input->nr*input->d, 16);
	if(input->residual_set==NULL) exit(-1);
	
	input->residual_codebook = (float*) _mm_malloc(sizeof(float)*input->k*input->d, 16);
	if(input->residual_codebook==NULL) exit(-1);

	input->qc = (float*) _mm_malloc(sizeof(float)*input->kc*input->d, 16);
	if(input->qc==NULL) exit(-1);

	//inizializza vettore
	input->qc_indexes = (int*) _mm_malloc(sizeof(int)*input->nr,16);
	if(input->qc_indexes==NULL) exit(-1);

	input->pq = (int*) _mm_malloc(input->nr*input->m*sizeof(int), 16);
	if(input->pq==NULL) exit(-1);


	printf("--2--\n");
	// Settagio parametri k-means
	data->source = input->ds;
	data->dest = input->qc;
	data->index = input->qc_indexes;
	data->index_colums=1;
	data->index_rows = input->nr;
	data->n_centroidi = input->kc;
	
	// printf("\n p:%p, val= %d\n",&(data->index[0]),data->index[0]);
	
	printf("--3--\n");	
	kmeans_from(input, data, 0, input->d); //calcolo dei q. grossolani memorizzati messi in codebook
	
	printf("--sono fuori dal kmeans #1 --\n");
	
	calcola_residui(input);

	printf("--eseguito calcola_residui --\n");
	// Settagio parametri k-means
	data->source = input->residual_set;
	data->dest = input->residual_codebook;
	data->index = input->pq;
	data->index_colums=input->m;
	data->index_rows = input->nr;
	data->n_centroidi = input->k;
	printf("struct data impostata m=%d \n", input->m);
	// calcolo dei quantizzatori prodotto
	for(i=0;i<input->m;i++){
		kmeans_from(input, data, i*dStar, (i+1)*dStar);
	}

	// for(i=0;i<input->k;i++){
	// 	for(int j=1;j<input->m;j+=2){
	// 		for(int k=0;k<dStar/4;k++){
	// 			printf(" %+3.2f ",input->residual_codebook[i*input->d+j*dStar+k]);
	// 		}
	// 		printf(" # ");
	// 	}
	// 	printf("\n");
	// }

	// for(i=0;i<input->nr;i++){
	// 	for(int j=0;j<input->m;j++){
	// 		printf("%2d ",input->pq[i*input->m+j]);
	// 	}
	// 	printf("\n");
	// }
	printf("Inizio: \"inizializzaSecLiv\" ");
   	inizializzaSecLiv(input);
	printf("Fine: \"inizializzaSecLiv\" ");
	
}

void pqnn_search_non_esaustiva(params* input){

}

void pqnn_index_esaustiva(params* input){
	int i, dStar;
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 16); 
	dStar=input->d/input->m;
	input->codebook = alloc_matrix(input->k, input->n); // row-major-order?
    if(input->codebook==NULL) exit(-1);
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
		//printf("break0\n");
		input->query_pq=(int*)_mm_malloc(input->nq*input->m*sizeof(int), 16);
		if(input->query_pq==NULL) exit(-1);
		c=input->d/input->m;
		//printf("break0.1\n");
		for(i=0; i<input->nq; i++){
			for(j=0; j<input->m; j++){
				input->query_pq[i*input->m+j]=calcolaQueryPQ(input, i, j*c, (j+1)*c);
			}
		}
	}
	//printf("break1\n");
	for(i=0; i<input->nq; i++){
		calcolaNN(input, i);
	}
	//printf("break2\n");
	_mm_free(input->codebook);
	_mm_free(input->pq);
	if(input->symmetric==1){
		_mm_free(input->distanze_simmetriche);
	}else{
		_mm_free(input->query_pq);
		_mm_free(input->distanze_asimmetriche);
	}
}

/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {
	// TODO: Gestire liberazione della memoria.
	if(input->exaustive==1){
		printf("ricerca esaustiva disattivata");
		//pqnn_index_esaustiva(input);
	}else{
		pqnn_index_non_esaustiva(input);
	}
    
    //pqnn32_index(input); // Chiamata funzione assembly

    // -------------------------------------------------
}


/*
 *	pqnn_search
 * 	===========
 */
void pqnn_search(params* input) {
	int i, j;
	if(input->exaustive==1){
		pqnn_search_esaustiva(input);
	}else{
		pqnn_search_non_esaustiva(input);
	}

    //pqnn32_search(input); // Chiamata funzione assembly

	// Restituisce il risultato come una matrice di nq * knn
	// identificatori associati agli ANN approssimati delle nq query.
	// La matrice è memorizzata per righe
    // -------------------------------------------------

}


int main(int argc, char** argv) {
	
	char fname[256];
	int i, j;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->filename = NULL;
	input->exaustive = 1;
	input->symmetric = 1;
	input->knn = 1;
	input->m = 8;
	input->k = 256;
	input->kc = 8192;
	input->w = 16;
	input->eps = 0.01;
	input->tmin = 10;
	input->tmax = 100;
	input->silent = 0;
	input->display = 0;

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while(par < argc) {
		if (par == 1) {
			input->filename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-knn") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing knn value!\n");
				exit(1);
			}
			input->knn = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-m") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing m value!\n");
				exit(1);
			}
			input->m = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-kc") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing kc value!\n");
				exit(1);
			}
			input->kc = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-w") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing w value!\n");
				exit(1);
			}
			input->w = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-nr") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing nr value!\n");
				exit(1);
			}
			input->nr = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-eps") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing eps value!\n");
				exit(1);
			}
			input->eps = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-tmin") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing tmin value!\n");
				exit(1);
			}
			input->tmin = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-tmax") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing tmax value!\n");
				exit(1);
			}
			input->tmax = atoi(argv[par]);
			par++;
 		} else if (strcmp(argv[par],"-exaustive") == 0) {
 			input->exaustive = 1;
 			par++;
 		} else if (strcmp(argv[par],"-noexaustive") == 0) {
 			input->exaustive = 0;
 			par++;
 		} else if (strcmp(argv[par],"-sdc") == 0) {
 			input->symmetric = 1;
 			par++;
 		} else if (strcmp(argv[par],"-adc") == 0) {
 			input->symmetric = 0;
 			par++;
		} else
			par++;
	}
	
	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if (!input->silent) {
		printf("Usage: %s <data_name> [-d][-s][-exaustive|-noexaustive][-sdc|-adc][...]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\t-d : display ANNs\n");
		printf("\t-s : silent\n");
		printf("\t-m: PQ groups\n");
		printf("\t-k: PQ centroids\n");
		printf("\t-kc: coarse VQ centroids\n");
		printf("\t-w: coarse VQ centroids to be selected\n");
		printf("\t-nr: residual sample size\n");
		printf("\t-eps: k-means termination threshold\n");
		printf("\t-tmin: min k-means iterations\n");
		printf("\t-tmax: max k-means iterations\n");
		printf("\n");
	}
	
	//
	// Legge il data set ed il query set
	//
	
	if (input->filename == NULL || strlen(input->filename) == 0) {
		printf("Missing input file name!\n");
		exit(1);
	}
	
	sprintf(fname, "%s.ds", input->filename);
	input->ds = load_data(fname, &input->n, &input->d);
	
	input->nr = input->n/20;

	sprintf(fname, "%s.qs", input->filename);
	input->qs = load_data(fname, &input->nq, &input->d);

	//
	// Visualizza il valore dei parametri
	//
	
	if (!input->silent) {
		printf("Input file name: '%s'\n", input->filename);
		printf("Data set size [n]: %d\n", input->n);
		printf("Number of dimensions [d]: %d\n", input->d);
		printf("Query set size [nq]: %d\n", input->nq);
		printf("Number of ANN [knn]: %d\n", input->knn);
		printf("PQ groups [m]: %d\n", input->m);
		printf("PQ centroids [k]: %d\n", input->k);
		if (!input->exaustive) {
			printf("Coarse VQ centroids [kc]: %d\n", input->kc);
			printf("Coarse VQ centroids to be selected [w]: %d\n", input->w);
			printf("Number of residuals used to determine PQ centroids [nr]: %d\n", input->nr);
		}
		printf("K-means parameters: eps = %.4f, tmin = %d, tmax = %d\n", input->eps, input->tmin, input->tmax);
	}
	
	//
	// Costruisce i quantizzatori
	//
	
	clock_t t = clock();
	pqnn_index(input);
	t = clock() - t;
	
	if (!input->silent)
		printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

	//
	// Determina gli ANN
	//
	
	input->ANN = calloc(input->nq*input->knn,sizeof(int));

	t = clock();
	pqnn_search(input);
	t = clock() - t;
	
	if (!input->silent)
		printf("\nSearching time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);
	
	//
	// Salva gli ANN
	//
	
 	if (input->ANN != NULL)
 	{
 		if (!input->silent && input->display) {
 			printf("\nANN:\n");
 			for (i = 0; i < input->nq; i++) {
				printf("query #%d:", i);
				for (j = 0; j < input->knn; j++)
					printf(" %d", input->ANN[i*input->knn+j]);
				printf("\n");
 			}
 		}
 		save_ANN(input->filename, input->ANN, input->nq, input->knn);
	}
	
	if (!input->silent)
		printf("\nDone.\n");

	return 0;
}
