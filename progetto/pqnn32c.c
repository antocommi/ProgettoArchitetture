/**************************************************************************************
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
 * DISATTIVARE ASSERT
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
#include <assert.h>

#define	MATRIX		float*
#define	VECTOR		float*
#define DATASET		0
#define QUERYSET	1
#define OFFSET 		0


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
	int* ANN; //dimensione: nq*knn
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
	int *qc_indexes;

	// Coarse q. dim: input.kc x input.d
	MATRIX qc;

	// Residual product quantizators in non exhaustive search
	// matrix type: Row-major-order
	MATRIX residual_codebook;

	// Learning set nella ricerca non esastiva. Contiene i residui dei primi nr vettori
	MATRIX residual_set;

	// Lista di liste (secondo livello dell'inverted index)
	struct entry *v; 

	float *zero;

	int* celle_voronoi;
	int* index_voronoi;
} params;

//Entry della s.d. multilivello
struct entry{
	int index;

	// Residuo q di dimensione d.
	VECTOR q;
	
	// VECTOR indexes;
	
	//temporaneo
	//Serve per gestire liste a dimensione sconosciuta. 
	struct entry * next;
};

typedef struct kmeans_data{

	// Sorgente da cui si impara il codebook
	float* source;

	int dim_source;
	
	// Per ogni vettore contiene il centroide di appartenenza
	int* index; 
	
	// Per ogni riga contiene il centroide per intero
	float* dest;

	// Dimensione degli index
	int index_rows,index_columns;

	// Numero di centroidi da calcolare
	int n_centroidi;

	int d;
	
} kmeans_data;

struct entry_heap{
	
	float dist;
	
	int index;
};

struct Heap{	
    struct entry_heap *arr;	
    int count;	
    int capacity;
};	

typedef struct Heap Heap;	

 // Metodi su heap 	
Heap* CreateHeap(int capacity);	

void insert(Heap *h, float key, int index);	

void addToVoronoi(int *celleVoronoi, int* posizioni, int* offset, int p, int k);

void heapify_bottom_top(Heap *h,int index);

void heapify_top_bottom(Heap *h, int parent_node);

Heap* CreateHeap(int capacity){	
    Heap *h = (Heap * ) _mm_malloc(sizeof(Heap),16); //one is number of heap	

     //check if memory allocation is fails	
    if(h == NULL) exit(-1);
    h->count=0;	
    h->capacity = capacity;	
    h->arr = _mm_malloc(capacity*sizeof(struct entry_heap),16); //size in bytes	

     //check if allocation succeed	
    if ( h->arr == NULL) exit(-1);	
    return h;	
}	

void insert(Heap *h, float key, int qc_index){	
    if( h->count < h->capacity){	
        (h->arr[h->count]).dist  = key;	
		(h->arr[h->count]).index = qc_index;
        heapify_bottom_top(h, h->count);	
        h->count++;	
    }else if( key < (h->arr[0]).dist ){
		(h->arr[0]).dist = key;
		(h->arr[0]).index = qc_index;
		heapify_top_bottom(h, 0);
	}
}

void heapify_bottom_top(Heap *h, int index){	
    int parent_node = (index-1)/2;
	float temp_dist;
	int temp_index;	

     if((h->arr[parent_node]).dist < (h->arr[index]).dist){	
        //swap and recursive call	
        temp_dist = (h->arr[parent_node]).dist;
		temp_index = (h->arr[parent_node]).index;

        // h->arr[parent_node] = h->arr[index];	
        (h->arr[parent_node]).dist = (h->arr[index]).dist;
		(h->arr[parent_node]).index = (h->arr[index]).index;
		
		// h->arr[index] = temp;	
        (h->arr[index]).dist = temp_dist;
		(h->arr[index]).index =temp_index;

		heapify_bottom_top(h,parent_node);	
    }	
}

void heapify_top_bottom(Heap *h, int parent_node){
    int left = parent_node*2+1;
    int right = parent_node*2+2;
    int min;
	float temp_dist;
	int temp_index;

	//forse si possono cacciare questi minori di zero. 
    if(left >= h->count || left <0)
        left = -1;
    if(right >= h->count || right <0)
        right = -1;

    if(left != -1 && (h->arr[left]).dist > (h->arr[parent_node]).dist )
        min=left;
    else
        min =parent_node;
    if(right != -1 && (h->arr[right]).dist > (h->arr[min]).dist )
        min = right;

    if(min != parent_node){
        // temp = h->arr[min];
        // h->arr[min] = h->arr[parent_node];
        // h->arr[parent_node] = temp;

        temp_dist = (h->arr[min]).dist;
		temp_index = (h->arr[min]).index;

        (h->arr[min]).dist = (h->arr[parent_node]).dist;
		(h->arr[min]).index = (h->arr[parent_node]).index;
		
        (h->arr[parent_node]).dist = temp_dist;
		(h->arr[parent_node]).index =temp_index;

        // recursive  call
        heapify_top_bottom(h, min);
    }
}


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

float pow2(float f, float s){
	//non modificata rispetto a file pqnn32_opt1.c
	return f*f;
}

void distanza(float* punto1, float* punto2, int dimensione, float* r){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	float* ind=punto1;
	float* ind2=punto2;
	for(i=0; i<dimensione; i++){
		ret+=pow2(*ind++ - *ind2++, 2.0);
	}
	*r=ret;
}

void creaMatricedistanze(params* input, float* codebook){
	// MODIFICATA SOLO CHIAMATA A FUNZIONE dist_simmetricaI(...) con aggiunta 
	// puntatore alla src dei centroidi
	int i, j, k;
	int dStar=input->d/input->m;
	int d=input->d;
	float temp;
	float *ind1, *ind2;
	input->nDist=input->k*(input->k+1)/2;
	
	input->distanze_simmetriche = alloc_matrix(input->m, input->nDist);
	if(input->distanze_simmetriche==NULL) exit(-1);
	for(i=1; i<input->k; i++){
		for(j=0; j<i; j++){
			ind1=codebook+i*d;
			ind2=codebook+j*d;
			for(k=0; k<input->m; k++){
				distanza(ind1, ind2, dStar, &temp);
				// dist_simmetricaI(input, input->residual_codebook, i, j, k*dStar, (k+1)*dStar, &temp);
				input->distanze_simmetriche[k+calcolaIndice(i, j)*input->m]=temp;
				ind1+=dStar;
				ind2+=dStar;
			}
		}
	}
}


void calcolaPQ(kmeans_data* data, int start, int end){
	// Dei primi input->k ed i primi input->kc già si conosce l'index
	// Per cui si può evitare di calcolare il più vicino.  
	int i, j;
	int m=data->index_columns;
	float min;
	float temp;
	float *ind1, *ind2;
	int* ind=data->index+start/((data->d)/(data->index_columns));
	ind1=data->source+start;
	for(i=0; i<data->dim_source; i++){// per ogni vettore del dataset
		min=FLT_MAX; //modificato
		ind2=data->dest+start;// destinazione
		for(j=0; j<data->n_centroidi; j++){// cerca il centroide + vicino
			distanza(ind1, ind2, end-start, &temp);//calcolando la distanza
			// printf("--2--\n");			
			if(temp<min){ 
				min=temp;
				*ind=j;
			}
			ind2+=data->d;
		}
		// printf("distanza minima trovata: %f\n",min );
		ind+=m;
		ind1+=data->d;
	}
}

void calcolaPQVoronoi(kmeans_data* data, int start, int end){
	// Dei primi input->k ed i primi input->kc già si conosce l'index
	// Per cui si può evitare di calcolare il più vicino.  
	int i, j;
	int m=data->index_columns;
	float min;
	float temp;
	float *ind1, *ind2;
	int* ind=data->index+start/((data->d)/(data->index_columns));
	ind1=data->source+start;
	for(i=0; i<data->dim_source; i++){// per ogni vettore del dataset
		min=FLT_MAX; //modificato
		ind2=data->dest+start;// destinazione
		for(j=0; j<data->n_centroidi; j++){// cerca il centroide + vicino
			distanza(ind1, ind2, end-start, &temp);//calcolando la distanza
			if(temp<min){ 
				min=temp;
				*ind=j;
			}
			ind2+=data->d;
		}
		// printf("distanza minima trovata: %f\n",min );
		ind+=m;
		ind1+=data->d;
	}
}

float absf(float f){
	if(f>0){
		return f;
	}
	return -f;
}

void stampaCentroidiGrossolani(params* input){
	int i;
	for(i=0;i<input->nr;i++){
		printf("%3d ",input->qc_indexes[i]);
		if(i%8==0 && i!=0) printf("\n");
		if(i==input->kc) printf("\n-------------------------\n");
	}
}

void kmeans(params* input, kmeans_data* data, int start, int end){
	// estremi start incluso ed end escluso
	int i, j, k, t;
	int count;
	float fob1, fob2;
	float temp;
	float *ind, *ind2, *ci, *distanze;
	int* ind3;
	int incr, incr2;
	int index_col=data->index_columns;
	int ipart=start/(input->d/input->m);
	calcolaPQ(data, start, end);
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && absf(fob1-fob2) > input->eps); t++){
		ci=data->dest+start;
		for(i=0; i<data->n_centroidi; i++){
			count=0;
			ind=ci;
			memset(ci, 0, (end-start)*sizeof(float));
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			ind3=data->index+ipart;
			for(j=0; j<data->dim_source; j++){
				if(*ind3==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					ind=ci;
					for(k=start; k<end; k++){
						*ind+=data->source[j*input->d+k];
						ind++;

					}
				}
				ind3+=index_col;
			}
			ind=ci;
			for(j=start; j<end; j++){
				if(count!=0){ 
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					*ind=*ind/count;
				}
				ind++;
			}
			
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
			ci+=input->d;
		}
		calcolaPQ(data, start, end);
		fob1=fob2;
		fob2=0;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		ind=data->dest+start;
		ind2=data->source+start;
		
		for(i=0; i<data->dim_source; i++){
			distanza(ind+data->index[i*index_col+ipart]*input->d, ind2, end-start, &temp);
			fob2+=pow2(temp, 2.0);
			ind2+=input->d;
		}
	}
}

// Ritorna il quantizzatore prodotto completo (con d dimensioni) del residuo r
VECTOR qp_of_r(params* input, int r){
	int qp_index, dStar, partial_index, initial_offset;
	float* res;
	dStar = input->d/input->m;
	
	res = _mm_malloc(sizeof(float)*input->d, 16);
	if(res==NULL) exit(-1);
	
	for(int i=0;i<input->m;i++){
		qp_index = input->pq[r*input->m+i];
		partial_index = qp_index*input->d+i*dStar;
		initial_offset = i*dStar;
		for(int j=0;j<dStar;j++){
			res[initial_offset+j] = input->residual_codebook[partial_index+j];
		}
	}
	return res;
}

// Aggiunge a input.v la entry new alla posizione i-esima
void add(struct entry * new, int i, params* input){
	struct entry* vett;
	vett=input->v;
	if(vett[i].next==NULL){
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
	int qc_i, y;//qc_i è l'indice del quantizzatore grossolano associato ad y
	struct entry *res; //res è il residuo
	
	input->v = _mm_malloc(sizeof(struct entry)*input->kc,16);
	if(input->v==NULL) exit(-1);

	for(y=0;y<input->nr;y++){
		
		res = _mm_malloc(sizeof(struct entry),16);
		if(res==NULL) exit(-1);

		qc_i = input->qc_indexes[y];
		assert(qc_i>=0 && qc_i<input->kc);
		res->index=y;
		res->next=NULL;
		res->q = qp_of_r(input, y);// ritorna il quantizzatore del residuo
		add(res, qc_i, input); 
		assert(res->index>=0 && res->index<input->nr);
	}
}


// Calcola il centroide grossolano associato ad y.
int qc_index(params* input, int y){ 
	return input->qc_indexes[y];
}

// extern void compute_residual_opt(params* input, float* res, int qc_i, int y,float* src);
void compute_residual(params* input, float* res, int qc_i, int y,float* src){
	// qc_i : corrisponde all' indice del quantizzatore grossolano nel codebook in input
	// y 	: indice del punto y appartenente al dataset ds in input
	// -----------------------------------------
	// ritorna un puntatore al residuo r(y)
	int i;
	for(i=0; i<input->d;i++){
		res[i]=src[y*input->d+i] - input->qc[qc_i*input->d+i]; // r(y) = y - qc(y)
	}
}

// Calcola tutti i residui dei vettori appartenenti al learning set
void calcola_residui(params* input){
	int qc_i; 
	float *ry;
	// int *indexes=input->qc_indexes; // puntatore al residuo corrente nel residual_codebook
	for(int y=0;y<input->n;y++){ // Per ogni y in Nr (learning-set):
		// qc_i = input->qc_indexes[y]; // Calcola il suo quantizzatore grossolano qc(y)
		ry = input->residual_set+y*input->d;
		
		compute_residual_opt(input, ry, *(input->qc_indexes+y), y, input->ds); // calcolo del residuo r(y) = y - qc(y)
	}
}

void pqnn_index_non_esaustiva(params* input){
	int i, j, d, m, knn, dStar, *offset;
	float* tmp;
	struct kmeans_data* data;
	// unsigned short *voronoi_cells, *count, *centroidi_assegnati;


	d = input->d;
	knn = input->knn;
	m = input->m;
	
	offset = _mm_malloc(sizeof(int)*input->k,16);
	if(offset==NULL) exit(-1);

	input->celle_voronoi = _mm_malloc(sizeof(int)*input->n*input->m,16);
	if(input->celle_voronoi==NULL) exit(-1);
	
	input->index_voronoi = _mm_malloc(sizeof(int)*input->k*input->m,16);
	if(input->index_voronoi==NULL) exit(-1);

	data = _mm_malloc(sizeof(struct kmeans_data),16);
	dStar = input->d/input->m;

	//AL momento sceglie i primi nr come elementi del learning set. 
	input->residual_set = (float*) _mm_malloc(sizeof(float)*input->n*input->d, 16);
	if(input->residual_set==NULL) exit(-1);
	
	input->residual_codebook = (float*) _mm_malloc(sizeof(float)*input->k*input->d, 16);
	if(input->residual_codebook==NULL) exit(-1);

	input->qc = (float*) _mm_malloc(sizeof(float)*input->kc*input->d, 16);
	if(input->qc==NULL) exit(-1);

	//inizializza vettore
	input->qc_indexes = (int*) _mm_malloc(sizeof(int)*input->n,16);
	if(input->qc_indexes==NULL) exit(-1);

	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 16);
	if(input->pq==NULL) exit(-1);
	
	//nuova aggiunta
	memcpy(input->qc, input->ds, input->kc*input->d*sizeof(float));
	
	// Settagio parametri k-means
	data->source = input->ds;
	data->dest = input->qc;
	data->index = input->qc_indexes;
	data->index_columns=1;
	data->index_rows = input->nr;
	data->n_centroidi = input->kc;
	data->d=input->d; 
	data->dim_source = input->nr;

	kmeans(input, data, 0, input->d); //calcolo dei q. grossolani memorizzati messi in codebook
	
	// Settagio parametri k-means
	data->source = & input->ds[(input->nr)*input->d];
	data->dest = input->qc;
	data->index = &input->qc_indexes[input->nr];
	data->index_columns=1;
	data->index_rows = input->n-input->nr;
	data->n_centroidi = input->kc;
	data->dim_source = input->n-input->nr;
	
	calcolaPQ(data, 0, input->d);
	
	calcola_residui(input);
	//nuova aggiunta
	memcpy(input->residual_codebook, input->residual_set, input->k*input->d*sizeof(float));
	data->dim_source=input->nr;
	// Settagio parametri k-means
	data->source = input->residual_set;
	data->dest = input->residual_codebook;
	data->index = input->pq;
	data->index_columns=input->m;
	data->index_rows =input->nr;
	data->n_centroidi = input->k;

	// calcolo dei quantizzatori prodotto
	for(i=0;i<input->m;i++){
		kmeans(input, data, i*dStar, (i+1)*dStar);
	}

	inizializzaSecLiv(input);

	// // Aggiunta degli n-nr
	data->source = input->ds+input->nr*input->d;
	data->dest = input->residual_codebook;
	data->dim_source = input->n-input->nr;
	data->index = input->pq+input->nr*input->m;
	data->index_columns=input->m;
	data->index_rows = input->n-input->nr;
	data->n_centroidi = input->k;
	
	for(i=0;i<input->m;i++){
		calcolaPQ(data, i*dStar, (i+1)*dStar);
	}
	
	for(j=0;j<m;j++){
		pq = input->pq + j;

		for(i=0;i<input->n;i++){
			pq += 
			pq += input->m;
		}
		calcola_posizioni(input->index_voronoi+j*input->k);

		memset(offset,0,input->k*sizeof(int));
	}
	

	addToVoronoi();
	_mm_free(data);
}



void pqnn_search_non_esaustiva(params* input){
	int i, q, p, h, s;
	int curr_qc;
	struct entry* curr_pq;
	Heap* qc_heap, *qp_heap;
	struct kmeans_data* data;
	float dist;
	struct entry_heap* arr;
	float * residuo, *q_x; 
	float somma=0, temp;
	int dS=((input->d)/(input->m)), *qp_query;
	qp_query = _mm_malloc(sizeof(int)*input->m,16);
	if(qp_query==NULL) exit(-1);
	residuo= _mm_malloc(sizeof(float)*input->d,16);
	if(residuo==NULL) exit(-1);
	data = _mm_malloc(sizeof(struct kmeans_data),16);
	if(data==NULL) exit(-1);
	data->source=input->qs;
	data->dest=input->qc;

	if(input->symmetric==1){
		creaMatricedistanze(input, input->residual_codebook);
	}

	for(int query=0; query<input->nq; query++){
		
		q_x = &input->qs[query*input->d]; //prende l indirizzo del vettore di query
		qc_heap = CreateHeap(input->w); //Creazione MAX-HEAP

		//potrei aggiungere un metodo restore su heap?
		for(int i=0;i<input->kc;i++){
			distanza(q_x, &input->qc[i*input->d], input->d, &dist); //distanza tra la query e il centroide grossolano
			insert(qc_heap, dist, i);
		}
		// for(int i=0;i<input->w;i++){
		// 	printf("distanza %.2f da %d \n", qc_heap->arr[i].dist,qc_heap->arr[i].index);
		// }

		arr = qc_heap->arr;

		qp_heap = CreateHeap(input->knn);
		//Ora in qc_heap ci sono i w centroidi grossolani più vicini. 
		
		for(int i=0; i<qc_heap->count; i++){
			curr_qc = (qc_heap->arr)[i].index;
			curr_pq = ((input->v)[curr_qc]).next;

			assert(curr_qc>-1 && curr_qc<input->kc);
			compute_residual(input, residuo, curr_qc, 0, input->qs);

			if(input->symmetric==1){
				data->source=residuo;
				data->dim_source=1;
				data->dest=input->residual_codebook;
				data->n_centroidi=input->k;
				data->index = qp_query;
				data->d = input->d;
				data->index_columns = input->m;
				data->index_rows = 1;
				for(int l=0;l<input->m;l++){
					calcolaPQ(data, l*dS, (l+1)*dS);
				}
			}	
			while(curr_pq!=NULL){
				somma=0;
				for (int j=0; j<input->m; j++){
					if(input->symmetric==1){
						if(qp_query[j] > input->pq[input->m*curr_pq->index+j]){
							h = qp_query[j];
							p = input->pq[input->m*curr_pq->index+j];
						}else{
							h = input->pq[input->m*curr_pq->index+j];
							p = qp_query[j];
						}
						if(h!=p){
							assert(j+calcolaIndice(h,p)*input->m>-1 && j+calcolaIndice(h,p)*input->m<input->m*input->k*(input->k-1)/2);
							temp = input->distanze_simmetriche[j+calcolaIndice(h,p)*input->m];
						}
						else
							temp=0;
					}
					else{
						// Si può aggiungere caching delle dist tra pq e query? Usando pq ed un vettore di bit
						// Sarebbe una matrice k*m*sizeof(float) più vettore/matrice k*m*sizeof(bit)  
						// dist_asimmetrica_ne(input, residuo, curr_pq->q, dS*j, dS*(j+1), &temp);
						distanza(residuo+dS*j,curr_pq->q+dS*j,dS,&temp);
					} 
					assert(temp>-1.0);
					somma += (temp*temp); 
				}
				insert(qp_heap, somma, curr_pq->index);
				curr_pq = curr_pq->next;
			}
		}

		
		
		//A questo punto i knn vicini sono in qp_heap->arr
		
		for(s=0;s<1;s++){
			input->ANN[query*input->knn+s] = (qp_heap->arr)[s].index;
		}
		

		_mm_free(qp_heap->arr);
		_mm_free(qc_heap->arr);
		_mm_free(qp_heap);
		_mm_free(qc_heap);
	}
	_mm_free(data);
}



/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {
	// TODO: Gestire liberazione della memoria.
	if(input->exaustive==1){
		printf("ricerca esaustiva disattivata");
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
	// int i, j;
	if(input->exaustive==0){
		// printf("simmetrica:%d\n",input->symmetric);
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
	input->nr = input->n/20;
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
