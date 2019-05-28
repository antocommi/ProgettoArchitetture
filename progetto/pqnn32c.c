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
sudo apt-get install libgtk2.0-dev

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

	// Coarse q. dim: input.kc x input.d
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
	// Residuo q di dimensione d.
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

void heapify_bottom_top(Heap *h,int index);

void heapify_top_bottom(Heap *h, int parent_node);

Heap* CreateHeap(int capacity){	
    Heap *h = (Heap * ) malloc(sizeof(Heap)); //one is number of heap	

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

void dist_eI(float* punto1, float* punto2, int dimensione, float* r){
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

// VECCHIA DIST
// float dist_eI(params* input, struct kmeans_data* data, int x, int y, int start, int end){
// 	// estremi start incluso ed end escluso
// 	int i;
// 	float ret=0, pow, tmp;
// 	for(i=start; i<end; i++){
// 		tmp = data->source[x*input->d+i]-data->dest[y*input->d+i];
// 		pow=tmp*tmp;
// 		ret =ret + pow;
// 	}
// 	return ret;
// }

void dist_asimmetrica_ne(params* input, VECTOR set1, VECTOR set2, int start, int end, float* r){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int i;
	float pow=0, ret=0, sott=0;
	float *ind=set1+start;
	float *ind2=set2+start;
	for(i=start; i<end; i++){
		sott=*ind - *ind2;
		pow=(sott)*(sott);
		ret += pow;
		ind++;  
		ind2++;
	}
	*r=ret;
}

float pow2(float f, float s){
	//non modificata rispetto a file pqnn32_opt1.c
	return f*f;
}

void dist_simmetricaI(params* input, float* src, int centroide1, int centroide2, int start, int end, float* r){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	float* ind=src+centroide1*input->d+start;
	float* ind2=src+centroide2*input->d+start;
	for(i=start; i<end; i++){
		ret+=pow2(*ind++ - *ind2++, 2);
	}
	*r=ret;
}

float* calcola_q(params* input, int query){
	// Ritorna un puntatore al quantizzatore esteso di q
	int dStar = input->d/input->m;
	struct kmeans_data* data;
	for(int i=0;i<input->m;i++){
		for(int c=0;c<input->k;c++){
			// todo
			// 
			// dist = dist_eI(input, data, x, y, i*dStar, (i+1)*dStar);
		}
	}
		
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
				dist_eI(ind1, ind2, dStar, &temp);
				// dist_simmetricaI(input, input->residual_codebook, i, j, k*dStar, (k+1)*dStar, &temp);
				input->distanze_simmetriche[k+calcolaIndice(i, j)*input->m]=temp;
				ind1+=dStar;
				ind2+=dStar;
			}
		}
	}
}

void dist_simmetrica_ne(params* input){
	
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


int PQ_non_esaustiva(params* input, int x, int start, int end, struct kmeans_data* data){
	// 	Estremi -> [start,end)
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    float min=FLT_MAX;
    int imin=-1, offset = x*input->d, dStar=input->d/input->m;
    float dist;
	float *src,*dest;
	src = data->source+offset+start;
    for(i=0; i<data->n_centroidi; i++){
        dist_eI(src,data->dest+i*input->d+start,dStar,&dist);
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
	float fob1, fob2, tmp_pow;
	int* index, dStar;
	float pow=0;
	dStar = input->d/input->m;
	srand( (time(NULL)) );
	//
	// Inizializzazione del codebook
	//		-Scelta dei k vettori casuali
	//
	// printf("\t --y--start=[%d]&end=[%d]\n",start,end);
    for(i=0; i<data->n_centroidi; i++){
		k = rand()% input->nr;
		//SCELTA DEI PRIMI n_centroidi
		// memcpy(&data->dest[i*input->d+start], &data->source[i*input->d+j], sizeof(float)*(end-start));
		for(j=start; j<end; j++){
			data->dest[i*input->d+j]=data->source[k*input->d+j];
		}
    }
	// Assegnazione dei vettori ai centroidi casuali individuati
	// stampa_matrice_flt(input->qc_indexes, input->nr, 1);
    for(i=0; i<data->index_rows; i++){
		data->index[i*data->index_colums+start/dStar] = PQ_non_esaustiva(input, i, start, end, data);
    } 
	// stampa_matrice_flt(input->qc_indexes, input->nr, 1);
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	
	for(t=0; !( input->tmin<=t && ( input->tmax<t || fob1-fob2 <= input->eps )); t++){
		for(i=0; i<data->n_centroidi; i++){
			count=0; 
			memset(&data->dest[i*input->d+start], 0, (end-start)*sizeof(float));
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//

			for(j=0; j<input->nr; j++){
				if(data->index[j*data->index_colums+start/dStar]==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
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
			data->index[i*data->index_colums+start/dStar]=PQ_non_esaustiva(input, i, start, end, data);
		}
		fob1=fob2;
		fob2=0;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(int y=0; y<input->nr; y++){
			c = data->index[y*data->index_colums+start/dStar];
			// tmp_pow = dist_eI(input, data, y, c, start, end);
			dist_eI(data->source+y*input->d+start,data->dest+c*input->d+start,end-start,&tmp_pow);
			pow = tmp_pow * tmp_pow;
			fob2+=pow;
		}
		
		// printf("------\nfob1=%.2f \nfob2=%.2f \n------\n", fob1, fob2 );
	}
	printf("# di iterazioni: %d\n",t);
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
	int qc_i, y;//qc_i è l'indice del quantizzatore grossolano associato ad y
	struct entry *res; //res è il residuo
	
	input->v = malloc(sizeof(struct entry)*input->kc);
	if(input->v==NULL) exit(-1);

	for(y=0;y<input->nr;y++){
		
		res = malloc(sizeof(struct entry));
		if(res==NULL) exit(-1);

		qc_i = input->qc_indexes[y];
		assert(qc_i>=0 && qc_i<input->kc);
		res->index=y;
		res->q = qp_of_r(input, y);// ritorna il quantizzatore del residuo
		add(res, qc_i, input); 
		assert(res->index>=0 && res->index<input->nr);
	}
	
	//
	// STAMPA DEL VETTORE DI LISTE V
	//
	// printf("---Inizio stampa lista a due livelli\n");
	// for(int i=0;i<input->kc;i++){
	// 	res = &(*(input->v+i));
	// 	while(res!=NULL){
	// 		printf("%d, ", res->index);
	// 		res=res->next;
	// 	}
	// 	printf("\n ");
	// }
	// printf("fine");
}

float dist_coarse_and_residual(params* input, int qc, int y){
	// qc 		: indice del quantizzatore grossolano nel codebook in input
	// y	: puntatore al vettore residuo pari a r(y)=y-qc(y)
	//	
	//	<-------------------------------------------------------------->
	//	
	//	return -> distanza euclidea tra qc e residual, entrambi vettori a d coordinate
	int i, offset1, offset2; 
	float sum=0, pow=0, tmp;
	offset2=y*input->d;
	offset1=qc*input->d;
	for(i=0; i<input->m; i++){
		tmp = input->codebook[offset1+i] - input->ds[offset2+i];
		pow=tmp*tmp;
		sum+=pow;
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
	
	// Settagio parametri k-means
	data->source = input->ds;
	data->dest = input->qc;
	data->index = input->qc_indexes;
	data->index_colums=1;
	data->index_rows = input->nr;
	data->n_centroidi = input->kc;
	kmeans_from(input, data, 0, input->d); //calcolo dei q. grossolani memorizzati messi in codebook
	calcola_residui(input);
	
	// Settagio parametri k-means
	data->source = input->residual_set;
	data->dest = input->residual_codebook;
	data->index = input->pq;
	data->index_colums=input->m;
	data->index_rows = input->nr;
	data->n_centroidi = input->k;
	
	// calcolo dei quantizzatori prodotto
	for(i=0;i<input->m;i++){
		kmeans_from(input, data, i*dStar, (i+1)*dStar);
	}

	inizializzaSecLiv(input);
	_mm_free(data);
	
}

void pqnn_search_non_esaustiva(params* input){
	int i, q;
	int curr_qc;
	struct entry* curr_pq;
	Heap* qc_heap, *qp_heap;
	struct kmeans_data* data;
	float dist;
	struct entry_heap* arr;
	float * residuo, *q_x; 
	float somma=0, temp;
	int dS=((input->d)/(input->m));

	residuo=_mm_malloc(sizeof(float)*input->d,16);
	if(residuo==NULL) exit(-1);

	data = _mm_malloc(sizeof(struct kmeans_data),16);
	if(data==NULL) exit(-1);
	
	data->source=input->qs;
	data->dest=input->qc;
	
	if(input->symmetric==0){
		creaMatricedistanze(input, input->residual_codebook);
	}
	

	//RIMETTERE IL VALORE INPUT->NQ È SOLO PER PROVA 
	for(int query=0; query<input->nq; query++){
		
		q_x = calcola_q(input, query);


		qc_heap = CreateHeap(input->w); //Creazione MAX-HEAP
		//potrei aggiungere un metodo restore su heap?
		


		for(int i=0;i<input->kc;i++){
			// dist = dist_eI(input, data, query, i, 0, input->d);
			dist_eI(input->qs+query*input->d, input->qc+i*input->d, input->d, &dist);
			insert(qc_heap,dist,i);
		}

		
		
		arr = qc_heap->arr;
		qp_heap = CreateHeap(input->knn);
		//Ora in qc_heap ci sono i w centroidi grossolani più vicini. 
		for(int i=0; i<input->w; i++){
			curr_qc = (qc_heap->arr)[i].index;
			curr_pq = ((input->v)[curr_qc]).next;
			// Calcolo r(x) rispetto al i-esimo centroide grossolano
			for(int j=0; j<input->d;j++){
					residuo[j]=input->qs[query*input->d+j] - input->qc[curr_qc*input->d+j]; // r(x) = y - qc(x)
			}
			while(curr_pq!=NULL){
				somma=0;
				for (int j=0; j<input->m; j++){
					if(input->symmetric==0){
						printf("simmetrica non implementata");
						// dist_simmetrica_ne();
					}
					else
						// Si può aggiungere caching delle dist tra pq e query? Usando pq ed un vettore di bit
						// Sarebbe una matrice k*m*sizeof(float) più vettore/matrice k*m*sizeof(bit)  
						dist_asimmetrica_ne(input, residuo, curr_pq->q, dS*j, dS*(j+1), &temp); 
					assert(temp>-1.0);
					somma += (temp*temp); 
				}
				insert(qp_heap, somma, curr_pq->index);
				curr_pq = curr_pq->next;
			}
		}
		// for(int k=0;k<input->knn;k++){
		// 		printf("vicino-%d-[%d]:%d\n",k+1,query,qp_heap->arr[k].index);
		// }
		// printf("fine0\n");
		_mm_free(qp_heap->arr);
		_mm_free(qc_heap->arr);
		_mm_free(qp_heap);
		_mm_free(qc_heap);
		// printf("fine1\n");
	}
	// printf("fine2\n");
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
		// pqnn_search_non_esaustiva(input);
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
	if(input->exaustive==0)
		pqnn_search_non_esaustiva(input);

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
