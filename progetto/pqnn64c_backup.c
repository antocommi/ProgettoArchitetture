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

#define	MATRIX		float*
#define	VECTOR		float*

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
	// struct entry *v; 

	float *zero;

	int* index_entry;
	int* celle_entry;
} params;

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

void heapify_bottom_top(Heap *h,int index);

void heapify_top_bottom(Heap *h, int parent_node);

Heap* CreateHeap(int capacity){	
    Heap *h = (Heap * ) _mm_malloc(sizeof(Heap),32); //one is number of heap	

     //check if memory allocation is fails	
    if(h == NULL) exit(-1);
    h->count=0;	
    h->capacity = capacity;	
    h->arr = _mm_malloc(capacity*sizeof(struct entry_heap),32); //size in bytes	

     //check if allocation succeed	
    if ( h->arr == NULL) exit(-1);	
    return h;	
}	

int PopMaxIndex(Heap *h){
    int pop_index;
	float pop_dist;
    if(h->count==0){
        printf("\n__Heap is Empty__\n");
        return -1;
    }
    // replace first node by last and delete last
    pop_index = h->arr[0].index;
	// pop_dist = h->arr[0].dist;
    h->arr[0].index = h->arr[h->count-1].index;
    h->arr[0].dist = h->arr[h->count-1].dist;
    h->count--;
    heapify_top_bottom(h, 0);
    return pop_index;
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



// Ritorna il quantizzatore prodotto completo (con d dimensioni) del residuo r
VECTOR qp_of_r(params* input, int r){
	int qp_index, dStar, partial_index, initial_offset;
	float* res;
	dStar = input->d/input->m;
	
	res = _mm_malloc(sizeof(float)*input->d, 32);
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


extern void compute_residual_opt(params* input, float* res, int qc_i, int y,float* src);

// Calcola tutti i residui dei vettori appartenenti al learning set
void calcola_residui(params* input){
	int *qc_i; 
	float *ry;
	ry = input->residual_set;
	qc_i=input->qc_indexes;
	// int *indexes=input->qc_indexes; // puntatore al residuo corrente nel residual_codebook
	for(int y=0;y<input->n;y++){ // Per ogni y in Nr (learning-set):
		// qc_i = input->qc_indexes[y]; // Calcola il suo quantizzatore grossolano qc(y)
		// OTTIMIZZABILE
		compute_residual_opt(input, ry, *qc_i++, y, input->ds); // calcolo del residuo r(y) = y - qc(y)
		ry += input->d;
	}

	
}

void pqnn_index_non_esaustiva(params* input){
	int i, j, l, d, m, knn, dStar, *offset, *offset2, n;
	int c, k, x, jk, *index, tmp;
	struct kmeans_data* data;

	k = input->k;
	d = input->d;
	knn = input->knn;
	m = input->m;
	n = input->n;

	input->zero = _mm_malloc(sizeof(float),32);
	if(input->zero==NULL) exit(-1);
	*(input->zero) = 0; 

	input->index_entry = _mm_malloc(sizeof(int)*input->kc,32);
	if(input->index_entry==NULL) exit(-1);
	memset(input->index_entry,0,input->kc*sizeof(int));

	input->celle_entry = _mm_malloc(sizeof(int)*input->n,32);
	if(input->celle_entry==NULL) exit(-1);

	data = _mm_malloc(sizeof(struct kmeans_data),32);
	dStar = input->d/input->m;

	//AL momento sceglie i primi nr come elementi del learning set. 
	input->residual_set = (float*) _mm_malloc(sizeof(float)*input->n*input->d, 32);
	if(input->residual_set==NULL) exit(-1);
	
	input->residual_codebook = (float*) _mm_malloc(sizeof(float)*input->k*input->d, 32);
	if(input->residual_codebook==NULL) exit(-1);

	input->qc = (float*) _mm_malloc(sizeof(float)*input->kc*input->d, 32);
	if(input->qc==NULL) exit(-1);

	//inizializza vettore
	input->qc_indexes = (int*) _mm_malloc(sizeof(int)*input->n,32);
	if(input->qc_indexes==NULL) exit(-1);

	// short ??
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 32);
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

	calcolaPQ(data,0, 0, input->d);
	
	// for(i=0;i<input->nr;i++){
	// 	printf("%d %d\n",i,input->qc_indexes[i]);
	// }

	calcola_residui(input);

	memcpy(input->residual_codebook, input->residual_set, input->k*input->d*sizeof(float));
	data->dim_source=input->nr;
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

	// // Aggiunta degli n-nr
	data->source = input->ds+input->nr*input->d;
	data->dest = input->residual_codebook;
	data->dim_source = input->n-input->nr;
	data->index = input->pq+input->nr*input->m;
	data->index_columns=input->m;
	data->index_rows = input->n-input->nr;
	data->n_centroidi = input->k;
	
	for(i=0;i<input->m;i++){
		calcolaPQ(data, i, i*dStar, (i+1)*dStar);
	}

	// TOLTO inizializzaSecLiv
	// Aggiunta dei punti r(y) del residual_set in celle_entry 
	// secondo index_entry[i] che indica l'inizio della i-esima cella 
	// sul vettore celle_entry.
	offset = _mm_malloc(sizeof(int)*input->kc,32);
	if(offset==NULL) exit(-1);
	memset(offset,0,input->kc*sizeof(int));

	for(i=0;i<n;i++){
		c = input->qc_indexes[i];
		input->index_entry[c]++;
	}

	x = input->index_entry[0];
	input->index_entry[0] = 0;
	for(l=1;l<input->kc;l++){
		tmp = input->index_entry[l];
		input->index_entry[l] = input->index_entry[l-1] + x;
		x = tmp;
	}

	for(i=0;i<n;i++){
		c = input->qc_indexes[i];
		l = input->index_entry[c] + offset[c]++;
		input->celle_entry[l] = i;
	}

	// printf("%d \n",input->index_entry[4]-input->index_entry[3] );
	// for(j=input->index_entry[3];j<input->index_entry[4];j++){
	// 	printf(" %d", input->celle_entry[j]);
	// 	if(j%8==0 && j!=0) printf("\n");
	// }

	_mm_free(input->residual_set);
	_mm_free(offset);
	_mm_free(data);
}

void creaMatricedistanzeAsimmetriche(params* input, float* residuo){
	int i, j, dStar;
	float *rx, *ci, *result, tmp;
	dStar = input->d/input->m;
	result = input->distanze_asimmetriche;
	rx = residuo;
	for(j=0;j<input->m;j++){
		ci = j*dStar + input->residual_codebook;
		for(i=0;i<input->k;i++){
			distanza(rx, ci, dStar, result);
			result++;
			ci += input->d;
		}
		rx += dStar;
	}

}

void pqnn_search_non_esaustiva(params* input){
	int i, q, query, j, p, h, s, residui_da_visitare, curr_residual, *ind_centroide;
	int curr_qc, indice_curr_pq, *pq_residuo, ci, cj;
	Heap* qc_heap, *qp_heap;
	struct kmeans_data* data;
	float dist;
	struct entry_heap* arr;
	float *residuo, *q_x, *dista; 
	float somma=0, temp;
	int dS=input->d/input->m;

	residuo= _mm_malloc(sizeof(float)*input->d,32);
	if(residuo==NULL) exit(-1);

	pq_residuo = _mm_malloc(sizeof(int)*input->m,32);
	if(pq_residuo==NULL) exit(-1);
	
	data = _mm_malloc(sizeof(struct kmeans_data),32);
	if(data==NULL) exit(-1);

	if(input->symmetric==1){
		creaMatricedistanze(input, input->residual_codebook);
		printf("\nSimmetrica\n");
	}else{
		input->distanze_asimmetriche = _mm_malloc(sizeof(float)*input->k*input->m,32);
		if(input->distanze_asimmetriche==NULL) exit(-1);
		printf("\nAsimmetrica\n");

	}

	// for(s=0;s<input->m;s++){
	// 	for(i=0;i<input->k;i++){
	// 		for(j=0;j<i;j++){
	// 			printf("(%d,%d)=%f ",i,j,input->distanze_simmetriche[s+calcolaIndice(i, j)*input->m]);
	// 		}
	// 	}
	// }
	
	// exit(-1);

	for(query=0; query<input->nq; query++){
		
		q_x = input->qs+query*input->d; //prende l indirizzo del vettore di query
		qc_heap = CreateHeap(input->w); //Creazione MAX-HEAP

		for(i=0;i<input->kc;i++){
			distanza(q_x, input->qc + i*input->d, input->d, &dist); //distanza tra la query e il centroide grossolano
			insert(qc_heap, dist, i);
		}

		arr = qc_heap->arr;		

		qp_heap = CreateHeap(input->knn);
		//Ora in qc_heap ci sono i w centroidi grossolani più vicini. 
		
		for(i=0; i<input->w; i++){
			curr_qc = arr[i].index;
			// curr_qc = PopMaxIndex(qc_heap); 
			indice_curr_pq = input->index_entry[curr_qc];
			compute_residual_opt(input, residuo, curr_qc, 0, q_x);
		
			if(input->symmetric==0){
				creaMatricedistanzeAsimmetriche(input,residuo);
				// for(h=0;h<input->k;h+=1){
				// 	printf("%.1f, ", input->distanze_asimmetriche[h]);
				// 	if(h != 0 && h%8==0) printf("\n");
				// }
			}else{
				data->source = residuo;
				data->d = input->d;
				data->dest = input->residual_codebook;
				data->index = pq_residuo;
				data->index_columns=input->m;
				data->index_rows = 1;
				data->n_centroidi = input->k;
				data->dim_source = 1;
				for(s=0;s<input->m;s++){
					calcolaPQ(data,s,s*dS,(s+1)*dS);
				}
			}

			if(curr_qc==input->kc-1) 
				residui_da_visitare=input->n;
			else 
				residui_da_visitare = input->index_entry[curr_qc+1];

			while(indice_curr_pq<residui_da_visitare){
				curr_residual = input->celle_entry[indice_curr_pq++];
				ind_centroide = input->pq+curr_residual*input->m;
				for(s=0;s<input->m;s++){
					if(input->symmetric==0){
						somma += input->distanze_asimmetriche[s*input->k+(*ind_centroide)];
						ind_centroide++;
					}else{
						ci = *ind_centroide++;
						cj = pq_residuo[s];
						somma += *(dist_matrix(input,ci,cj,s));
					}
				}
				insert(qp_heap, somma, curr_residual);
				somma=0;
			}
		}
		arr = qp_heap->arr;
		for(s=input->knn-1;s>=0;s--){
			// input->ANN[query*input->knn+s] = arr[s].index;
			// printf("%.2f ",sqrtf(qp_heap->arr[0].dist));
			input->ANN[query*input->knn+s] = PopMaxIndex(qp_heap);
			// printf("%.2d ", input->ANN[query*input->knn+s]);
		}
		// printf("\n");

		_mm_free(qp_heap->arr);
		_mm_free(qp_heap);
		_mm_free(qc_heap->arr);
		_mm_free(qc_heap);
	}
	_mm_free(residuo);
	_mm_free(data);
}

