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
	return _mm_malloc(elements*size,32);
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

float pow2(float f, float s){
	return f*f;
	//return pow(f, 2);
}

extern int calcolaIndice(int i, int j);
// int calcolaIndice(int i, int j){
// 	//funzione che calcola l'indice per la matrice delle distanze_simmetriche
// 	return i*(i-1)/2+j;
// }

extern float distanza(float* punto1, float* punto2, int dimensione);
// float distanza(float* punto1, float* punto2, int dimensione){
// 	int i;
// 	float ret=0;
// 	float* ind=punto1;
// 	float* ind2=punto2;
// 	for(i=0; i<dimensione; i++){
// 		ret+=pow2(*ind++ - *ind2++, 2.0);
// 	}
// 	return ret;
// }

float dist_asimmetrica(params* input, int punto2){
	int i;
	float sum=0;
	int* c2=input->pq+punto2*input->m;
	for(i=0; i<input->m; i++){
		//sum+=distanza(ind1, ind2+(*c2)*input->d, dStar);
		sum+=input->distanze_asimmetriche[(*c2)*input->m+i];
		c2++;
	}
	return sum;
}

extern float* dist_matrix(params* input, int centroide1, int centroide2, int ipart);
// float* dist_matrix(params* input, int centroide1, int centroide2, int ipart){
// 	// estremi start incluso ed end escluso
// 	if(centroide1==centroide2){
// 		return input->zero;
// 	}else{
// 		//column major order-------------------------------------------
// 		if(centroide1<centroide2){
// 			return input->distanze_simmetriche+ipart+calcolaIndice(centroide2, centroide1)*input->m;
// 		}else{
// 			return input->distanze_simmetriche+ipart+calcolaIndice(centroide1, centroide2)*input->m;
// 		}
// 		//-------------------------------------------
// 	}
// }

extern float dist(params* input, int* quantizer, int punto1, int punto2);
// float dist(params* input, int* quantizer, int punto1, int punto2){
// 	int i;
// 	float sum=0;
// 	int* c1=quantizer+punto1*input->m;
// 	int* c2=input->pq+punto2*input->m;
// 	for(i=0; i<input->m; i++){
// 		sum+=*dist_matrix(input, *c1, *c2, i);
// 		c2++;
// 		c1++;
// 	}
// 	return sum;
// }

extern void calcolaPQ(kmeans_data* data, int partition, int start, int end);
// void calcolaPQ(kmeans_data* data, int partition, int start, int end){
// 	int i, j;
// 	int m=data->index_columns;
// 	float min;
// 	float temp;
// 	float *ind1, *ind2;
// 	int* ind=data->index+partition;
// 	ind1=data->source+start;
// 	for(i=0; i<data->dim_source; i++){
// 		min=1.79E+308;
// 		ind2=data->dest+start;
// 		for(j=0; j<data->n_centroidi; j++){
// 			if(start>0)	printf("calcolapq %d %d\n", i, j);
// 			temp=distanza(ind1, ind2, end-start);
// 			if(start>0)	printf("calcolapq %d %d\n", i, j);
// 			if(temp<min){
// 				min=temp;
// 				*ind=j;
// 			}
// 			ind2+=data->d;
// 		}
// 		ind+=m;
// 		ind1+=data->d;
// 	}
// }

extern float calcolaFob(params* input, kmeans_data* data, int ipart, int start, int end);
// float calcolaFob(params* input, kmeans_data* data, int ipart, int start, int end){
// 	int i;
// 	float* ind=data->dest+start;
// 	float* ind2=data->source+start;
// 	float ret=0;
// 	for(i=0; i<data->dim_source; i++){
// 		ret+=distanza(ind+data->index[i*input->m+ipart]*data->d, ind2, end-start);
// 		ind2+=data->d;
// 	}
// 	return ret;
// }

extern void somma(float* source, float* dest, int dim);
// void somma(float* source, float* dest, int dim){
// 	for(int i=0; i<dim; i++){
// 		*dest+=*source;
// 		dest++;
// 		source++;
// 	}
// }

void kmeans(params* input, kmeans_data* data, int start, int end){
	// estremi start incluso ed end escluso
	int i, j, k, t;
	int count;
	float fob1, fob2;
	float fob11, fob22;
	float temp;
	float *ind, *ind2, *ci;
	int* ind3;
	int incr, incr2;
	int m=data->index_columns;
	int ipart=start/(input->d/m);
	//printf("prima calcolapq %ld %d %d %d\n", (long)data, ipart, start, end);
	//printf("%ld %ld\n", (long)data->source, (long)data->dest);
	calcolaPQ(data, ipart, start, end);
	//printf("dopo calcolapq\n");

	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && fabs(fob1-fob2)/fob1 > input->eps); t++){
		//printf("kmeans interation %d\n", t);
		ci=data->dest+start;
		for(i=0; i<data->n_centroidi; i++){
			count=0;
			ind=ci;
			memset(ci, 0, (end-start)*sizeof(float));
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			//printf("it:%d centr:%d 0\n", t, i);
			ind3=data->index+ipart;
			//printf("it:%d centr:%d 0.0\n", t, i);
			ind=data->source+start;
			//printf("it:%d centr:%d 0.1\n", t, i);
			for(j=0; j<data->dim_source; j++){
				//if(start>0)
				//printf("it:%d centr:%d 0.1.1 %d\n", t, i, j);
				if(*ind3==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					//printf("prima somma\n");
					somma(ind, ci, end-start);
				}
				//if(start>0)
				//printf("it:%d centr:%d 0.2 %d\n", t, i, j);
				ind3+=m;
				ind+=input->d;
			}
			//printf("it:%d centr:%d 1\n", t, i);
			ind=ci;
			for(j=start; j<end; j++){
				if(count!=0){
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					*ind=*ind/count;
				}
				ind++;
			}
			//printf("it:%d centr:%d 2\n", t, i);

			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
			ci+=input->d;
		}
		//printf("breakpoint\n");
		calcolaPQ(data, ipart, start, end);

		fob1=fob2;
		//fob11=fob22;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		fob2=calcolaFob(input, data, ipart, start, end);
		//printf("fob1:%f fob2:%f d:%f dn:%f", fob1, fob2, fabs(fob1-fob2), fabs(fob1-fob2)/fob1);
		//getchar();
	}
	//printf("%d\n", t);
}

extern void creaMatriceDistanze(params* input, float* codebook);
// void creaMatriceDistanze(params* input, float* codebook){
// 	// MODIFICATA SOLO CHIAMATA A FUNZIONE dist_simmetricaI(...) con aggiunta
// 	// puntatore alla src dei centroidi
// 	int i, j, k;
// 	int dStar=input->d/input->m;
// 	int d=input->d;
// 	float *ind1, *ind2;
// 	int count=0;
// 	for(i=1; i<input->k; i++){
// 		ind2=codebook;
// 		for(j=0; j<i; j++){
// 			ind1=codebook+i*d;
// 			for(k=0; k<input->m; k++){
// 				input->distanze_simmetriche[count]=distanza(ind1, ind2, dStar);
// 				ind1+=dStar;
// 				ind2+=dStar;
// 				count++;
// 			}
// 		}
// 	}
// }

extern void calcolaSimmetriche(params* input, float* ind, int query);
// void calcolaSimmetriche(params* input, float* ind, int query){
// 	for(int i=0; i<input->n; i++){
// 		*ind=dist(input, input->query_pq, query, i);
// 		ind++;
// 	}
// }

void calcolaNN(params* input, int query, VECTOR m){
	int i, j, k;
	VECTOR distanze=alloc_matrix(input->n, 1);
	int* di;
	float* ind=distanze;
	int* ind2;
	float *ind3, *ind4, *ind5;
	int dStar=input->d/input->m;

	if(input->symmetric==0){ 
		input->distanze_asimmetriche=alloc_matrix(input->k, input->m);
		ind3=input->distanze_asimmetriche;
		ind4=input->codebook;
		for(i=0; i<input->k; i++){
			ind5=input->qs+query*input->d;
			for(j=0; j<input->m; j++){
				*ind3=distanza(ind4, ind5, dStar);
				ind3++;
				ind4+=dStar;
				ind5+=dStar;
			}
		}
		for(i=0; i<input->n; i++){
			*ind++=dist_asimmetrica(input, i);
		}
		_mm_free(input->distanze_asimmetriche);
	}else{
		calcolaSimmetriche(input, distanze, query);
	}

	if(input->knn>1){
		//knn>1
		ind=m;
		for(i=0; i<input->knn; i++){
			*ind++=1.79E+308;
		}
		ind3=distanze;
		for(i=0; i<input->n; i++){
			ind=m;
			ind2=input->ANN+query*input->knn;
			for(j=0; j<input->knn; j++){
				if(*ind3<*ind){
					for(k=input->knn-1; k>j; k--){
						if(m[k-1!=-1]) break;
					}
					for(k; k>j; k--){
						input->ANN[query*input->knn+k]=input->ANN[query*input->knn+k-1];
						m[k]=m[k-1];
					}
					*ind2=i;
					*ind=*ind3;
					break;
				}
				ind++;
				ind2++;
			}
			ind3++;
		}
	}else{
		//knn=1
		float min=1.79E+308;
		ind3=distanze;
		ind2=input->ANN+query;
		for(i=0; i<input->n; i++){
			if(*ind3<min){
				min=*ind3;
				*ind2=i;
			}
			ind3++;
		}
		//printf("%.2f\n", sqrtf(min));
	}
	dealloc_matrix(distanze);
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

extern void compute_residual(params* input, float* res, int qc_i, int y,float* src);

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
		compute_residual(input, ry, *qc_i++, y, input->ds); // calcolo del residuo r(y) = y - qc(y)
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
			*result=distanza(rx, ci, dStar);
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
		creaMatriceDistanze(input, input->residual_codebook);
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
			dist=distanza(q_x, input->qc + i*input->d, input->d); //distanza tra la query e il centroide grossolano
			insert(qc_heap, dist, i);
		}

		arr = qc_heap->arr;		

		qp_heap = CreateHeap(input->knn);
		//Ora in qc_heap ci sono i w centroidi grossolani più vicini. 
		
		for(i=0; i<input->w; i++){
			curr_qc = arr[i].index;
			// curr_qc = PopMaxIndex(qc_heap); 
			indice_curr_pq = input->index_entry[curr_qc];
			compute_residual(input, residuo, curr_qc, 0, q_x);
		
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

void pqnn_index_esaustiva(params* input){
	int i, j, dStar;
	int d2=0;
	input->zero=_mm_malloc(sizeof(float), 32);
	*input->zero=0;
	float *ind1, *ind2;
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 32);
	if(input->pq==NULL) exit(-1);
	dStar=input->d/input->m;
	input->codebook = alloc_matrix(input->k, input->d); // row-major-order?
	if(input->codebook==NULL) exit(-1);
	ind1=input->codebook;
	ind2=input->ds;
	memcpy(input->codebook, input->ds, input->k*input->d*sizeof(float));
	kmeans_data* data=_mm_malloc(sizeof(kmeans_data), 32);
	if(data==NULL) exit(-1);
	data->source=input->ds;
	data->dim_source=input->n;
	data->index=input->pq;
	data->dest=input->codebook;
	data->index_rows=input->n;
	data->index_columns=input->m;
	data->n_centroidi=input->k;
	data->d=input->d;
	//printf("prima kmeans\n");
	for(i=0; i<input->m; i++){
		//printf("prima kmeans %d\n", i);
		kmeans(input, data, d2, d2+dStar);
		//printf("dopo kmeans %d\n", i);
		d2+=dStar;
	}
	//printf("fine kmeans\n");

	if(input->symmetric==1){
		input->nDist=input->k*(input->k+1)/2;
		input->distanze_simmetriche = (float*) alloc_matrix(input->m, input->nDist);
		if(input->distanze_simmetriche==NULL) exit(-1);
		creaMatriceDistanze(input, input->codebook);
	}
	_mm_free(data);
}

void pqnn_search_esaustiva(params* input){
	int i, j, c, part;
	int *ipq, *ind;
	if(input->symmetric==1){
		input->query_pq=(int*)_mm_malloc(input->nq*input->m*sizeof(int), 32);
		if(input->query_pq==NULL) exit(-1);
		c=input->d/input->m;
		kmeans_data* data=_mm_malloc(sizeof(kmeans_data), 32);
		data->source=input->qs;
		data->dim_source=input->nq;
		data->index=input->query_pq;
		data->dest=input->codebook;
		data->index_rows=input->nq;
		data->index_columns=input->m;
		data->n_centroidi=input->k;
		data->d=input->d;
		part=0;
		for(j=0; j<input->m; j++){
			calcolaPQ(data, j, part, part+c);
			part+=c;
		}
		_mm_free(data);
	}
	VECTOR m=(VECTOR) _mm_malloc(input->knn*sizeof(float),32);
	for(i=0; i<input->nq; i++){
		calcolaNN(input, i, m);
	}
	//printf("prova 1\n");
	_mm_free(m);
	_mm_free(input->codebook);
	_mm_free(input->pq);
	if(input->symmetric==1){
		_mm_free(input->distanze_simmetriche);
	}else{
		_mm_free(input->query_pq);
	}
}

/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {
	// TODO: Gestire liberazione della memoria.
	if(input->exaustive==1){
		pqnn_index_esaustiva(input);
	}else{
		pqnn_index_non_esaustiva(input);
	}

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
	input->nr = 0;

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

	if (input->nr == 0){
		input->nr = input->n/20;
	}
    	
	

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