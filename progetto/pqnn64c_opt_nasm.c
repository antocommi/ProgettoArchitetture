/**************************************************************************************
<<<<<<< HEAD
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

=======
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
 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
 Software necessario per l'esecuzione:

     NASM (www.nasm.us)
     GCC (gcc.gnu.org)

<<<<<<< HEAD
 entrambi sono disponibili come pacchetti software
 installabili mediante il packaging tool del sistema
=======
 entrambi sono disponibili come pacchetti software 
 installabili mediante il packaging tool del sistema 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
 operativo; per esempio, su Ubuntu, mediante i comandi:

     sudo apt-get install nasm
     sudo apt-get install gcc

 potrebbe essere necessario installare le seguenti librerie:

     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
     sudo apt-get install libc6-dev-i386

 Per generare il file eseguibile:

 nasm -f elf32 pqnn32.nasm && gcc -O0 -m32 -msse pqnn32.o pqnn32c.c -o pqnn32c && ./pqnn32c
<<<<<<< HEAD

 oppure

=======
 
 oppure
 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
 ./runpqnn32

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <limits.h>
#include <stddef.h>

#define	MATRIX		float*
#define	VECTOR		float*

#define DATASET		0
#define QUERYSET	1


typedef struct {
	char* filename; //
<<<<<<< HEAD
	MATRIX ds; // data set
=======
	MATRIX ds; // data set 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD
	float eps; //
=======
	float eps; // 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	MATRIX distanze_simmetriche;
	int nDist;
	MATRIX distanze_asimmetriche;
	// ---------------------------------------
	// Strutture ad-hoc ricerca non esaustiva

	// Vettore contenente alla posizione i l'indice di qc(Y_i) in codebook
	int * qc_indexes;

	MATRIX qc;

	// Residual product quantizators in non exhaustive search
	// matrix type: Row-major-order
	MATRIX residual_codebook;

	// Learning set nella ricerca non esastiva. Contiene i residui dei primi nr vettori
	MATRIX residual_set;

	// Lista di liste (secondo livello dell'inverted index)
<<<<<<< HEAD
	struct entry* v;
=======
	struct entry* v; 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61

	float* zero;
} params;

//Entry della s.d. multilivello
struct entry{
	int index;
	VECTOR q;
	//temporaneo
<<<<<<< HEAD
	//Serve per gestire liste a dimensione sconosciuta.
=======
	//Serve per gestire liste a dimensione sconosciuta. 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	struct entry * next;
};

typedef struct{

	// Sorgente da cui si impara il codebook
	float* source;

	int dim_source;
<<<<<<< HEAD

	// Per ogni vettore contiene il centroide di appartenenza
	int* index;

=======
	
	// Per ogni vettore contiene il centroide di appartenenza
	int* index; 
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	// Per ogni riga contiene il centroide per intero
	float* dest;

	// Dimensione degli index
	int index_rows,index_columns;

	// Numero di centroidi da calcolare
	int n_centroidi;

	int d;
} kmeans_data;



/*
<<<<<<< HEAD
 *
 *	Le funzioni sono state scritte assumento che le matrici siano memorizzate
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere
 * 	memorizzate mediante array di array (float**).
 *
=======
 * 
 *	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere 
 * 	memorizzate mediante array di array (float**).
 * 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
 * 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
 * 	matrici per righe (row-major order) o per colonne (column major-order).
 *
 * 	L'assunzione corrente è che le matrici siano in row-major order.
<<<<<<< HEAD
 *
 */


void* get_block(int size, int elements) {
	return _mm_malloc(elements*size,32);
}


void free_block(void* p) {
=======
 * 
 */


void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}


void free_block(void* p) { 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	_mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(float),rows*cols);
}


void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}


/*
<<<<<<< HEAD
 *
 * 	load_data
 * 	=========
 *
 *	Legge da file una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 *
=======
 * 
 * 	load_data
 * 	=========
 * 
 *	Legge da file una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 * 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
 * 	Codifica del file:
 * 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
 * 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
 * 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
<<<<<<< HEAD
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

=======
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
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	if (fp == NULL) {
		printf("'%s' : bad data file name!\n", filename);
		exit(0);
	}
<<<<<<< HEAD

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);

	*n = rows;
	*d = cols;

=======
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
		
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*d = cols;
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	return data;
}


<<<<<<< HEAD
void save_ANN(char* filename, int* ANN, int nq, int knn) {
	FILE* fp;
	int i, j;
	char fpath[256];

=======
void save_ANN(char* filename, int* ANN, int nq, int knn) {	
	FILE* fp;
	int i, j;
	char fpath[256];
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	sprintf(fpath, "%s.ann", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < nq; i++) {
		for (j = 0; j < knn; j++)
			fprintf(fp, "%d ", ANN[i*knn+j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}


//extern void pqnn32_index(params* input);
//extern int* pqnn32_search(params* input);

//funzioni fatte da noi

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

float dist_asimmetrica(params* input, MATRIX set, int punto1, int punto2){
	int i;
	float sum=0;
	float *ind1, *ind2;
	int dStar=input->d/input->m;
	int* c2=input->pq+punto2*input->m;
	ind1=set+punto1*input->d;
	ind2=input->codebook;
	for(i=0; i<input->m; i++){
		sum+=pow2(distanza(ind1, ind2+(*c2)*input->d, dStar), 2);
		ind1+=dStar;
		ind2+=dStar;
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
// 		sum+=pow2(*dist_matrix(input, *c1, *c2, i), 2);
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
<<<<<<< HEAD
// 			if(start>0)	printf("calcolapq %d %d\n", i, j);
// 			temp=distanza(ind1, ind2, end-start);
// 			if(start>0)	printf("calcolapq %d %d\n", i, j);
// 			if(temp<min){
=======
// 			temp=distanza(ind1, ind2, end-start);
// 			if(temp<min){ 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD
// 		ret+=pow2(distanza(ind+data->index[i*input->m+ipart]*data->d, ind2, end-start), 2.0);
=======
// 		ret+=pow2(distanza1(ind+data->index[i*input->m+ipart]*data->d, ind2, end-start), 2.0);
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
// 		ind2+=data->d;
// 	}
// 	return ret;
// }

extern void somma(float* source, float* dest, int dim);
<<<<<<< HEAD
void sommaC(float* source, float* dest, int dim){
	for(int i=0; i<dim; i++){
		*dest+=*source;
		dest++;
		source++;
	}
}
=======
// void somma(float* source, float* dest, int dim){
// 	for(int i=0; i<dim; i++){
// 		*dest+=*source;
// 		dest++;
// 		source++;
// 	}
// }
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61

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
	int m=input->m;
<<<<<<< HEAD
	int ipart=start/(input->d/m);
	//printf("prima calcolapq %ld %d %d %d\n", (long)data, ipart, start, end);
	//printf("%ld %ld\n", (long)data->source, (long)data->dest);
	calcolaPQ(data, ipart, start, end);
	//printf("dopo calcolapq\n");

	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && fabs(fob1-fob2)/fob1 > input->eps); t++){
		//printf("kmeans interation %d\n", t);
=======
	int ipart=start/(input->d/input->m);

	calcolaPQ(data, ipart, start, end);
	
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && fabs(fob2-fob1) > input->eps); t++){
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
		ci=data->dest+start;
		for(i=0; i<data->n_centroidi; i++){
			count=0;
			ind=ci;
			memset(ci, 0, (end-start)*sizeof(float));
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
<<<<<<< HEAD
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
					sommaC(ind, ci, end-start);
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
=======
			ind3=data->index+ipart;
			ind=data->source+start;
			for(j=0; j<data->dim_source; j++){
				if(*ind3==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					somma(ind, ci, end-start);
				}
				ind3+=m;
				ind+=input->d;
			}

			ind=ci;
			for(j=start; j<end; j++){
				if(count!=0){ 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					*ind=*ind/count;
				}
				ind++;
			}
<<<<<<< HEAD
			//printf("it:%d centr:%d 2\n", t, i);

=======
			
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
			ci+=input->d;
		}
<<<<<<< HEAD
		//printf("breakpoint\n");
		calcolaPQ(data, ipart, start, end);

		fob1=fob2;
		//fob11=fob22;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		fob2=calcolaFob(input, data, ipart, start, end);
		//printf("fob1:%f fob2:%f d:%f dn:%f", fob1, fob2, fabs(fob1-fob2), fabs(fob1-fob2)/fob1);
		//getchar();
	}
	printf("%d\n", t);
=======
		
		calcolaPQ(data, ipart, start, end);
		
		fob1=fob2;
		fob11=fob22;
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		fob2=calcolaFob(input, data, ipart, start, end);
	}
	//printf("%d\n", t);
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
}

extern void creaMatriceDistanze(params* input, float* codebook);
// void creaMatriceDistanze(params* input, float* codebook){
<<<<<<< HEAD
// 	// MODIFICATA SOLO CHIAMATA A FUNZIONE dist_simmetricaI(...) con aggiunta
=======
// 	// MODIFICATA SOLO CHIAMATA A FUNZIONE dist_simmetricaI(...) con aggiunta 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
// 	// puntatore alla src dei centroidi
// 	int i, j, k;
// 	int dStar=input->d/input->m;
// 	int d=input->d;
// 	float *ind1, *ind2;
// 	int count=0;
//
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

<<<<<<< HEAD
void bubbleSort(VECTOR arr, int* arr2, int n, int nit){
	int i, j, t1;
	float t2;
	int scambi=1;
	for (i = 0; i < nit && scambi==1; i++){
		scambi=0;
    	for (j = n-2; j > i-1; j--)
=======
void bubbleSort(VECTOR arr, int* arr2, int n, int nit){ 
	int i, j, t1;
	float t2;
	int scambi=1; 
	for (i = 0; i < nit && scambi==1; i++){
		scambi=0;
    	for (j = n-2; j > i-1; j--)  
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
	VECTOR a=(VECTOR) _mm_malloc(k-i+1 ,32);
	int* b=(int*) _mm_malloc(k-i+1 ,32);
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

typedef struct{
	float dist;
	int ind;
} vec;

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
	float* ind3;

	if(input->knn<4500){
		if(input->symmetric==0){
			for(i=0; i<input->n; i++){
				*ind++=dist_asimmetrica(input, input->qs, query, i);
			}
		}else{
			calcolaSimmetriche(input, distanze, query);
		}

		ind=m;
<<<<<<< HEAD

		for(i=0; i<input->knn; i++){
			*ind++=1.79E+308;
		}

=======
		
		for(i=0; i<input->knn; i++){
			*ind++=1.79E+308;
		}
		
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
		ind3=distanze;
		for(i=0; i<input->n; i++){
			ind=m;
			ind2=input->ANN+query*input->knn;
			for(j=0; j<input->knn; j++){
				if(*ind3<*ind){
					for(k=input->knn-1; k>j; k--){
						if(m[k-1]!=-1){
							input->ANN[query*input->knn+k]=input->ANN[query*input->knn+k-1];
							m[k]=m[k-1];
						}
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
		di=(int*) _mm_malloc(input->n*sizeof(int), 32);
		// //printf("breakpoint 1\n");
		// if(input->symmetric==0){
		// 	for(i=0; i<input->n; i++){
		// 		distanze[i]=dist_asimmetrica(input, input->qs, query, i);
		// 		di[i]=i;
		// 	}
		// }else{
		// 	for(i=0; i<input->n; i++){
		// 		distanze[i]=dist(input, input->query_pq, query, i);
		// 		di[i]=i;
		// 	}
		// }
		// bubbleSort(distanze, di, input->n, input->knn);
		// //mergesort(distanze, di, 0, input->n);
		// for(i=0; i<input->knn; i++){
		// 	input->ANN[query*input->knn+i]=di[i];
		// }
		_mm_free(di);
	}
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	dealloc_matrix(distanze);
}

void pqnn_index_esaustiva(params* input){
	int i, j, dStar;
	int d2=0;
	input->zero=_mm_malloc(sizeof(float), 32);
	*input->zero=0;
	float *ind1, *ind2;
<<<<<<< HEAD
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 32);
	if(input->pq==NULL) exit(-1);
=======
	input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 32); 
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	dStar=input->d/input->m;
	input->codebook = alloc_matrix(input->k, input->d); // row-major-order?
	if(input->codebook==NULL) exit(-1);
	ind1=input->codebook;
	ind2=input->ds;
	memcpy(input->codebook, input->ds, input->k*input->d*sizeof(float));
	kmeans_data* data=_mm_malloc(sizeof(kmeans_data), 32);
<<<<<<< HEAD
	if(data==NULL) exit(-1);
=======
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	data->source=input->ds;
	data->dim_source=input->n;
	data->index=input->pq;
	data->dest=input->codebook;
	data->index_rows=input->n;
	data->index_columns=input->m;
	data->n_centroidi=input->k;
	data->d=input->d;
<<<<<<< HEAD
	//printf("prima kmeans\n");
	for(i=0; i<input->m; i++){
		//printf("prima kmeans %d\n", i);
		kmeans(input, data, d2, d2+dStar);
		//printf("dopo kmeans %d\n", i);
		d2+=dStar;
	}
	//printf("fine kmeans\n");

=======
	for(i=0; i<input->m; i++){
		//printf("d2:%d d22:%d dStar:%d\n", d2, d2+dStar, dStar);
		kmeans(input, data, d2, d2+dStar);
		d2+=dStar;
	}
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61

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
	_mm_free(m);
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
		pqnn_index_esaustiva(input);
	}else{
		//pqnn_index_non_esaustiva(input);
	}
<<<<<<< HEAD

=======
    
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
		//pqnn_search_non_esaustiva(input);
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
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD

	//
	// Legge il data set ed il query set
	//

=======
	
	//
	// Legge il data set ed il query set
	//
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	if (input->filename == NULL || strlen(input->filename) == 0) {
		printf("Missing input file name!\n");
		exit(1);
	}
<<<<<<< HEAD

	sprintf(fname, "%s.ds", input->filename);
	input->ds = load_data(fname, &input->n, &input->d);

=======
	
	sprintf(fname, "%s.ds", input->filename);
	input->ds = load_data(fname, &input->n, &input->d);
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	input->nr = input->n/20;

	sprintf(fname, "%s.qs", input->filename);
	input->qs = load_data(fname, &input->nq, &input->d);

	//
	// Visualizza il valore dei parametri
	//
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	//
	// Costruisce i quantizzatori
	//
	clock_t t = clock();
	pqnn_index(input);
	t = clock() - t;
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	if (!input->silent)
		printf("\nIndexing time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

	//
	// Determina gli ANN
	//
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	input->ANN = calloc(input->nq*input->knn,sizeof(int));

	t = clock();
	pqnn_search(input);
	t = clock() - t;
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	if (!input->silent)
		printf("\nSearching time = %.3f secs\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);
<<<<<<< HEAD

	//
	// Salva gli ANN
	//

=======
	
	//
	// Salva gli ANN
	//
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
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
<<<<<<< HEAD

=======
	
>>>>>>> 7409d2b8a25360b6831cc8a233aae675d7771b61
	if (!input->silent)
		printf("\nDone.\n");

	return 0;
}