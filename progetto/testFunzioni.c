
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
} params;

typedef struct{

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

//funzioni necessarie--------------------
float pow2(float f, float e){
	return f*f;
}
//---------------------------------------

extern int calcolaIndice(int i, int j);
int calcolaIndice2(int i, int j){
	return i*(i-1)/2+j;
}
void calcolaIndice_Test(){
	int a, b;
	for(int i=1; i<10; i++){
		for(int j=0; j<i; j++){
			a=calcolaIndice(i, j);
			b=calcolaIndice2(i, j);
			printf("%d %d\nnasm:%d\n   c:%d\n", i, j, a, b);
		}
	}
}

// extern void dist_eI(params* input, MATRIX set, int punto1, int punto2, int start, int end, float* r);
// void dist_eI2(params* input, MATRIX set, int punto1, int punto2, int start, int end, float* r){
// 	// estremi start incluso ed end escluso
// 	int i;
// 	float ret=0;
// 	float* ind=set+punto1*input->d+start;
// 	float* ind2=input->ds+punto2*input->d+start;
// 	for(i=start; i<end; i++){
// 		//printf("%f\n", pow2(*ind - *ind2, 2.0));
// 		ret+=pow2(*ind++ - *ind2++, 2.0);
// 		//ret+=(*ind++ - *ind2++);
// 	}
// 	*r=ret;
// }
// void dist_eI_Test(){
// 	int i;
// 	params* input=_mm_malloc(sizeof(params), 16);
// 	input->d=16;
// 	MATRIX set=_mm_malloc(2*(input->d)*sizeof(float), 16);
// 	srand((long)clock());
// 	for(i=0; i<input->d*2; i++){
// 		set[i]=100*(float)rand()/(float)RAND_MAX;
// 		printf("%d %f\n", i, set[i]);
// 	}

// 	input->ds=set;
// 	int punto1=0;
// 	int punto2=1;
// 	int start=0;
// 	int end=input->d;
// 	float r;
// 	float r2;
// 	//printf("Inizio c\n");
// 	dist_eI2(input, set, punto1, punto2, start, end, &r2);
// 	//printf("Fine c\n");
// 	printf("Inizio nasm\n");
// 	dist_eI(input, set, punto1, punto2, start, end, &r);
// 	printf("Fine nasm\n");
// 	printf("nasm:%f\n   c:%f\n", r, r2);
// 	_mm_free(set);
// 	_mm_free(input);
// }

//extern void dist_simmetricaI(params* input, int centroide1, int centroide2, int start, int end, float* r); 
// void dist_simmetricaI2(params* input, int centroide1, int centroide2, int start, int end, float* r){
// 	// estremi start incluso ed end escluso
// 	int i;
// 	float ret=0;
// 	float* ind=input->codebook+centroide1*input->d+start;
// 	float* ind2=input->codebook+centroide2*input->d+start;
// 	for(i=start; i<end; i++){
// 		ret+=pow2(*ind++ - *ind2++, 2);
// 	}
// 	*r=ret;
// }
// void dist_simmetricaI_Test(){
// 	params* input=_mm_malloc(sizeof(params), 16);
// 	input->d=8;
// 	input->codebook=_mm_malloc(2*input->d*sizeof(float), 16);
// 	for(int i=0; i<input->d*2; i++){
// 		input->codebook[i]=100*(float)rand()/(float)RAND_MAX;
// 		//printf("%f\n", set[i]);
// 	}
// 	int centroide1=0;
// 	int centroide2=1;
// 	int start=0;
// 	int end=input->d;
// 	float r;
// 	float r2;
// 	printf("Inizio c\n");
// 	dist_simmetricaI2(input, centroide1, centroide2, start, end, &r2);
// 	printf("Fine c\n");
// 	printf("Inizio nasm\n");
// 	//dist_simmetricaI(input, centroide1, centroide2, start, end, &r);
// 	printf("Fine nasm\n");
// 	printf("nasm:%f\n   c:%f\n", r, r2);
// 	_mm_free(input->codebook);
// 	_mm_free(input);
// }
// void printReg(){
// 	register int i asm("esi");
// 	printf("%d\n", i);
// 	register float f asm("xmm0");
// 	printf("%f\n", f);
// }

//extern void distanza(float* punto1, float* punto2, int dimensione, float* r);
//extern void calcolaFob(params* input, kmeans_data* data, int ipart, int start, int end, float* r);
// void calcolaFob1(params* input, kmeans_data* data, int ipart, int start, int end, float* r){
// 	int i;
// 	float* ind=data->dest+start;
// 	float* ind2=data->source+start;
// 	float ret=0;
// 	float temp;
// 	for(i=0; i<data->dim_source; i++){
// 		//distanza(input, input->codebook, input->pq[i*m+ipart], i, start, end, &temp);
// 		distanza(ind+data->index[i*input->m+ipart]*data->d, ind2, end-start, &temp);
// 		//printf("%f\n", temp);
// 		printf("%f\n", temp);
// 		ret+=pow2(temp, 2.0);
// 		//printf("%f\n", fob2);
// 		ind2+=data->d;
// 	}
// 	*r=ret;
// }
// void calcolaFob_Test(){
// 	params* input=_mm_malloc(sizeof(params), 16);
// 	kmeans_data* data=_mm_malloc(sizeof(kmeans_data), 16);
// 	data->d=16;
// 	input->m=4;
// 	int sx=10;
// 	int sy=data->d;
// 	int dx=10;
// 	int dy=data->d;
// 	int ix=sx;
// 	int iy=input->m;
// 	data->source=_mm_malloc(sx*sy*sizeof(double), 16);
// 	data->dest=_mm_malloc(dx*dy*sizeof(double), 16);
// 	data->dim_source=sx;
// 	data->index=_mm_malloc(ix*iy*sizeof(int), 16);
// 	printf("source\n");
// 	for(int i=0; i<sx; i++){
// 		for(int j=0; j<sy; j++){
// 			data->source[i*sy+j]=(float)rand()/(float)RAND_MAX;
// 			printf("%f ", data->source[i*sy+j]);
// 		}
// 		printf("\n");
// 	}
// 	printf("\ndest\n");
// 	for(int i=0; i<dx; i++){
// 		for(int j=0; j<dy; j++){
// 			data->dest[i*dy+j]=(float)rand()/(float)RAND_MAX;
// 			printf("%f ", data->dest[i*dy+j]);
// 		}
// 		printf("\n");
// 	}
// 	printf("\nindex\n");
// 	for(int i=0; i<ix; i++){
// 		for(int j=0; j<iy; j++){
// 			data->index[i*iy+j]=i;
// 			printf("%d ", data->index[i*iy+j]);
// 		}
// 		printf("\n");
// 	}
// 	int ipart=1;
// 	int dm=data->d/input->m;
// 	float t1;
// 	float t2;
// 	printf("prima c\n");
// 	calcolaFob1(input, data, ipart, ipart*dm, (ipart+1)*dm, &t1);
// 	printf("dopo c\nprima assembly\n");
// 	calcolaFob(input, data, ipart, ipart*dm, (ipart+1)*dm, &t2);
// 	printf("dopo assembly\n");
// 	printf("C:\t\t%f\nAssembly:\t%f\n", t1, t2);
// }

//extern void testdiv(int n);
int main(int argc, char** argv) {
	calcolaIndice_Test();
	//dist_eI_Test();
	//dist_simmetricaI_Test();
	//testdiv(2);
	//calcolaFob_Test();
	return 0;
}