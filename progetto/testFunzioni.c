
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

//funzioni necessarie--------------------
float pow2(float f, float e){
	return f*f;
}
//---------------------------------------

// extern int calcolaIndice(int i, int j);
int calcolaIndice2(int i, int j){
	return i*(i-1)/2+j;
}
void calcolaIndice_Test(){
	int a, b;
	for(int i=1; i<10; i++){
		for(int j=0; j<i; j++){
			//a=calcolaIndice(i, j);
			b=calcolaIndice2(i, j);
			printf("%d %d\nnasm:%d\n   c:%d\n", i, j, a, b);
		}
	}
}

extern void dist_eI(params* input, MATRIX set, int punto1, int punto2, int start, int end, float* r);
void dist_eI2(params* input, MATRIX set, int punto1, int punto2, int start, int end, float* r){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	float* ind=set+punto1*input->d+start;
	float* ind2=input->ds+punto2*input->d+start;
	for(i=start; i<end; i++){
		//printf("%f\n", pow2(*ind - *ind2, 2.0));
		ret+=pow2(*ind++ - *ind2++, 2.0);
		//ret+=(*ind++ - *ind2++);
	}
	*r=ret;
}
void dist_eI_Test(){
	int i;
	params* input=_mm_malloc(sizeof(params), 16);
	input->d=16;
	MATRIX set=_mm_malloc(2*(input->d)*sizeof(float), 16);
	srand((long)clock());
	for(i=0; i<input->d*2; i++){
		set[i]=100*(float)rand()/(float)RAND_MAX;
		printf("%d %f\n", i, set[i]);
	}

	input->ds=set;
	int punto1=0;
	int punto2=1;
	int start=0;
	int end=input->d;
	float r;
	float r2;
	//printf("Inizio c\n");
	dist_eI2(input, set, punto1, punto2, start, end, &r2);
	//printf("Fine c\n");
	printf("Inizio nasm\n");
	dist_eI(input, set, punto1, punto2, start, end, &r);
	printf("Fine nasm\n");
	printf("nasm:%f\n   c:%f\n", r, r2);
	_mm_free(set);
	_mm_free(input);
}

//extern void dist_simmetricaI(params* input, int centroide1, int centroide2, int start, int end, float* r); 
void dist_simmetricaI2(params* input, int centroide1, int centroide2, int start, int end, float* r){
	// estremi start incluso ed end escluso
	int i;
	float ret=0;
	float* ind=input->codebook+centroide1*input->d+start;
	float* ind2=input->codebook+centroide2*input->d+start;
	for(i=start; i<end; i++){
		ret+=pow2(*ind++ - *ind2++, 2);
	}
	*r=ret;
}
void dist_simmetricaI_Test(){
	params* input=_mm_malloc(sizeof(params), 16);
	input->d=8;
	input->codebook=_mm_malloc(2*input->d*sizeof(float), 16);
	for(int i=0; i<input->d*2; i++){
		input->codebook[i]=100*(float)rand()/(float)RAND_MAX;
		//printf("%f\n", set[i]);
	}
	int centroide1=0;
	int centroide2=1;
	int start=0;
	int end=input->d;
	float r;
	float r2;
	printf("Inizio c\n");
	dist_simmetricaI2(input, centroide1, centroide2, start, end, &r2);
	printf("Fine c\n");
	printf("Inizio nasm\n");
	//dist_simmetricaI(input, centroide1, centroide2, start, end, &r);
	printf("Fine nasm\n");
	printf("nasm:%f\n   c:%f\n", r, r2);
	_mm_free(input->codebook);
	_mm_free(input);
}

void printReg(){
	register int i asm("esi");
	printf("%d\n", i);
	register float f asm("xmm0");
	printf("%f\n", f);
}

int main(int argc, char** argv) {
	//calcolaIndice_Test();
	dist_eI_Test();
	//dist_simmetricaI_Test();
	return 0;
}