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

int main(int argc, char** argv) {
    printf("filename index %d\n", offsetof(params, filename));
    printf("ds index %d\n", offsetof(params, ds));
    printf("qs index %d\n", offsetof(params, qs));
    printf("n index %d\n", offsetof(params, n));
    printf("d index %d\n", offsetof(params, d));
    printf("nq index %d\n", offsetof(params, nq));
    printf("knn index %d\n", offsetof(params, knn));
    printf("m index %d\n", offsetof(params, m));
    printf("k index %d\n", offsetof(params, k));
    printf("kc index %d\n", offsetof(params, kc));
    printf("w index %d\n", offsetof(params, w));
    printf("nr index %d\n", offsetof(params, nr));
    printf("eps index %d\n", offsetof(params, eps));
    printf("tmin index %d\n", offsetof(params, tmin));
    printf("tmax index %d\n", offsetof(params, tmax));
    printf("exaustive index %d\n", offsetof(params, exaustive));
    printf("symmetric index %d\n", offsetof(params, symmetric));
    printf("silent index %d\n", offsetof(params, silent));
    printf("display index %d\n", offsetof(params, display));
    printf("ANN index %d\n", offsetof(params, ANN));
    printf("vq index %d\n", offsetof(params, vq));
    printf("pq index %d\n", offsetof(params, pq));
    printf("query_pq index %d\n", offsetof(params, query_pq));
    printf("codebook index %d\n", offsetof(params, codebook));
    printf("distanze_simmetriche index %d\n", offsetof(params, distanze_simmetriche));
    printf("nDist index %d\n", offsetof(params, nDist));
    printf("distanze_asimmetriche index %d\n", offsetof(params, distanze_asimmetriche));
    printf("qc_indexes index %d\n", offsetof(params, qc_indexes));
    printf("residual_codebook index %d\n", offsetof(params, residual_codebook));
    printf("residual_set index %d\n", offsetof(params, residual_set));
}