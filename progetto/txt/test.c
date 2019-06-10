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
	int * qc_indexes;

	MATRIX qc;

	// Residual product quantizators in non exhaustive search
	// matrix type: Row-major-order
	MATRIX residual_codebook;

	// Learning set nella ricerca non esastiva. Contiene i residui dei primi nr vettori
	MATRIX residual_set;

	// Lista di liste (secondo livello dell'inverted index)
	struct entry* v; 

	float* zero;
} params;

struct entry{
	int index;
	VECTOR q;
	//temporaneo
	//Serve per gestire liste a dimensione sconosciuta. 
	struct entry * next;
};

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

int main(int argc, char** argv) {
    printf("filename equ %ld\n", offsetof(params, filename));
    printf("ds equ %ld\n", offsetof(params, ds));
    printf("qs equ %ld\n", offsetof(params, qs));
    printf("n equ %ld\n", offsetof(params, n));
    printf("d equ %ld\n", offsetof(params, d));
    printf("nq equ %ld\n", offsetof(params, nq));
    printf("knn equ %ld\n", offsetof(params, knn));
    printf("m equ %ld\n", offsetof(params, m));
    printf("k equ %ld\n", offsetof(params, k));
    printf("kc equ %ld\n", offsetof(params, kc));
    printf("w equ %ld\n", offsetof(params, w));
    printf("nr equ %ld\n", offsetof(params, nr));
    printf("eps equ %ld\n", offsetof(params, eps));
    printf("tmin equ %ld\n", offsetof(params, tmin));
    printf("tmax equ %ld\n", offsetof(params, tmax));
    printf("exaustive equ %ld\n", offsetof(params, exaustive));
    printf("symmetric equ %ld\n", offsetof(params, symmetric));
    printf("silent equ %ld\n", offsetof(params, silent));
    printf("display equ %ld\n", offsetof(params, display));
    printf("ANN equ %ld\n", offsetof(params, ANN));
    printf("vq equ %ld\n", offsetof(params, vq));
    printf("pq equ %ld\n", offsetof(params, pq));
    printf("query_pq equ %ld\n", offsetof(params, query_pq));
    printf("codebook equ %ld\n", offsetof(params, codebook));
    printf("distanze_simmetriche equ %ld\n", offsetof(params, distanze_simmetriche));
    printf("nDist equ %ld\n", offsetof(params, nDist));
    printf("distanze_asimmetriche equ %ld\n", offsetof(params, distanze_asimmetriche));
    printf("qc_indexes equ %ld\n", offsetof(params, qc_indexes));
    printf("qc equ %ld\n", offsetof(params, qc));
    printf("residual_codebook equ %ld\n", offsetof(params, residual_codebook));
    printf("residual_set equ %ld\n", offsetof(params, residual_set));
    printf("v equ %ld\n", offsetof(params, v));
    printf("zero equ %ld\n", offsetof(params, zero));


	printf("source equ %ld\n", offsetof(kmeans_data, source));
    printf("dim_source equ %ld\n", offsetof(kmeans_data, dim_source));
    printf("index %ld\n", offsetof(kmeans_data, index));
    printf("dest %ld\n", offsetof(kmeans_data, dest));
    printf("index_rows %ld\n", offsetof(kmeans_data, index_rows));
    printf("index_columns %ld\n", offsetof(kmeans_data, index_columns));
    printf("n_centroidi %ld\n", offsetof(kmeans_data, n_centroidi));
    printf("d %ld\n", offsetof(kmeans_data, d));
    
}