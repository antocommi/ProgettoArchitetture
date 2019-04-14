/**************************************************************************************
 * 
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2017/18
 * 
 * Progetto dell'algoritmo di PageRank
 * in linguaggio assembly x86-32 + SSE
 * 
 * Fabrizio Angiulli, 23 novembre 2017
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

 nasm -f elf32 pagerank32.nasm && gcc -O0 -m32 -msse pagerank32.o pagerank32c.c -o pagerank32c && ./pagerank32c
 
 oppure
 
 ./runpagerank32

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>


#define	MATRIX		double*
#define	VECTOR		double*
#define	GRAPH		double*	// DECIDERE LA RAPPRESENTAZIONE IN MEMORIA (dev'essere un puntatore)


typedef struct {
	int x;
	int y;
} location;


typedef struct {
	char* file_name;
	MATRIX P; // codifica dense
	GRAPH G; // codifica full
	int N; // numero di nodi
	int M; // numero di archi
	double c; // default=0.85
	double eps; // default=1e-5
	int format; // 0=sparse, 1=full
	int prec; // 0=single, 1=double
	int opt; // 0=nopt, 1=opt
	double* pagerank;
	int silent;
	int display;
} params;


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
	return (MATRIX) get_block(sizeof(double),rows*cols);
}


void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}


/*
 * 
 * 	load_dense
 * 	===========
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
MATRIX load_dense(char* filename, int *n, int *m) {	
	FILE* fp;
	int rows, cols, status, i;
	char fpath[256];
	
	sprintf(fpath, "%s.matrix", filename);
	fp = fopen(fpath, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad matrix file name!\n", fpath);
		exit(0);
	}
	
	status = fread(&rows, sizeof(int), 1, fp);
	status = fread(&cols, sizeof(int), 1, fp);
		
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(double), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*m = rows*cols;
	
	return data;
}

/*
 * 
 * 	load_sparse
 * 	===========
 * 
 *	Legge da file il grafo
 * 
 * 	Codifica del file:
 * 	primi 4 byte: numero di nodi (N) --> numero intero
 * 	successivi 4 byte: numero di archi (M) --> numero intero
 * 	successivi M*2*4 byte: M archi rappresentati come coppie (i,j) di interi a 32 bit 
 * 
 */
GRAPH load_sparse(char* filename, int *n, int *m) {	
	FILE* fp;
	int nodes, arcs, status, i;
	char fpath[256];
	int sorg, dest;
	
	sprintf(fpath, "%s.graph", filename);
	fp = fopen(fpath, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad graph file name!\n", fpath);
		exit(0);
	}
	
	status = fread(&nodes, sizeof(int), 1, fp);
	status = fread(&arcs, sizeof(int), 1, fp);

	GRAPH g; // = ...; // alloca la struttura dati contenente il grafo
	
	for (i = 0; i < arcs; i++) {
		status = fread(&sorg, sizeof(int), 1, fp);
		status = fread(&dest, sizeof(int), 1, fp);
		// aggiungi l'arco (sorg,dest) a g
	}
	fclose(fp);
	
	*n = arcs;
	*m = nodes;
	
	// return g;
	return NULL;
}


void save_pageranks(char* filename, int n, VECTOR pagerank) {	
	FILE* fp;
	int i;
	char fpath[256];
	
	sprintf(fpath, "%s_pageranks.txt", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < n; i++)
		fprintf(fp, "%.14g\n", pagerank[i]);
	fclose(fp);
}


extern void pagerank32(params* input);


/*
 *	pagerank
 * 	====
 * 
 *	img contiene l'immagine codificato come una matrice di N righe
 * 	ed M colonne memorizzata in un array lineare in row-major order
 * 
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	dell'immagine.
 * 
 * 
 */
void pagerank(params* input) {
	
    // -------------------------------------------------
    // Codificare qui l'algoritmo risolutivo
    // -------------------------------------------------
    
    pagerank32(input); // Esempio di chiamata di funzione assembly

    // -------------------------------------------------

}


#define	SPARSE	0
#define	DENSE	1
#define	SINGLE	0
#define	DOUBLE	1
#define	NOPT	0
#define	OPT		1



int main(int argc, char** argv) {
	
	params* input = malloc(sizeof(params));

	input->file_name = NULL;
	input->P = NULL; // dense format
	input->G = NULL; // sparse format
	input->N = 0; // number of nodes
	input->M = 0; // number of arcs
	input->c = 0.85;
	input->eps = 1e-5;
	input->format = SPARSE; // 0=sparse, 1=dense
	input->prec = SINGLE; // 0=single, 1=double
	input->opt = NOPT; // 0=nopt, 1=opt
	input->silent = 0; // 0=false,<>0=true
	input->display = 0; // 0=false <>0=true

	int i, j;

	int par = 1;
	while (par < argc) {
		if (par == 1) {
			input->file_name = argv[par];
			par++;
		} else if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-c") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing c value!\n");
				exit(1);
			}
			input->c = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-eps") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing eps value!\n");
				exit(1);
			}
			input->eps = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sparse") == 0) {
			input->format = SPARSE;
			par++;
		} else if (strcmp(argv[par],"-dense") == 0) {
			input->format = DENSE;
			par++;
		} else if (strcmp(argv[par],"-single") == 0) {
			input->prec = SINGLE;
			par++;
		} else if (strcmp(argv[par],"-double") == 0) {
			input->prec = DOUBLE;
			par++;
		} else if (strcmp(argv[par],"-nopt") == 0) {
			input->opt = NOPT;
			par++;
		} else if (strcmp(argv[par],"-opt") == 0) {
			input->opt = OPT;
			par++;
		} else
			par++;
	}
	
	if (!input->silent) {
		printf("Usage: %s <input_file_name> [-d][-s][-sparse|-dense][-single|-double][-nopt|-opt][-c <value>][-eps <value>]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\t-d : display input and output\n");
		printf("\t-s : silent\n");
		printf("\t-sparse/-full: input format (sparse=list of arcs,full=matrix)\n");
		printf("\t-single/-double: floating-point precision (only sparse format)\n");
		printf("\t-nopt/-opt: disable/enable optimizations\n");
		printf("\t-c <value> : 1-teleportation_probability (default 0.85)\n");
		printf("\t-eps <value> : termination error (default 1e-5)\n");
		printf("\n");
	}
	
	if (input->file_name == NULL || strlen(input->file_name) == 0) {
		printf("Missing input file name!\n");
		exit(1);
	}
	
	if (input->format == 0)
		input->G = load_sparse(input->file_name, &input->N, &input->M);
	else
		input->P = load_dense(input->file_name, &input->N, &input->M);
		
	if (!input->silent) {
		printf("Input file name: '%s'\n", input->file_name);
		printf("Number of nodes: %d\n", input->N);
		printf("Number of arcs: %d\n", input->M);
		printf("Parameter c: %f\n", input->c);
		printf("Parameter eps: %f\n", input->eps);
	}
	
	clock_t t = clock();
	pagerank(input);
	t = clock() - t;
	
	if (!input->silent)
		printf("\nExecution time = %.3f seconds\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);
			
	if (input->pagerank != NULL)
	{
		if (!input->silent && input->display) {
			printf("\nPageRanks:\n");
			for (i = 0; i < input->N; i++) {
				printf("%d %.14g\n", i+1, input->pagerank[i]);
			}
		}
		save_pageranks(input->file_name, input->N, input->pagerank);
	}
	
	return 0;
}
