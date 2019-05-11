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

#define	MATRIX		double*
#define	VECTOR		double*

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
	//
	// Inserire qui i campi necessari a memorizzare i Quantizzatori
	//
	VECTOR vq;
	int* pq;
	int* query_pq;
	MATRIX codebook;
	MATRIX distanze_simmetriche;
	int nDist;
	MATRIX distanze_asimmetriche;
	// ...
	// ...
	// ...
	struct entry* v;
	//
} params;

//Entry della s.d. multilivello
typedef struct {
	int index;
	VECTOR q;

	//temporaneo
	//Serve per gestire liste a dimensione sconosciuta. 
	entry next;
	
} entry;

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
	status = fread(data, sizeof(double), rows*cols, fp);
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

// Fuzioni fatte da noi

int calcolaIndice(int i, int j){
	//funzione che calcola l'indice per la matrice delle distanze_simmetriche
	return i*(i-1)/2+j;
}

int dist_eI(params* input,int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	int i;
	int ret=0;
	for(i=start; i<end; i++){
		ret += pow(input->ds[punto1*input->d+i]-input->ds[punto2*input->d+i], 2.0);
	}
	return ret;
}

int dist_e(params* input,int punto1, int punto2){
	int i;
	double sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_eI(input, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

int calcolaPQ(params* input, int x, int start, int end){
	// estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    double min=1.79E+308;
    int imin=-1;
    double temp;
    for(i=0; i<input->k; i++){
        temp=dist_eI(input, x, i, start, end);
        if(temp<min){ 
            min=temp;
            imin=i;
        }
    }
    return imin;
}

double dist_simmetricaI(params* input, int centroide1, int centroide2, int start, int end){
	// estremi start incluso ed end escluso
	int i;
	double ret=0;
	for(i=start; i<end; i++){
		ret += pow( input->codebook[centroide1*input->d+i] - input->codebook[centroide2*input->d+i] , 2);
	}
	return ret;
}

double dist_simmetrica(params* input, int centroide1, int centroide2){
	int i;
	double sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_simmetricaI(input, centroide1, centroide2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

double dist_asimmetricaI(params* input, MATRIX set, int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int i, c;
	double ret=0;
	c=input->pq[punto2];
	for(i=start; i<end; i++){
		ret += pow( set[punto1*input->d+i] - input->codebook[c*input->d+i] , 2);
	}
	return ret;
}

double dist_asimmetrica(params* input, MATRIX set, int punto1, int punto2){
	int i;
	double sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(dist_asimmetricaI(input, set, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

double distI(params* input, int* quantizer, int punto1, int punto2, int start, int end){
	// estremi start incluso ed end escluso
	// punto2 è un punto del dataset
	//
	// punto1 può essere del dataset o del query set, quindi in set si passa
	// la constante DATASET o QUERYSET
	int c1, c2;
	if(punto1==punto2){
		return 0;
	}else{
		c1=quantizer[punto1*input->m+(start/input->m)];
		c2=input->pq[punto2*input->m+(start/input->m)];
		if(c1<c2){
			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(c2, c1)];
		}else{
			return input->distanze_simmetriche[(input->nDist*start/input->m)+calcolaIndice(c1, c2)];
		}
	}
}
double dist(params* input, int* quantizer, int punto1, int punto2){
	int i;
	double sum=0;
	for(i=0; i<input->m; i++){
		sum+=pow(distI(input, quantizer, punto1, punto2, i*input->m, (i+1)*input->m), 2);
	}
	return sum;
}

int calcolaQueryPQ(params* input, int x, int start, int end){
	// estremi start incluso ed end escluso
    //
    //	INPUT: 	Punto x di dimensione d.
    //	OUTPUT: indice del centroide c più vicino ad x. 
    //
    int i;
    double min=1.79E+308;
    int imin=-1;
    double temp;
	if(input->symmetric==1){
		for(i=0; i<input->k; i++){
			temp=distI(input, input->query_pq, x, i, start, end);
			if(temp<min){ 
				min=temp;
				imin=i;
			}
		}
	}else{
		for(i=0; i<input->k; i++){
			temp=dist_asimmetricaI(input, input->qs, x, i, start, end);
			if(temp<min){ 
				min=temp;
				imin=i;
			}
		}
	}
	return imin;
}

void kmeans(params* input, int start, int end, int n_centroidi){
	// estremi start incluso ed end escluso
	int i, j, k, t;
	int count;
	double fob1, fob2;
	double* codebook;

	codebook = alloc_matrix(n_centroidi, input->n); // row-major-order?
    if(codebook==NULL) exit(-1);
	
	//
	// Inizializzazione del codebook
	//		-Scelta dei k vettori casuali
	//
	
    for(i=0; i<n_centroidi; i++){
		k=rand()%input->n;
		for(j=start; j<end; j++){
			codebook[i*input->d+j]=input->ds[k*input->d+j];
		}
    }
    
    for(i=0; i<input->n; i++){
        input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
    }
    
	fob1=0; //Valori della funzione obiettivo
	fob2=0;
	for(t=0; t<input->tmin || (t<input->tmax && (fob2-fob1) > input->eps); t++){
		for(i=0; i<n_centroidi; i++){
			count=0;
			for(j=start; j<end; j++){
				codebook[i*input->d+j]=0; // con calloc forse è più veloce. 
			}
			
			//
			// INIZIO: RICALCOLO NUOVI CENTROIDI
			//
			
			for(j=0; j<input->n; j++){
				if(input->pq[j*input->m+(start/input->m)]==i){ // se q(Yj)==Ci -- se Yj appartiene alla cella di Voronoi di Ci
					count++;
					for(k=start; k<end; k++){
						codebook[i*input->d+k]+=input->ds[j*input->d+k];
					}
				}
			}
			
			for(j=start; j<end; j++){
				if(count!=0){ 
					// Alcune partizioni potrebbero essere vuote
					// Specie se ci sono degli outliers
					codebook[i*input->d+j]=codebook[i*input->d+j]/count;
				}
			}
			
			//
			// FINE: RICALCOLO NUOVI CENTROIDI
			//
		}
		
		for(i=0; i<input->n; i++){
			input->pq[i*input->m+(start/input->m)]=calcolaPQ(input, i, start, end);
		}
		
		fob1=fob2;
		fob2=0;
		
		//CALCOLO NUOVO VALORE DELLA FUNZIONE OBIETTIVO
		for(i=0; i<input->n; i++){
			fob2+=pow(dist_eI(input, i, input->pq[i*input->m+(start/input->m)], start, end), 2.0);
		}
	}

	input->codebook=codebook;
}

void creaMatricedistanze(params* input){
	int i, j, k;
	MATRIX distanze_simmetriche;
	input->nDist=input->k*(input->k+1)/2;
	distanze_simmetriche = alloc_matrix(input->m, input->nDist);
	if(distanze_simmetriche==NULL) exit(-1);
	for(k=0; k<input->m; k++){
		for(i=1; i<input->k; i++){
			for(j=0; j<i; j++){
				distanze_simmetriche[k*input->nDist+calcolaIndice(i, j)] = dist_simmetricaI(input, i, j, k*input->m, (k+1)*input->m);
				// verificare se qui va usata la distanza simmetrica o no
			}
		}
	}
	input->distanze_simmetriche=distanze_simmetriche;
}

void bubbleSort(double* arr, int* arr2, int n){ 
   int i, j, temp; 
   for (i = 0; i < n-1; i++)    
       for (j = 0; j < n-i-1; j++)  
           if (arr[j] > arr[j+1]){
			   temp=arr[j];
			   arr[j]=arr[j+1];
			   arr[j+1]=temp;
			   temp=arr2[j];
			   arr2[j]=arr2[j+1];
			   arr2[j+1]=temp;
		   }
} 

void calcolaNN(params* input, int query){
	int i;
	VECTOR distanze=alloc_matrix(input->n, 1);
	int* d2=(int*) _mm_malloc(input->n*sizeof(int),16);
	if(input->symmetric==0){
		for(i=0; i<input->n; i++){
			distanze[i]=dist_asimmetrica(input, input->qs, query, i);
			d2[i]=i;
		}
	}else{
		for(i=0; i<input->n; i++){
			distanze[i]=dist(input, input->query_pq, query, i);
			d2[i]=i;
		}
	}
	bubbleSort(distanze, d2, input->n);
	dealloc_matrix(distanze);
	for(i=0; i<input->knn; i++){
		input->ANN[query*input->knn+i]=d2[i];
	}
	_mm_free(d2);
}

void inizializza_learning_set(params* input){
	//to do
}

void inizializzaSecLiv(params* input){

	input->v = _mm_malloc(sizeof(entry)*input->kc,16);
	if(input->v==NULL) return;
	for(int i=0;i<input->kc;i++){
		input->v[i].next=NULL;

	}	
}

void add (entry e,int i,params* input){
	if(input->v==)
	input->v=
	
}

/*
 *	pqnn_index
 * 	==========
 */
void pqnn_index(params* input) {
	int i, dStar;
	// TODO: Gestire liberazione della memoria.
	if(input->exaustive==1){
		input->pq = (int*) _mm_malloc(input->n*input->m*sizeof(int), 16); 
		dStar=input->d/input->m;
		for(i=0; i<input->m; i++){
			kmeans(input, i*dStar, (i+1)*dStar, input->k);
		}
		// controllare caso in cui d non sia multiplo di m
		if(input->symmetric==1){
			creaMatricedistanze(input);
		}
	}
	else{
		//
		// RICERCA NON ESAUSTIVA
		//
		inizializza_learning_set(input);//selezionati i primi nr del dataset
		kmeans(input, 0, input->d, input->kc);//calcolo i grossolani
		inizializzaSecLiv(input);
	}
    
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
		//ricerca esaustiva
		input->query_pq=(int*)_mm_malloc(input->nq*input->m*sizeof(int), 16);
		for(i=0; i<input->nq; i++){
			for(j=0; j<input->m; j++){
				input->query_pq[i*input->nq+j]=calcolaPQ(input, i, i*input->m, (i+1)*input->m);
			}
		}
		input->ANN=(int*)_mm_malloc(input->nq*input->knn*sizeof(int),16);
		for(i=0; i<input->nq; i++){
			calcolaNN(input, i);
		}
	}else{
		//ricerca non esaustiva
	}
    // -------------------------------------------------
    // Codificare qui l'algoritmo di interrogazione
    // -------------------------------------------------
    
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
