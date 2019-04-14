/**************************************************************************************
 * 
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2016/17
 * 
 * Progetto di un algoritmo di Corner Detection
 * in linguaggio assembly x86-32 + SSE
 * 
 * Fabrizio Angiulli, 20 aprile 2016
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

 nasm -f elf64 corner64.nasm && gcc -O0 -m64 -msse corner64.o corner64c.c -o corner64c && ./corner64c
 
 oppure
 
 ./runcorner64

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>


#define	MATRIX		float*
#define	VECTOR		float*
#define	IMAGE		float*


typedef struct {
	int x;
	int y;
} location;


typedef struct {
	int method; //
	char* img_name; //
	IMAGE img; //
	int N; //
	int M; //
// corner detection parameters
	float sigma; // 
	float theta; // 
	int t; // 
// filtering parameters
	char* filt_name; // 
	IMAGE filt; // 
	int N_filt; // 
	int M_filt; // 
// generic parameters
	int silent; // 
	int display; // 
// filtering output
	IMAGE filtered_img; // 
// corner detection output
	IMAGE R; // 
	int n_corners; // 
	location* corners; // 
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
 * 	load_input
 * 	===========
 * 
 *	Legge da file l'immagine codificata come una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 * 
 * 	Codifica del file:
 * 	primi 4 byte: numero di colonne (M) --> numero intero
 * 	4 byte successivi: numero di righe (N) --> numero intero
 * 	N*M*4 byte successivi: imag data in row-major order --> numeri floating-point a precisione singola 
 * 
 *****************************************************************************
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	dell'immagine. 
 *****************************************************************************
 * 
 */
IMAGE load_img(char* filename, int *n, int *m) {	
	FILE* fp;
	int rows, cols, status, i;
	char fpath[256];
	
	sprintf(fpath, "%s.img", filename);
	fp = fopen(fpath, "rb");
	
	if (fp == NULL) {
		printf("'%s' : bad image file name!\n", fpath);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	IMAGE data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(float), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*m = cols;
	
	return data;
}


void save_img(char* filename, IMAGE img, int n, int m) {
	FILE* fp;
	char fpath[256];
	
	sprintf(fpath, "%s_output.img", filename);
	fp = fopen(fpath, "wb");
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	if (img != NULL)
		fwrite(img, sizeof(float), n*m, fp);
	fclose(fp);
}


void save_corners(char* filename, int n_corners, location* corners) {	
	FILE* fp;
	int i;
	char fpath[256];
	
	sprintf(fpath, "%s_corners.txt", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < n_corners; i++)
		fprintf(fp, "%d %d\n", corners[i].x, corners[i].y);
	fclose(fp);
}


extern void corner64(params* input);


/*
 *	corner
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
void corner(params* input) {
	
    // -------------------------------------------------
    // Codificare qui l'algoritmo risolutivo
    // -------------------------------------------------
    
    corner64(input); // Esempio di chiamata di funzione assembly

    // -------------------------------------------------

}


#define	CORNER	((int)0)
#define	FILTERING	((int)1)


//char* method_name[2] = {"Corner detection", "Filtering"};



int main(int argc, char** argv) {
	
	params* input = malloc(sizeof(params));

	input->method = CORNER;
	input->img_name = "";
	input->img = NULL;
	input->N = 0;
	input->M = 0;
	input->sigma = 1.5f;
	input->theta = 0.2f;
	input->t = 1;
	input->filt_name = "";
	input->filt = NULL;
	input->N_filt = 0;
	input->M_filt = 0;
	input->silent = 0;
	input->display = 0;
	input->filtered_img = NULL;
	input->R = NULL;
	input->n_corners = 0;
	input->corners = NULL;

	int i, j;

	int par = 1;
	while (par < argc) {
		if (par == 1) {
			input->img_name = argv[par];
			par++;
		} else if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-sigma") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing sigma value!\n");
				exit(1);
			}
			input->sigma = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-theta") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing theta value!\n");
				exit(1);
			}
			input->theta = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-t") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing t value!\n");
				exit(1);
			}
			input->t = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-corner") == 0) {
			input->method = CORNER;
			par++;
		} else if (strcmp(argv[par],"-filtering") == 0) {
			input->method = FILTERING;
			par++;
			if (par >= argc) {
				printf("Missing filter file name!\n");
				exit(1);
			}
			input->filt_name = argv[par];
			par++;
		} else
			par++;
	}
	
	if (!input->silent) {
		printf("Usage: %s <img_file_name> [-d][-s] [[-corner][-sigma <value>][-theta <value>][-t <value>]] [-filtering <filter_file_name>]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\t-corner: perform corner detection (default)\n");
		printf("\t-d : display input and output\n");
		printf("\t-s : silent\n");
		printf("\t-sigma <value> : st.dev. of Gaussian filter (default 1.5)\n");
		printf("\t-theta <value> : thresholding parameter (default 0.2)\n");
		printf("\t-t <value> : non-maxima suppression parameter (default 1)\n");
		printf("\t-filtering <filter_file_name> : perform image filtering\n");
		printf("\n");
	}
	
	if (strlen(input->img_name) == 0) {
		printf("Missing image file name!\n");
		exit(1);
	}
	
	input->img = load_img(input->img_name, &input->N, &input->M);
	
	if (input->method == FILTERING)
		input->filt = load_img(input->filt_name, &input->N_filt, &input->M_filt);

	if (!input->silent) {
		printf("Image file name: '%s'\n", input->img_name);
		printf("Image rows: %d\n", input->N);
		printf("Image columns: %d\n", input->M);
		if (input->method == CORNER) {
			printf("Parameter sigma: %f\n", input->sigma);
			printf("Parameter theta: %f\n", input->theta);
			printf("Parameter t: %d\n", input->t);
		} else if (input->method == FILTERING) {
			printf("Filter file name: '%s'\n", input->filt_name);
			printf("Filter rows: %d\n", input->N_filt);
			printf("Filter columns: %d\n", input->M_filt);
		}
	}
	
	clock_t t = clock();
	corner(input);
	t = clock() - t;

	
	if (!input->silent)
		printf("\nExecution time = %.3f seconds\n", ((float)t)/CLOCKS_PER_SEC);
	else
		printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);
		
	if (!input->silent && input->display) {
		printf("\nCorner locations:\n");
		for (i = 0; i < input->n_corners; i++) {
			printf("(%d,%d)\n", input->corners[i].x, input->corners[i].y);
		}
	}
	
	if (input->method == CORNER) {
		save_corners(input->img_name, input->n_corners, input->corners);
		save_img(input->img_name, input->R, input->N, input->M);
	} else if (input->method == FILTERING) {
		save_img(input->img_name, input->filtered_img, input->N, input->M);
	}

	return 0;
}
