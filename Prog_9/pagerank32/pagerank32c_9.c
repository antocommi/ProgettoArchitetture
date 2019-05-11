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
#define	VECTORD		double*
#define VECTORF     float*
#define GRAPH       location*
#define	SPARSE	0
#define	DENSE	1
#define	SINGLE	0
#define	DOUBLE	1
#define	NOPT	0
#define	OPT		1

typedef struct {
	int x;
	int y;
} location;

typedef struct {
    double *dd;
    float *df;
    int *col;
    int *row;
} csc;


typedef struct {
	char* file_name;
	MATRIX P; // codifica dense
	csc G; // codifica full
	int N; // numero di nodi
        int restoN;
	int M; // numero di archi
        int restoM;
        GRAPH L;
	double c; // default=0.85
	double eps; // default=1e-5
	int format; // 0=sparse, 1=full
	int prec; // 0=single, 1=double
	int opt; // 0=nopt, 1=opt
	VECTORD pagerank;
	int silent;
	int display;
} params;


void* get_block(int size, int elements) {
	return _mm_malloc(elements*size,16);
}

double* get_blockD(int size, int elements) {
	double* punt= _mm_malloc(elements*size,16);
        return punt;
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

MATRIX load_dense(char* filename, int *n, int *m, int *nr, int *mr) {
	FILE* fp;
	int rows, cols, status, i,j,restoN,restoM;
	char fpath[256];

	sprintf(fpath, "%s.matrix", filename);
	fp = fopen(fpath, "rb");

	if (fp == NULL) {
		printf("'%s' : bad matrix file name!\n", fpath);
		exit(0);
	}

	status = fread(&rows, sizeof(int), 1, fp);
	status = fread(&cols, sizeof(int), 1, fp);
        restoN = 8-cols%8;
        if(restoN==8) restoN=0;
        restoM = 16-cols%16;
        if(restoM==16) restoM=0;

	MATRIX data = alloc_matrix(rows,cols+restoN);
        for(i=0;i<rows;i++){
           for(j=0;j<cols;j++){
              status = fread(data+i+j*(rows+restoN), sizeof(double), 1, fp);  
           }
        }
    
        for(i=rows;i<rows+restoN;i++){
           for(j=0;j<cols;j++) data[i+j*(rows+restoN)]=0.0;
        }
	fclose(fp);

        *mr=restoM;
        *nr = restoN;
	*n = rows;
	*m = rows*cols;
    
	return data;
}

extern void initOutF(GRAPH L, int* d, VECTORF df,int m);
extern void initOutD(GRAPH L, int* d, VECTORD df,int m);
extern void initCsc(GRAPH L,int* col,int* row,int M);

csc creaGF(location* L, double c, int N, int M, int restoM){
        int i,j;
	float* df;
	int* col;
	int* row;
	csc graph;
	int sorg,dest;

	df = (float *) _mm_malloc ((N+restoM)*sizeof(float),16);
	col = (int *) _mm_malloc ((N+1)*sizeof(int),16);
	row = (int *) _mm_malloc (M*sizeof(int),16);
	for(i=0;i<N;i++){
	   col[i]=0;
	   df[i]=0.0;

        }
        for(i=N;i<N+restoM;i++){
	    df[i]=0.0;

        }
	col[N]=0;

	graph.col=col;
	graph.row=row;
	graph.df=df;
	/*for (i = 0; i < M; i++) {
		sorg = L[i].x;
		dest = L[i].y;
		col[dest+1]= col[dest+1]+1;
		df[sorg]= df[sorg]+1.0;
	}*/
        initOutF(L,col,df,M);

	for(i=1; i<N; i++){
		col[i]= col[i-1]+col[i];
	}

	/*for(i=0; i<M; i++){
		j=col[L[i].y];
		col[L[i].y]++;
		row[j]=L[i].x;
	}*/
        initCsc(L,col,row,M);
	for(i=0;i<N;i++){
		   if(df[i]==0.0)df[i]=-1.0;
		   else df[i]=c/df[i];
	}
         _mm_free(L);
	return graph;
}

csc creaGD(location* L, double c, int N, int M, int restoM){
	double* dd;
	int* col;
	int* row;
	csc graph;
	int sorg,dest,i,j;

	dd = (double *) _mm_malloc ((N+restoM)*sizeof(double),16);
	col = (int *) _mm_malloc ((N+1)*sizeof(int),16);
	row = (int *) _mm_malloc (M*sizeof(int),16);
	for(i=0;i<N;i++){
	   col[i]=0;
	   dd[i]=0.0;

   }
   for(i=N;i<N+restoM;i++){
	 dd[i]=0.0;

   }
	col[N]=0;

	graph.col=col;
	graph.row=row;
	graph.dd=dd;
	/*for (i = 0; i < M; i++) {
		sorg = L[i].x;
		dest = L[i].y;
		col[dest+1]= col[dest+1]+1;
		dd[sorg]= dd[sorg]+1.0;
	}*/

        initOutD(L,col,dd,M);

	for(i=1; i<N; i++){
		col[i]= col[i-1]+col[i];
	}

	/*for(i=0; i<M; i++){
		j=col[L[i].y];
		col[L[i].y]++;
		row[j]=L[i].x;
	}*/
        initCsc(L,col,row,M);
	for(i=0;i<N;i++){
		   if(dd[i]==0.0)dd[i]=-1.0;
		   else dd[i]=c/dd[i];
	}
         _mm_free(L);
	return graph;
}
GRAPH load_sparse(char* filename, int *n, int *m,double c,int prec, int* rN, int* rM) {
	FILE* fp;
	int nodes, arcs, status, i,j,restoN,restoM;
	char fpath[256];
	int sorg, dest;
	location* g;



 	sprintf(fpath, "%s.graph", filename);
	fp = fopen(fpath, "rb");

	if (fp == NULL) {
		printf("'%s' : bad graph file name!\n", fpath);
		exit(0);
	}
    
	status = fread(&nodes, sizeof(int), 1, fp);
	status = fread(&arcs, sizeof(int), 1, fp);
	if(prec==SINGLE){
		restoN= 16-nodes%16;
		if(restoN==16)restoN=0;
		restoM= 32-nodes%32;
		if(restoM==32) restoM=0;
	}
	else{
		restoN= 8-nodes%8;
		if(restoN==8)restoN=0;
		restoM= 16-nodes%16;
		if(restoM==16) restoM=0;
	}
    g =(location *) _mm_malloc(arcs*sizeof(location),16);
	for (i = 0; i < arcs; i++) {
		status = fread(&sorg, sizeof(int), 1, fp);
		status = fread(&dest, sizeof(int), 1, fp);
		g[i].x= sorg-1;
		g[i].y= dest-1;

	}
	fclose(fp);
      
        
        *rN= restoN;
        *rM= restoM;
	*m = arcs;
	*n = nodes;
        
	return g;
}


void save_pageranks(char* filename, int n, VECTORD pagerank) {
	FILE* fp;
	int i;
	char fpath[256];

	sprintf(fpath, "%s_pageranks.txt", filename);
	fp = fopen(fpath, "w");
	for (i = 0; i < n; i++)
		fprintf(fp, "%.14g\n", pagerank[i]);
	fclose(fp);
}

//extern void norma1d( VECTORD v, int dim, double* norma);
/*void norma1d2( VECTORD v, int dim, double* norma) {
	int i;
	double norm=0.0;

	for(i=0;i<dim;i++){
		norm += fabs(v[i]);
	}

	*norma= norm;
}*/
//extern void norma1f( VECTORF v, int dim, float* norma);
/*void norma1f2( VECTORF v, int dim, float* d) {
    int i;
    float norm=0.0;
    
    for(i=0;i<dim;i++){
        norm += fabs(v[i]);
    }
    
    *d=norm;
}*/

extern void sommaf(VECTORF v, int dim,float* s);
/*void sommaf2( VECTORF v, int dim,float* s) {
    int i;
    float norm=0.0;
    
    for(i=0;i<dim;i++){
        norm += v[i];
    }
    
    *s=norm;
}*/

extern double sommad( VECTORD v, int dim, double* s);

/*void sommad2( VECTORD v, int dim, double* s) {
    int i;
    double norm=0.0;
    
    for(i=0;i<dim;i++){
        norm += v[i];
    }
    
    *s=norm;
}*/

extern void differenceNormd(VECTORD v1, VECTORD v2,int dim,double* delta);
/*void differenceNormd2(VECTORD v1, VECTORD v2,int dim,double* delta){
	int i;
        double norm;
	for(i=0;i<dim;i++){
		v1[i] = v1[i]-v2[i];
	}
        norma1d2(v1,dim,&norm);
        *delta=norm;
}*/

extern void differenceNormf(VECTORF v1, VECTORF v2,int dim, float* delta);
/*float differenceNormf2(VECTORF v1, VECTORF v2,int dim, float* delta){
    int i;
    float norm;
    for(i=0;i<dim;i++){
        v1[i] = v1[i]-v2[i];
    }
    norma1f2(v1,dim,&norm);
    *delta=norm;
}*/

extern void prod(MATRIX P, VECTORD v,int dim,int n,VECTORD temp);
/*void prod2(MATRIX P, VECTORD v,int dim,int n,VECTORD temp){
	int i,j;
	double somma;

	for(i=0;i<n;i++){
		somma=0.0;
		for(j=0;j<dim;j++){
			somma = somma + P[j+i*dim]*v[j];
                        
		}
		temp[i]=somma;
	}
}*/

extern void prodSparseF(int* row, int* col, VECTORF vf, int N, VECTORF xk1, VECTORF xk,VECTORF d,int M);
/*void prodSparseF2(int* row, int* col, VECTORF vf, int N, VECTORF xk1, VECTORF xk,float*d,int M){
    int i,j,colel,currow=0;
    float s;
    for(i=0;i<M;i++){
      vf[i]= (xk1[i])*(d[i]);
    }
    colel=col[0];
    for(i=0;i<N;i++){
        s=0.0;
        for(j=0; j<colel; j++){
            s += vf[row[currow]];
            currow++;
        }
        xk[i]=s;
        colel=col[i+1]-col[i];
        
    }
    
}*/

/*void prodSparseD2(int* row, int* col, VECTORD vf, int N, VECTORD xk1, VECTORD xk,double*d, int M){
    int i,j,rowel,curcol=0;
    double s;
    for(i=0;i<M;i++){
      vf[i]= (xk1[i])*(d[i]);
    }
    rowel=col[0];
    for(i=0;i<N;i++){
        s=0.0;
        for(j=0; j<rowel; j++){
            s+= vf[row[curcol]];
            curcol++;

        }
        xk[i]=s;
        rowel=col[i+1]-col[i];
        //printf("%d %f\n",i,xk[i]);
        
    }
    
}*/
/*void completaSolD2(VECTORD xk,int dim,double val){
    int i;
    for(i=0;i<dim;i++){
        xk[i]=xk[i]+val;
    }
}*/

/*void completaSolF2(VECTORF xk,int dim,float val){
    int i;
    for(i=0;i<dim;i++){
        xk[i]=xk[i]+val;
    }
}*/
extern void completaSolD(VECTORD xk,int dim,double val);
extern void completaSolF(VECTORF xk,int dim,float val);
extern void prodSparseD(int* row, int* col, VECTORD vd, int N, VECTORD xk1, VECTORD xk,VECTORD d,int M);
extern void normalizzazione(VECTORD v, int dim,int dim2);
/*void normalizzazione2(VECTORD v, int dim,int dim2){
	double norma;
	int i;
	norma1d2(v,dim,&norma);
	for(i=0;i<dim2;i++){
		v[i] = v[i]/norma;
	}
}*/
extern void assegnaSoluzioneF(float* xk1,float* xk,int n);
/*void assegnaSoluzioneF2(float* xk1,float* xk,int n){
    int i;
    for(i=0;i<n;i++){
        xk1[i] = xk[i];
        xk[i]=0.0;
    }
}*/
extern void assegnaSoluzioneD(double* pagerank, double* xk, int n);
/*void assegnaSoluzioneD2(double* pagerank, double* xk, int n){
    int i;
    for(i=0;i<n;i++){
        pagerank[i] = xk[i];
        xk[i]=0.0;
    }
}*/

extern void assegnaD(double* v1, double* v2,int n);
/*void assegnaD2(double* v1, double* v2,int n){
    int i;
    for(i=0;i<n;i++){
        v1[i] = v2[i];
    }
}*/ 

/*void assegnaF2(float* v1, float* v2,int n){
    int i;
    for(i=0;i<n;i++){
        v1[i] = v2[i];
    }
}*/

void converti(double* pagerank,float* xk1,int n){
    int i;
    for(i=0;i<n;i++){
        pagerank[i]=xk1[i];
    }
}

double power(double c,int d){
    if(d==0) return 1.0;
    if(d==1) return c;
    if(d%2==0) return power(c*c, d/2);
    else return c* power(c*c,d/2); 
}



void pagerankfullNoOpt(params* input){
	int i;
	double delta;
	double* temp = _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    for(i=0;i<input->N;i++){
        input->pagerank[i] = 1.0/input->N;
    }
    for(i=input->N;i<input->N+input->restoM;i++){
        input->pagerank[i]=0.0;
        temp[i]=0.0;
    }
	delta = 1.0;
	while(delta>=input->eps){
		prod(input->P,input->pagerank,input->N+input->restoN,input->N,temp);
	        differenceNormd(input->pagerank,temp,input->N+input->restoN,&delta);
                assegnaD(input->pagerank,temp,input->N+input->restoM);
	}
      normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
      _mm_free(temp);
      _mm_free(input->P);
    
}

void pageranksparseNoOptFloat(params* input){
    int i,j;
    float delta,s;
    float* xk= _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
    float* xk1= _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
    float* v= _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
   
    for(i=0;i<input->N;i++){
        xk1[i] = 1.0/input->N;
        xk[i]=0.0;
        v[i]=0.0;

    }
    for(i=input->N;i<input->N+input->restoM;i++){
        xk[i]=0.0;
        v[i]=0.0;
    }
    
    delta = 1.0;
    while(delta>=input->eps){
        prodSparseF(input->G.row,input->G.col,v,input->N,xk1,xk,input->G.df,input->N+input->restoM);
        sommaf(xk,input->N+input->restoM,&s);
        completaSolF(xk,input->N+input->restoM,(1.0-s)/input->N);
        for(i=input->N;i<input->N+input->restoM;i++){
             xk[i]=0.0;
        }
        differenceNormf(xk1,xk,input->N+input->restoN,&delta);
        assegnaSoluzioneF(xk1,xk,input->N+input->restoM);
    }
    converti(input->pagerank,xk1,input->N);
    normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
    _mm_free(xk);
    _mm_free(xk1);
    _mm_free(v);
}

void pageranksparseNoOptDouble(params* input){
    int i,j;
    double delta,somma;
    double* xk= _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    double* v= _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    for(i=0;i<input->N;i++){
        input->pagerank[i] = 1.0/input->N;
        xk[i]=0.0;
        v[i]=0.0;
    }
    for(i=input->N;i<input->N+input->restoM;i++){
        xk[i]=0.0;
        v[i]=0.0;
    }
    delta = 1.0;
    while(delta>=input->eps){
       	  prodSparseD(input->G.row,input->G.col,v,input->N,input->pagerank,xk,input->G.dd,input->N+input->restoM);
          sommad(xk,input->N+input->restoM,&somma);
          completaSolD(xk,input->N+input->restoM,(1.0-somma)/input->N);
          for(i=input->N;i<input->N+input->restoM;i++){
              xk[i]=0.0;
          }
          differenceNormd(input->pagerank,xk,input->N+input->restoN,&delta);
          assegnaSoluzioneD(input->pagerank,xk,input->N+input->restoM);
          
    }
    normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
     _mm_free(xk);
    _mm_free(v);
}

extern void powerExtrd(VECTORD l, VECTORD temp,double p,int n);
/*void powerExtrd2(VECTORD l, VECTORD temp,double p,int n){
    int i;
    for(i=0;i<n;i++){
        temp[i]=(temp[i]-p*l[i])/(1.0-p);
    }
}*/
extern void powerExtrf(VECTORF l, VECTORF temp,float p,int n);
/*void powerExtrf2(VECTORF l, VECTORF temp,float p,int n){
    int i;
    for(i=0;i<n;i++){
        temp[i]=(temp[i]-p*l[i])/(1.0-p);
    }
}*/

void pagerankfullOpt(params* input){
    int i,k,d=6;
    double delta,p;
    double* temp = _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    double* l = _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    
    for(i=0;i<input->N;i++){
        input->pagerank[i] = 1.0/input->N;
        temp[i]=0.0;
        l[i]=0.0;
    }
    for(i=input->N;i<input->N+input->restoM;i++){
        input->pagerank[i]=0.0;
        temp[i]=0.0;
        l[i]=0.0;
    }
    delta = 1.0;
    k=1;
    while(delta>=input->eps){
        prod(input->P,input->pagerank,input->N+input->restoN,input->N,temp);
        if(k==2){
            assegnaD(l,temp,input->N+input->restoM);
        }
        if(k==d+2){ p=power(input->c,d); powerExtrd(l,temp,p,input->N+input->restoN);}
        differenceNormd(input->pagerank,temp,input->N+input->restoN,&delta);
        assegnaD(input->pagerank,temp,input->N+input->restoM);
        k=k+1;
    }
    normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
    _mm_free(temp);
    _mm_free(l);
    _mm_free(input->P);
}


void pageranksparseOptDouble(params* input){
    int i,rowel,curcol=0,j,k,d=6;
    double delta,p,somma;
    double* v= _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    double* xk = _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    double* l = _mm_malloc(sizeof(double)*(input->N+input->restoM),16);
    
    for(i=0;i<input->N;i++){
        input->pagerank[i] = 1.0/input->N;
        xk[i]=0.0;
        l[i]=0.0;
        v[i]=0.0;
    }
    for(i=input->N;i<input->N+input->restoM;i++){
        input->pagerank[i]=0.0;
        xk[i]=0.0;
        l[i]=0.0;
        v[i]=0.0;
    }
    delta = 1.0;
    k=1;
    while(delta>=input->eps){
        prodSparseD(input->G.row,input->G.col,v,input->N,input->pagerank,xk,input->G.dd,input->N+input->restoM);
        sommad(xk,input->N+input->restoM,&somma);
        completaSolD(xk,input->N+input->restoM,(1.0-somma)/input->N);
        for(i=input->N;i<input->N+input->restoM;i++){
            xk[i]=0.0;
        }
        if(k==2){
            assegnaD(l,xk,input->N+input->restoM);
        }
        if(k==d+2){ p=power(input->c,d);powerExtrd(l,xk,p,input->N+input->restoN);}
        differenceNormd(input->pagerank,xk,input->N+input->restoN,&delta);
        assegnaSoluzioneD(input->pagerank,xk,input->N+input->restoM);
        k=k+1;
    }
    normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
    _mm_free(l);
    _mm_free(xk);
    _mm_free(v);
}

extern void assegnaF(float* v1, float* v2, int n);

void pageranksparseOptFloat(params* input){
    int i,rowel,curcol=0,j,k,d=6;
    float delta,p,s;
    float* xk = _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
    float* xk1 = _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
    float* l = _mm_malloc(sizeof(float)*(input->N+input->restoM),16);
    float* v= _mm_malloc(sizeof(float)*(input->N+input->restoM),16);

    for(i=0;i<input->N;i++){
        xk1[i] = 1.0/input->N;
        xk[i]=0.0;
        l[i]=0.0;
        v[i]=0.0;
    }
    for(i=input->N;i<input->N+input->restoM;i++){
        input->pagerank[i]=0.0;
        xk[i]=0.0;
        l[i]=0.0;
        v[i]=0.0;
    }
    delta = 1.0;
    k=1;
    while(delta>=input->eps){
        prodSparseF(input->G.row,input->G.col,v,input->N,xk1,xk,input->G.df,input->N+input->restoM);
        sommaf(xk,input->N+input->restoM,&s);
        completaSolF(xk,input->N+input->restoM,(1.0-s)/input->N);
        for(i=input->N;i<input->N+input->restoM;i++){
             xk[i]=0.0;
        }
        if(k==2){
            assegnaF(l,xk,input->N+input->restoM);
        }
        if(k==d+2){ p=power(input->c,d);powerExtrf(l,xk,p,input->N+input->restoN);}
        differenceNormf(xk1,xk,input->N+input->restoN,&delta);
        assegnaSoluzioneF(xk1,xk,input->N+input->restoM);
        k=k+1;
    }
    converti(input->pagerank,xk1,input->N);
    normalizzazione(input->pagerank,input->N+input->restoN,input->N+input->restoM);
    _mm_free(l);
    _mm_free(xk1);
    _mm_free(v);
    _mm_free(xk);
}

void pagerank(params* input) {
    int i;
    input->pagerank=(double*) _mm_malloc((input->N+input->restoM)*sizeof(double),16);
    for(i=input->N;i<input->N+input->restoM;i++) input->pagerank[i]=0.0;

	if (input->format==SPARSE ){
		if(input->prec == SINGLE){
                    input->G=creaGF(input->L,input->c,input->N,input->M,input->restoM);
			if(input->opt == NOPT)
			   pageranksparseNoOptFloat(input);
            else
            	pageranksparseOptFloat(input);
		}
		else{
                    input->G=creaGD(input->L,input->c,input->N,input->M,input->restoM);
		    if(input->opt == NOPT)
				pageranksparseNoOptDouble(input);
			else
				pageranksparseOptDouble(input);
		}
	}
    else{
        if(input->opt== NOPT)
        	pagerankfullNoOpt(input);
        else
        	pagerankfullOpt(input);
    }
}


int main(int argc, char** argv) {
        
	params* input = malloc(sizeof(params));

	input->file_name = NULL;
	input->P = NULL; // dense format
	input->L=NULL;
	input->G.row = NULL; // sparse format
	input->G.col = NULL;
        input->G.df = NULL;
        input->G.dd = NULL;
	input->N = 0; // number of nodes
        input->restoN=0;
	input->M = 0; // number of arcs
        input->restoM=0;
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

    if (input->format == 0){
        if(input->prec == SINGLE)
          input->L= load_sparse(input->file_name, &input->N, &input->M,input->c,input->prec,&input->restoN, &input->restoM);
        else
          input->L= load_sparse(input->file_name, &input->N, &input->M,input->c,input->prec,&input->restoN, &input->restoM);
    }
    else
		input->P = load_dense(input->file_name, &input->N, &input->M, &input->restoN, &input->restoM);

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
		printf("\nExecution time = %.8f seconds\n", ((float)t)/CLOCKS_PER_SEC);
    else{
		printf("%.8f\n", ((float)t)/CLOCKS_PER_SEC);
        
    }

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
