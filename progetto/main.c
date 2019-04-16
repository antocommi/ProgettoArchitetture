#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int boolean;
#define true 1
#define false 0

double dist(int* p1, int* p2, int d){
    int i;
    double dist=0;

    for(i=0; i<d; i++){
        dist+=(p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(dist);
}

int calcolaQ(int* x, int** codebook, int K, int d){
    int i;
    double min=0.0;
    int imin=-1;
    double temp;

    for(i=0; i<K; i++){
        temp=dist(x, codebook[i], d);
        if(temp<min){
            min=temp;
            imin=i;
        }
    }
    return imin;
}

int** kmeans(int** Y, int n, int d, int K){
    int i;
    int** codebook=(int**) calloc(K, sizeof(int*));
    
    if(codebook==NULL) exit(-1);
    for(i=0; i<K; i++){
        codebook[i]=Y[rand()%n];
    }
    int* q=(int*) calloc(n, sizeof(int));
    for(i=0; i<n; i++){
        q[i]=calcolaQ(Y[i], codebook, K, d);
    }
    //da completare
    return codebook;
}

int main (int argc, char *argv[]){
    int K=1;
    int m=8;
    int kstar=256;
    int kc=8192;
    int w=16;
    int n=10;
    int d=5;
    int nr=n/20;  
    double soglia=0.01;
    int tmin=10;
    int tmax=100; 
    boolean esaustiva=true;
    boolean simmetrica=true;
    int i=0;
    int j=0;
    int** codebook;
    int** Y;

    Y=(int**) calloc(n, sizeof(int*));
    if(Y==NULL) exit(-1);
    for(i=0; i<n; i++){
        Y[i]=(int*) calloc(d, sizeof(int));
        if(Y[i]==NULL) exit(-1);
    }

    Y[0][0]=1;
    Y[0][1]=2;
    Y[0][2]=3;
    Y[0][3]=2;
    Y[0][4]=3;

    Y[1][0]=1;
    Y[1][1]=8;
    Y[1][2]=7;
    Y[1][3]=9;
    Y[1][4]=3;

    Y[2][0]=4;
    Y[2][1]=2;
    Y[2][2]=5;
    Y[2][3]=2;
    Y[2][4]=3;

    Y[3][0]=6;
    Y[3][1]=7;
    Y[3][2]=3;
    Y[3][3]=2;
    Y[3][4]=3;

    Y[4][0]=1;
    Y[4][1]=9;
    Y[4][2]=3;
    Y[4][3]=8;
    Y[4][4]=3;

    Y[5][0]=7;
    Y[5][1]=2;
    Y[5][2]=3;
    Y[5][3]=2;
    Y[5][4]=3;

    Y[6][0]=6;
    Y[6][1]=6;
    Y[6][2]=3;
    Y[6][3]=6;
    Y[6][4]=3;

    Y[7][0]=7;
    Y[7][1]=2;
    Y[7][2]=3;
    Y[7][3]=7;
    Y[7][4]=3;

    Y[8][0]=9;
    Y[8][1]=2;
    Y[8][2]=7;
    Y[8][3]=2;
    Y[8][4]=8;

    Y[9][0]=11;
    Y[9][1]=3;
    Y[9][2]=33;
    Y[9][3]=22;
    Y[9][4]=3;

    while(i<argc){
        if(strcmp(argv[i++], "-exaustive"))
            continue;
        else if(strcmp(argv[i++], "-noexaustive"))
            esaustiva=false;
        else if(strcmp(argv[i++], "-sdc"))
            continue;
        else if(strcmp(argv[i++], "-adc"))
            simmetrica=false;
        else if(strcmp(argv[i++], "-knn"))
            K=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-m"))
            m=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-k"))
            kstar=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-kc"))
            kc=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-w"))
            w=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-nr"))
            nr=atoi(argv[i++]);
        else if(strcmp(argv[i++], "-kmeans")){
            soglia=atoi(argv[i++]);
            if(i<argc)
                tmin=atoi(argv[i++]);
            if(i<argc)
                tmax=atoi(argv[i++]);
        }
    }

    codebook=kmeans(Y, n, d, K);

    //da continuare
}