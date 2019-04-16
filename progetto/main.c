#include<stdio.h>
#include<string.h>

typedef int boolean;
#define true 1
#define false 0

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

    int** Y;
    Y=(int**)calloc(sizeof(int*)*n);
    for(i=0; i<n; i++){
        Y[i]=(int*)calloc(sizeof(int)*d);
    }

    Y[0][0]=1;
    Y[0][1]=2;
    Y[0][2]=3;
    Y[0][3]=2;
    Y[0][4]=3;

    Y[1][0]=1;
    Y[1][1]=2;
    Y[1][2]=3;
    Y[1][3]=2;
    Y[1][4]=3;

    Y[2][0]=1;
    Y[2][1]=2;
    Y[2][2]=3;
    Y[2][3]=2;
    Y[2][4]=3;

    Y[3][0]=1;
    Y[3][1]=2;
    Y[3][2]=3;
    Y[3][3]=2;
    Y[3][4]=3;

    Y[4][0]=1;
    Y[4][1]=2;
    Y[4][2]=3;
    Y[4][3]=2;
    Y[4][4]=3;

    Y[5][0]=1;
    Y[5][1]=2;
    Y[5][2]=3;
    Y[5][3]=2;
    Y[5][4]=3;

    Y[6][0]=1;
    Y[6][1]=2;
    Y[6][2]=3;
    Y[6][3]=2;
    Y[6][4]=3;

    Y[7][0]=1;
    Y[7][1]=2;
    Y[7][2]=3;
    Y[7][3]=2;
    Y[7][4]=3;

    Y[8][0]=1;
    Y[8][1]=2;
    Y[8][2]=3;
    Y[8][3]=2;
    Y[8][4]=3;

    Y[9][0]=1;
    Y[9][1]=2;
    Y[9][2]=3;
    Y[9][3]=2;
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
    int** codebook;
    codebook=kmeans(Y, d, K);
}

int** kmeans(int** Y, int d, int K){
    int i;
    int** codebook=(int**)calloc(sizeof(int*)*K);
    for(i=0; i<K, i++){
        codebook[i]=(int*)calloc(sizeof(int)*d);
    }

    return codebook;
}
