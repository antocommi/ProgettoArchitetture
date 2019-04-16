#include<stdio.h>
#include<string.h>
typedef int boolean;
#define true 1
#define false 0
//prova
int main (int argc, char *argv[]){
    int K=1;
    int m=8;
    int kstar=256;
    int kc=8192;
    int w=16;
    int n;
    int nr=n/20;  
    double soglia=0.01;
    int tmin=10;
    int tmax=100; 
    boolean esaustiva=true;
    boolean simmetrica=true;
    int i=0;

    while(i<argc)
    {
        if(argv[i++]=="-exaustive")
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
    
}
