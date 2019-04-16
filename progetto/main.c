#include <stdio.h>
typedef int boolean;
#define true 1
#define false 0

main (int argc, char *argv[]){
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

    while(i<argc){
        if(argv[i++]=="-exaustive")
            continue;
        else if(argv[i++]=="-noexaustive")
            esaustiva=false;
        else if(argv[i++]=="-sdc")
            continue;
        else if(argv[i++]=="-adc")
            simmetrica=false;
        else if(argv[i++]=="-knn")
            K=atoi(argv[i++]);
        else if(argv[i++]=="-m")
            m=atoi(argv[i++]);
        else if(argv[i++]=="-k")
            kstar=atoi(argv[i++]);
        else if(argv[i++]=="-kc")
            kc=atoi(argv[i++]);
        else if(argv[i++]=="-w")
            w=atoi(argv[i++]);
        else if(argv[i++]=="-nr")
            nr=atoi(argv[i++]);
        else if(argv[i++]=="-kmeans"){
            soglia=atoi(argv[i++]);
            if(i<argc)
                tmin=atoi(argv[i++]);
            if(i<argc)
                tmax=atoi(argv[i++]);
        }
    }
    
}
