#include <stdio.h>

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
    int esaustiva=1;
    int simmetrica=1;
    int i=0;

    while(i<argc){
        if(argv[i]=="-exaustive"){
            esaustiva=1;
            i++;
        }else if(argv[i]=="-noexaustive"){
            esaustiva=0;
            i++;
        }else if(argv[i]=="-sdc"){
            simmetrica=1;
            i++;
        }else if(argv[i]=="-adc"){
            simmetrica=0;
            i++;
        }else if(argv[i]=="-knn"){
            K=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-m"){
            m=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-k"){
            kstar=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-kc"){
            kc=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-w"){
            w=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-nr"){
            nr=atoi(argv[i+1]);
            i+=2;
        }else if(argv[i]=="-kmeans"){
            soglia=atoi(argv[i+1]);
            i+=2
            if(i<argc){
                tmin=atoi(argv[i]);
                i++;
            }
            if(i<argc){
                tmax=atoi(argv[i]);
                i++;
            }
        }
    }

    
}