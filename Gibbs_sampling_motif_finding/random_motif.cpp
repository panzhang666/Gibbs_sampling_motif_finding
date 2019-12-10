#include "random_motif.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>

using namespace std;

float producep(float ICPC){
    float p;

    switch (ICPC) {
        case 1.0:
            p = 0.8105;
            break;
        case 1.5:
            p = 0.9245;
            break;
        case 2.0:
            p = 1;
            break;
    }
    return p;

    if(ICPC==1.0){
        p=0.8105;
    }
    else{
        if(ICPC==1.5){
            p=0.9245;
        }
        else{
            p=1.0;
        }
    }
    return p;
}
int* generateoneline(float p ,int n){
    int v1;
    static int a[4];
    a[0]=0;
    a[1]=0;
    a[2]=0;
    a[3]=0;
    int v2;
    p=p*10000;
    v2=rand()%4+1 ;
    for(int i=0;i<n;i++){
        v1=rand()%10000+1 ;
        if(v1<p){
            v1=1;
        }
        if(v1>=p&&v1<p+(1.0/3.0)*(10000-p)){
            v1=2;
        }
        if(v1>=p+(1.0/3.0)*(10000-p)&&v1<p+(2.0/3.0)*(10000-p)){
            v1=3;
        }
        if(v1>=p+(2.0/3.0)*(10000-p)){
            v1=4;
        }
        if(v2==1){
            switch(v1){
                case 1: a[0]++;break;
                case 2: a[1]++;break;
                case 3: a[2]++;break;
                case 4: a[3]++;break;
            }
        }
        if(v2==2){
            switch(v1){
                case 1: a[1]++;break;
                case 2: a[0]++;break;
                case 3: a[2]++;break;
                case 4: a[3]++;break;
            }
        }
        if(v2==3){
            switch(v1){
                case 1: a[2]++;break;
                case 2: a[1]++;break;
                case 3: a[0]++;break;
                case 4: a[3]++;break;
            }
        }
        if(v2==4){
            switch(v1){
                case 1: a[3]++;break;
                case 2: a[1]++;break;
                case 3: a[2]++;break;
                case 4: a[0]++;break;
            }
        }
    }
    return a;
}

string generatebindsite(int numA,int numC ,int numG ,int numT,int SC){
    string bindsite;
    int v1=rand()%SC+1 ;
    if(v1<=numA){
        bindsite.append("A");
    }
    if(v1>numA&&v1<=(numA+numC)){
        bindsite.append("C");
    }
    if(v1>(numA+numC)&&v1<=(numA+numC+numG)){
        bindsite.append("G");
    }
    if(v1>(numA+numC+numG)&&v1<=(numA+numC+numG+numT)){
        bindsite.append("T");
    }
    return bindsite;
}

vector<string> random_motif(float ICPC,int ML,int SC,int n){
    //create p value
    float p = producep(ICPC);
//    int n=20;
    int* a;
    int b[ML][4];
    //generate motif
    string bindsite;
    for(int i =0;i<ML;i++){
        a = generateoneline(p,SC);
        //cout<<*a<<"     "<<*(a+1)<<"     "<<*(a+2)<<"    "<<*(a+3)<<endl;
        //generate motif matrix
        for(int j=0;j<4;j++){
            b[i][j]=*(a+j);
        }
    }
    //gererate bindsite according to motif matrix
    vector<string> bindsite_vec;
    string motiffinal;
    for(int k=0;k<SC;k++){
        bindsite = "";
        for(int i=0;i<ML;i++){
            int numA=b[i][0];
            int numC=b[i][1];
            int numG=b[i][2];
            int numT=b[i][3];
            string t = generatebindsite(numA, numC ,numG ,numT,SC);
            bindsite = bindsite.append(t);
            if(i==ML-1){
                //cout<<bindsite<<endl;
            }
        }
        bindsite_vec.push_back(bindsite);
    }
    //output the motif
    string str = to_string(ML);
    motiffinal.append(">MOTIF"+to_string(n)+"    "+str);
    for(int i =0;i<ML;i++){
        motiffinal.append("\n");
        for(int j=0;j<4;j++){
            motiffinal.append(to_string(b[i][j])+"    ");
        }
    }
    motiffinal.append("\n");
    motiffinal.append("<");
    cout<<motiffinal<<endl;
    int mm=n;
    //write down the motif
    string route1 = "motif_finding/data_set copy ";
    route1.append(to_string(n));
    route1.append("/motif.txt");
    std::ofstream out(route1);
    out<<motiffinal;
    //write down the motif length
    string route2 = "motif_finding/data_set copy ";
    route2.append(to_string(mm));
    route2.append("/motiflength.txt");
    std::ofstream out2(route2);
    out2<<to_string(ML);
    
    return bindsite_vec;
}

vector<string> random_motif_better_filename(float ICPC,int ML,int SL, int SC, int vers){

    string identifier = "";
    identifier.append(to_string(ICPC));
    identifier.append("_");
    identifier.append(to_string(ML));
    identifier.append("_");
    identifier.append(to_string(SL));
    identifier.append("_");
    identifier.append(to_string(SC));
    identifier.append("_");
    identifier.append(to_string(vers));

    //create p value
    float p = producep(ICPC);
//    int n=20;
    int* a;
    int b[ML][4];
    //generate motif
    string bindsite;
    for(int i =0;i<ML;i++){
        a = generateoneline(p,SC);
        //cout<<*a<<"     "<<*(a+1)<<"     "<<*(a+2)<<"    "<<*(a+3)<<endl;
        //generate motif matrix
        for(int j=0;j<4;j++){
            b[i][j]=*(a+j);
        }
    }
    //gererate bindsite according to motif matrix
    vector<string> bindsite_vec;
    string motiffinal;
    for(int k=0;k<SC;k++){
        bindsite = "";
        for(int i=0;i<ML;i++){
            int numA=b[i][0];
            int numC=b[i][1];
            int numG=b[i][2];
            int numT=b[i][3];
            string t = generatebindsite(numA, numC ,numG ,numT,SC);
            bindsite = bindsite.append(t);
            if(i==ML-1){
                //cout<<bindsite<<endl;
            }
        }
        bindsite_vec.push_back(bindsite);
    }
    //output the motif
    string str = to_string(ML);
    motiffinal.append(">MOTIF"+identifier+"    "+str);
    for(int i =0;i<ML;i++){
        motiffinal.append("\n");
        for(int j=0;j<4;j++){
            motiffinal.append(to_string(b[i][j])+"    ");
        }
    }
    motiffinal.append("\n");
    motiffinal.append("<");
    //cout<<motiffinal<<endl;
    //write down the motif
    string directory = "motif_finding/data_set_";
    directory.append(identifier);

    string route1 = directory;
    route1.append("/motif.txt");
    std::ofstream out(route1);
    out<<motiffinal;
    //write down the motif length
    string route2 = directory;
    route2.append("/motiflength.txt");
    std::ofstream out2(route2);
    out2<<to_string(ML);
    
    return bindsite_vec;
}

