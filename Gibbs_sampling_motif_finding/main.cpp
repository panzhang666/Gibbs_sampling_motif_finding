#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include "SC_random.hpp"
#include "random_motif.hpp"
#include <vector>
using namespace std;

int main(int argc, const char * argv[]) {
    srand (10);
    string SC_rand ;
    float ICPC;
    int ML;
    int SL;
    int SC;
//    cout<<"please input information content per column(ICPC)";
//    cin>>ICPC;
//    cout<<"please input motif length(ML)";
//    cin>>ML;
//    cout<<"please input sequence length(SL)";
//    cin>>SL;
//    cout<<"please input sequence count(SC)";
//    cin>>SC;
    ICPC=2;
    ML=8;
    SL=500;
    SC=10;
    for(int n=1;n<11;n++){
    vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
    vector<string> sequence_ini ;
    vector<string> sequence_planted;
    for(int i=0;i<SC;i++){
        sequence_ini.push_back(SC_random(SL));
    }
    sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
//    Let the default parameter combination be “ICPC = 2, ML = 8, SL = 500, SC = 10”
//    a) set ICPC = 1, 1.5, 2, while other parameters are at default values.
//        b) Set ML = 6,7,8, while other parameters are at default values.
//            c) Set SC = 5, 10, 20, while other parameters are at default values.
    ICPC=1;
    for(int n=11;n<21;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
    ICPC=1.5;
    for(int n=21;n<31;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
    ICPC=2;
    ML=7;
    for(int n=31;n<41;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
    ML=6;
    for(int n=41;n<51;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
    ML=8;
    SC=20;
    for(int n=51;n<61;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }
    SC=5;
    for(int n=61;n<71;n++){
        vector<string> bindsite_vec = random_motif(ICPC, ML,SC,n);
        vector<string> sequence_ini ;
        vector<string> sequence_planted;
        for(int i=0;i<SC;i++){
            sequence_ini.push_back(SC_random(SL));
        }
        sequence_planted = plant_site(bindsite_vec,sequence_ini,SC,ML,SL,n);
    }

    ICPC=0;
    ML=8;
    SL=500;
    SC=10;
    for (ICPC = 0.1; ICPC<2.1; ICPC+=0.1){
        for (int vers = 1; vers < 11; vers++) {
            vector<string> bindsite_vec = random_motif_better_filename(ICPC, ML, SL, SC, vers);
            vector<string> sequence_ini ;
            vector<string> sequence_planted;
            for(int i=0;i<SC;i++){
                sequence_ini.push_back(SC_random(SL));
            }
            sequence_planted = plant_site_better_filename(bindsite_vec,sequence_ini,ICPC, SC,ML,SL, vers);
        }
    }

    ICPC=2;
    ML=8;
    SL=500;
    SC=10;
    for (SC = 5; SC < 21; SC+=1){
        for (int vers = 1; vers < 11; vers++) {
            vector<string> bindsite_vec = random_motif_better_filename(ICPC, ML, SL, SC, vers);
            vector<string> sequence_ini ;
            vector<string> sequence_planted;
            for(int i=0;i<SC;i++){
                sequence_ini.push_back(SC_random(SL));
            }
            sequence_planted = plant_site_better_filename(bindsite_vec,sequence_ini,ICPC, SC,ML,SL, vers);
        }
    }

}
