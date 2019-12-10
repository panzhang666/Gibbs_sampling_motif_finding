#include "SC_random.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>

using namespace std;

string SC_random(int L) {
    string S1;
    int a[L];
    // insert code here...
    for(int i=0;i<L;i++){
        a[i]=rand()%100+1 ;
        if(a[i]<26){
            S1.append("A");
        }
        if(a[i]>25&&a[i]<51){
            S1.append("T");
        }
        if(a[i]>50&&a[i]<76){
            S1.append("C");
        }
        if(a[i]>75&&a[i]<101){
            S1.append("G");
        }
    }
 // cout << S1<<endl;
    return S1;
}

vector<string> plant_site(vector<string> bindsite_vec,vector<string> sequence_ini, int SC,int ML,int SL,int n){
    vector<string> sequence_planted;
    string finalsequence;
    string location;
    for(int i=0;i<SC;i++){
        int p =rand()%(SL-ML);
        string sequenceold=sequence_ini.at(i);
        string part1;
        string part2;
        string part3;
        part1 = sequenceold.substr(0,p);
        part2 = bindsite_vec.at(i);
        part3 = sequenceold.substr(p+ML,SL);
        string sequencenew="";
        sequencenew.append(part1);
        sequencenew.append(part2);
        sequencenew.append(part3);
        sequence_planted.push_back(sequencenew);
        //cout<<sequencenew<<endl;
        //convert int to string;
        string str = to_string(i);
        string loc = to_string(p);
        //create finalsequence data sets
        finalsequence.append("<sequence");
        finalsequence.append(str);
        finalsequence.append("\n");
        finalsequence.append(sequencenew);
        finalsequence.append("\n");
        //create location of data sets
        location.append(loc);
        location.append("\n");
    }
    int mm=n;
    cout<<finalsequence<<endl;
    string route1 = "motif_finding/data_set copy ";
    route1.append(to_string(n));
    route1.append("/sequences.fa");
    std::ofstream out(route1);
    out<<finalsequence;
    cout<<location<<endl;
    string route2 = "motif_finding/data_set copy ";
    route2.append(to_string(mm));
    route2.append("/sites.txt");
    std::ofstream out2(route2);
    out2<<location;
    
    return sequence_planted;
}

vector<string> plant_site_better_filename(vector<string> bindsite_vec,vector<string> sequence_ini, float ICPC, int SC,int ML,int SL, int vers){
    vector<string> sequence_planted;
    string finalsequence;
    string location;
    for(int i=0;i<SC;i++){
        int p =rand()%(SL-ML);
        string sequenceold=sequence_ini.at(i);
        string part1;
        string part2;
        string part3;
        part1 = sequenceold.substr(0,p);
        part2 = bindsite_vec.at(i);
        part3 = sequenceold.substr(p+ML,SL);
        string sequencenew="";
        sequencenew.append(part1);
        sequencenew.append(part2);
        sequencenew.append(part3);
        sequence_planted.push_back(sequencenew);
        //cout<<sequencenew<<endl;
        //convert int to string;
        string str = to_string(i);
        string loc = to_string(p);
        //create finalsequence data sets
        finalsequence.append("<sequence");
        finalsequence.append(str);
        finalsequence.append("\n");
        finalsequence.append(sequencenew);
        finalsequence.append("\n");
        //create location of data sets
        location.append(loc);
        location.append("\n");
    }
    //cout<<finalsequence<<endl;
    string directory = "motif_finding/data_set_";
    directory.append(to_string(ICPC));
    directory.append("_");
    directory.append(to_string(ML));
    directory.append("_");
    directory.append(to_string(SL));
    directory.append("_");
    directory.append(to_string(SC));
    directory.append("_");
    directory.append(to_string(vers));
    cout<<directory<<endl;
    string route1 = directory;
    route1.append("/sequences.fa");
    std::ofstream out(route1);
    out<<finalsequence;
    out.close();
    cout<<location<<endl;
    string route2 = directory;
    route2.append("/sites.txt");
    std::ofstream out2(route2);
    out2<<location;
    out2.close();

    int rand_seed = rand();
    string rand_file = directory;
    rand_file.append("/rand.txt");
    std::ofstream rout(rand_file);
    rout << rand_seed;
    rout.close();
    
    return sequence_planted;
}

//  cout<<S1<<endl;
//    string route1;
//    string route2;
//    cout<<"please input the first route and file name of the sequence"<<endl;
//   cin>>route1;
//    std::ofstream out(route1);
//    out << S2;
//    cout<<"please input the second route and file name of the sequence"<<endl;
//    cin>>route2;
//    std::ofstream out2(route2);
//    out2 << S3;
