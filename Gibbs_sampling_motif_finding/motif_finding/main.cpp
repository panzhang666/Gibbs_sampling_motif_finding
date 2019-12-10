#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

using namespace std;
bool myfunction (int i,int j) { return (i<j); }

struct myclass {
    bool operator() (int i,int j) { return (i<j);}
} myobject;

string find_dir_data(float ICPC, int ML, int SL, int SC, int vers) {
    string directory = "data_set_";
    directory.append(to_string(ICPC));
    directory.append("_");
    directory.append(to_string(ML));
    directory.append("_");
    directory.append(to_string(SL));
    directory.append("_");
    directory.append(to_string(SC));
    directory.append("_");
    directory.append(to_string(vers));

    return directory;
}

string find_dir(float ICPC, int ML, int SL, int SC, int vers) {
    string directory = "predicted_data_set_";
    directory.append(to_string(ICPC));
    directory.append("_");
    directory.append(to_string(ML));
    directory.append("_");
    directory.append(to_string(SL));
    directory.append("_");
    directory.append(to_string(SC));
    directory.append("_");
    directory.append(to_string(vers));

    return directory;
}

float caculate(int p , int SC){
    float result;
    float x = (p+0.5)/(SC+1);
    result = x;
    return result;
}
//P[SL-ML]
int pchoose(double p[],double sum,int count,int SL,int ML){
    double sum2 = sum*100000;
    int sum3 = (int)sum2;
    long k = rand()%sum3+1 ;
    double temp = k/100000.0;
    //cout<<temp<<endl;
    for(int i=1;i<count;i++){
        p[i]=p[i]+p[i-1];
       // cout<<p[i]<<endl;
    }
    for(int i =1;i<count;i++){
        if(temp>=p[i]&&temp<p[i+1]){
           // cout<<i+1<<endl;
            return (i+1);
        }
    }
    return 0;
}
//int pchoose(double p[],double sum,int count,int SL,int ML){
//    std::vector<double> myvector (p, p+(SL-ML));
//    std::sort (myvector.begin(), myvector.end(), myobject);
//    double max = myvector.back();
//    for(int i=0;i<count;i++){
//        if (p[i]==max){
//            return i;
//        }
//    }
//    return 0;
//}



int readMLlength(string filename){
    ifstream in;
    string str1;
    vector<string> sequence;
    string file_contents;
    in.open(filename);
    
    while(in) {
        while (std::getline(in, str1)) {
            file_contents.append(str1);
        }
    }
    int result = atoi(file_contents.substr(0,1).c_str());
    return result;
}

vector<string> str(string filename){
    ifstream in;
    string str1;
    vector<string> sequence;
    string file_contents;
    in.open(filename);
    
    while(in) {
        while (std::getline(in, str1)) {
            str1.erase(std::remove(str1.begin(), str1.end(), ' '), str1.end());
            file_contents.append(str1);
            //file_contents.append("  ");
        }
        
    }
    file_contents.append("  ");
    //cout<<file_contents<<endl;
    while(file_contents.at(0)=='<'){
        string temp ="ATCG";
        string temp2;
        std::size_t found_start = file_contents.find_first_of(temp);
        file_contents.erase(0,found_start);
        std::size_t found_end = file_contents.find_first_not_of(temp);
        sequence.push_back(file_contents.substr(0,found_end));
        file_contents.erase(0,found_end);
    }
    
    return sequence;
}
void predictedmotif(vector<string> bindsitesequence,int SC ,int ML,int n){
    int a[SC][4];
    for(int i=0;i<ML;i++){
        for(int j=0;j<4;j++){
            a[i][j]=0;
        }
    }
    for(int i=0;i<SC;i++){
        for(int j=0;j<ML;j++){
            switch(bindsitesequence.at(i).at(j)){
                case 'A': a[j][0]++;break;
                case 'C': a[j][1]++;break;
                case 'G': a[j][2]++;break;
                case 'T': a[j][3]++;break;
            }
        }
    }
    for(int i=0;i<ML;i++){
        //cout<<endl;
        for(int j=0;j<4;j++){
            //cout<<a[i][j]<<"   ";
        }
    }
    string motiffinal;
    string str = to_string(ML);
    motiffinal.append(">MOTIF"+to_string(n)+"    "+str);
    for(int i =0;i<ML;i++){
        motiffinal.append("\n");
        for(int j=0;j<4;j++){
            motiffinal.append(to_string(a[i][j])+"    ");
        }
    }
    motiffinal.append("\n");
    motiffinal.append("<");
    //write down the motif
    //string route1 = "predicted/predictedmotif.txt";
    string route1 = "predicted copy ";
    route1.append(to_string(n));
    route1.append("/predictedmotif.txt");
    std::ofstream out(route1);
    out<<motiffinal;
}
void predictedmotif_better_filename(vector<string> bindsitesequence,float ICPC, int ML, int SL, int SC, int vers){

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

    int a[SC][4];
    for(int i=0;i<ML;i++){
        for(int j=0;j<4;j++){
            a[i][j]=0;
        }
    }
    for(int i=0;i<SC;i++){
        for(int j=0;j<ML;j++){
            switch(bindsitesequence.at(i).at(j)){
                case 'A': a[j][0]++;break;
                case 'C': a[j][1]++;break;
                case 'G': a[j][2]++;break;
                case 'T': a[j][3]++;break;
            }
        }
    }
    for(int i=0;i<ML;i++){
        //cout<<endl;
        for(int j=0;j<4;j++){
            //cout<<a[i][j]<<"   ";
        }
    }
    string motiffinal;
    string str = to_string(ML);
    motiffinal.append(">MOTIF"+identifier+"    "+str);
    for(int i =0;i<ML;i++){
        motiffinal.append("\n");
        for(int j=0;j<4;j++){
            motiffinal.append(to_string(a[i][j])+"    ");
        }
    }
    motiffinal.append("\n");
    motiffinal.append("<");
    //write down the motif
    //string route1 = "predicted/predictedmotif.txt";
    string directory = find_dir(ICPC, ML, SL, SC, vers);
    string route1 = directory;
    route1.append("/predictedmotif.txt");
    std::ofstream out(route1);
    out<<motiffinal;
}
string motif_finding(vector<string> sequence,int ML,int inc){
    int qq=0;
    string result;
    int SL = sequence.at(0).length();
    int SC=sequence.size();
    float a[4][ML];
    float ap[4][ML];
    double q[SL-ML];
    double sum[SC];
    int postion[SC][2];
    int index=0;
    for(int i=0;i<SC;i++){
        for(int tt =0;tt<2;tt++){
        postion[i][tt]=0;
        }
    }
    for(int i=0;i<SL-ML;i++){
        q[i]=1;
    }
    for(int i=0;i<4;i++){
        for(int j=0;j<ML;j++){
            a[i][j]=0;
            ap[i][j]=0;
        }
    }
    vector<string> bindsitesequence;
    for(int i=0;i<SC;i++){
        int p =rand()%(SL-ML);
        string bindsite = sequence.at(i).substr(p,ML);
        //cout<<bindsite<<endl;
        bindsitesequence.push_back(bindsite);
    }
    //start iterate
    while(qq==0&&index<20000){
    for(int k=0;k<SC;k++){
        sum[k]=0;
        qq=0;
        for(int i=0;i<SL-ML;i++){
            q[i]=1;
        }
        for(int i=0;i<4;i++){
            for(int j=0;j<ML;j++){
                a[i][j]=0;
                ap[i][j]=0;
            }
        }
        for(int j=0;j<SC;j++){
            for(int i=0;i<SL-ML;i++){
                q[i]=1;
            }
            if(j!=k){
                for(int t=0;t<ML;t++){
                    switch (bindsitesequence.at(j).at(t))
                    {
                        case 'A': a[0][t]++;
                            break;
                        case 'C': a[1][t]++;
                            break;
                        case 'G': a[2][t]++;
                            break;
                        case 'T': a[3][t]++;
                            break;
                    }
                }
            }
            for(int i=0;i<4;i++)
                for(int j=0;j<ML;j++){
                    ap[i][j]= caculate(a[i][j],SC);
                }
        }
        for(int m1=0;m1<(SL-ML);m1++){
            for(int m2=0;m2<ML;m2++){
                switch (sequence.at(k).at(m1+m2))
                {
                    case 'A': q[m1]=q[m1]*ap[0][m2]*4;
                        break;
                    case 'C': q[m1]=q[m1]*ap[1][m2]*4;
                        break;
                    case 'G': q[m1]=q[m1]*ap[2][m2]*4;
                        break;
                    case 'T': q[m1]=q[m1]*ap[3][m2]*4;
                        break;
                }
            }
            sum[k]=sum[k]+q[m1];
        }
        int positionchoice = pchoose(q,sum[k],(SL-ML),SL,ML);
  //      cout<<positionchoice<<endl;
        //cout<<sequence.at(k).substr(positionchoice,ML)<<endl;
        vector<string> temp = bindsitesequence;
        bindsitesequence.clear();
        for(int ii=0;ii<k;ii++){
            bindsitesequence.push_back(temp.at(ii));
        }
        bindsitesequence.push_back(sequence.at(k).substr(positionchoice,ML));
        for(int ii=k;ii<SC;ii++){
            bindsitesequence.push_back(temp.at(ii));
        }
        postion[k][0]=postion[k][1];
        postion[k][1]=positionchoice;
    }
        int ss=0;
        for(int r=0;r<SC;r++){
            if(postion[r][0]!=postion[r][1]){
                
                ss++;
            }
        
        }
        if(ss==0){
            qq=1;
        }
        if(SC!=5){
        int a[SC][4];
        for(int i=0;i<ML;i++){
            for(int j=0;j<4;j++){
                a[i][j]=0;
            }
        }
        for(int i=0;i<SC;i++){
            for(int j=0;j<ML;j++){
                switch(bindsitesequence.at(i).at(j)){
                    case 'A': a[j][0]++;break;
                    case 'C': a[j][1]++;break;
                    case 'G': a[j][2]++;break;
                    case 'T': a[j][3]++;break;
                }
            }
        }
        int aa=0;
                for(int i=0;i<ML;i++){
                    for(int j=0;j<4;j++){
                        if(a[i][j]>=(SC-2)){
                            aa++;
                           // cout<<a[i][j]<<endl;
                        }
                    }
                }
        if(aa!=ML){
            qq=0;
        }
        }
        index++;
    }
    
    
    
    //the next step
    int b[SC];
    string location;
    for(int i=0;i<SC;i++){
        //cout<<postion[i][1]<<endl;
        b[i]=postion[i][1];
        string loc = to_string(b[i]);
        //cout<<bindsitesequence.at(i)<<endl;
        //create location of data sets
        location.append(loc);
        location.append("\n");
    }
    //string route2 = "predicted/predictedsites.txt";
    string route2 = "predicted copy ";
    route2.append(to_string(inc));
    route2.append("/predictedsites.txt");
    std::ofstream out2(route2);
    out2<<location;
    predictedmotif(bindsitesequence,SC,ML,inc);
    return result;
}

string motif_finding_better_filename(int randSeed, vector<string> sequence, float ICPC, int ML, int vers){
    int qq=0;
    string result;
    int SL = sequence.at(0).length();
    int SC=sequence.size();
    float a[4][ML];
    float ap[4][ML];
    double q[SL-ML];
    double sum[SC];
    int postion[SC][2];

    srand(randSeed);

    int index=0;
    for(int i=0;i<SC;i++){
        for(int tt =0;tt<2;tt++){
        postion[i][tt]=0;
        }
    }
    for(int i=0;i<SL-ML;i++){
        q[i]=1;
    }
    for(int i=0;i<4;i++){
        for(int j=0;j<ML;j++){
            a[i][j]=0;
            ap[i][j]=0; }
    }
    vector<string> bindsitesequence;
    for(int i=0;i<SC;i++){
        int p =rand()%(SL-ML);
        string bindsite = sequence.at(i).substr(p,ML);
        //cout<<bindsite<<endl;
        bindsitesequence.push_back(bindsite);
    }
    //start iterate
    while(qq==0&&index<20000){
    for(int k=0;k<SC;k++){
        sum[k]=0;
        qq=0;
        for(int i=0;i<SL-ML;i++){
            q[i]=1;
        }
        for(int i=0;i<4;i++){
            for(int j=0;j<ML;j++){
                a[i][j]=0;
                ap[i][j]=0;
            }
        }
        for(int j=0;j<SC;j++){
            for(int i=0;i<SL-ML;i++){
                q[i]=1;
            }
            if(j!=k){
                for(int t=0;t<ML;t++){
                    switch (bindsitesequence.at(j).at(t))
                    {
                        case 'A': a[0][t]++;
                            break;
                        case 'C': a[1][t]++;
                            break;
                        case 'G': a[2][t]++;
                            break;
                        case 'T': a[3][t]++;
                            break;
                    }
                }
            }
            for(int i=0;i<4;i++)
                for(int j=0;j<ML;j++){
                    ap[i][j]= caculate(a[i][j],SC);
                }
        }
        for(int m1=0;m1<(SL-ML);m1++){
            for(int m2=0;m2<ML;m2++){
                switch (sequence.at(k).at(m1+m2))
                {
                    case 'A': q[m1]=q[m1]*ap[0][m2]*4;
                        break;
                    case 'C': q[m1]=q[m1]*ap[1][m2]*4;
                        break;
                    case 'G': q[m1]=q[m1]*ap[2][m2]*4;
                        break;
                    case 'T': q[m1]=q[m1]*ap[3][m2]*4;
                        break;
                }
            }
            sum[k]=sum[k]+q[m1];
        }
        int positionchoice = pchoose(q,sum[k],(SL-ML),SL,ML);
  //      cout<<positionchoice<<endl;
        //cout<<sequence.at(k).substr(positionchoice,ML)<<endl;
        vector<string> temp = bindsitesequence;
        bindsitesequence.clear();
        for(int ii=0;ii<k;ii++){
            bindsitesequence.push_back(temp.at(ii));
        }
        bindsitesequence.push_back(sequence.at(k).substr(positionchoice,ML));
        for(int ii=k;ii<SC;ii++){
            bindsitesequence.push_back(temp.at(ii));
        }
        postion[k][0]=postion[k][1];
        postion[k][1]=positionchoice;
    }
        int ss=0;
        for(int r=0;r<SC;r++){
            if(postion[r][0]!=postion[r][1]){
                
                ss++;
            }
        
        }
        if(ss==0){
            qq=1;
        }
        if(SC!=5){
        int a[SC][4];
        for(int i=0;i<ML;i++){
            for(int j=0;j<4;j++){
                a[i][j]=0;
            }
        }
        for(int i=0;i<SC;i++){
            for(int j=0;j<ML;j++){
                switch(bindsitesequence.at(i).at(j)){
                    case 'A': a[j][0]++;break;
                    case 'C': a[j][1]++;break;
                    case 'G': a[j][2]++;break;
                    case 'T': a[j][3]++;break;
                }
            }
        }
        int aa=0;
                for(int i=0;i<ML;i++){
                    for(int j=0;j<4;j++){
                        if(a[i][j]>=(SC-2)){
                            aa++;
                           // cout<<a[i][j]<<endl;
                        }
                    }
                }
        if(aa!=ML){
            qq=0;
        }
        }
        index++;
    }
    
    
    
    //the next step
    int b[SC];
    string location;
    for(int i=0;i<SC;i++){
        //cout<<postion[i][1]<<endl;
        b[i]=postion[i][1];
        string loc = to_string(b[i]);
        //cout<<bindsitesequence.at(i)<<endl;
        //create location of data sets
        location.append(loc);
        location.append("\n");
    }
    //string route2 = "predicted/predictedsites.txt";
    string directory = find_dir(ICPC, ML, SL, SC, vers);
    string route2 = directory;
    route2.append("/predictedsites.txt");
    std::ofstream out2(route2);
    out2<<location;
    predictedmotif_better_filename(bindsitesequence,ICPC, ML, SL, SC, vers);
    return result;
}


    

int main(int argc, const char * argv[]) {
    vector<string> sequence;

    struct timespec start, finish;
    double elapsed;

    srand(100);

    if (argc == 1) {
        for(int n=1;n<71;n++){
            string directory = "data_set copy ";
            directory.append(to_string(n));
            int ML = readMLlength(directory + "/motiflength.txt");
            sequence = str(directory + "/sequences.fa");
            //sequence = str("data_set/sequences.fa");
            //    for(int i=0;i<sequence.size();i++){
            //        cout<<sequence.at(i)<<endl;}
            //cout<<ML<<endl;


            clock_gettime(CLOCK_MONOTONIC, &start);
            string result = motif_finding(sequence,ML,n);
            clock_gettime(CLOCK_MONOTONIC, &finish);
            elapsed = finish.tv_sec - start.tv_sec;
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            std::cout<<n<<": "<<elapsed<<std::endl;

            string directory2 = "predicted copy ";
            directory2.append(to_string(n));
            std::ofstream out(directory2 + "/time.txt");
            out << elapsed;
            out.close();
        }
    }

    //argv = icpc, ml, sl, sc
    
    if (argc != 9) return -1;

    if (argv[1] != 0) {
        float ICPC=0;
        int ML=8;
        int SL=500;
        int SC=10;
        for (ICPC = stod(argv[1]); ICPC<stod(argv[2]); ICPC+=0.1){
            for (int vers = 1; vers < 11; vers++) {
                string directory = find_dir_data(ICPC, ML, SL, SC, vers);
                int ML = readMLlength(directory + "/motiflength.txt");
                sequence = str(directory + "/sequences.fa");

                string randFile = directory + "/rand.txt";

                ifstream randStream(randFile);
                string randString;
                std::getline(randStream, randString);
                int randSeed = stoi(randString, nullptr);

                clock_gettime(CLOCK_MONOTONIC, &start);
                string result = motif_finding_better_filename(randSeed, sequence,ICPC, ML, vers);
                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = finish.tv_sec - start.tv_sec;
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                std::cout<<vers<<": "<<elapsed<<std::endl;

                string dir_target = find_dir(ICPC, ML, SL, SC, vers);
                std::ofstream out(directory + "/time.txt");
                out << elapsed;
                out.close();
            }
        }
    }

    if (argv[5] != 0) {
        float ICPC=2;
        int ML=8;
        int SL=500;
        int SC=10;
        for (SC = stoi(argv[5]); SC<stoi(argv[6]); SC+=1){
            for (int vers = 1; vers < 11; vers++) {
                string directory = find_dir_data(ICPC, ML, SL, SC, vers);
                int ML = readMLlength(directory + "/motiflength.txt");
                sequence = str(directory + "/sequences.fa");

                string randFile = directory + "/rand.txt";

                ifstream randStream(randFile);
                string randString;
                std::getline(randStream, randString);
                int randSeed = stoi(randString, nullptr);

                clock_gettime(CLOCK_MONOTONIC, &start);
                string result = motif_finding_better_filename(randSeed, sequence,ICPC, ML, vers);
                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = finish.tv_sec - start.tv_sec;
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                std::cout<<vers<<": "<<elapsed<<std::endl;

                string dir_target = find_dir(ICPC, ML, SL, SC, vers);
                std::ofstream out(directory + "/time.txt");
                out << elapsed;
                out.close();
            }
        }
    }
    return 0;
}
