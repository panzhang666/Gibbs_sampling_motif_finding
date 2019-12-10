#ifndef random_motif_hpp
#define random_motif_hpp

#include <stdio.h>
#include <vector>
using namespace std;
#include <stdio.h>
#include <iostream>
#include <string>
vector<string> random_motif(float ICPC,int ML, int SC, int n);
vector<string> random_motif_better_filename(float ICPC,int ML,int SL, int SC, int vers);
float producep(float ICPC);
int *generateoneline(float p ,int n);
string generatebindsite(int numA,int numC ,int numG ,int numT);

#endif /* random_motif_hpp */
