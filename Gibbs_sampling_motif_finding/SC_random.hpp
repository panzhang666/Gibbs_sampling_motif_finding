#ifndef SC_random_hpp
#define SC_random_hpp
using namespace std;
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
std::string SC_random(int L);
vector<string> plant_site(vector<string> bindsite_vec,vector<string> sequence_ini,int SC,int ML,int SL,int n);
vector<string> plant_site_better_filename(vector<string> bindsite_vec,vector<string> sequence_ini, float ICPC, int SC,int ML,int SL, int vers);
#endif /* SC_random_hpp */
