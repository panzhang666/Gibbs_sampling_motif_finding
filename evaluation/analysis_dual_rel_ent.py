from __future__ import division
import math

import numpy as np
from matplotlib import pyplot as plt

def read_motif(filename, nLines):
    with open(filename, 'r') as f:
        f.readline()
        motif = {}

        motif['nLines'] = nLines
        vals = f.readline().strip().split()
        motif['0_A'] = int(vals[0])
        motif['0_C'] = int(vals[1])
        motif['0_G'] = int(vals[2])
        motif['0_T'] = int(vals[3])
        motif['total'] = motif['0_A']+motif['0_C']+motif['0_G']+motif['0_T']
        motif['0_A'] /= motif['total']
        motif['0_C'] /= motif['total']
        motif['0_G'] /= motif['total']
        motif['0_T'] /= motif['total']

        for i in range(1, nLines):
            vals = f.readline().strip().split()
            motif[str(i)+'_A'] = int(vals[0]) / motif['total']
            motif[str(i)+'_C'] = int(vals[1]) / motif['total']
            motif[str(i)+'_G'] = int(vals[2]) / motif['total']
            motif[str(i)+'_T'] = int(vals[3]) / motif['total']
        
    return motif
    

def open_files(folder_num):
    dir_ds = "../CS466_mini/motif_finding/data_set copy " + str(folder_num) + "/"
    dir_pred = "../CS466_mini/motif_finding/predicted copy " + str(folder_num) + "/"
    
    with open(dir_ds+"motiflength.txt", 'r') as f:
        nLines = int(f.readline().strip())
        
    motif = read_motif(dir_ds+"motif.txt", nLines)
    pred_motif = read_motif(dir_pred+"predictedmotif.txt", nLines)
    
    with open(dir_ds+"sites.txt", 'r') as f:
        sites = [int(i) for i in f.readlines()]
       
    with open(dir_pred+"predictedsites.txt", 'r') as f:
        pred_sites = [int(i) for i in f.readlines()]
        
    with open(dir_pred+"time.txt", 'r') as f:
        time = [float(i) for i in f.readlines()][0]
    return motif, pred_motif, sites, pred_sites, time
        
def open_files_better_filename(icpc, ml, sl, sc, vers):
    dir_ds  = "../CS466_mini/motif_finding/data_set_"
    dir_pred = "../CS466_mini/motif_finding/predicted_data_set_"
    dirpath = ''
    dirpath+="{0:.6f}".format(icpc)
    dirpath+="_"
    dirpath+=str(ml)
    dirpath+="_"
    dirpath+=str(sl)
    dirpath+="_"
    dirpath+=str(sc)
    dirpath+="_"
    dirpath+=str(vers)
    dir_ds += dirpath
    dir_pred += dirpath
    
    with open(dir_ds+"/motiflength.txt", 'r') as f:
        nLines = int(f.readline().strip())
        
    motif = read_motif(dir_ds+"/motif.txt", nLines)
    pred_motif = read_motif(dir_pred+"/predictedmotif.txt", nLines)
    
    with open(dir_ds+"/sites.txt", 'r') as f:
        sites = [int(i) for i in f.readlines()]
       
    with open(dir_pred+"/predictedsites.txt", 'r') as f:
        pred_sites = [int(i) for i in f.readlines()]
        
    with open(dir_ds+"/time.txt", 'r') as f:
        time = [float(i) for i in f.readlines()][0]
    return motif, pred_motif, sites, pred_sites, time
    
def rel_entropy(motif, pred_motif):
    """
    Calculate relative entropy of pred_motif relative to motif
    """
    
    to_ret = 0
    missing = 0
    for i in range(motif['nLines']):
        if motif[str(i)+'_A'] == 0:
            missing += pred_motif[str(i)+'_A']
        if (motif[str(i)+'_T'] == 0):
            missing += pred_motif[str(i)+'_T']
        if motif[str(i)+'_C'] == 0:
            missing += pred_motif[str(i)+'_C']
        if (motif[str(i)+'_G'] == 0):
            missing += pred_motif[str(i)+'_G']
        if motif[str(i)+'_A'] != 0 and pred_motif[str(i)+'_A'] != 0:
            to_ret += pred_motif[str(i)+'_A'] * \
                                 math.log(pred_motif[str(i)+'_A'] / motif[str(i)+'_A'])
        if motif[str(i)+'_T'] != 0 and pred_motif[str(i)+'_T'] != 0:
            to_ret += pred_motif[str(i)+'_T'] * \
                                 math.log(pred_motif[str(i)+'_T'] / motif[str(i)+'_T'])
        if motif[str(i)+'_C'] != 0 and pred_motif[str(i)+'_C'] != 0:
            to_ret += pred_motif[str(i)+'_C'] * \
                                 math.log(pred_motif[str(i)+'_C'] / motif[str(i)+'_C'])
        if motif[str(i)+'_G'] != 0 and pred_motif[str(i)+'_G'] != 0:
            to_ret += pred_motif[str(i)+'_G'] * \
                                 math.log(pred_motif[str(i)+'_G'] / motif[str(i)+'_G'])
    return abs(to_ret), missing
    
def site_overlap(sites, pred_sites, ml):
    """
    return amount of overlap
    """
    
    to_ret = 0
    
    for i in range(len(sites)):
        diff = (min(sites[i], pred_sites[i]) + ml) - max(sites[i], pred_sites[i]) 
    
        if diff > 0:
            to_ret += diff/ml
    return to_ret / len(sites)
    
class ParamSet:
    def __init__(self, icpc, ml, sl, sc):
        self.icpc = icpc
        self.ml = ml
        self.sl = sl
        self.sc = sc
        self.rel_ent = []
        self.rel_ent_missing = []
        self.sites_overlap = []
        self.time = []

    def add_rel_ent(self, rel_ent_tuple):
        """ add relative entropy value """
        rel_ent = rel_ent_tuple[0]
        rel_ent_miss = rel_ent_tuple[1]
        self.rel_ent.append(rel_ent)
        self.rel_ent_missing.append(rel_ent_miss)
    def min_rel_ent(self):
        """ calculate minimum relative entropy of all sets included """
        return min(self.rel_ent), min(self.rel_ent_missing)
    def max_rel_ent(self):
        """ calculate maximum relative entropy of all sets included """
        return max(self.rel_ent), max(self.rel_ent_missing)
    def avg_rel_ent(self):
        """ calculate average relative entropy of all sets included """
        return sum(self.rel_ent) / len(self.rel_ent), sum(self.rel_ent_missing) / len(self.rel_ent_missing)
    def std_rel_ent(self):
        """ calculate standard deviation of relative entropy of all sets included """
        return np.std(self.rel_ent), np.std(self.rel_ent_missing)
    
    def add_sites(self, sites, pred_sites):
        """ calculate and add the overlap amount """
        self.sites_overlap.append(site_overlap(sites, pred_sites, self.ml))
    def min_sites(self):
        """ calculate minimum site overlap of all sets included """
        return min(self.sites_overlap)
    def max_sites(self):
        """ calculate maximum site overlap of all sets included """
        return max(self.sites_overlap)
    def avg_sites(self):
        """ calculate average site overlaps of all sets included """
        return sum(self.sites_overlap) / len(self.sites_overlap)
    def std_sites(self):
        """ calculate standard deviation of site overlaps of all sets included """
        return np.std(self.sites_overlap)
        
    def add_time(self, time):
        """ add the time """
        self.time.append(time)
    def min_time(self):
        """ calculate minimum time of all sets included """
        return min(self.time)
    def max_time(self):
        """ calculate maximum time of all sets included """
        return max(self.time)
    def avg_time(self):
        """ calculate average time of all sets included """
        return sum(self.time) / len(self.time)
    def std_time(self):
        """ calculate standard deviation of time of all sets included """
        return np.std(self.time)
    
 
ps_2_8_500_10 = ParamSet(2, 8, 500, 10)
for folder_num in range(1,11):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_2_8_500_10.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_2_8_500_10.add_sites(sites, pred_sites)
    ps_2_8_500_10.add_time(time)
ps_1_8_500_10 = ParamSet(1, 8, 500, 10)
for folder_num in range(11,21):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_1_8_500_10.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_1_8_500_10.add_sites(sites, pred_sites)
    ps_1_8_500_10.add_time(time)
ps_15_8_500_10 = ParamSet(1.5, 8, 500, 10)
for folder_num in range(21,31):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_15_8_500_10.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_15_8_500_10.add_sites(sites, pred_sites)
    ps_15_8_500_10.add_time(time)
ps_2_7_500_10 = ParamSet(2, 7, 500, 10)
for folder_num in range(31,41):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_2_7_500_10.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_2_7_500_10.add_sites(sites, pred_sites)
    ps_2_7_500_10.add_time(time)
ps_2_6_500_10 = ParamSet(2, 6, 500, 10)
for folder_num in range(41,51):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_2_6_500_10.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_2_6_500_10.add_sites(sites, pred_sites)
    ps_2_6_500_10.add_time(time)
ps_2_8_500_20 = ParamSet(2, 8, 500, 20)
for folder_num in range(51,61):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_2_8_500_20.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_2_8_500_20.add_sites(sites, pred_sites)
    ps_2_8_500_20.add_time(time)
ps_2_8_500_5 = ParamSet(2, 8, 500, 5)
for folder_num in range(61,71):
    motif, pred_motif, sites, pred_sites, time = open_files(folder_num)
    ps_2_8_500_5.add_rel_ent(rel_entropy(motif, pred_motif))
    ps_2_8_500_5.add_sites(sites, pred_sites)
    ps_2_8_500_5.add_time(time)
        
sc_data = []
for i in range(6, 21, 1):
    sc_data.append(ParamSet(2, 8, 500, i))
    for vers in range(1, 11, 1):
        motif, pred_motif, sites, pred_sites, time = open_files_better_filename(2, 8, 500, i, vers)
        sc_data[-1].add_rel_ent(rel_entropy(motif, pred_motif))
        sc_data[-1].add_sites(sites, pred_sites)
        sc_data[-1].add_time(time)
    
def rel_icpc():
    
    x = [1, 1.5, 2]
    d_icpc = [ps_1_8_500_10, ps_15_8_500_10, ps_2_8_500_10]
    y1 = [i.avg_rel_ent()[0] for i in d_icpc]
    y2 = [i.avg_rel_ent()[1] for i in d_icpc]
    yerr_lower = [i.avg_rel_ent()[0] - i.min_rel_ent()[0] for i in d_icpc]
    yerr_upper = [i.max_rel_ent()[0] - i.avg_rel_ent()[0] for i in d_icpc]
    yerr_std = [i.std_rel_ent()[0] for i in d_icpc]
    yerr_std2 = [i.std_rel_ent()[1] for i in d_icpc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('ICPC')
    ax.set_ylabel('Relative Entropy')
    plt.errorbar(x, y1, yerr=yerr_std, ecolor='r')
    plt.errorbar(x, y2, yerr=yerr_std2, ecolor= 'y')
    plt.title("Relative Entropy changes with changes in ICPC")
    plt.show()
    
def rel_ml():
    
    x = [6, 7, 8]
    d_ml = [ps_2_6_500_10, ps_2_7_500_10, ps_2_8_500_10]
    y1 = [i.avg_rel_ent()[0] for i in d_ml]
    y2 = [i.avg_rel_ent()[1] for i in d_ml]
    yerr_lower = [i.avg_rel_ent()[0] - i.min_rel_ent([0]) for i in d_ml]
    yerr_upper = [i.max_rel_ent()[0] - i.avg_rel_ent()[0] for i in d_ml]
    yerr_std = [i.std_rel_ent()[0] for i in d_ml]
    yerr_std2 = [i.std_rel_ent()[1] for i in d_ml]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Motif Length')
    ax.set_ylabel('Relative Entropy')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.errorbar(x, y2, yerr=yerr_std2, ecolor= 'y')
    plt.title("Relative Entropy changes with changes in motif length")
    plt.show()
    
def rel_sc():
    
    x = [5, 10, 20]
    d_sc = [ps_2_8_500_5, ps_2_8_500_10, ps_2_8_500_20]
    y1 = [i.avg_rel_ent()[0] for i in d_sc]
    y2 = [i.avg_rel_ent()[1] for i in d_sc]
    yerr_lower = [i.avg_rel_ent()[0] - i.min_rel_ent()[0] for i in d_sc]
    yerr_upper = [i.max_rel_ent()[0] - i.avg_rel_ent()[0] for i in d_sc]
    yerr_std = [i.std_rel_ent()[0] for i in d_sc]
    yerr_std2 = [i.std_rel_ent()[1] for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Relative Entropy')
    plt.errorbar(x, y1, yerr=yerr_std, ecolor='r')
    plt.errorbar(x, y2, yerr=yerr_std2, ecolor= 'y')
    plt.title("Relative Entropy changes with changes in sequence count")
    plt.show()
    
def sites_icpc():
    
    x = [1, 1.5, 2]
    d_icpc = [ps_1_8_500_10, ps_15_8_500_10, ps_2_8_500_10]
    y = [i.avg_sites() for i in d_icpc]
    yerr_lower = [i.avg_sites() - i.min_sites() for i in d_icpc]
    yerr_upper = [i.max_sites() - i.avg_sites() for i in d_icpc]
    yerr_std = [i.std_sites() for i in d_icpc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('ICPC')
    ax.set_ylabel('Sites Overlap')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Sites overlap with changes in ICPC")
    plt.show()

def sites_ml():
    
    x = [6, 7, 8]
    d_ml = [ps_2_6_500_10, ps_2_7_500_10, ps_2_8_500_10]
    y = [i.avg_sites() for i in d_ml]
    yerr_lower = [i.avg_sites() - i.min_sites() for i in d_ml]
    yerr_upper = [i.max_sites() - i.avg_sites() for i in d_ml]
    yerr_std = [i.std_sites() for i in d_ml]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Motif Length')
    ax.set_ylabel('Sites Overlap')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Sites overlap with changes in motif length")
    plt.show()
    
def sites_sc():
    
    x = [5, 10, 20]
    d_sc = [ps_2_8_500_5, ps_2_8_500_10, ps_2_8_500_20]
    y = [i.avg_sites() for i in d_sc]
    yerr_lower = [i.avg_sites() - i.min_sites() for i in d_sc]
    yerr_upper = [i.max_sites() - i.avg_sites() for i in d_sc]
    yerr_std = [i.std_sites() for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Sites Overlap')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Sites overlap with changes in sequence count")
    plt.show()
    
def time_icpc():
    
    x = [1, 1.5, 2]
    d_icpc = [ps_1_8_500_10, ps_15_8_500_10, ps_2_8_500_10]
    y = [i.avg_time() for i in d_icpc]
    yerr_lower = [i.avg_time() - i.min_time() for i in d_icpc]
    yerr_upper = [i.max_time() - i.avg_time() for i in d_icpc]
    yerr_std = [i.std_time() for i in d_icpc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('ICPC')
    ax.set_ylabel('Time taken')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Time taken with changes in ICPC")
    plt.show()
    
def time_ml():
    
    x = [6, 7, 8]
    d_ml = [ps_2_6_500_10, ps_2_7_500_10, ps_2_8_500_10]
    y = [i.avg_time() for i in d_ml]
    yerr_lower = [i.avg_time() - i.min_time() for i in d_ml]
    yerr_upper = [i.max_time() - i.avg_time() for i in d_ml]
    yerr_std = [i.std_time() for i in d_ml]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Motif Length')
    ax.set_ylabel('Time taken')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Time taken with changes in motif length")
    plt.show()
    
def time_sc():
    
    x = [5, 10, 20]
    d_sc = [ps_2_8_500_5, ps_2_8_500_10, ps_2_8_500_20]
    y = [i.avg_time() for i in d_sc]
    yerr_lower = [i.avg_time() - i.min_time() for i in d_sc]
    yerr_upper = [i.max_time() - i.avg_time() for i in d_sc]
    yerr_std = [i.std_time() for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Time taken')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Time taken with changes in sequence count")
    plt.show()

def rel_sc_better_filename():
    
    x = [i for i in range(6, 21, 1)]
    d_sc = sc_data
    y1 = [i.avg_rel_ent()[0] for i in d_sc]
    y2 = [i.avg_rel_ent()[1] for i in d_sc]
    yerr_lower = [i.avg_rel_ent()[0] - i.min_rel_ent()[0] for i in d_sc]
    yerr_upper = [i.max_rel_ent()[0] - i.avg_rel_ent()[0] for i in d_sc]
    yerr_std = [i.std_rel_ent()[0] for i in d_sc]
    yerr_std2 = [i.std_rel_ent()[1] for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Relative Entropy')
    plt.errorbar(x, y1, yerr=yerr_std, ecolor='r')
    plt.errorbar(x, y2, yerr=yerr_std2, ecolor= 'y')
    plt.title("Relative Entropy changes with changes in sequence count")
    plt.show()
    
def sites_sc_better_filename():
    
    x = [i for i in range(6, 21, 1)]
    d_sc = sc_data
    y = [i.avg_sites() for i in d_sc]
    yerr_lower = [i.avg_sites() - i.min_sites() for i in d_sc]
    yerr_upper = [i.max_sites() - i.avg_sites() for i in d_sc]
    yerr_std = [i.std_sites() for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Sites Overlap')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Sites overlap with changes in sequence count")
    plt.show()

def time_sc_better_filename():
    
    x = [i for i in range(6, 21, 1)]
    d_sc = sc_data
    y = [i.avg_time() for i in d_sc]
    yerr_lower = [i.avg_time() - i.min_time() for i in d_sc]
    yerr_upper = [i.max_time() - i.avg_time() for i in d_sc]
    yerr_std = [i.std_time() for i in d_sc]
                  
    plt.figure()
    ax = plt.gca()
    ax.margins(0.05)
    ax.set_xlabel('Sequence Count')
    ax.set_ylabel('Time taken')
    plt.errorbar(x, y, yerr=yerr_std, ecolor='r')
    plt.title("Time taken with changes in sequence count")
    plt.show()
