import os

def create_dir(icpc, ml, sl, sc, i):
    dirpath = "motif_finding/data_set_"
    dirpath+="{0:.6f}".format(icpc)
    dirpath+="_"
    dirpath+=str(ml)
    dirpath+="_"
    dirpath+=str(sl)
    dirpath+="_"
    dirpath+=str(sc)
    dirpath+="_"
    dirpath+=str(i)
    return dirpath

def create_dir_tar(icpc, ml, sl, sc, i):
    dirpath = "motif_finding/predicted_data_set_"
    dirpath+="{0:.6f}".format(icpc)
    dirpath+="_"
    dirpath+=str(ml)
    dirpath+="_"
    dirpath+=str(sl)
    dirpath+="_"
    dirpath+=str(sc)
    dirpath+="_"
    dirpath+=str(i)
    return dirpath

#motif_finding/data_set_1.300000_8_500_10
for icpc in [x / 10.0 for x in range(1, 21, 1)]:
    for j in range(1, 11):
        dirpath = create_dir(icpc, 8, 500, 10, j)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        print dirpath

        dirpath_tar = create_dir_tar(icpc, 8, 500, 10, j)
        if not os.path.exists(dirpath_tar):
            os.makedirs(dirpath_tar)
        print dirpath_tar

for sc in range(5, 21):
    for j in range(1, 11):
        dirpath = create_dir(2, 8, 500, sc, j)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        print dirpath

        dirpath_tar = create_dir_tar(2, 8, 500, sc, j)
        if not os.path.exists(dirpath_tar):
            os.makedirs(dirpath_tar)
        print dirpath_tar
