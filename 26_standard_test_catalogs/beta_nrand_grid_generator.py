#!/lustre/apps/anaconda/3.7/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3.7 module.

import numpy as np
import subprocess as subp

beta  = np.linspace(1,2,6)
nrand = np.linspace(1,2,6)

counter = 0
for b in beta:
    for n in nrand:
        #print(b, n, counter)      
        process_string = "make main_sequence FILENAME=v300.dat FILENUM={} BETA={} NRAND={}".format(counter, b, n)
        subp.run(process_string , shell=True, check=True)
        counter +=1

