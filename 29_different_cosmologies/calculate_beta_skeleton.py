#!/lustre/apps/anaconda/3.7/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3.7 module.

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import os
import subprocess as subp
import timeit
import datetime

import argparse


#############################################################
#                                                           #
#                        DOCUMENTATION                      #                       
#                                                           #
#############################################################

description = "This script generate a Random Catalog, the Full Catalog and calls the Xiao-Dong Li's Fortran Beta-Skeleton Calculator"
epilog = "At the end, the script stores a BetaSkeleton (.bsk) file."
parser = argparse.ArgumentParser(description=description, epilog=epilog)

parser.add_argument('ocfilein', type=str)

parser.add_argument('-b', '--beta', type=float,
                    default=1.0,
                    help='Beta Skeleton Value, a float value "b>=1". Default Value = 1.0')
parser.add_argument('-n', '--nrand', type=float,
                    default=1.0,
                    help='The ratio between Number of Random Points and Number of Observational Points (nrand= N_random/N_obs)')
parser.add_argument('-A', '--ALGORITHM', type=str,
                    default="XDL",
                    help='Algorithm used to calculate .BSKIndex, Options: "NGL", "XDL"')
parser.add_argument('-T', '--TEST', 
                    action='store_true',
                    default=False,
                    help='Tests filenames and folders generating empty files, does not runs the hard calculations.')

arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
ALGORITHM   = arg.ALGORITHM
OC_FILE_IN  = arg.ocfilein
BSK_FILE_IN = OC_FILE_IN
    
#############################################################
#############################################################
##                                                         ##
##                                                         ##
##                Begins the Main Routine                  ##                       
##                                                         ##
##                                                         ##
#############################################################
#############################################################
    
now = datetime.datetime.now()
toc = timeit.default_timer()

print("toc {}".format(toc))


#############################################################
#                                                           #
#                 Load Observed Catalogs                    #                       
#                                                           #
#############################################################

OC_path = "./observed_catalogs/"
RC_path = "./random_catalogs/"
FC_path = "./full_catalogs/"
BS_path = "./xdl_beta_skeleton/"
ML_path = "./masterlists/"
FG_PATH = "./figures/"

OC_filename = "{}.cat".format(OC_FILE_IN)
RC_filename = "{}.cat".format(OC_FILE_IN)
FC_filename = "{}.cat".format(OC_FILE_IN)
BS_filename = "{}.BSKIndex".format(OC_FILE_IN)
ML_filename = "{}.mls".format(OC_FILE_IN)
VE_filename = "{}.vae".format(OC_FILE_IN)

word_count =  subp.getoutput("wc -l " + OC_path + OC_filename).split()
try:
    OBS_CAT_SIZE = int( subp.getoutput("wc -l " + OC_path + OC_filename).split()[0] )
except:
    print("\n\t*** ERROR ***\n\n\t  --filein error: file '{}' not found. Are you sure that the filein file is placed inside the './observed_catalogs/' folder? \n".format(OC_filename))

prog = "progress.txt"

progress_string = "echo  Init time: {} >> {}".format(now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
    
#############################################################
#                                                           #
#                 Generate FULL Catalogs                    #                       
#                                                           #
#############################################################    

### LOAD Random Catalog
RC = np.loadtxt(RC_path + RC_filename)
OC = np.loadtxt(OC_path + RC_filename)

### CREATE Full Catalog stacking RC and OC
FC = np.vstack([RC, OC])
 
now = datetime.datetime.now()
progress_string = "echo  Generating Full Catalog: {} at time {} >> {}".format(RC_filename, now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
np.savetxt( FC_path + FC_filename, FC)


#############################################################
#                                                           #
#     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
#                                                           #
#############################################################

subp.run("LSS_BSK_calc -input  " + FC_path + FC_filename + 
         " -output " + str(BSK_FILE_IN) + 
         " -beta " + str(BETA) + 
         " -printinfo True -numNNB 300"
         , shell=True, check=True)

now = datetime.datetime.now()
progress_string = "echo  Generating {} Beta Skeleton: {} >> {}".format(BS_filename, now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
tic = timeit.default_timer()


print("Time Elapsed: ", tic - toc)

