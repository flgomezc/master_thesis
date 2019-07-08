#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

import argparse
import subprocess as subp
import random


#############################################################
#                                                           #
#                        DOCUMENTATION                      #                       
#                                                           #
#############################################################

description = ("This script generates a range of numbers between N and M"
               " (not included). An 'i' index runs over the set, and generates"
               " a Random Catalog of data points with the same shape thab the"
               " Illustris corner of radius 300 Mpc. Then calls the Xiao-Dong"
               " Li's Fortran Beta-Skeleton Calculator and performs the whole"
               " analysis.\n\n"
               "Defaul FILEIN is 'corner_data_cut_0.dat'\n\n"
               "The idea is to run in paralell using different N,M pairs on "
               "each processor.")
epilog = ("At the end we have the analysis of each random pair of Beta and "
          "Nrandom parameters.")
parser = argparse.ArgumentParser(description=description, epilog=epilog)
parser.add_argument('N', type=int,
                    help='N, a filenumber to start the range(N,M) to run in parallel.')
parser.add_argument('M', type=int,
                    help='M, a filenumber to finish the range(N,M) to run in parallel.')

#############################################################
#                                                           #
#                Initialize, load variables                 #                       
#                                                           #
#############################################################


arg = parser.parse_args()
N = arg.N
M = arg.M


FILENAME = "corner_data_cut_0.dat"

py00_BETA_SKEL = "calculate_beta_skeleton.py"
py01_VOID_ID_R = "void_identifier-recursion.py"
py02_ELLIPS_AN = "ellipsoid_analysis.py"
py03_MORPHO_AN = "morphological_analysis.py"

random.seed(N)

def rand_par(l=1.0, h=3.0):
    rand_l = l
    rand_h = h
    rand_del = rand_h - rand_l
    r = random.random() * rand_del + rand_l
    return r

N_RAND = []
BETA = []

for i in range(N,M):
    N_RAND.append(rand_par(0.7, 3.0)) # rand_range(0.7-3.0)
    BETA.append(rand_par(1.0, 3.0))   # rand_range(1.0-3.0)


#############################################################
#                                                           #
#            Do the main sequence:                          #
#   RandomCat, BetaSkel, FindVoids, Analyze Voids.          #
#                                                           #
#############################################################

    
with open("N_RAND_and_BETA_list.txt", "w") as f:
    f.write("n_rand\tbeta\n")
    for i in range(N,M):
        j = i - N
        n_rand = N_RAND[j]
        beta   = BETA[j]
        FILENUM = i
        f.write("{}\t{}\t{}\n".format(i,n_rand, beta))

for i in range(N,M):
    j = i - N
    n_rand = N_RAND[j]
    beta   = BETA[j]
    FILENUM = i
    
    message = "\n\n########################################################\n"
    message += 'n_rand = {}\n'.format(n_rand)
    message += 'beta = {}\n'.format(beta)
    message+= '\n>>> Run {}\n\n'.format(py00_BETA_SKEL)
    print(message)
    subp.run("python {} {} {} --beta {} --nrand {}".format(py00_BETA_SKEL, FILENAME, FILENUM, beta, n_rand), shell=True)

    message = "\n\n-------------------------------------------------------\n"
    message += 'n_rand = {}\n'.format(n_rand)
    message += 'beta = {}\n'.format(beta)
    message+= '\n>>> Run {}\n\n'.format(py01_VOID_ID_R)
    print(message)
    subp.run("python {} {} {} --beta {} --nrand {}".format(py01_VOID_ID_R, FILENAME, FILENUM, beta, n_rand), shell=True)
    
    
    message = "\n\n-------------------------------------------------------\n"
    message += 'n_rand = {}\n'.format(n_rand)
    message += 'beta = {}\n'.format(beta)
    message+= '\n>>> Run {}\n\n'.format(py02_ELLIPS_AN)
    print(message)
    subp.run("python {} {} {} --beta {} --nrand {}".format(py02_ELLIPS_AN, FILENAME, FILENUM, beta, n_rand), shell=True)

    
    message = "\n\n-------------------------------------------------------\n"
    message += 'n_rand = {}\n'.format(n_rand)
    message += 'beta = {}\n'.format(beta)
    message+= '\n>>> Run {}\n\n'.format(py03_MORPHO_AN)
    print(message)
    subp.run("python {} {} {} --beta {} --nrand {}".format(py03_MORPHO_AN, FILENAME, FILENUM, beta, n_rand), shell=True)

