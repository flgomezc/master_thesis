#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

TIMES = 1
DO_SKELETON = True

import numpy as np
import subprocess as subp
import timeit

toc = timeit.default_timer()

print("toc {}".format(toc))


#############################################################
#                                                           #
#                     Create Folders                        #                       
#                                                           #
#############################################################

fig_path = "./figures/"
oc_path  = "./observed_catalogs/"
rc_path  = "./random_catalogs/"
fc_path  = "./full_catalogs/"

subp.run("mkdir -p " + fig_path, shell=True, check=True)
subp.run("mkdir -p " + oc_path, shell=True, check=True)
subp.run("mkdir -p " + rc_path, shell=True, check=True)
subp.run("mkdir -p " + fc_path, shell=True, check=True)
    
oc_filename = "SDSS_data_Planck15.txt"
rc_filename = "SDSS_data_random_Planck15.txt"

subp.run( "echo # Count, Beta, n_rand > progress.txt", shell=True, check=True)


#############################################################
#                                                           #
#               Load and create catalogs                    #                       
#                                                           #
#############################################################
     
## Read the file
OC = np.loadtxt(oc_path + oc_filename)

## Number of observations.
N_obs = OC.shape[0]

## Beta and random_points_catalog_size parameters
beta, n_rand = 1, 1

## File Identifier:
## Will be useful to have the same name here, there and everywhere.
fileId = "SDSS_data_Plank15"

### LOAD Random Catalog
RC = np.loadtxt(rc_path + rc_filename)

### CREATE Full Catalog stacking RC and OC
FC = np.vstack([RC, OC])
# Store FullCatalog
FC_filename = fc_path + "FC_" + fileId + ".cat"
np.savetxt(FC_filename, FC)



#############################################################
#                                                           #
#     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
#                                                           #
#############################################################
maxiCounter = 0

if(DO_SKELETON):
    subp.run("LSS_BSK_calc -input  " + FC_filename + 
             " -output " + "BetaSkeleton_" + fileId + 
             " -beta " + str(beta) + 
             " -printinfo True -numNNB 300"
             , shell=True, check=True)

t = timeit.default_timer()
subp.run( "echo {}, {}, {}, {} >> progress.txt".format(maxiCounter, beta, n_rand, t), shell=True, check=True)

print( "Step {} of {} done".format(maxiCounter+1, TIMES) )



tic = timeit.default_timer()

print("tic {}".format(tic))

print("Time Elapsed: ", tic - toc)
