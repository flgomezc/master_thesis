#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

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
parser.add_argument('filein', type=str,
                    help='Name of the Observed Catalog file, must be stored in the folder OC_PATH="./observed_catalgos/"')
parser.add_argument('filenumber', type=int,
                     help='The consecutive number of the file, this is necessary to generate the Beta-Skeleton file.')
parser.add_argument('-b', '--beta', type=float,
                    default=1.0,
                    help='Beta Skeleton Value, a float value "b>=1". Default Value = 1.0')
parser.add_argument('-n', '--nrand', type=float,
                    default=1.0,
                    help='The ratio between Number of Random Points and Number of Observational Points (nrand= N_random/N_obs)')
parser.add_argument('-T', '--TEST', 
                    action='store_true',
                    default=False,
                    help='Tests filenames and folders generating empty files, does not runs the hard calculations.')

arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
OC_FILE_IN = arg.filein
FILENUM = arg.filenumber
TEST = arg.TEST

if(TEST):
    DO_SKELETON = False
else:
    DO_SKELETON = True


    
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

oc_path  = "./observed_catalogs/"
rc_path  = "./random_catalogs/"
fc_path  = "./full_catalogs/"
fig_path = "./figures/"

oc_filename = OC_FILE_IN
rc_filename = str(FILENUM) + ".cat"
fc_filename = str(FILENUM) + ".cat"

word_count =  subp.getoutput("wc -l " + oc_path + oc_filename).split()
try:
    OBS_CAT_SIZE = int( subp.getoutput("wc -l " + oc_path + oc_filename).split()[0] )
except:
    print("\n\t*** ERROR ***\n\n\t  --filein error: file '{}' not found. Are you sure that the filein file is placed inside the './observed_catalogs/' folder? \n".format(oc_filename))

prog = "progress.txt"

progress_string = "echo  Init time: {} >> {}".format(now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)


#############################################################
#                                                           #
#                      Plot Catalogs                        #                       
#                                                           #
############################################################# 
    
def PLOT_CATALOG(CAT, n, figname):

    fig = plt.figure(figsize=(10,10))

    ax1 = fig.add_subplot(221, adjustable='box', aspect=1)

    ax1.scatter(CAT[:,0], CAT[:,1], s= 0.01)
    ax1.set_xlabel("x [Mpc]")
    ax1.set_ylabel("y [Mpc]")


    ax2 = fig.add_subplot(223, adjustable='box', aspect=1)
    ax2.scatter(CAT[:,0], CAT[:,2], s= 0.01)
    ax2.set_xlabel("x [Mpc]")
    ax2.set_ylabel("z [Mpc]")
    ax2.set_ylim(-10,250)

    ax3 = fig.add_subplot(224, adjustable='box', aspect=1)
    ax3.scatter(CAT[:,1], CAT[:,2], s= 0.01)
    ax3.set_xlabel("y [Mpc]")
    ax3.set_ylabel("z [Mpc]")
    ax3.set_ylim(-10,250)

    ax4 = fig.add_subplot(222, adjustable='box', aspect=1)
    ax4.scatter(1,1, s=0.1)
    ax4.text(1, 1, figname + n, fontsize=20, 
         horizontalalignment='center',
         verticalalignment='center')
    
    plt.tight_layout()

    fig_name = "{}_{}.png".format(figname, n)
    plt.savefig(fig_path + fig_name)
    plt.close()

    
    
#############################################################
#                                                           #
#               Generate Random Catalogs                    #                       
#                                                           #
#############################################################
    
def getPoint2(r_0, r_1, t_0, t_1, p_0, p_1):

    #t's & p's must be in degrees

    
    t_0 = t_0 * np.pi / 180
    t_1 = t_1 * np.pi / 180
    
    u = np.random.rand()    
    theta = u * (t_1 - t_0 ) + t_0
    
    v_0 = np.cos((90-p_0)*np.pi/180)
    v_1 = np.cos((90-p_1)*np.pi/180)
    v = np.random.rand()
    phi = np.arccos( (v_0 - v_1)* v + v_1  )
    
    R = (1 - r_0/ r_1) * np.random.rand() + r_0/ r_1 
    r = r_1 * R **(1/3.)

    sinTheta = np.sin(theta)
    cosTheta = np.cos(theta)
    sinPhi = np.sin(phi)
    cosPhi = np.cos(phi)

    x = r * sinPhi * cosTheta
    y = r * sinPhi * sinTheta
    z = r * cosPhi

    return [x, y, z]



    
now = datetime.datetime.now()
progress_string = "echo  Generating {} Random Catalog: {} >> {}".format(FILENUM, now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
np.random.seed(FILENUM)
points = []
R = 300

N_rand = int(OBS_CAT_SIZE * nrand) # Number Points in Random Catalog
for i in range(N_rand):      
    points.append(getPoint2(0,R,0,90,45,0))
points = np.array(points)
np.savetxt( rc_path + rc_filename , points )
PLOT_CATALOG(points, FILENUM, "rc")
   
    
#############################################################
#                                                           #
#                 Generate FULL Catalogs                    #                       
#                                                           #
#############################################################    

### LOAD Random Catalog
RC = np.loadtxt(rc_path + rc_filename)
OC = np.loadtxt(oc_path + oc_filename)

### Create Figures
PLOT_CATALOG(OC, FILENUM, "oc")

### CREATE Full Catalog stacking RC and OC
FC = np.vstack([RC, OC])
 
now = datetime.datetime.now()
progress_string = "echo  Generating {} Full Catalog: {} >> {}".format(FILENUM, now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
np.savetxt( fc_path + fc_filename, FC)


#############################################################
#                                                           #
#     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
#                                                           #
#############################################################

if(DO_SKELETON):
    subp.run("LSS_BSK_calc -input  " + fc_path + fc_filename + 
             " -output " + str(FILENUM) + 
             " -beta " + str(BETA) + 
             " -printinfo True -numNNB 300"
             , shell=True, check=True)

now = datetime.datetime.now()
progress_string = "echo  Generating {} Beta Skeleton: {} >> {}".format(FILENUM, now.isoformat(), prog)
subp.run( progress_string, shell=True, check=True)
tic = timeit.default_timer()


print("Time Elapsed: ", tic - toc)
