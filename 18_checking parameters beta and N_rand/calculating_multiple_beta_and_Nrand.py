#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

TIMES = 500
DO_SKELETON = True

import numpy as np
import subprocess as subp
import timeit

toc = timeit.default_timer()

print("toc {}".format(toc))


def sph_random_point(Radius):
    u = np.random.rand()
    x1 = np.random.normal()
    x2 = np.random.normal()
    x3 = np.random.normal()
    
    norm = np.sqrt( x1**2 + x2**2 + x3**2)
    x1 /= norm
    x2 /= norm
    x3 /= norm
    
    r = Radius * u ** (1/3)
    
    return  [r*x1,r*x2,r*x3] 

def generateParams():
    beta = 1 + 4 * np.random.rand()
    N = 1  + 4 * np.random.rand()

    return beta, N


    
cut = 0
radius = 100

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
    
filename = "sphere_data_cut_{}.dat".format(cut)

subp.run( "echo # Count, Beta, n_rand > progress.txt", shell=True, check=True)
     
## Read the file
OC = np.loadtxt(oc_path + filename)

## Number of observations.
N_obs = OC.shape[0]

## Center of the dataset cut=0
center = np.array([-400,-400,-400])

## Move to Center of Data Coord. Syst. 
r_list = []
r_max  = 0
for i in range(N_obs):
    OC[i] = OC[i] - center

    r = 0
    for j in range(3):
        r += OC[i,j]**2
        r = r**0.5
        r_list.append(r)
        if r > r_max:
            r_max = r

r_list = np.array(r_list)



np.random.seed(27)
B, N = [], []


for maxiCounter in range(TIMES):
    beta, n_rand = generateParams()
    B.append(beta)
    N.append(n_rand)

    #############################################################
    #                                                           #
    #     Generate 2N Random Catalog (R = 0.95 r_max)           #
    #                                                           #
    #############################################################


    ### Using the maxiCounter as seed, because "cut" is constant.
    np.random.seed(maxiCounter)

    ### A sigthly smaller radius to avoid surface percolation.
    R_rc = radius * 0.95


    ### File Identifier:
    ### Will be useful to have the same name here, there and everywhere.
    fileId = "cat_R_{}_cut_{}_nrand_{}_Beta_{}_".format(R_rc,cut,n_rand,beta)

    ### Create Random Catalog

    # N_rand: Number of Random must be integer. 
    N_rand = round(n_rand * N_obs)

    RC = np.zeros([N_rand, 3])
    for x in RC:
        x += sph_random_point(R_rc)
    np.savetxt(rc_path + "rnd_sph_" + fileId + ".cat", RC)

    # Create the full catalog (to calculate Beta Skeleton Graph)
    FC = np.vstack([RC,OC])

    FC_filename = fc_path + "FC_" + fileId + ".cat"

    # Store FullCatalog
    np.savetxt(FC_filename, FC)


    #############################################################
    #                                                           #
    #     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
    #                                                           #
    #############################################################


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
