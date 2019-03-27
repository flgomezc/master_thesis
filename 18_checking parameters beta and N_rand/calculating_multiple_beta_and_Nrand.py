#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

import numpy as np
import matplotlib.pyplot as plt
import subprocess as subp


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

np.random.seed(27)
B, N = [], []

for i in range(300):
    beta, N_rand = generateParams()
    B.append(beta)
    N.append(N_rand)
    
# plt.scatter(B,N, s=1)

cut = 0
radius = 100

# Create figures folder
fig_path = "./figures/"
subp.run("mkdir -p " + fig_path, shell=True, check=True)

# Create Observed Catalogs folder
oc_path = "./observed_catalogs/"
subp.run("mkdir -p " + oc_path, shell=True, check=True)

# Create Random Catalog folder                                                                                                                                            
rc_path = "./random_catalogs/"
subp.run("mkdir -p " + rc_path, shell=True, check=True)

# Create Full Catalog folder                                                                                                                                              
fc_path = "./full_catalogs/"
subp.run("mkdir -p " + fc_path, shell=True, check=True)


filename = "sphere_data_cut_{}.dat".format(cut)


## Read the file
OC = np.loadtxt(oc_path + filename)

## Number of observations.
N_obs = OC.shape[0]

## Center of Data
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


# Checking outer shell.
fig, ax = plt.subplots()
index = np.where(r_list> 0.95*r_max )[0]
plt.scatter( OC[index,1], OC[index,0], s=0.1)
ax.set_aspect(1.0)
plt.tight_layout()

plt.savefig(fig_path+"outer_shell_cut_{}.pdf".format(cut))              

 
#############################################################
#                                                           #
#     Generate 2N Random Catalog (R = 0.95 r_max)           #
#                                                           #
#############################################################

### Using the cut as seed
np.random.seed(cut)

### A sigthly smaller radius to avoid surface percolation.
R_rc = radius * 0.95

### Create Random Catalog
RC = np.zeros([2*N_obs,3])
for x in RC:
    x += sph_random_point(R_rc)
np.savetxt(rc_path + "rnd_sph_cat_R_{}_cut_{}_Nrand_{}_Beta_{}_.cat".format(R_rc,cut,N_rand,beta), RC)

# Create the full catalog (to calculate Beta Skeleton Graph)
FC = np.vstack([RC,OC])
FC_filename = fc_path + "FC_CUT_{}_Nrand_{}_Beta_{}_.cat".format(cut, N_rand, beta)

# Store FullCatalog
np.savetxt(FC_filename, FC)



#############################################################
#                                                           #
#     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
#                                                           #
#############################################################

"""
subp.run("LSS_BSK_calc -input  " + FC_filename + 
         " -output " + "BetaSkeleton_CUT_{}_Nrand_{}_Beta_{}" + 
         " -beta " + str(beta) + 
         " -printinfo True -numNNB 300"
         , shell=True, check=True)

subp.run( "echo " + str(cut) + " >> progress.txt", shell=True, check=True)

cut += 1
"""