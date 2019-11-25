import numpy as np
import subprocess as subp
import os
import timeit
import argparse

#############################################################
#                                                           #
#                        DOCUMENTATION                      # 
#                                                           # 
############################################################# 


description = "This script uses a recursive method to identify the voids from the Beta-Skeleton File, and the Random, Observed and Full Catalog."
epilog = "At the end, the script stores a Masterlist (.mls) file."
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


arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
ALGORITHM   = arg.ALGORITHM
OC_FILE_IN = arg.ocfilein


#############################################################
#############################################################
##                                                         ## 
##                                                         ##
##                Begins the Main Routine                  ## 
##                                                         ## 
##                                                         ##
#############################################################
############################################################# 


prog = "progress_{}.txt".format(OC_FILE_IN)

### Paths
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


#################################################################
#                                                               #
#                        Load Catalogs                          #
#                              and                              #
#                    Corrects BetaSkeleton Format               #
#                                                               #
#################################################################


beta = BETA
n_rand = nrand
toc = timeit.default_timer()
subp.run( "echo  Finding voids in {} >> {}".format(OC_FILE_IN, prog), shell=True, check=True)

InitialMessage  = "\n\n\n #########################################################"
InitialMessage += "\n\n Running with:\n"
InitialMessage += "\n\t n_rand = {}".format(n_rand)
InitialMessage += "\n\t beta   = {}".format(beta)

print(InitialMessage)
print(" Using this files:\n\t{}\n\t{}\n\t{}\n\t{}".format(RC_filename,FC_filename,BS_filename,OC_filename,ML_filename))

RC    = np.loadtxt(RC_path + RC_filename)
OC    = np.loadtxt(OC_path + OC_filename)
FC    = np.loadtxt(FC_path + FC_filename)
BS    = np.loadtxt(BS_path + BS_filename)

N_rnd = RC.shape[0]
N_obs = OC.shape[0]

print(" Reading Full Catalogs and Beta-Skeleton from:\n\t{}\n\t{}".format(FC_filename,BS_filename))
print(" Previous BetaSkeleton Shape before Stacking: ", BS.shape)



### Transforms Xiao-Dong Li's Beta Skeleton Index to long list
try:
    if (ALGORITHM == "XDL"):
        a = BS[:,0].astype(int)
        a = list(a)
        b = BS[:,1].astype(int)
        b = list(b)

        c = []
        c.extend(a)
        c.extend(b)
        d = []
        d.extend(b)
        d.extend(a)

        c = np.array(c, dtype=int)
        d = np.array(d, dtype=int)

        fcBSkel = np.vstack((c,d)).T
        print("Next BetaSkeleton Shape after Stacking: ", fcBSkel.shape)
    elif (ALGORITHM == "NGL"):
        fcBSkel = BS
except ValueError:
    print("Invalid ALGORITHM, must be 'XDL' or 'NGL'")
        



#################################################################
#                                                               #
#              Identify TRUE VOID random POINTS                 #
#                        and                                    #
#                Their CONNECTIONS in Beta-Skeleton             #
#                                                               #
#################################################################

### Search for the first N_rnd points in the FC.

# Find RANDOM POINTS in the fcBeta-Skeleton Graph.
first_filter_index = np.where(fcBSkel[:,0] < N_rnd)  
    
# Store the partial Beta-Skeleton Graph of Random Points and
# its connections. They may have connections with Obs. points
# and other Random points.
first_filter_BSkel = np.array(fcBSkel[first_filter_index]).astype(int)

# Find the Random Points connected only to Random Points.

# To do this, first we find those points whom are connected to 
# observational points.
second_filter_index = np.where( first_filter_BSkel[:,1] >= N_rnd )[0]
    
# They are going to be dropped.
particle_ID_to_drop = first_filter_BSkel[second_filter_index,0]
particle_ID_to_drop.sort()
# A set of the Random Points connected to Observational points
# is created, there are not repeated items.
droplist = set(particle_ID_to_drop)

print( "First filter shape:", first_filter_BSkel.shape, 
       "\nHow many of them have direct connections"+
       " with galaxies (i.e. droplist length)", 
       len(droplist),
       "\nThen, must survive", len(set(first_filter_BSkel[:,0])) -len(droplist), 
       "trueVoidPoints")

# We have the Random points set:
# Maybe not all Random Particles are connected to the Skeleton. (large Beta)
# Because of this, we doesn't take into account something like
# Points_in_Skeleton = range(0,N_rnd).


print('Checking Random Points in the Beta Skeleton')
Points_in_Skeleton = set(first_filter_BSkel[:,0])
print('Random Points in Beta Skeleton Checked') 

# and the droplist. The complement(difference) is the
# pure void points set.

print('')
trueVoidPointsIndex = Points_in_Skeleton.difference(droplist)
# This set is converted to list, it will be used as an index to find 
# True Voids.
trueVoidPointsIndex = list(trueVoidPointsIndex)
trueVoidPointsIndex.sort()

# This is the first definition of TRUE VOID POINTS.
# Catalog of particles in voids
void_cat = FC[trueVoidPointsIndex]

### True Voids have been foud. #########################################
########################################################################


### Looking for the connections of the TrueVoidPoints
    
index=[]
for k in trueVoidPointsIndex:
    #index.extend( list( np.where( (fcBSkel[:,0] == k) &
    #                             (fcBSkel[:,1] > k))[0].astype(int) ) )
    
    index.extend( list( np.where( (fcBSkel[:,0] == k))[0].astype(int) ) )
    
    
    index = list(set(index ) )
    index.sort()
#DEBUG
#print(fcBSkel[index])
#print(fcBSkel[index].shape)

# Beta-Skeleton of TrueVoidPoints (includes Frontier Points)
VoidsBS = np.array(fcBSkel[index]).astype(int)
trueVoidPointsIndex.sort()
    
print(" Void BetaSkeleton Shape: ", VoidsBS.shape)
print(" The len of trueVoidPointsIndex", len(trueVoidPointsIndex))

### The list of All (random) Void Points
#
# Each TrueVoidPoint_index is checked, and their connections.
# They are random points connected to galaxies and other random
# points.

VoidPoints = []

for truepoint in trueVoidPointsIndex:
    aux = [truepoint]
    index = np.where(VoidsBS[:,0] == truepoint)[0]
    aux.extend(list(VoidsBS[index,1]))
    aux.sort()
    VoidPoints.extend(aux)

VoidPoints = list(set(VoidPoints))
VoidPoints.sort()
VoidPoints = np.array(VoidPoints).astype(int)
trueVoidPointsIndex = np.array(trueVoidPointsIndex).astype(int)

#################################################################
#                                                               #
#      RECURSIVE FRIEND OF FRIENDS VOID FINDER FUNCTION         #
#                                                               #
#################################################################

# In order to avoid infinite loops, an array is created to check
# if the True Random Point has been checked before. 
iscounted = np.zeros([len(trueVoidPointsIndex)], dtype=int)

def neighbours(ID):
    """
    Returns the list of neighbours of the ID point (included)
    IF ID is a TrueVoidPoint and has not been counted yet.
    """
    
    # IDin: Index of the void particle #ID in the list of VoidPoints
    IDin = np.where(trueVoidPointsIndex == ID)[0]
    # Empty list of of neighbours
    neigh = []
    
    # Base case: The Random Point is not a TRUE VOID POINT
    if IDin.size == 0:
        iscounted[IDin] += 1
        return neigh

    # Base case: The TrueVoidPOint has been counted. Skip
    if iscounted[IDin]:
        return neigh
    
    else:
        # Avoid infinite loops increasing the counter.
        # Append the True Random Point to the list.
        # Find RANDOM POINT neighbours in the Beta Skeleton.
        # They are stored in the "friends" list.
        iscounted[IDin] += 1               
        neigh = [ID]
        index = np.where(VoidsBS[:,0] == ID)[0]
        friends = list(VoidsBS[index,1])
        friends.sort()     

        # Now, the magic happens:
        for friend in friends:      
            neigh.append(friend)
            # In order to understand recursion,
            # first you need to understand recursion.
            neigh.extend(neighbours(friend))
        
        neigh = list(set(neigh))
        neigh.sort()
        return neigh

#                                                               #
#################################################################



#################################################################
#                                                               #
#                GENERATE THE Void Masterlist                   #
#                                                               #
#                 calling the recursive function                #
#                                                               #
#################################################################

VOIDS = []
VPlen = len(trueVoidPointsIndex)

for ID in trueVoidPointsIndex:
    index = np.where(trueVoidPointsIndex==ID)[0]
    if(ID // 100 == 0):
        print("Progress: {:04.2f}%".format(100.0*index[0]/VPlen))
    candidate = neighbours(ID)
    candidate = list(set(candidate))
    candidate.sort()
    if len(candidate) > 0:
        VOIDS.append(candidate)

#################################################################
#                                                               #
#                Store the Void Masterlist                      #
#                                                               #
#################################################################
MasterList = VOIDS

X = RC[:,0]
Y = RC[:,1]
Z = RC[:,2]

with open(ML_path + ML_filename, 'w') as file:

    for k in range(len(MasterList)):   
        for particle in MasterList[k]:
            line  = str( k) + " , " 
            line += str(X[particle]) + "," 
            line += str(Y[particle]) + ", " 
            line += str(Z[particle] ) + "\n"
            
            file.write(line)

    tic = timeit.default_timer()

FinalMessage = "\n\n We have fihished the process of "
FinalMessage += "\n Finding Voids in the file with "
FinalMessage += "\n beta: \t{}".format(beta)
FinalMessage += "\n n_rnd:\t{}".format(n_rand)
FinalMessage += "\n\n\tOutput written in \n" + ML_path + ML_filename
FinalMessage += "\n\n\t Time elapsed:{} seconds".format(tic - toc)
FinalMessage += "\n ###############################################"

print(FinalMessage)

subp.run( "echo MasterList {}, beta {}, n_rand {}, voids found in {} seconds. >> {}".format(ML_filename, beta, n_rand, tic-toc, prog), shell=True, check=True)
