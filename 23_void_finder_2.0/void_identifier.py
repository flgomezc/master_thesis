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

description = "This script identifies the voids from the Beta-Skeleton File, and the Random,Observed and Full Catalog."
epilog = "At the end, the script stores a Masterlist (.mls) file."
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

arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
OC_FILE_IN = arg.filein
FILENUM = arg.filenumber



#############################################################
#############################################################
##                                                         ## 
##                                                         ##
##                Begins the Main Routine                  ## 
##                                                         ## 
##                                                         ##
#############################################################
############################################################# 


prog = "progress.txt"

### Paths
RC_path = "./random_catalogs/"
FC_path = "./full_catalogs/"
OC_path = "./observed_catalogs/"
BS_path = "./xdl_beta_skeleton/"
ML_path = "./masterlists/"

OC_filename = OC_FILE_IN
FC_filename = "{}.cat".format(FILENUM)
RC_filename = "{}.cat".format(FILENUM)
BS_filename = "{}.BSKIndex".format(FILENUM)
ML_filename = "{}.mls".format(FILENUM)


#################################################################
#                                                               #
#                    Void Finder Main Loop                      #
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
       "\nThen, must survive", first_filter_BSkel.shape[0]-len(droplist), 
       "trueVoidPoints")

# We have the Random points set:
# Maybe not all Random Particles are connected to the Skeleton. (large Beta)
# Because of this, we doesn't take into account something like
# Points_in_Skeleton = range(0,N_rnd).
Points_in_Skeleton = set(first_filter_BSkel[:,0])
# and the droplist. The complement(difference) is the
# pure void points set.
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
    index.extend( list( np.where( fcBSkel[:,0] == k)[0].astype(int) ) )
    
    index = list(set(index ) )
    index.sort()

# Beta-Skeleton of TrueVoidPoints (includes Frontier Points)
VoidsBS = np.array(fcBSkel[index]).astype(int)
trueVoidPointsIndex.sort()
    
print(" Void BetaSkeleton Shape: ", VoidsBS.shape)
print(" The len of trueVoidPointsIndex", len(trueVoidPointsIndex))

# This is the MasterList of Voids.

# Each TrueVoidPoint_index is checked.
# If doesn't exists, is identified as a new void (a new sublist is created)
# If it exists already, the point and its connections are added to the existing void sublist.
    
print(" Initialize MasterList.")


MasterList = []

for search in trueVoidPointsIndex:
    # Does the TrueVoidPoint belongs to any existing Void?
    is_in_master = any( search in sublist for sublist in MasterList)
    
    # Create a new void.
    # If is the first time it appears on the MasterList
    if not is_in_master:
        my_list = []

        # Find B-Skeleton connections of the TrueVoidPoint.
        index = np.where(VoidsBS[:,0] == search)
        # Append them to the auxiliar list.
        my_list.append(search)
        my_list.extend( list(VoidsBS[index,1][0]) )
        my_list.sort()
        # Append this auxiliar list to the MasterList.
        # a new void has been appended. :)
        MasterList.append(my_list)
        
    # If the TrueVoidPoint already exists in the Masterlist
    if is_in_master:
        repetitions = []
        
        # Find how many times the TrueVoidPoint has appeared before
        # in the MasterList. (Store the search)
        for k in range(len(MasterList)):
            if(search in MasterList[k]):
                sublist = MasterList[k]
                #print("Si est'a en la sublista", k, sublist)
                repetitions.append(k)
                
                #print(search, "appears in sublists:" , repetitions)

        # If it appears only in one time in the sublitst, append the 
        # TrueVoidPont and its B-Skeleton connections to the existing void.
        if ( len(repetitions) == 1 ):
            j = repetitions[0]
            index = np.where(VoidsBS[:,0] == search)
            my_list = list(VoidsBS[index,1][0])
            my_list.sort()
            MasterList[j].extend(my_list)


            # If the TrueVoidPoint appears more than one time, the lists
            # will merge into a new one.
            # Old lists will be empty. They will be removed after this
            # cicle ends.
        elif (len(repetitions) > 1):            
        
            # print("Friend of many friends, n=", len(repetitions))
            my_list = []
            for j in repetitions:
                my_list.extend(MasterList[j]) # collect data create a new list
                MasterList[j]=[]              # Empty the merged lists.
                my_list.sort()
                MasterList.append(my_list)

# Emtpy lists removed. 
while( [] in MasterList):
    MasterList.remove([])


# Sort each void list at the MasterList
for j in range(len(MasterList)):
    MasterList[j] = list( set(MasterList[j]))
    MasterList[j].sort()

print("\nRaw Masterlist created\n")           




# Some lists may share elements (FrontierPoints).
# Those FrontierPoints may be in a BottleNeck
# 
#                       
#     T       F    Halo     F     T   
#                                     T
#       T       F <----> F      T
#   T                                 T
#        T     F    Halo    F   T
#
#   Void_1                       Void_2
#
#

# This is a list of voids and the common particles.
to_merge = []
for i in range(len(MasterList)):              # For each VoidList
    for j in range(len(MasterList)):          # compare against other VoidLists
        if (j > i):
            # Check if two or more lists have common elements.
            aux = [x for x in MasterList[i] if x in MasterList[j]]
            # Store the common elements in the aux list.
            
            if (len(aux)>0):              # If the list is not empty
                # Void[i], Void[j], common particles.
                # print( i, j, aux)
                to_merge.append([i,j])   # Store the two void indices.
                
                ## NOTE: Many voids can be concatenated.


# This is a recount to concatenate those voids that share
# particles. Runs over the previous "to merge" list.

to_merge2 = []
for i in range( len(to_merge)):
    x = to_merge[i][0]
    y = to_merge[i][1]
    
    is_in_list1 = any( x in sublist for sublist in to_merge2 )
    is_in_list2 = any( y in sublist for sublist in to_merge2 )

    if( (is_in_list1 == False) & (is_in_list2 == False) ):
        to_merge2.append(to_merge[i])
        
    elif( (is_in_list1 == False) & (is_in_list2 == True) ):
        aux = []
        for j in range(len(to_merge2)):
            if( y in to_merge2[j]):
                to_merge2[j].append(x)

    elif( (is_in_list1 == True) & (is_in_list2 == False) ):
        aux = []
        for j in range(len(to_merge2)):
            if( x in to_merge2[j]):
                to_merge2[j].append(y)


for x in to_merge2:
    x.sort()

print( " ---> Control: Length of the list 'to_merge2': " + str( len(to_merge2) ) )



for sublist in to_merge2:
    aux = []
    for x in sublist:
        aux.extend(MasterList[x])
        MasterList[x] = []
        aux.sort()
        MasterList.append(aux)
        
while( [] in MasterList):
    MasterList.remove([])

    
print(" Total number of Void Particles\n\t", len(trueVoidPointsIndex))

aux = 0
for Void in MasterList:
    #print( len(Void))
    aux += len(Void)
    
    if len(Void)<2:
        print( Void)

print("Total number of particles in Voids and close to filaments:", aux , "(Void + Frontier particles)")





#################################################################
#                                                               #
#                Store the Void Masterlist                      #
#                                                               #
#################################################################

X = RC[:,0]
Y = RC[:,1]
Z = RC[:,2]

with open(ML_path + ML_filename, 'w') as file:

    for k in range(len(MasterList)):   
        for particle in MasterList[k]:
            line  = str( k *1.0) + " " 
            line += str(X[particle]) + " " 
            line += str(Y[particle]) + " " 
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

subp.run( "echo Filenumber {}, beta {}, n_rand {}, voids found in {} seconds. >> {}".format(FILENUM, beta, n_rand, tic-toc, prog), shell=True, check=True)
