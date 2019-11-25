#! /usr/bin/python
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

description = "This script does the volume and excentricity calculation of voids from the set of points in Masterlists."
epilog = "At the end, the script stores a VolumeAndExcentricity (.vae) file."
parser = argparse.ArgumentParser(description=description, epilog=epilog)

parser.add_argument('ocfilein', type=str)

parser.add_argument('-b', '--beta', type=float,
                    default=1.0,
                    help='Beta Skeleton Value, a float value "b>=1". Default Value = 1.0')
parser.add_argument('-n', '--nrand', type=float,
                    default=1.0,
                    help='The ratio between Number of Random Points and Number of Observational Points (nrand= N_random/N_obs)')

arg = parser.parse_args()

BETA  = arg.beta
nrand = arg.nrand
OC_FILE_IN = arg.ocfilein
prog = "progress.txt"

#############################################################
#############################################################
##                                                         ## 
##                                                         ##
##                Begins the Main Routine                  ## 
##                                                         ## 
##                                                         ##
#############################################################
############################################################# 

### Paths
OC_path = "./observed_catalogs/"
RC_path = "./random_catalogs/"
FC_path = "./full_catalogs/"
BS_path = "./xdl_beta_skeleton/"
ML_path = "./masterlists/"
FG_PATH = "./figures/"
VE_path = "./volume_and_excentricity/"

OC_filename = "{}.cat".format(OC_FILE_IN)
RC_filename = "{}.cat".format(OC_FILE_IN)
FC_filename = "{}.cat".format(OC_FILE_IN)
BS_filename = "{}.BSKIndex".format(OC_FILE_IN)
ML_filename = "{}.mls".format(OC_FILE_IN)
VE_filename = "{}.vae".format(OC_FILE_IN)

beta = BETA
n_rand = nrand

#################################################################
#                                                               #
#                    Void Finder Main Loop                      #
#                                                               #
#################################################################
subp.run( "echo Analysis of Ellipticity >> {}".format(prog), shell=True, check=True)

toc = timeit.default_timer()


print("Loading MasterList into 'Void Particle Cat'")
VoidParticleCat = np.loadtxt(ML_path + ML_filename, delimiter=',', dtype=[('ID', 'i8'),('x', 'f8'), ('y', 'f8'),('z', 'f8')])
print("Void Particle Cat")
print("VoidID, x, y, z")
print("Shape", VoidParticleCat.shape)
    
    
Void_N_particles_list = VoidParticleCat["ID"]
print("Len Void_N_particles_list", len(Void_N_particles_list)) 
    
M = Void_N_particles_list.max()

print('Maximum number of particles in a void: {}'.format(M))
   
if M > 0:
    print("\n\tNormal operation: Beta Skeleton file seems OK.             :) \n")
    
    # Void ID, Particles in this void.
    Void_N_Particles = []
    
    for i in range( M ):
        index = np.where(Void_N_particles_list == i)[0]
        Void_N_Particles.append([i, index.shape[0]])

    print("len Void_N_Particles", len(Void_N_Particles))

    Void_N_Particles = np.array(Void_N_Particles).astype("int")
    
    print("\tVoid number of Particles", Void_N_Particles.shape)


    #################################################################
    #                                                               #
    #                   Calcule and Store                           #
    #               Volume and Excentricity                         #
    #                                                               #
    #################################################################

    with open(VE_path + VE_filename, "w") as file:
        file.write("# Void_ID, X_c, Y_c, Z_c, N_particles, R=(abc)^(-1/3), a, b, c, vector_a, vector_b, vector_c \n\n")
        
        for ID in Void_N_Particles[:,0]:
            
            index = np.where( VoidParticleCat["ID"] == ID)[0]
            Void_n = VoidParticleCat[["x", "y", "z"]][index]
            particles_void_N = Void_n.shape[0]
            
            # Avoid 2 and three particle Voids
            if particles_void_N > 4:
                ###########################################
                #       Move to center of mass
                #           frame of reference.
                x = Void_n["x"]
                y = Void_n["y"]
                z = Void_n["z"]
                
                x_c = x.mean()
                y_c = y.mean()
                z_c = z.mean()
                
                x -= x.mean()
                y -= y.mean()
                z -= z.mean()
                #
                ###########################################
                

                ##########################################
                #        Calcule Inertia Tensor
                #              for Void Particles.
                I_11 = 0
                I_22 = 0
                I_33 = 0
                
                I_12 = 0
                I_13 = 0
                I_23 = 0

                for k in range( particles_void_N ):
                    I_11 += y[k]**2 + z[k]**2
                    I_22 += z[k]**2 + x[k]**2
                    I_33 += x[k]**2 + y[k]**2
                    
                    I_12 += - x[k] * y[k]
                    I_13 += - z[k] * x[k]
                    I_23 += - y[k] * z[k]

                I_21 = I_12
                I_31 = I_13
                I_32 = I_23

                I = np.array([
                    [I_11, I_12, I_13],
                    [I_21, I_22, I_23],
                    [I_31, I_32, I_33]])

                I = I/particles_void_N  # mass = \sum_1^N 1 / N = 1
                #
                ###################################################
                
                
                ###################################################
                #   Eigenvalues (w[i])
                #          and
                #   EigenVector  (v[i])
                #
                w, v = np.linalg.eig(I)

                w_max_index = np.where(w == w.max())[0]
                w_min_index = np.where(w == w.min())[0]
                w_med_index = np.where( (w != w.max()) & (w != w.min()))[0]
                
                
                # Sort semiaxes by length.
                # a - major semi-axis
                # b - med. semi-axis
                # c - minor semi-axis
                
                I_a = w[ w_min_index[0] ]
                I_b = w[ w_med_index[0] ]
                I_c = w[ w_max_index[0] ]
                
                a = ((5.0/2.0)*( -I_a + I_b + I_c))**0.5
                b = ((5.0/2.0)*( +I_a - I_b + I_c))**0.5
                c = ((5.0/2.0)*( +I_a + I_b - I_c))**0.5
                
                axis_a = v[:, w_min_index]
                axis_b = v[:, w_med_index]
                axis_c = v[:, w_max_index]
                #
                ##################################################

                vectors = [axis_a, axis_b, axis_c]
                
                str_vectors = ""
                
                for i in range(3):
                    for j in range(3):
                        str_vectors += str(vectors[i][j][0]) + ", "
                str_vectors = str_vectors[:-2]

                r = (a*b*c)**(1/3.0)
                
                stats = "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format( ID,   #1 
                                                                   x_c, y_c, z_c,    #2,3,4
                                                                   particles_void_N, #5
                                                                   r,
                                                                   a, b, c,          #7,8,9
                                                                   str_vectors)      
                file.write(stats)
else:
    print("\n\n\t%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("\n  \t% Beta Skeleton Failure  %")
    print("\n  \t%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("\n\t Wrong operation: Beta Skeleton file seems bad             :( \n")


tic = timeit.default_timer()

FinalMessage = "\n We have fihished the process of "
FinalMessage += "\n Finding Voids excentricity and Volume in the file with "
FinalMessage += "\n beta: \t{}".format(beta)
FinalMessage += "\n n_rand:\t{}".format(n_rand)
FinalMessage += "\n\n\tOutput written in \n" + VE_path + VE_filename
FinalMessage += "\n\n\t Time elapsed:{} seconds".format(tic - toc)
FinalMessage += "\n\n ########################"

print(FinalMessage)