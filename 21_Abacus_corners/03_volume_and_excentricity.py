#! /usr/bin/python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import subprocess as subp
import os
import timeit

TESTING = False

### Paths
RC_path = "random_catalogs/"
FC_path = "full_catalogs/"
OC_path = "observed_catalogs/"
BS_path = "xdl_beta_skeleton/"
ML_path = "masterlists/"
FG_path = "figures/"
VE_path = "volume_and_excentricity/"

List_Random_Catalogs   = os.listdir(RC_path)
List_Full_Catalogs     = os.listdir(FC_path)
List_Observed_Catalogs = os.listdir(OC_path)
List_BetaSkeleton      = os.listdir(BS_path)
List_MasterLists       = os.listdir(ML_path)

List_Random_Catalogs.sort()
List_Full_Catalogs.sort()
List_Observed_Catalogs.sort()
List_BetaSkeleton.sort()
List_MasterLists.sort()

List_Figures = []
List_VolAndExc = []

Times = len(List_Random_Catalogs)

for n in range(Times):
    List_Figures.append("VoidVolumes_for_cut_{}.png".format(n))
    List_VolAndExc.append("VandE_cut_{}.vae".format(n))
List_Figures.sort()
List_VolAndExc.sort()

beta = np.ones(Times)
n_rand = np.ones(Times)

#################################################################
#                                                               #
#                    Void Finder Main Loop                      #
#                                                               #
#################################################################
subp.run( "echo # Count, Beta, n_rand > 03_Excentricity_progress.txt", shell=True, check=True)

toc = timeit.default_timer()

for maxiCounter in range(Times):


    InitialMessage  = "\n\n ########################"
    InitialMessage += "\n Running over cut {} with:\n".format(maxiCounter)
    InitialMessage += "\n\t n_rand = {}".format(n_rand[maxiCounter])
    InitialMessage += "\n\t beta   = {}".format(beta[maxiCounter])
    InitialMessage += "\n"

    print(InitialMessage)


    RC_filename = List_Random_Catalogs[maxiCounter]
    FC_filename = List_Full_Catalogs[maxiCounter]
    BS_filename = List_BetaSkeleton[maxiCounter]
    OC_filename = List_Observed_Catalogs[maxiCounter]
    ML_filename = List_MasterLists[maxiCounter]
    FG_filename = List_Figures[maxiCounter]
    VE_filename = List_VolAndExc[maxiCounter]
    
    print(ML_filename)


    RC = np.loadtxt(RC_path + RC_filename)
    N_rnd = RC.shape[0]

    OC = np.loadtxt(OC_path + OC_filename)
    N_obs = OC.shape[0]

    FC = np.loadtxt(FC_path + FC_filename)
    BS = np.loadtxt(BS_path + BS_filename)

    print("Previous BetaSkeleton Shape before Stacking: ", BS.shape)

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
    a = b = c = d = 0
    print("Next BetaSkeleton Shape after Stacking: ", fcBSkel.shape)

    print("Loading MasterList into 'Void Particle Cat'")
    VoidParticleCat = np.loadtxt(ML_path + ML_filename)
    print("Void Particle Cat")
    print("VoidID, x, y, z")
    print("Shape", VoidParticleCat.shape)


    Void_N_particles_list = VoidParticleCat[:,0].astype(int)
    print("Len Void_N_particles_list", len(Void_N_particles_list)) 

    # Void ID, Particles in this void.
    Void_N_Particles = []
    
    print("len Void_N_Particles", len(Void_N_Particles))

    M = Void_N_particles_list.max()
    
    if M > 0:
        print("\n\tNormal operation: Beta Skeleton file seems OK.             :) \n")

        for i in range( Void_N_particles_list.max()):
            index = np.where(Void_N_particles_list == i)[0]
            Void_N_Particles.append([i, index.shape[0]])

        Void_N_Particles = np.array(Void_N_Particles).astype("int")

        print("\tVoid number of Particles", Void_N_Particles.shape)

        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter( Void_N_Particles[:,0], Void_N_Particles[:,1])
        ax.set_yscale("log")
        ax.set_ylabel("Particles per Void")
        ax.set_xlabel("Void ID")
        plt.savefig=( FG_path + FG_filename)
        plt.close()


        #################################################################
        #                                                               #
        #                   Calcule and Store                           #
        #               Volume and Excentricity                         #
        #                                                               #
        #################################################################

        with open(VE_path + VE_filename, "w") as file:
            file.write("# Void_ID, Void_Volume, r = (a*b*c)**(1/3.0), a, b, c, vector_a, vector_b, vector_c \n\n")

            for ID in Void_N_Particles[:,0]:

                index = np.where( VoidParticleCat[:,0] == ID)[0]
                Void_n = VoidParticleCat[index, 1:4]
                particles_void_N = Void_n.shape[0]

                # Avoid 2 and three particle Voids
                if particles_void_N > 4:
                    ###########################################
                    #       Move to center of mass
                    #           frame of reference.
                    x = Void_n[:,0]
                    y = Void_n[:,1]
                    z = Void_n[:,2]

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

                    V = (4.0 / 3.0) * np.pi * a * b * c   
                    r = (a*b*c)**(1/3.0)

                    # debug
                    # print(r, V)
                    stats = "{}, {}, {}, {}, {}, {}, {}, {} \n".format(ID, particles_void_N, V, r, a, b, c, str_vectors)
                    file.write(stats)
    else:
        print("\n\n\t%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("\n  \t% Beta Skeleton Failure  %")
        print("\n  \t%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("\n\t Wrong operation: Beta Skeleton file seems bad             :( \n")


    tic = timeit.default_timer()

    FinalMessage = "\n We have fihished the process of "
    FinalMessage += "\n Finding Voids excentricity and Volume in the file with "
    FinalMessage += "\n Step {} of {}".format(maxiCounter+1, Times)
    FinalMessage += "\n beta: \t{}".format(beta[maxiCounter])
    FinalMessage += "\n n_rnd:\t{}".format(n_rand[maxiCounter])
    FinalMessage += "\n\n\tOutput written in \n" + VE_path + VE_filename
    FinalMessage += "\n\n\t Time elapsed:{} seconds".format(tic - toc)
    FinalMessage += "\n\n ########################"

    print(FinalMessage)
