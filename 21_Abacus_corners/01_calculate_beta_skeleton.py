#!/lustre/apps/anaconda/3/bin/python3
# Runing on the Magnus Cluster.
# Be sure to load the anaconda/python3 module.

BETA = 1
TEST = False
DO_SKELETON = True

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import os
import subprocess as subp
import timeit
import datetime

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

oc_filename = "SDSS_data_Planck15.txt"
rc_filename = "SDSS_data_random_Planck15.txt"

progress_string = "echo  Init time: {} > 01_beta_skeleton.progress".format(now.isoformat())
subp.run( progress_string, shell=True, check=True)




List_Observed_Catalogs = os.listdir(oc_path)
List_Observed_Catalogs.sort()

N_OBS_CATS = len(List_Observed_Catalogs)

OBS_CAT_SIZE = []


"""
for i in range( N_OBS_CATS):
    filename = oc_path + List_Observed_Catalogs[i]    
    OBS_CAT_SIZE.append(int( subp.getoutput("wc -l " + filename).split()[0]))
"""


"""    
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
    ax4.text(1, 1, "holi", fontsize=20, 
         horizontalalignment='center',
         verticalalignment='center')
    
    plt.tight_layout()

    fig_name = "{}_{}.png".format(figname, n)
    plt.savefig(fig_path + fig_name)
    plt.close()
""" 
    

    
"""     
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
""" 

    

if(TEST == True):
    N_OBS_CATS = 3
    
    
"""     
for n in range(N_OBS_CATS):
    now = datetime.datetime.now()
    progress_string = "echo  Generating {} Random Catalog: {} >> 01_beta_skeleton.progress".format(n, now.isoformat())
    subp.run( progress_string, shell=True, check=True)

    np.random.seed(n)
    points = []
    R = 300

    for i in range(OBS_CAT_SIZE[n]):      
        points.append(getPoint2(0,R,0,90,45,0))
    points = np.array(points)

    filename = "random_cut_{}.txt".format(n)   

    np.savetxt( rc_path + filename , points )
    
    PLOT_CATALOG(points, n, "Random_cut")
    
    
    
    
    
#############################################################
#                                                           #
#                 Generate FULL Catalogs                    #                       
#                                                           #
#############################################################    

for n in range(N_OBS_CATS):
    print(n)
    
    ### LOAD Random Catalog
    rc_filename = "random_cut_{}.txt".format(n)
    RC = np.loadtxt(rc_path + rc_filename)

    oc_filename = "corner_data_cut_{}.dat".format(n)
    OC = np.loadtxt(oc_path + oc_filename)

    ### Create Figures
    PLOT_CATALOG(OC, n, "Abacus_cut")
    
    ### CREATE Full Catalog stacking RC and OC
    FC = np.vstack([RC, OC])
 

    now = datetime.datetime.now()
    progress_string = "echo  Generating {} Full Catalog: {} >> 01_beta_skeleton.progress".format(n, now.isoformat())
    subp.run( progress_string, shell=True, check=True)

    fc_filename = "FC_{}.cat".format(n)
    np.savetxt( fc_path + fc_filename, FC)

""" 
    
List_Full_Catalogs = os.listdir(fc_path)
List_Full_Catalogs.sort()






#############################################################
#                                                           #
#     Call BetaSkeletonCalc from Xiao-Dong Li's Lib         #
#                                                           #
#############################################################

for n in range(N_OBS_CATS):

    fc_filename = List_Full_Catalogs[n]


    print(fc_filename)


    if(DO_SKELETON):
        subp.run("LSS_BSK_calc -input  " + fc_path + fc_filename + 
                 " -output " + "BetaSkeleton_" + str(n) + 
                 " -beta " + str(BETA) + 
                 " -printinfo True -numNNB 300"
                 , shell=True, check=True)

    now = datetime.datetime.now()
    progress_string = "echo  Generating {} Beta Skeleton: {} >> 01_beta_skeleton.progress".format(n, now.isoformat())
    subp.run( progress_string, shell=True, check=True)
    print( "Step {} of {} done".format(n + 1, N_OBS_CATS) )
    tic = timeit.default_timer()
    print("tic {}".format(tic))

print("Time Elapsed: ", tic - toc)
