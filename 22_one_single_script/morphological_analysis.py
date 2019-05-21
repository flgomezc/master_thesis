#! /usr/bin/python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import subprocess as subp
import os
import timeit

import scipy
from scipy.integrate import quad, dblquad, tplquad
from numpy import *
import scipy.stats

import argparse

#############################################################
#                                                           #
#                        DOCUMENTATION                      #      
#                                                           #
#############################################################

description = "This script uses the Masterlist (.mls) file to do the morphological analysis: volume, effective radius, ellipse semiaxes (a,b,c)."
epilog = "At the end, the script stores the info in the INFO (.inf) file and some graphs."
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


def Plot_Control(VAE_catalog, fig_name, n ):
    n = int(n)
    
    FG_format = "png"
    ID     = VAE_catalog[:,0]   # Void ID
    N      = ID.shape[0]        # Number of Voids in Catalog
    N_part = VAE_catalog[:,1]   # Number of Random Particles in Void
    V      = VAE_catalog[:,2]   # Void Volume 
    r      = VAE_catalog[:,3]   # r_eff
    a      = VAE_catalog[:,4]   # Major semi-axe
    b      = VAE_catalog[:,5]   # Medium semi-axe?
    c      = VAE_catalog[:,6]   # Minor semi-axe

    
    plt.scatter(ID,N_part)
    plt.xlabel("Void ID")
    plt.ylabel("Number of Particles")
    fg_filename = "n_particles_per_void_{}.{}".format(fig_name,FG_format) 
    plt.savefig(FG_path + fg_filename)
    plt.close()
    
    plt.scatter(ID, r)
    plt.xlabel("Void ID")
    plt.ylabel("Radius (Mpc)")
    fg_filename = "radius_per_void_{}.{}".format(fig_name,FG_format)
    plt.savefig(FG_path + fg_filename)
    plt.close()
    
    plt.scatter(N_part,r)
    plt.xlabel("Number of Particles")
    plt.ylabel("Radius (Mpc)")
    fg_filename = "radius_vs_number_of_particles_{}.{}".format(fig_name,FG_format)
    plt.savefig(FG_path + fg_filename)
    plt.close()


def PLOT_ANALYSIS(n, ID, N_part, r):
    plt.scatter(ID,N_part)
    plt.xlabel("Void ID")
    plt.ylabel("Number of Particles")
    plt.savefig(fig_path + "{}_paricles per Void.pdf".format(filein))
    plt.close()
    
    plt.scatter(ID, r)
    plt.xlabel("Void ID")
    plt.ylabel("Radius (Mpc)")
    plt.savefig(fig_path + "Planck15 Radius per Void.pdf")
    plt.close()
    
    plt.scatter(N_part,r)
    plt.xlabel("Number of Particles")
    plt.ylabel("Radius (Mpc)")
    plt.savefig(fig_path + "Planck15 Radius vs Particle Number.pdf")
    plt.close()    


def diff_volume(p,t,r):
    return r**2*sin(p)


def GET_VOLUME(r1, r2, t1, t2, p1, p2):
    """ 
    VOLUME # In Mpc**3
    Arguments: Radius (low, high), Theta, Phi
    """
    ## limits for radius
    # r1 = 0.
    # r2 = 300.
    ## limits for theta
    t1 = t1 * np.pi / 180  # 0
    t2 = t2 * np.pi / 180  # 2*pi
    # limits for phi
    p1 = (90 - p1 ) * np.pi / 180  # 0
    p2 = (90 - p2 ) * np.pi / 180                 # pi
    volume = tplquad(diff_volume, r1, r2, lambda r:   t1, lambda r:   t2,
                                      lambda r,t: p1, lambda r,t: p2)[0]
    return volume


## To do the prolate/oblate scatter plot.
def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z
    

RC_path = "random_catalogs/"
FC_path = "full_catalogs/"
OC_path = "observed_catalogs/"
BS_path = "xdl_beta_skeleton/"
ML_path = "masterlists/"
FG_path = "figures/"
VE_path = "volume_and_excentricity/"
AN_path = "analysis/"

OC_filename = "corner_data_cut_{}.dat".format(FILENUM)
RC_filename = "{}.cat".format(FILENUM)
FC_filename = "{}.cat".format(FILENUM)
BS_filename = "{}.bsk".format(FILENUM)
ML_filename = "{}.mls".format(FILENUM)
FG_filename = "{}".format(FILENUM)
VE_filename = "{}.vae".format(FILENUM)
AN_filename = "{}.info".format(FILENUM)

FG_format = "png"

OC = np.loadtxt( OC_path + OC_filename )
RC = np.loadtxt( RC_path + RC_filename )
FC = np.loadtxt( FC_path + FC_filename )
VE = np.genfromtxt( VE_path + VE_filename, delimiter=', ')[:,:-1]

volume = GET_VOLUME(0, 300, 0, 90, 45, 0)
# 0-300 Mpc, 0-90 deg, 0-45 deg

print("Volume = ", volume)

ID = VE[:,0]
N = ID.shape[0]
N_part = VE[:,1]
V = VE[:,2]
r = VE[:,3]
a = VE[:,4]
b = VE[:,5]
c = VE[:,6]


###########################################################
#                                                         #
#                  Excentricity Function                  #
#                                                         #
###########################################################

x = 1-c/a

fig = plt.figure(figsize=(4,4))
plt.hist(x, bins=20, density=True, histtype="step")
plt.ylabel(r"$f(N)dN$")
plt.xlabel(r"$\epsilon = 1 - c/a$")
plt.xlim(0,1)
plt.ylim(0,4)
plt.tight_layout()
plt.savefig(FG_path + "void_ellipticity_{}.{}".format(FILENUM, FG_format))
plt.close()



###########################################################
#                                                         #
#                   Prolate / Oblate                      #
#                                                         #
###########################################################

m1, m2 = b/a, c/b
xmin, xmax = 0, 1
ymin, ymax = 0, 1

X, Y, Z = density_estimation(m1, m2)

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111)            
# Show density 
ax.imshow(np.rot90(Z), cmap=plt.cm.terrain_r, extent=[xmin, xmax, ymin, ymax])
# Add contour lines
plt.contour(X, Y, Z, cmap="terrain")                                                                           
ax.plot(m1, m2, 'k.', markersize=2)    
ax.set_xlim([xmin, xmax])                                                                           
ax.set_ylim([ymin, ymax])                                                                           

unitary = np.linspace(0,1)
ax.plot(unitary, unitary, c="k")
ax.text(0.9, 0.1,  "Oblate", horizontalalignment="right", verticalalignment="center" )
ax.text(0.1, 0.9, "Prolate", horizontalalignment="left", verticalalignment="center" )

plt.xlabel("b / a")
plt.ylabel("c / b")
plt.tight_layout()
plt.savefig(FG_path + "void_two_axis_ratios_{}.{}".format(FILENUM, FG_format))
plt.close()


###########################################################
#                                                         #
#                Volume Density Function                  #
#                                                         #
###########################################################

# Observational data
ron_64  = np.loadtxt("data/ronconi_2019_z0_b64.dat" , delimiter="\t")
ron_128 = np.loadtxt("data/ronconi_2019_z0_b128.dat", delimiter="\t")
ron_256 = np.loadtxt("data/ronconi_2019_z0_b256.dat", delimiter="\t")
ron_500 = np.loadtxt("data/ronconi_2019_z0_b500.dat", delimiter="\t")
ade_PDF = np.loadtxt("data/adermann_2018_PDF.dat")

h = 0.6774
V_h = volume * (h**3)

Bins = np.linspace(0.1,1.2, 20)
deltaBins = Bins[1]-Bins[0]

y, Bins = np.histogram( np.log10(r*h), bins=Bins, density=False )
x = (Bins[:-1] + Bins[1:])/2

# Plot Volume Density Function
#
fig = plt.figure(figsize=[8,6])
fs = 20
plt.title("Using Uniform-Random Points", fontsize=fs)
# Adermann data
plt.plot(( ade_PDF[:,0] * ( 3 / (4 * pi)) ) ** (1 / 3.0), ade_PDF[:,1] , label="Adermann et al.")
# Ronconi Data
plt.plot( ron_64[:,0], ron_64[:,1],   label="Ronconi BZ 64")
plt.plot(ron_128[:,0], ron_128[:,1], label="Ronconi BZ 128")
plt.plot(ron_256[:,0], ron_256[:,1], label="Ronconi BZ 256")
plt.plot(ron_500[:,0], ron_500[:,1], label="Ronconi BZ 500")
# This Work!
plt.scatter(10**x, y/( deltaBins * V_h), label="This Work - SDSS-Planck15")

# Labels, legends, scale.
plt.xscale("log")
plt.yscale("log")
plt.ylim(0.00000001,0.1)
plt.xlim(0.4,30)
plt.xlabel(r"$r_{eff} \mathrm{ [Mpc/h]}$", fontsize=fs)
plt.ylabel(r"$ dn / d \ln r (\mathrm{ [h/Mpc]}^3)$", fontsize=fs)
plt.legend(loc=3, fontsize=int(fs*0.7))
plt.savefig(FG_path + "volume_density_function_{}.{}".format(FILENUM, FG_format))
plt.close()


## Store histogram data.
data = vstack([x,y / (deltaBins * V_h)]).T
np.savetxt(AN_path + "volume_pdf_" + AN_filename,data)


###########################################################
#                                                         #
#                  Galaxy/Halo Density                    #
#                                                         #
###########################################################

R_xyz = (OC[:,0]**2 + OC[:,1]**2 + OC[:,2]**2 )**0.5

NBINS = 20
Y, BINS = np.histogram(R_xyz, bins = NBINS)
X = []
BINVOLUME = []
for i in range(len(BINS)-1):
    X.append( (BINS[i+1] + BINS[i]) / 2)

X = np.array(X)
Y = np.array(Y)

BINVOLUME = []
for i in range(len(BINS)-1):
    r1 = BINS[i]
    r2 = BINS[i+1]
    # limits for theta

    ## limits for theta
    t1 = 0 * np.pi / 180  # 0
    t2 = 90 * np.pi / 180  # 2*pi
    # limits for phi
    p1 = (90 - 45 ) * np.pi / 180  # 0
    p2 = (90 - 0 ) * np.pi / 180                 # pi
    
    
    def diff_volume(p,t,r):
        return r**2*sin(p)

    BinVolume = np.abs(tplquad(diff_volume, r1, r2, lambda r:   t1, lambda r:   t2,
                                          lambda r,t: p1, lambda r,t: p2)[0])
    print(BinVolume)
    BINVOLUME.append(BinVolume)
    
BINVOLUME = np.array(BINVOLUME)

fig = plt.figure()
plt.scatter(X, Y / BINVOLUME)
plt.yscale("log")
plt.ylim(0.001,1)
plt.xlabel("R (Mpc)")
plt.ylabel("dN / dV")
plt.title("Random Points Density Number (per shell)")
plt.savefig(FG_path + "density_per_shell_{}.{}".format(FILENUM, FG_format))
