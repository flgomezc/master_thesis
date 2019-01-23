
# coding: utf-8

# # Generator of a Random Field Catalog with one Spherical DIRTY Void
# 
# Spatial distribution of points: uniform
# 
# Box Length: 100 Mpc/h
# 
# Spherical Void radius: 30 Mpc/h
# 
# Void position: center of the Box
# 
# 
# ### Abacus Cosmos 720 box Plank
# 
# Number of Halos: 8'714.934
# 
# Volume = ( 720 Mpc/h )**3
# 
# Halo density: 0.02334891010802469 halo/ (Mpc/h)**3

# In[1]:


print( 8714934 / (720.0)**3 )


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import random as rand

# Halos per (Mpc/h)**3
D_halo = 0.02334891010802469

# Cubic box size (Mpc/h)
BoxLength = 100 

# Void radius (Mpc/h)
R = 30

# Number of particles
N = int( D_halo * (BoxLength)**3 ) 

print("Number of Particles:",N)
# Initialize random seed
rand.seed(1)


# In[3]:


def IsInsideVoid(v):
    # Position vector components
    x = v[0]
    y = v[1]
    z = v[2]
    
    # Void sphere center
    x0 = 0.50 * BoxLength
    y0 = 0.50 * BoxLength
    z0 = 0.50 * BoxLength

    # Distance to the void center
    d = ( (x0 - x)**2 + (y0 - y)**2 + (z0 - z)**2)**0.5 
    
    if ( d < R ):
        return 1;
    else:
        return 0;
    

def rand_pos():
    x = rand.random() * BoxLength
    y = rand.random() * BoxLength
    z = rand.random() * BoxLength   
    return ([x,y,z])


# In[4]:


positions = []

for i in range(N): 
    r = rand_pos();
    
    if ( IsInsideVoid(r) ):
        r = rand_pos()
    
    positions.append(r)  

positions = np.array(positions)

positions.shape


# In[5]:


np.savetxt("../data/Testing_Data/Spherical_Void_100Mpc_Box_Dirty.csv", positions)


# In[6]:


# Plot a slice centered at Z_c, slice tickness s_t

Z_c = 50

s_t = 10

positions_slice = positions[ np.where( abs(positions[:,2]- Z_c )< s_t )]

X = positions_slice[:,0]
Y = positions_slice[:,1]

fig = plt.figure(figsize=[5,5])
plt.scatter(X,Y, s=0.4)
plt.show()

