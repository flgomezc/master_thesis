
# coding: utf-8

# # Finding voids by comparison of 0.90-1.00-Skeleton and 0.99-1.00-Skeleton

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def distance(p1,p2):
    r = 0
    for i in range(3):
        r += (p1[i]-p2[i])**2
    return ( r )**0.5


BoxLength = 60

filename = "EllipVoid_R_[20, 10, 10, 10]_BL_60_clean.cat"

catalog = np.loadtxt(filename)
Slice = 2

X = catalog[:,0]
Y = catalog[:,1]
Z = catalog[:,2]

print("Catalog Shape:", catalog.shape)


# In[2]:


index = np.where( abs(catalog[:,2] - BoxLength / 2) < Slice )


fig = plt.figure(figsize=[5,5])
plt.scatter(catalog[index,0], catalog[index,1], s=0.2)


# In[3]:


# Beta_list = [ "0.7", "0.8", "0.90", "0.99", "1.00", "1.3"]

# Beta_list = [  "1.00", "1.3"]

Beta_list = ["0.7"]

"""
for m in range(len(Beta_list)):
    beta_x = np.loadtxt("IRR_V_R20_BL60_B_" + Beta_list[m] + ".bsk")

    fig = plt.figure(figsize=(6,6))
    for k in beta_x:
        i = int(k[0])
        j = int(k[1])
        if ( (abs(Z[i] - BoxLength / 2) < Slice ) or (abs(Z[j] - BoxLength /2 ) < Slice ) ) :
            plt.plot( (X[i],X[j]), (Y[i],Y[j]) )
        # print(k)
    plt.xlim(0,BoxLength)
    plt.ylim(0,BoxLength)

    plt.title("Irregular Void (R = 20Mpc/h)\n" + r"$\beta =$ " + Beta_list[m])

    plt.xlabel("X (Mpc/h)")
    plt.ylabel("Y (Mpc/h)")


    plt.savefig( "beta_" + Beta_list[m] + ".pdf" )
    
    # plt.close()
"""


# # Where Connection Length (@ $\beta = 0.8$) is above Mean Length (@ $\beta = 1.3$)

# In[4]:


a_label = "0.8"
b_label = "1.3"



a = np.loadtxt("IRR_V_R20_BL60_B_" + a_label +".bsk")
b = np.loadtxt("IRR_V_R20_BL60_B_" + b_label +".bsk")

a = a.astype(int)
b = b.astype(int)

a_connections = []
b_connections = []

for i in range( catalog.shape[0] ):
    a_connections.append(np.where(a[:,0]==i)[0].shape[0])
    b_connections.append(np.where(b[:,0]==i)[0].shape[0])


#  # Connection Length

# In[5]:


a_dist = []
b_dist = []

for k in a:
    a_dist.append( distance( catalog[k[0]], catalog[k[1]]) )
    
for k in b:
    b_dist.append( distance( catalog[k[0]], catalog[k[1]]) )    
    
a_dist = np.array(a_dist)
b_dist = np.array(b_dist)

R_cut = max(b_dist)

bins = np.linspace(0,60,25)

fig = plt.figure(figsize=(7,5))
plt.hist( a_dist, bins, log=True, alpha=0.5, label=r"$\beta=$" + a_label)
plt.hist( b_dist, bins, log=True, alpha=0.9, label=r"$\beta=$" + b_label, histtype="step")
plt.xlabel("Mpc/h")
plt.legend()
plt.title("Connection Length")
plt.axvline(R_cut)


# In[6]:


index = np.where( (a_dist[:] > R_cut) )
sphere = catalog[a[index[0],0]]
print(sphere.shape)

x = sphere[:,0]
y = sphere[:,1]
z = sphere[:,2]

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.scatter(x, y, z, s=1)

plt.show()


# # Where Connection Length (@ $\beta = 0.90$) is above Mean Length (@ $\beta = 1.0$)

# In[7]:


a_label = "0.7"
b_label = "1.3"

R_cut = 13

a = np.loadtxt("IRR_V_R20_BL60_B_" + a_label +".bsk")
b = np.loadtxt("IRR_V_R20_BL60_B_" + b_label +".bsk")

a = a.astype(int)
b = b.astype(int)

a_connections = []
b_connections = []

for i in range( catalog.shape[0] ):
    a_connections.append(np.where(a[:,0]==i)[0].shape[0])
    b_connections.append(np.where(b[:,0]==i)[0].shape[0])


#  # Connection Length

# In[8]:


a_dist = []
b_dist = []

for k in a:
    a_dist.append( distance( catalog[k[0]], catalog[k[1]]) )
    
for k in b:
    b_dist.append( distance( catalog[k[0]], catalog[k[1]]) )    
    
a_dist = np.array(a_dist)
b_dist = np.array(b_dist)

R_cut = max(b_dist)

bins = np.linspace(0,60,25)

fig = plt.figure(figsize=(7,5))
plt.hist( a_dist, bins, log=True, alpha=0.5, label=r"$\beta=$" + a_label)
plt.hist( b_dist, bins, log=True, alpha=0.9, label=r"$\beta=$" + b_label, histtype="step")
plt.xlabel("Mpc/h")
plt.legend()
plt.title("Connection Length")
plt.axvline(R_cut)


# In[9]:


index = np.where( (a_dist[:] > R_cut) )
sphere = catalog[a[index[0],0]]
print(sphere.shape)

x = sphere[:,0]
y = sphere[:,1]
z = sphere[:,2]

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.scatter(x, y, z, s=0.3)

plt.show()


# In[10]:


index[0].shape

