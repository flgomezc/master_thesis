
# coding: utf-8

# In[50]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# In[2]:


BoxLength = 60
R = 20
filename = "A123456789B.cat"

catalog = np.loadtxt(filename)

X = catalog[:,0]
Y = catalog[:,1]
Z = catalog[:,2]


# In[3]:


N = len(X)


# In[4]:


rndcat = np.zeros([N,3])

for r in rndcat:
    a = np.random.random()
    b = np.random.random()
    c = np.random.random()
    
    r += BoxLength * np.array([a,b,c])


# In[5]:


rndcat


# In[6]:


indexC = np.where( abs(catalog[:,2] - BoxLength / 2) < R / 2 )


fig = plt.figure(figsize=[5,5])
plt.scatter(catalog[indexC,0], catalog[indexC,1], s=4)

indexR = np.where( abs(catalog[:,2] - BoxLength / 2) < R / 2 )

plt.scatter(rndcat[indexR,0], rndcat[indexR,1], s=4)


# In[7]:


fullcat = np.vstack((rndcat, catalog))


# In[8]:


fullcat.shape


# In[9]:


filename = "fullCat.cat"

with open( filename, "w") as file:
    for v in fullcat:
        file.write(str(v[0])+" "+str(v[1])+" "+str(v[2])+"\n")

print("Data saved to '" + filename + "'.") 


# # Now, go to terminal and run Beta Skeleton over the FullCatalog.

# ![title](a_few_moments_later.jpg)
# 
# 
# 

# In[10]:


filename = "fullCat_BETA_1.0.bsk"

# Full Catalog (cat + rndm) Beta Skeleton
fcBSkel = np.loadtxt(filename)


# In[11]:


fcBSkel.shape


# In[12]:


index = np.where(fcBSkel[:,0] < N)


# In[13]:


first_filter = np.array(fcBSkel[index])


# In[14]:


first_filter.shape


# In[15]:


index = np.where(first_filter[:,1] > N)


# In[16]:


index[0].shape


# In[17]:


droplist_raw = first_filter[index[0],0]


# In[18]:


droplist_raw.shape


# In[19]:


droplist_raw = droplist_raw.astype(int)


# In[20]:


set(droplist_raw)


# In[21]:


droplist = list(set(droplist_raw))


# In[22]:


len(droplist)


# The idea is to do something like this
# 
# index = np.where( first_filter[:,0] != droplist[])

# In[45]:


first_filter = first_filter.astype(int)


# In[46]:


void_points = []

for  x in list(set(first_filter[:,0])):
    if( x not in droplist):
        void_points.append(x)


# In[47]:


void_points


# In[51]:


sphere = rndcat[void_points]
print(sphere.shape)

x = sphere[:,0]
y = sphere[:,1]
z = sphere[:,2]

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.scatter(x, y, z, s=5)
ax.view_init(30, 30)

plt.title("Spherical Void Recognition\nUsing RandomPoints & 1-Skeleton")
plt.tight_layout()
#plt.savefig("IrregularVoidSurface_for_" + a_label + "-" + b_label + "-skeletons.pdf")
#plt.close()

plt.show()

