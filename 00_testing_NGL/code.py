import numpy as np
import matplotlib.pyplot as plt


from mpl_toolkits.mplot3d import Axes3D


# In[2]:


filename = "test_in"

data = np.loadtxt(filename)

data


# In[3]:


x = data[:,0]
y = data[:,1]
z = data[:,2]


# In[4]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z)


# In[5]:


filename = "test_out"
beta_skeleton = np.loadtxt(filename)

beta_skeleton.shape


# In[6]:


data.shape


# In[11]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

BS_size = beta_skeleton.shape[0]

ax.scatter(x, y, z)
ax.set_aspect('equal')

for n in range(BS_size):
    i, j = beta_skeleton[n]
    
    i = int(i)
    j = int(j)
    
    ax.plot([x[i], x[j]], [y[i],y[j]], zs=[z[i],z[j]])
    
plt.show()
