#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import scipy.signal as signal
from scipy.optimize import curve_fit


# In[2]:


df = pd.read_csv("dados.csv", sep = ";", header=0)
df.head()
df = pd.read_csv("dados.csv", sep = ";", header=0, usecols=["Distancia km","Velocidade km/s"])
df.head()


# In[3]:


dist = df["Distancia km"].to_numpy()
vel = df["Velocidade km/s"].to_numpy()
print('Distancia = ',dist,'\n\n Velocidade',vel)


# In[4]:


##### Definir os minimos quadrados #######
ln_d = []
ln_v = []

for i in range(len(dist)):
    vel[i] = (10**3)*vel[i]
    dist[i] = (10**3)*dist[i]
print('Distancia = ',dist,'\n\n Velocidade',vel)

for i in range(len(dist)):
    ln_d.append(np.log(dist[i]))
    ln_v.append(2*np.log(vel[i]))
    
print('\n\n Log(d)',ln_d,'\n\n Log(v)',ln_v)

#### Plotar logs ###############################
fig, ax = plt.subplots(figsize=(10,10))
ax.scatter(ln_d, ln_v, s=25, color='red', label='Dados')
ax.set_title("Dispersão $ln(v)\ x\ ln(r)$")
ax.set_xlabel("ln(r)")
ax.set_ylabel("ln(v) ")
ax.grid(color='black')
plt.legend(loc='best')
plt.savefig('ln.png', format='png')
plt.show()


# In[5]:


###### Definição dos parametros para  a regressão mininmos quadrados ######
G = 6.6743e-11
x_quad =0
xy = 0
x = 0
y = 0
for i in range(len(vel)):
    x_quad += ln_d[i]**2
    xy+= ln_d[i]*ln_v[i]
    x +=ln_d[i]
    y+=ln_v[i]

print(x_quad,'\n\n',xy,'\n\n',x,'\n\n',y)
####Regressão primeira#####

b1 = (x_quad*y - xy*x)/(len(vel)*x_quad - x**2)
a1 = (y - len(vel)*b1)/x

print('\n\n A:',a1,'\n\n B:',b1)
    
#### Plotar logs ###############################
def func1(x1):
    return a1*x1 + b1

xx = np.linspace(24,30,20)

fig, ax = plt.subplots(figsize=(10,10))
ax.scatter(ln_d, ln_v, s=25, color='red', label='Dados')
ax.plot(xx,func1(xx),linestyle = '--', linewidth = 1.2, color='black', label='Ajuste' )
ax.set_title("Interpolação $ln(v)\ x\ ln(r)$")
ax.set_xlabel("ln(r)")
ax.set_ylabel("ln(v) ")
ax.grid(color = 'black')
plt.legend(loc='best')

plt.show()

sq_r=0
sq_t=0
for i in range(len(vel)):
    sq_r += (func1(ln_d[i]) - np.mean(ln_v))**2
    sq_t += (func1(ln_d[i]) - np.mean(ln_v))**2
    #print(a1*ln_d[i]+ln_d[i])
    print('\n\n M:', (np.exp(b1)*x**((a1+2)/2))/G)

#R = np.sqrt(1)


# In[6]:


x_quad =0
xy = 0
x = 0
y = 0
for i in range(len(vel)):
    x_quad += ln_d[i]**2
    xy+= ln_d[i]*ln_v[i]
    x +=ln_d[i]
    y+=ln_v[i]
    
#### Regressão segunda considerando a = -1 ###
G = 6.6743e-11
b = (y+x)/len(vel)
print('B:',b)

def func(x1):
    return -x1 + b

###### Incerteza do b ######

b_i = []

for i in range(len(vel)):
    b_i.append(np.abs((2*ln_v[i] + ln_d[i]) - b))

print("\n\n $\Delta$b:", b_i)
xx = np.linspace(24,30,20)

fig, ax = plt.subplots(figsize=(10,10))
ax.scatter(ln_d, ln_v, s=25, color='red', label='Dados')
ax.plot(xx,func(xx),linestyle = '--', linewidth = 1.2, color='black', label='Ajuste' )
ax.set_title("Interpolação $ln(v)\ x\ ln(r)$")
ax.set_xlabel("ln(r)")
ax.set_ylabel("ln(v) ")
ax.grid(color = 'black')
plt.legend(loc='best')
plt.savefig('linear_fit.png', format='png')
plt.show()

M = np.exp(b)/G

print('\n\n M = ',M)


##### Incerteza de M ########

delta_M = (np.exp(b)/G)*np.sqrt(G*np.mean(b_i)**2 )
print('\n\n Incerteza M:', delta_M)


# In[7]:


def expression(x2):
    return np.sqrt((G*M)/x2)

xx = np.linspace(0, 7e12, 200)
yy = expression(xx)

fig, ax = plt.subplots(figsize=(10,10))
#ax = fig.add_subplot()

ax.scatter(dist, vel, s=25, color='red', label = 'Dados')
ax.plot(xx, expression(xx), linestyle = '--', linewidth = 1.2, color='black', label ='Ajuste')
ax.grid(color = 'black')
ax.set_title("Interpolação $V\ x\ r$")
ax.set_xlabel("Distância (m)")
ax.set_ylabel("Velocidade ($\dfrac{m}{s}$)")
plt.legend(loc = 'best')
plt.savefig('poli_fit.png', format='png')
plt.show()


# In[8]:


fig, ax = plt.subplots(figsize=(10,10))
ax.scatter(dist, vel, s=25, color='red', label = 'Dados')
ax.grid(color = 'black')
ax.set_title("Dispersão $V\ x\ r$")
ax.set_xlabel("Distância (m)")
ax.set_ylabel("Velocidade ($\dfrac{m}{s}$)")
plt.legend(loc = 'best')
plt.savefig('scatter.png', format='png')
plt.show()


# In[ ]:




