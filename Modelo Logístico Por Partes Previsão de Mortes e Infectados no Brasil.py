#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#####IMPORTAR  TRANSFORMAR EM ARRAY OS DADOS ##################################################################
df = pd.read_csv("data_15_03_2021_Brazil.csv", sep = ";", header=0)
df.head()
df = pd.read_csv("data_15_03_2021_Brazil.csv", sep = ";", header=0, usecols=["casosAcumulado","obitosAcumulado"])
df.head()
infec = df["casosAcumulado"].to_numpy()
M = df["obitosAcumulado"].to_numpy()
print(len(infec),len(M))
######################### TRANSFORMANDO EM ARRAY PARA TRABALHAR NUMERICAMENTE ################################
t = []
casos = []
obitos = []
for i in range(len(infec)):
    t.append(i+1)
    casos.append(infec[i])
    obitos.append(M[i])
print(len(t),len(casos),len(obitos))
print(t[1])
######################### MODELO QUE PODE SER AJUSTADO PARA INFECTADOS ########################################
####### BASTA TROCAR OS INTERVALOS ######################################## PARA SIMULAR PARA DIFERENTES DADOS


# metodo da diferença finita para infectados no Br
#para i = 1
df_infec = []
primeiro_infec = (casos[1] - casos[0])/(t[1]-t[0])
df_infec.append(primeiro_infec)
print(df_infec)

#para 1<i<n
for i in range(len(t) - 2):
    j = i+1
    k = i+2
    y_1 = ((casos[k] - casos[j])/(t[k]-t[j]))/casos[j]
    y_2 = ((casos[j] - casos[i])/(t[j]-t[i]))/casos[j]
    delta = (y_1+y_2)/2
    df_infec.append(delta)

#para i = n

n = len(casos)-1

ultimo_infec = ((casos[n] - casos[n-1])/(t[n]-t[n-1]))/casos[n]
df_infec.append(ultimo_infec)

#Minimos quadrados inicio a partir do primeiro dia de março:
quadinf = []
prodxy = []
casos1 = []
df_infec1 = []
for i in range(len(df_infec)-371):
    j = i+371
    cquad = casos[j]**2
    pxy = casos[j]*df_infec[j]
    quadinf.append(cquad)
    prodxy.append(pxy)
    casos1.append(casos[j])
    df_infec1.append(df_infec[j])
    
    
x_m = np.mean(casos1)
y_m = np.mean(df_infec1)
xq_m= np.mean(quadinf)
pxy_m=np.mean(prodxy)

    #determinação dos coeficientes

mc = (x_m*y_m - pxy_m)/(x_m**2 - xq_m)
bc = y_m - mc*x_m
        
    #fim da determinação dos coeficientes

    #inicio função
def minquad(x):
    return mc*x + bc
    #fim da função
    
    #calculo de R^2 sobre a variancia de y que explica x
s = 0
s1 = 0
for i in range(len(df_infec1)):
    s += (df_infec1[i] - minquad(casos1[i]))**2
    s1+= (df_infec1[i] - y_m)**2

R = 1 - s/s1

    #Fim do calculo R^2


    #Inicio da determinação do erro dos coeficientes
    
a = ((2*pxy_m*x_m - (x_m**2)*y_m - y_m*xq_m)/((x_m**2 - xq_m)**2))*np.var(x_m)
b = (x_m/(x_m**2 - xq_m))*np.var(y_m)
c = (1/(x_m**2 - xq_m))*np.var(pxy_m)
d = ((x_m*y_m - pxy_m)/((x_m**2 - xq_m)**2))*np.var(xq_m)

dmc = np.sqrt(a**2 + b**2 + c**2 + d**2)

e = np.var(y_m)
f = -x_m*dmc
g = -mc*np.var(x_m)

dbc = np.sqrt(e**2 + f**2 + g**2)

    #fim da determinação do erro dos coeficientes




#fim dos minimos quadrados



x = np.linspace(0,1.2e7)
y = minquad(x)

print(len(df_infec),len(df_infec1), len(quadinf))

fig, ax = plt.subplots(figsize=(8,8))
ax.scatter(casos, df_infec, s=10, color='red', label='Dados Brasil')
ax.scatter(x, y, s=10, color='black', label='Minimo Quadrado')
#ax.set_yscale('log')
#plt.ylim((-1, 1.25e7)) 
#ax.plot(tamanho, area, '-', color='black', label='Suavização')
ax.set_title('Diferença Finita')
ax.set_xlabel('Casos')
ax.set_ylabel('Diferença Finita')
ax.grid()
plt.legend(loc='best')
plt.show()

#print(df_infec)
#np.savetxt('difinf.txt', df_infec, newline='\n')
#np.savetxt('casos.txt', casos, newline='\n')
#print(df_infec[3])

print('\n\nCoeficientes:','\na: ',mc,'+\-',dmc,'\nb: ',bc,'+\-',dbc,'\n\nR²: ', round(R,4)*100,'%')


#definindo para todos os dias de março até o dia atual

k = -mc
y_M = bc/k
#definir a função y(t)

def equ(x):
    return y_M/(1+ ((y_M - casos[371])/casos[371])*np.exp(-k*y_M*(x - t[371])))


x = np.arange(0,450,10)
y = equ(x)


fig, ax = plt.subplots(figsize=(8,8))
#ax.scatter(x, y, s=10, color='red', label='Dados Brasil')
#ax.set_yscale('log')
#plt.ylim((0, 690)) 
ax.plot(x, y, '-', color='black', label='Suavização')
ax.set_title('Mês de Março - Casos')
ax.set_xlabel('Dias dia de Março')
ax.set_ylabel('Infectados')
ax.grid()
plt.legend(loc='best')
plt.show()

for i in range(len(df_infec1)):
    print('Certeza: ',round(equ(i+372)/casos[i+371],4)*100,'%','\nCasos: ',casos1[i],' Teórico: ', round(equ(i+372)))


print('Dia 0: ',round(equ(0)),'\nÚltimo dia teorico',round(equ(394)) )
print('Ultimo dia',casos[393])
print('Certeza: ',round(equ(394)/casos[393],4)*100,'%')

####################################################################################################################


####################################################################################################################


        ################## IMPLEMENTAÇÃO PARA ÓBITOS ###############################################
    
####################################################################################################################


####################################################################################################################


# metodo da diferença finita para obitos no Br
#para i = 1
df_obt = []
primeiro_obt = (obitos[1] - obitos[0])/(t[1]-t[0])
df_obt.append(primeiro_obt)
#print(df_obt)

#para 1<i<n
for i in range(len(t) - 2):
    j = i+1
    k = i+2
    if obitos[i]!=0:
        y_1 = ((obitos[k] - obitos[j])/(t[k]-t[j]))/obitos[j]
        y_2 = ((obitos[j] - obitos[i])/(t[j]-t[i]))/obitos[j]
        delta = (y_1+y_2)/2
        df_obt.append(delta)
    else:
        (obitos[k] - obitos[j])/(t[k]-t[j])
        (obitos[j] - obitos[i])/(t[j]-t[i])
        delta = (y_1+y_2)/2
        df_obt.append(delta)
#para i = n

n = len(obitos)-1

ultimo_obt = ((obitos[n] - obitos[n-1])/(t[n]-t[n-1]))/obitos[n]
df_obt.append(ultimo_obt)

#Minimos quadrados inicio *****************************************************************************************
quadinf = []
prodxy = []
obitos1 = []
df_obt1 = []
for i in range(len(df_obt)-390):
    j = i+390
    cquad = obitos[j]**2
    pxy = obitos[j]*df_obt[j]
    quadinf.append(cquad)
    prodxy.append(pxy)
    obitos1.append(obitos[j])
    df_obt1.append(df_obt[j])

x_m = np.mean(obitos1)
y_m = np.mean(df_obt1)
xq_m= np.mean(quadinf)
pxy_m=np.mean(prodxy)

    #determinação dos coeficientes

mo = (x_m*y_m - pxy_m)/(x_m**2 - xq_m)
bo = y_m - mo*x_m
    #fim da determinação dos coeficientes
    
   #inicio função
def minquad(x):
    return mo*x + bo
    #fim da função
    
     #calculo de R^2 sobre a variancia de y que explica x
s = 0
s1 = 0
for i in range(len(df_obt1)):
    s += (df_obt1[i] - minquad(obitos1[i]))**2
    s1+= (df_obt1[i] - y_m)**2

R = 1 - s/s1

    #Fim do calculo R^2
    
    
    #Inicio da determinação do erro dos coeficientes
    
a = ((2*pxy_m*x_m - (x_m**2)*y_m - y_m*xq_m)/((x_m**2 - xq_m)**2))*np.var(x_m)
b = (x_m/(x_m**2 - xq_m))*np.var(y_m)
c = (1/(x_m**2 - xq_m))*np.var(pxy_m)
d = ((x_m*y_m - pxy_m)/((x_m**2 - xq_m)**2))*np.var(xq_m)

dmo = np.sqrt(a**2 + b**2 + c**2 + d**2)

e = np.var(y_m)
f = -x_m*dmo
g = -mo*np.var(x_m)

dbo = np.sqrt(e**2 + f**2 + g**2)

    #fim da determinação do erro dos coeficientes

#fim dos minimos quadrados ****************************************************************************************

x1 = np.linspace(0,3e5)
y1 = minquad(x1)


print(len(df_obt), obitos[384],df_obt[384])

fig, ax = plt.subplots(figsize=(8,8))
ax.scatter(obitos, df_obt, s=10, color='red', label='Dados Brasil')
ax.scatter(x1, y1, s=10, color='black', label='Minimo Quadrado')
#ax.set_yscale('log')
#plt.ylim((-1, 1.25e7)) 
#ax.plot(tamanho, area, '-', color='black', label='Suavização')
ax.set_title('Diferença Finita')
ax.set_xlabel('Obitos')
ax.set_ylabel('Diferença Finita')
ax.grid()
plt.legend(loc='best')
plt.show()

#print(df_infec)
#print(df_obt)

#np.savetxt('difobt.txt', obitos, newline='\n')
print(df_obt[3])
print('\n\nCoeficientes:','\na: ',mo,'+\-',dmo,'\nb: ',bo,'+\-',dbo,'\n\nR²: ', round(R,3)*100,'%')

#definindo para semana atual de março para obitos ****************************************************************

k = -mo
y_M = bo/k
#definir a função y(t)

def equ(x):
    return y_M/(1+ ((y_M - obitos[390])/obitos[390])*np.exp(-k*y_M*(x - 391)))


x = np.arange(1,250,1)
y = equ(x)
fig, ax = plt.subplots(figsize=(8,8))
ax.scatter(x, y, s=10, color='red', label='Dados Brasil')
#ax.set_yscale('log')
#plt.ylim((0, 690)) 
ax.plot(x, y, '-', color='black', label='Suavização')
ax.set_title('Obitos')
ax.set_xlabel('Dias dia de Março')
ax.set_ylabel('Obitos')
ax.grid()
plt.legend(loc='best')
plt.show()

print('Teste de certeza: ',(round(equ(392),3)/obitos[391])*100,'%')

for i in range(len(df_obt1)):
    print('\n',19+i,'de Março de 2021','Obitos: ',obitos1[i],' - Teórico: ', round(equ(i+391)), 'Certeza: ',(round(equ(i+391)/obitos[i+390],3))*100,'%')
print('\n\n\nOBITOS 22/03/2021: ',obitos[394],'\nPREVISÃO Total de OBITOS 19->23/03/2021: ', round(equ(396)),'\nDiário: ',round(equ(396))-obitos[394] )
print('\nPrevisão de 14-22 de Março:', 296628,'Diário:',296628-obitos[394],'\nPrevisão de 17-22 de Março:',296720,'Diário',296720-obitos[394])


########################################################################################################################################################

########################################################################################################################################################

########################################################################################################################################################

########################################### FIM ##################################### FIM ##############################################################

