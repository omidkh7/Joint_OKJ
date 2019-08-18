#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun June 23 17:13:20 2019

@author: omidkhajehdehi
"""
import numpy as np
from scipy.fftpack import fft2, ifft2 ,fft,ifft
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft2, ifft2 ,fft,ifft
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 4, figsize=(12,8) )  

n = 1052* 2          # number of bonds in the 2D square lattic based on my assumption about the lattice
m = 1052          # size of lattice 
ngrid= 1051
sigma = 1.0
Mean = 0.0
No_of_sigma = 2.0 # number of sigma away from average

G = np.zeros((n,m))
for i in range(n):
    for j in range(m):
        x = Mean + No_of_sigma * sigma
        while abs(x - Mean) >=   No_of_sigma * sigma :
            x = np.random.normal(loc= Mean , scale =sigma )
        assert abs(x - Mean) <   No_of_sigma * sigma
        G[i][j] = x

sigma_rescale = np.std(G)    
       
a = 1.0 # Variance
b = 0.0  # mean
G =( a*G+b ) # Transfering G to a  mean and sqrt(1/2^5) variance. 

HIST_G=[]
for i in range(np.shape(G)[0]):
    for j in range(np.shape(G)[1]):
        HIST_G.append(G[i][j])
    
plt.sca( axes[0][0] ) 
plt.hist(HIST_G,128,color='r',label='G')

plt.sca( axes[1][0] )   
neg=plt.imshow(G, cmap='inferno', vmin=np.min(G), vmax=np.max(G),interpolation='none')
cbar =fig.colorbar(neg)
    

fft2G=fft2(G)                     # Fast Fourier Transform of uncorrelated generator array
fft2Gshift=np.fft.fftshift(fft2G) # invert array from peripheral to radial storage

x = np.linspace( -n/2,n/2-1,n)     # set up filter array
y = np.linspace( -(m+1)/2,(m+1)/2-1,m)
X,Y = np.meshgrid(y,x)
R= np.sqrt(X**2+Y**2) + 1

p= 1.2    #set filter power….p ~ 1 gives 1/k noise
                         
fft2Gshift_filt=fft2Gshift/(R**p) #apply filter array raised to exponent ‘p’
fft2G_filt=np.fft.fftshift(fft2Gshift_filt) # invert array from radial to peripheral storage
F= np.real(ifft2(fft2G_filt))      # take real part of inverse transform of filtered array

F = (F-(np.mean(F))) / (np.std(F)) #zero-mean/unit-variance spatially correlate array
F =( sigma_rescale * F + b )
F = np.exp(F) # making a log-normal distribution

# In this part, the most permeable site would be found and
# then using numpy roll the whole permeability field would be shifted
Index_max_F = np.where(F == np.amax(F))
Index_max_row = Index_max_F[0][0]
Index_max_column = Index_max_F[1][0]
row_center = int(n/2)
column_center = int(m/2)

if( Index_max_row > row_center ):
    diff = int(Index_max_row - row_center)
    F = np.roll(F, -diff, axis = 0) #up
        
if( Index_max_row < row_center ):
    diff = int(row_center -  Index_max_row )
    F = np.roll(F, diff, axis = 0) #down       
        
if( Index_max_column < column_center ):
    diff = int(column_center -  Index_max_column) 
    F = np.roll(F, diff, axis = 1) #right            
    
if( Index_max_column > column_center ):
    diff = int(Index_max_column  - column_center)
    F = np.roll(F, -diff, axis = 1) #left      
        
print(np.where(F == np.amax(F)))


print('Mean Permeability : ',np.mean(F))
print('Std  Permeability : ',np.std(F))
print('Max  Permeability : ',np.max(F))
print('Min  Permeability : ',np.min(F))
Real_max_F = np.max(F)


#Calculation of alpha 
# 0.35 maximum porosity we want, 0.01 minimum porosity we want
alpha = (1 / (0.35 - 0.1) ) * np.log(np.max(F) / np.min(F) )
print('Calculated Alpha : ',alpha)
#porosity should be in between 0.01 to 0.35
Porosity =( ( (1/alpha) * np.log( F / np.min(F) ) ) + 0.1 ) # alpha = 20.0 comes form empirical data
print('Mean Porosity: ',np.mean(Porosity))
print('Std Porosity : ',np.std(Porosity))
print('Max porosity : ',np.max(Porosity))
print('Min porosity : ',np.min(Porosity))

plt.sca( axes[0][1] )
plt.title('porosity')
HIST_Porosity=[]

for i in range(np.shape(Porosity)[0]):
    for j in range(np.shape(Porosity)[1]):
        HIST_Porosity.append(Porosity[i][j]) 
        
plt.hist(HIST_Porosity,128,color='b',label='F')

plt.sca( axes[1][1] ) 
neg=plt.imshow(Porosity, cmap='inferno', vmin=np.min(Porosity), vmax=np.max(Porosity))
cbar =fig.colorbar(neg)


plt.sca( axes[0][2] )
plt.title('Permeability')       
HIST_F=[]
Counter = 0 
for i in range(np.shape(F)[0]):
    for j in range(np.shape(F)[1]):
        if (F[i][j] > 20):
            Counter +=1
            #F[i][j] = np.random.uniform(low= 10.0,high= 20.0)
        HIST_F.append(F[i][j]) 

#print('percentile of sites bigger than sth:', (Counter * 100 /(m*n) ) )
#
#print('After putting upper limit')        
#print('Mean Permeability : ',np.mean(F))
#print('Std  Permeability : ',np.std(F))
#print('Max  Permeability : ',np.max(F))
#print('Min  Permeability : ',np.min(F))



plt.hist(HIST_F,128,color='b',label='F')

plt.sca( axes[1][2] ) 
neg=plt.imshow(F,cmap='inferno', vmin=np.min(F), vmax= 5.0 ,interpolation='none')
cbar =fig.colorbar(neg)


plt.sca( axes[0][3] )
HIST_Per_pos=[]

Final = np.zeros((n,m))
for i in range(np.shape(Final)[0]):
    for j in range(np.shape(Final)[1]):
        Final[i][j] = F[i][j] / Porosity[i][j]
        HIST_Per_pos.append( Final[i][j] ) 
        
plt.hist(HIST_Per_pos,128,color='b',label='F')

plt.sca( axes[1][3] ) 
neg=plt.imshow(Final, cmap='inferno', vmin=np.min(Final), vmax=40)#np.max(Final))
cbar =fig.colorbar(neg)

print('average perpor: ',np.mean(Final))
print('Std perpor : ',np.std(Final))
print('Max perpor : ',np.max(Final))
print('Min perpor : ',np.min(Final))

#calculating the optimum dt for solving diffusion equation
viscosity = 1e-4
compressibility = 1e-8
K_char = 1e-16 # characteristic Permeability
dt = 0.25 * ( (viscosity * compressibility * np.min(Porosity) )) * (1/(np.max(F) * K_char)) 
print('Optimum dt : ', dt )

sfile = open( 'K_x.txt', 'w' )
file = open( 'K_y.txt', 'w' )
phiX = open( 'Phi_x.txt', 'w' )
phiY = open( 'Phi_y.txt', 'w' )
for i in range( 2 * ngrid ):
    for j in range( ngrid ):
        if ( i % 2 == 0 ):
            file.write('%s\n'%  F[i][j] ) 
            phiX.write('%s\n'%  Porosity[i][j])
        else:    
            sfile.write('%s\n'% F[i][j])
            phiY.write('%s\n'%  Porosity[i][j])
sfile.close()
file.close()
phiX.close()
phiY.close()

ssfile = open( 'yield.txt', 'w' )
for i in range( ngrid * ngrid ):
    ssfile.write('%s\n'% np.random.uniform(low= 4.0,high= 5.0))
ssfile.close()

sssfile = open( 'shear.txt', 'w' )
for i in range( ngrid * ngrid ):
    sssfile.write('%s\n'%  np.random.uniform(low= 0.0,high= 1.0))
sssfile.close()