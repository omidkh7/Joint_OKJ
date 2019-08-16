#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 12:24:14 2019

@author: omidkhajehdehi
"""

import numpy as np
from scipy.fftpack import fft2, ifft2 ,fft,ifft
import matplotlib.pyplot as plt

#fig, axes = plt.subplots(2, 4, figsize=(12,8) )  

n = (4001+1)*2 #502* 2          # number of bonds in the 2D square lattic based on my assumption about the lattice 
m = 4001+1 #502          # size of lattice 
ngrid= 4001 #501
G = np.zeros((n,m))

sigma=1.0
for i in range(n):
#	for j in xrange(m):
#		x=3*sigma
#		while abs(x)>2*sigma:
#			x=np.random.normal()
#		assert abs(x)<=2*sigma
#		G[i,j]=x
#    G[i]=np.random.randn(m)
#    G[i]=np.random.randn(m)
	G[i]=[0.0 for i in xrange(m)] #kam #np.random.randn(m)
    
       
a = 1.0 # Variance
b = 0.0  # mean
G =( a*G+b ) # Transfering G to a  mean and sqrt(1/2^5) variance. 

HIST_G=[]
for i in range(np.shape(G)[0]):
    for j in range(np.shape(G)[1]):
        HIST_G.append(G[i][j])
    
#plt.sca( axes[0][0] ) 
#plt.hist(HIST_G,128,color='r',label='G')

#plt.sca( axes[1][0] )   
#neg=plt.imshow(G, cmap='inferno', vmin=np.min(G), vmax=np.max(G),interpolation='none')
#cbar =fig.colorbar(neg)
    
F = np.exp(G) # making a log-normal distribution
print('average F: ',np.mean(F))
print('Std F : ',np.std(F))
print('Max Permeability : ',np.max(F))
print('Min Permeability : ',np.min(F))
#porosity should be in between 0.01 to 0.35
Porosity =( ( (1/15) * np.log( F / np.min(F) ) ) + 0.01 ) # alpha = 20.0 comes form empirical data
Porosity = [[1.0 for j in range(m)]for i in range(n)] #kam
print('average Porosity: ',np.mean(Porosity))
print('Std Porosity : ',np.std(Porosity))
print('Max porosity : ',np.max(Porosity))
print('Min porosity : ',np.min(Porosity))

#plt.sca( axes[0][1] )
#plt.title('porosity')
HIST_Porosity=[]

for i in range(np.shape(Porosity)[0]):
    for j in range(np.shape(Porosity)[1]):
        HIST_Porosity.append(Porosity[i][j]) 
        
#plt.hist(HIST_Porosity,128,color='b',label='F')

#plt.sca( axes[1][1] ) 
#neg=plt.imshow(Porosity, cmap='inferno', vmin=np.min(Porosity), vmax=np.max(Porosity))
#cbar =fig.colorbar(neg)


#plt.sca( axes[0][2] )
#plt.title('Permeability')       
HIST_F=[]
beshmore = 0
for i in range(np.shape(F)[0]):
    for j in range(np.shape(F)[1]):
        if (F[i][j] > 10):
            beshmore += 1
            F[i][j] = np.random.uniform(low= 5.0,high= 10.0)
        HIST_F.append(F[i][j]) 


print( 'fraction of number bigger than sth: ' ,(beshmore / (n*m*2)) * 100 )
#plt.hist(HIST_F,128,color='b',label='F')

#plt.sca( axes[1][2] ) 
#neg=plt.imshow(F,cmap='inferno', vmin=np.min(F), vmax= 10 ,interpolation='none')
#cbar =fig.colorbar(neg)


#plt.sca( axes[0][3] )
HIST_Per_pos=[]

Final = np.zeros((n,m))
for i in range(np.shape(Final)[0]):
    for j in range(np.shape(Final)[1]):
        Final[i][j] = F[i][j] / Porosity[i][j]
        HIST_Per_pos.append( Final[i][j] ) 
        
#plt.hist(HIST_Per_pos,128,color='b',label='F')

#plt.sca( axes[1][3] ) 
#neg=plt.imshow(Final, cmap='inferno', vmin=np.min(Final), vmax=np.max(Final))
#cbar =fig.colorbar(neg)


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
    ssfile.write('%s\n'% np.random.uniform(low= 1.0,high= 1.5))
ssfile.close()

s3sfile = open( 'shear.txt', 'w' )
for i in range( ngrid * ngrid ):
    s3sfile.write('%s\n'%  np.random.uniform(low= 0.0,high= 1.0))
s3sfile.close()


#s1sfile = open( 'yfunc.txt', 'w' )
#for i in range( ngrid * ngrid ):
#    s1sfile.write('%s\n'%  - np.random.uniform(low= 0,high= 1.0))
#s1sfile.close()
#
#s2sfile = open( 'simga.txt', 'w' )
#for i in range( ngrid * ngrid ):
#    s2sfile.write('%s\n'%   np.random.uniform())
#s2sfile.close()
#
#file = open( 'Permeability.xyz', 'w' )
#particles = ngrid * ngrid
#file.write('%d\n'% particles )
#file.write('Permeability\n')
#for i in range( 2 * ngrid ):
#    for j in range( ngrid ):
#        if ( i % 2 == 0 ):
#            file.write('%d\t%d\t%f\n' % (i,j,F[i][j]))
#        else:    
#            file.write('%d\t%d\t%f\n' % (i,j,F[i][j]))
#
#file.close()   
