# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:35:30 2020

@author: gfloerchinger
"""

import numpy as np

T = 873
P = 101325
R = 8.3145

d_pore_an = 1E-6
eps_an = 0.26
tau_an = 3

d_pore_ca = 1E-6
eps_ca = 0.3
tau_ca = 3


#moleqular masses
MM_H2 = 2.016
MM_H2O = 18.015
MM_N2 = 28.0134
MM_O2 = 16.000

MM_f = [MM_H2,MM_H2O]
MM_a = [MM_O2,MM_N2]

X_H2 = 0.97
X_H2O = 0.03
X_O2 = 0.21
X_N2 = 0.79

X_f = [X_H2,X_H2O]
X_a = [X_O2,X_N2]


#diffusion volumes for fuller approx 
V_H2 = 6.12
V_H2O = 12.7
V_N2 = 17.9
V_O2 = 16.6

V_f = [V_H2,V_H2O]
V_a = [V_O2,V_N2]


#def difn_coeff(param,gas_props):
N = len(X_f)
    
D_m_ij = np.zeros((N,N))
D_kn_i = np.zeros((N,1))
D_m_k_an = np.zeros((N,1))


"anode"

for i, k in enumerate(D_m_ij):
    #knudsen Diffusion

    D_kn_i[i] = 2/3*d_pore_an * np.sqrt(8*R*T/np.pi/MM_f[i])
    print(D_kn_i)
    for j, l in enumerate(k):
         #Fuller Approximation for effective mixture diffusion   
        #print(i,j)
        MM_ij = 2*(1/MM_f[i] + 1/MM_f[j])**(-1)
        D_m_ij[i,j] = 0.00143*T**(1.75) / (P*MM_ij**0.5 * (V_f[i]**(1/3) + (V_f[j]**(1/3))))
        
     
    #mixture averaged diffusion coeff    
    D_m_k_an[i] = np.sum(x for l,x in enumerate(np.divide(X_f,D_m_ij[i,:])) if l!=i)
    D_m_k_an[i] = (1-X_f[i])/D_m_k_an[i]


#Effective diffusion coefficent takes knudsen and mixture difn in series
D_eff_k_an = eps_an/tau_an * 1./(1./D_kn_i + 1./D_m_k_an)
    


