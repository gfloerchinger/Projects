# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:35:30 2020

@author: gfloerchinger
"""

import numpy as np

T = 873
P = 101325
R = 8.3145
d_pore = 0.5E-6



D_m_ij = np.zeros((2,2))

#moleqular masses
MM_H2 = 2.016
MM_H2O = 18.015
MM_N2 = 28.0134
MM_O2 = 16.000

MM_f = [MM_H2,MM_H2O]
MM_a = [MM_O2,MM_N2]

#difusion volumes for fuller approx 
V_H2 = 7.07
V_H2O = 12.7
V_N2 = 17.9
V_O2 = 16.6

V_f = [V_H2,V_H2O]
V_a = [V_O2,V_N2]


for i, k in enumerate(D_m_ij):
    #knudsen Diffusion

    D_k_i = 2/3*d_pore* np.sqrt(8*R*T/np.pi/MM_f[i])
    
    for j, l in enumerate(k):
         #Fuller Approximation for effective mixture diffusion   
        #print(i,j)
        MM_ij = 2*(1/MM_f[i] + 1/MM_f[j])**(-1)
        D_m_ij[i,j] = 0.00143*T**(1.75) / (P*MM_ij**0.5 * (V_f[i]**(1/3) + (V_f[j]**(1/3))))

