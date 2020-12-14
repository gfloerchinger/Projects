# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 13:26:41 2020

@author: gusfl
"""
import numpy as np
#from SOFC_inputs import param, SV_0, ptr

R = 8.3145
F = 96485

###############################################


def gas_flux(state1, state2, param):
    
    N_k = np.zeros_like(state1['C_k'])
    
    f1 = state1['dY']/(state1['dY'] + state2['dY'])
    f2 = 1-f1

    C_int = f1*state1['C_k'] + f2*state2['C_k']

    X_k_1 = state1['C_k']/np.sum(state1['C_k'])
    X_k_2 = state2['C_k']/np.sum(state2['C_k'])
    #print(X_k_1, X_k_2)
    X_k_int = f1*X_k_1 + f2*X_k_2
    
    dY = (state1['dY'] + state2['dY'])/2
  

    # #Use Fuller Approximation to find Diffusion coefficents
    D_eff_k = anode_difn_coeffs(X_k_2, param)
    D_eff_k = np.reshape(D_eff_k,(1,2))
    #D_eff_k = np.array([5.48e-4, 6.13e-5])


    # #Convection velocity
    V_conv_k = 0 
    
    # #Diffusion Velocity
    V_difn_k = -D_eff_k * (X_k_2-X_k_1)/dY/X_k_int

    #total flux
    N_k = C_int*X_k_int*(V_conv_k + V_difn_k)

    return N_k


###############################################


def anode_difn_coeffs(X_f, param):
    
    N = len(X_f)
        
    D_m_ij = np.zeros((N,N))
    D_kn_i = np.zeros_like(X_f)
    D_m_k_an = np.zeros_like(X_f)
    
    
    #iterate over speceis, this will be the speceis in question k
    for i in range(N):
        #define matrix to store components of mixture averaged doeff's to be summed over
        D_m_i = np.zeros(N**2-N)
        
        #knudsen Diffusion
    
        D_kn_i[i] = 2/3*param.d_pore_an * np.sqrt(8*R*param.T/np.pi/param.MM_f[i])
        #print(D_kn_i)
        
        #iterate overall species with relation to species k, this is species l
        for j in range(N):
             #Fuller Approximation for effective mixture diffusion   
            
            
            MM_ij = 2*(1/param.MM_f[i] + 1/param.MM_f[j])**(-1)
            
            D_m_ij[i,j] = 0.00143*param.T**(1.75) / (param.P_an*MM_ij**0.5 * (param.V_f[i]**(1/3) + (param.V_f[j]**(1/3))))
               
            
            #check if i=j then calculate components of mixture averaged difn coeff
            if j!=i:
                D_m_i[j] = X_f[j]/D_m_ij[i,j]
            
        #mixture averaged diffusion coeff   
        D_m_k_an[i] = np.sum(D_m_i)

        D_m_k_an[i] = (1-X_f[i])/D_m_k_an[i]
    
    
    #Effective diffusion coefficent takes knudsen and mixture difn in series
    D_eff_k_an = param.eps_an/param.tau_an * 1./(1./D_kn_i + 1./D_m_k_an)
        
    return D_eff_k_an





###############################################

# SV = SV_0

# C_k_an_sub = SV[ptr.C_k_an_sub] 

# C_k_an_CL = SV[ptr.C_k_an_CL] 

# # State variables for node 1: substrate
# state1 = {'C_k':C_k_an_sub, 'T':param.T,'dY': param.th_sub,'mu': param.mu_g_an}
# # State variables for node 2: anode CL
# state2 = {'C_k':C_k_an_CL, 'T':param.T, 'dY': param.th_an, 'mu': param.mu_g_an}

# N_k_an = gas_flux(state1,state2,param)


###############################################