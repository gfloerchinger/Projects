# -*- coding: utf-8 -*-
"""
This function creates a resirual vector for  1-D IT-SOFC cell. The function 
uses Butler Volmer kinetics to resolve current densities and (eventually ) 
incorporates the Dusty gas model to model massdiffusion of species.


@author: gusfl
"""
import numpy as np
from math import exp
from SOFC_charge_transport import Butler_Volmer
from SOFC_gas_transport import  gas_flux

####constants and useful values######

R = 8.3145 #universal gas constant [J/mol-K]
F = 96485 #Faradays constant [C/mol]


#Solution vector takes the form: [phi_an_DL,phi_ca_DL,C_k_an,C_k_ca]

def residual (t,SV,param,ptr):

    dSV_dt = np.zeros_like(SV)
    
    
    
    #######################Charge transport############################
    #Fetch faradiac current from Butler Volmer equation
    i_far =  Butler_Volmer(SV, param, ptr)

    #print(i_far)
    #calculate anode double layer current
    i_dl_an = i_far[0] + param.i_ext
    
    #calculate anode double layer current
    i_dl_ca = i_far[1] + param.i_ext

    #governing eqns for charge transport
    dSV_dt[ptr.phi_dl_an] = -i_dl_an/param.C_dl_an
    dSV_dt[ptr.phi_dl_ca] = -i_dl_ca/param.C_dl_ca 

    
   ######################Species Transport################################ 
   
    C_k_an_sub = SV[ptr.C_k_an_sub] 
    C_k_an_CL = SV[ptr.C_k_an_CL] 
    
    # State variables for node 1: substrate
    state1 = {'C_k':C_k_an_sub, 'T':param.T,'dY': param.th_sub,'mu': param.mu_g_an}
    # State variables for node 2: anode CL
    state2 = {'C_k':C_k_an_CL, 'T':param.T, 'dY': param.th_an, 'mu': param.mu_g_an}
    
    sdot_k_an = i_far[0]*param.nu_k_an/param.n/F
    

    N_k_an = gas_flux(state1,state2,param)
    
    dSV_dt[ptr.C_k_an_CL] = (N_k_an + sdot_k_an*param.A_fac_an)/state2['dY']/param.eps_g_CL


 

    return dSV_dt








                    