# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:49:19 2020

@author: gfloerchinger
"""

import numpy as np
from math import exp

####constants and useful values######

R= 8.3145 #universal gas constant [J/mol-K]
F = 96485 #Faradays constant [C/mol]

def Butler_Volmer(SV, param, ptr):

    "Anode"    
      
    C_k_an_sub = SV[ptr.C_k_an_sub] 
    C_k_an_CL = SV[ptr.C_k_an_CL]
    
    X_k_an = C_k_an_CL/np.sum(C_k_an_CL)
    # print(X_k_an)
    # print(C_k_an_sub)
    

    
    #get partial pressure of hydrogen
    # P_H2 = C_k_an_CL[0]*R*param.T
    # P_H2O = C_k_an_CL[1]*R*param.T
    #print(P_H2,P_H2O)
    
    P_H2 = param.P_an*X_k_an[0]
    P_H2O = param.P_an*X_k_an[1]
    
    P_O2 = param.P_an*param.X_k_ca_0[0]
    
    #print(P_H2/P_H2O)
    
    delta_phi_an_eq = -param.delta_G_ca_0/(param.n*F) - R*param.T/(param.n*F)*np.log(np.product(np.power(param.X_k_an_0,param.nu_k_ca))) - np.log(param.P_an/param.P_amb)
    #print(delta_phi_an_eq)
    #delta_phi_ca_eq = -param.delta_G_ca_0/(param.n*F) - R*param.T/(param.n*F)*np.log(np.product(np.power(param.X_k_ca_0,param.nu_k_ca)))
    
    #print(P_H2)
    
    #Activation Overpotental
    eta_an = SV[ptr.phi_dl_an] - param.delta_phi_an_eq
    #print(eta_an)

    #The equation for anode exchange current density is given in Leah et al
    #i_o_an = (P_H2)**0.5 * param.K_a* (R*param.T)/(param.n*F)*exp(-param.E_a_an/(R*param.T))
    i_o_an = 15600*(P_H2/100000)**0.5*(P_H2O/100000)**0.5 * exp(-60000/R*(1/param.T - 1/873))

    
    #Butler-Volmer 
    i_far_an = i_o_an*(exp((-param.alpha_an_a*param.n*F*eta_an)/(R*param.T)) - 
                       exp((param.alpha_an_c*param.n*F*eta_an)/(R*param.T)))

 
    "Cathode"    
    
    #Activation Overpotental
    eta_ca = SV[ptr.phi_dl_ca] - param.delta_phi_ca_eq
    
    #The equation for anode exchange current density is given in Leah et al
    #i_o_ca = param.K_c* (R*param.T)/(param.n*F)*exp(-param.E_a_ca/(R*param.T))
    i_o_ca = 2470*(P_O2/100000)**0.2 * exp(-162000/R*(1/param.T-1/873))
    

    #Butler-Volmer 
    i_far_ca = i_o_ca*(exp((-param.alpha_ca_a*param.n*F*eta_ca)/(R*param.T)) - 
                       exp((param.alpha_ca_c*param.n*F*eta_ca)/(R*param.T)))

    ##########################################################
    
    i_far = np.array([i_far_an, i_far_ca])
    #print(eta_an, eta_ca)
    #print(i_far)
    #print(i_o_an,i_o_ca)
    
    
    
    
    
    return i_far