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
    
    #get partial pressure of hydrogen
    P_H2 = C_k_an_CL[0]*R*param.T/100000
    
    
    #Activation Overpotental
    eta_an = SV[ptr.phi_dl_an] - param.delta_phi_an_eq
    

    #The equation for anode exchange current density is given in Leah et al
    i_o_an = (P_H2)**0.5 * param.K_a* (R*param.T)/(param.n*F)*exp(-param.E_a_an/(R*param.T))


    
    #Butler-Volmer 
    i_far_an = i_o_an*(exp((param.alpha_an_a*param.n*F*eta_an)/(R*param.T)) - 
                       exp((-param.alpha_an_c*param.n*F*eta_an)/(R*param.T)))

 
    "Cathode"    
    
    #Activation Overpotental
    eta_ca = SV[ptr.phi_dl_ca] - param.delta_phi_ca_eq
    
    #The equation for anode exchange current density is given in Leah et al
    i_o_ca = param.K_c* (R*param.T)/(param.n*F)*exp(-param.E_a_ca/(R*param.T))
    
    #Butler-Volmer 
    i_far_ca = i_o_ca*(exp((param.alpha_ca_a*param.n*F*eta_ca)/(R*param.T)) - 
                       exp((-param.alpha_ca_c*param.n*F*eta_ca)/(R*param.T)))

    ##########################################################
    
    i_far = np.array([i_far_an, i_far_ca])


    return i_far