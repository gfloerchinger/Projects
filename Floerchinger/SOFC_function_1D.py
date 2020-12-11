# -*- coding: utf-8 -*-
"""
This function creates a resirual vector for  1-D IT-SOFC cell. The function 
uses Butler Volmer kinetics to resolve current densities and (eventually ) 
incorporates the Dusty gas model to model massdiffusion of species.


@author: gusfl
"""
import numpy as np
from math import exp


####constants and useful values######

R= 8.3145 #universal gas constant [J/mol-K]
F = 96485 #Faradays constant [C/mol]


#Solution vector takes the form: [phi_an_DL,phi_ca_DL,C_k_an,C_k_ca]

def residual (t,SV,param,ptr):
    dSV_dt = np.zeros_like(SV)
    
    
    ###########charge transport ###################3
    "Anode"    
    #Activation Overpotental
    eta_an = SV[ptr.phi_dl_an] - param.delta_phi_an_eq
    
    #The equation for anode exchange current density is given in Leah et al
    i_o_an = param.K_a* (R*param.T)/(param.n*F)*exp(-param.E_a_an/(R*param.T))
    
    #Butler-Volmer 
    i_far_an = i_o_an*(exp((param.alpha_an_a*param.n*F*eta_an)/(R*param.T)) - 
                       exp((-param.alpha_an_c*param.n*F*eta_an)/(R*param.T)))


    i_dl_an = i_far_an - param.i_ext

    #governing eqns
    dSV_dt[ptr.phi_dl_an] = -i_dl_an/param.C_dl_an 
    

    
    "Anode"    
    #Activation Overpotental
    eta_ca = SV[ptr.phi_dl_ca] - param.delta_phi_ca_eq
    
    
    #The equation for anode exchange current density is given in Leah et al
    i_o_ca = param.K_c* (R*param.T)/(param.n*F)*exp(-param.E_a_ca/(R*param.T))
    
    #Butler-Volmer 
    i_far_ca = i_o_ca*(exp((param.alpha_ca_a*param.n*F*eta_ca)/(R*param.T)) - 
                       exp((-param.alpha_ca_c*param.n*F*eta_ca)/(R*param.T)))


    i_dl_ca = i_far_ca - param.i_ext

    #governing eqns
    dSV_dt[ptr.phi_dl_ca] = -i_dl_ca/param.C_dl_ca 
    

    return dSV_dt




def species_flux:
    
















                    