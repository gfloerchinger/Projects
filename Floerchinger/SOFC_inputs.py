# -*- coding: utf-8 -*-
"""
This Function collects all user defined unputsand physical parameters 
for the model that can be passed into the solver.

@author: gusfl
"""

import numpy as np

R = 8.3145
F = 96485

n = 2   #charge transfer number

T = 873     #Temperature [K]


i_ext =  0       #applied current [A/m^2]

t_final =  100    # solve time [sec]


#Mol fractions of fuel in anode {H2, H2O}
X_k_an_0 = np.array([0.97, 0.03])
#Mol fractions of fuel in anode {O2, N2}
X_k_ca_0 = np.array([0.21, 0.79])


 

" THIS IS A GUESS :eventually it will depend on the gibbs gree energy and the activities"
delta_phi_an_eq =  0.6 #Anode overpotental [V]
delta_phi_ca_eq =  0.4 #Cathode overpotental [V]

#delta_phi_an_eq = delta_G_0_an - R*T/(n*F)*logproduct(np.power(X_an_0,nu_an))
#delta_phi_ca_eq = delta_G_0_ca - R*T/(n*F)logproduct(np.power(X_ca_0,nu_ca))


" THIS IS A GUESS :Charge transfer inputs "
C_dl_an = 1e2 # F/m2
C_dl_ca = 1e2 # F/m2

"THese are GUESSES"
phi_an_0 = 0
phi_elyte_0 = 0.6
phi_ca_0 = 1

"IT-SOFC metal supported parameters from expirimental data (Leah et al.)"
th_an = 1.5E-5 #anode thiskness[m]
th_ca = 5E-5 #cathode thickness [m]

th_sub = 3.0E-4 #substraight thickness [m]
th_elyte = 1.5E-5 #electrolyte thickness [m] 

sigma_anode = 8.0E4 #anode electrical conductivity [S/m]
sigma_cathode = 8.4E3 #cathode electrical conductivity [S/m]

E_a_ca = 1.309E5   #Cathode activation energy [J/mol]
E_a_an = 1.295E5   #Anode activation energy [J/mol]

K_a = 3.2E13 #anode pre-expoenntal factor [S/m^2]
K_c = 7.0E11 #cathode pre-expoenntal factor [S/m^2]


#From Chris Wendel's thesis

alpha_an_a = 0.4
alpha_an_c = 0.6

alpha_ca_a = 0.5
alpha_ca_c = 0.5


delta_phi_dl_an_0 = phi_elyte_0 - phi_an_0
delta_phi_dl_ca_0 = phi_ca_0 - phi_elyte_0

SV_0 = np.array([delta_phi_dl_an_0,delta_phi_dl_ca_0])



class param:
    
    time_span = np.array([0,t_final])
    
    n = n
    T = T

    i_ext = i_ext  #applied current [A/m^2]

    t_final = t_final    # solve time [sec]

    delta_phi_an_eq = delta_phi_an_eq  #Anode overpotental [V]
    delta_phi_ca_eq = delta_phi_ca_eq #Cathode overpotental [V]
    

    #IT-SOFC metal supported parameters from expirimental data (Leah et al.)
    E_a_ca = E_a_ca   #Cathode activation energy [J/mol]
    E_a_an = E_a_an   #Anode activation energy [J/mol]
    
    C_dl_an = C_dl_an
    C_dl_ca = C_dl_ca

    K_a = K_a
    K_c = K_c

    E_a_an = E_a_an
    E_a_ca = E_a_ca
    
    alpha_an_a = alpha_an_a
    alpha_an_c = alpha_an_c

    alpha_ca_a = alpha_ca_a
    alpha_ca_c = alpha_ca_c
    
    
"Need pointers"


class ptr:
    phi_dl_an = 0
    
    phi_dl_ca = 1


