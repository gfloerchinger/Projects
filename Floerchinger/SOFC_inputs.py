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
P_amb = 100000# ambiant presure [Pa]
P_an_0 = 100000  #Pressure [Pa]

i_ext =  10000       #applied current [A/m^2]

t_final =  1000000    # solve time [sec]


#Mol fractions of fuel in anode {H2, H2O}
X_k_an_0 = np.array([0.97, 0.03])
#Mol fractions of fuel in anode {O2, N2}
X_k_ca_0 = np.array([0.21, 0.79])


mu_g_an = 3.22E-5  # Dynamic viscsity, Pa-s

#from table 3 of Zerki et al.
eps_g_CL =  0.38 # Volume fraction of gasses in anode
eps_solid_CL = 0.62 # Volume fraction of solids (Ni+CGO) in anode

d_part_CL = 1E-6 #particle diamater of anode

#stoich coeff's for anode
nu_k_an = np.array([-1,1])
nu_k_ca = np.array([-0.5,0])
" THIS IS A GUESS :eventually it will depend on the gibbs gree energy and the activities"
delta_phi_an_eq =  -0.6 #Anode overpotental [V]
delta_phi_ca_eq =  0.4 #Cathode overpotental [V]

delta_G_an_0 = -228572 #Gibb free energy at standard state [J/mol]
delta_G_ca_0 = 0 #Gibb free energy at standard state [J/mol]

#Equilibrium potentals
delta_phi_an_eq = -delta_G_an_0/(n*F) - R*T/(n*F)*np.log(np.product(np.power(X_k_an_0,nu_k_an)))
delta_phi_ca_eq = -delta_G_ca_0/(n*F) - R*T/(n*F)*np.log(np.product(np.power(X_k_ca_0,nu_k_ca)))


" THIS IS A GUESS :Charge transfer inputs "
C_dl_an = 1e6 # F/m2
C_dl_ca = 1e6 # F/m2

"THese are GUESSES"
phi_an_0 =1.1
phi_elyte_0 = 0
phi_ca_0 = -0.2

"IT-SOFC metal supported parameters from expirimental data (Leah et al.)"
th_an = 1.5E-5 #anode thickness[m]
th_ca = 5E-5 #cathode thickness [m]

th_sub = 3.0E-4 #substraight thickness [m]
th_elyte = 1.5E-5 #electrolyte thickness [m] 

sigma_anode = 8.0E4 #anode electrical conductivity [S/m]
sigma_cathode = 8.4E3 #cathode electrical conductivity [S/m]

E_a_ca = 1.309E5   #Cathode activation energy [J/mol]
E_a_an = 1.295E5   #Anode activation energy [J/mol]

K_a = 3.2E13 #anode pre-expoenntal factor [S/m/bar^0.5]
K_c = 7.0E11 #cathode pre-expoenntal factor [S/m^2]


#From Chris Wendel's thesis

alpha_an_a = 0.4
alpha_an_c = 0.6

alpha_ca_a = 0.5
alpha_ca_c = 0.5


delta_phi_dl_an_0 = phi_an_0 - phi_elyte_0
delta_phi_dl_ca_0 = phi_ca_0 - phi_elyte_0


C_k_an_0 = P_an_0*X_k_an_0/R/T


SV_0 = np.hstack([delta_phi_dl_an_0, delta_phi_dl_ca_0, C_k_an_0, C_k_an_0])


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


#diffusion volumes for fuller approx 
V_H2 = 6.12
V_H2O = 12.7
V_N2 = 17.9
V_O2 = 16.6


class param:
    
    time_span = np.array([0,t_final])
    
    n = n
    T = T
    P_an = P_an_0
    P_amb = P_amb
    
    i_ext = i_ext  #applied current [A/m^2]

    t_final = t_final    # solve time [sec]

    delta_phi_an_eq = delta_phi_an_eq  #Anode overpotental [V]
    delta_phi_ca_eq = delta_phi_ca_eq #Cathode overpotental [V]
    
    delta_G_an_0 = delta_G_an_0 #Gibb free energy at standard state [J/mol]
    delta_G_ca_0 = delta_G_ca_0 #Gibb free energy at standard state [J/mol]
    
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
    
    th_an = th_an
    th_ca = th_ca

    th_sub = th_sub
    th_elyte = th_elyte
    
    mu_g_an = mu_g_an
    
    nu_k_an = nu_k_an
    nu_k_ca = nu_k_ca
    
    eps_g_CL = eps_g_CL

    d_part_CL = d_part_CL
    
    A_fac_dl = 3*th_an*eps_solid_CL/(d_part_CL/2)    
    
    d_pore_an = d_pore_an
    d_pore_ca = d_pore_ca
    
    eps_an = eps_an
    eps_ca = eps_ca
    
    tau_an = tau_an
    tau_ca = tau_ca
    
    MM_f = [MM_H2,MM_H2O]
    MM_a = [MM_O2,MM_N2]

    X_k_an_0 = X_k_an_0
    X_k_ca_0 = X_k_ca_0

    V_f = [V_H2,V_H2O]
    V_a = [V_O2,V_N2]

"Need pointers"


class ptr:
    phi_dl_an = 0
    
    phi_dl_ca = 1
    
    # C_k in anode substrate: starts just after phi_dl, is same size as X_k_an:
    C_k_an_sub = np.arange(phi_dl_ca+1, phi_dl_ca+1+X_k_an_0.shape[0])
    
    # C_k in anode CL: starts just after substrate, is same size as X_k_an:
    C_k_an_CL = np.arange(C_k_an_sub[-1]+1, C_k_an_sub[-1]+1+X_k_an_0.shape[0])


