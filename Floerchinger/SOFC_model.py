"""
    This file runs and executes a youyr model, calculating the cell/device/system properties of interest as a function of time. 

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Reads inputs and initializes the model
        
        2 - Calls the residual function and integrates over the user-defined    
            time span.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

# Import necessary modules:
from matplotlib import pyplot as plt
import numpy as np
from SOFC_inputs import param


    
    
def SOFC_model(i_ext = None,P = None):     
    from scipy.integrate import solve_ivp #integration function for ODE system.
    from SOFC_function_1D import residual # point the model to the residual function
    from SOFC_inputs import param, SV_0, ptr
   
    if i_ext:
        param.i_ext = i_ext
    if P:
        param.P_an = P
    
    solution = solve_ivp(lambda t, y: residual(t, y, param, ptr),param.time_span, SV_0, rtol=1e-9, atol=1e-7, method='BDF')
    
    return solution


# if __name__ == '__main__':
#     SOFC_model()
    

    
solution = SOFC_model(1000, 500000)
V_cell = solution.y[1,-1] - solution.y[0,-1]
print(V_cell)

# C_H2_CL = solution.y[4,-1]
# C_H2O_CL = solution.y[5,-1]

# print(solution.y[4,-1],solution.y[5,-1])


########################Current Study######################################    
    
# i_array = np.linspace(0.00000001,10000,50)

# V_cell = np.zeros_like(i_array)
# #Dphi_an = np.zeros_like(i_array)
# #Dphi_ca = np.zeros_like(i_array)
# C_H2_CL = np.zeros_like(i_array)
# C_H2O_CL = np.zeros_like(i_array)

# for j, current in enumerate(i_array):
#     #print(current)
#     solution = SOFC_model(current)
#     V_cell[j] = solution.y[1,-1] - solution.y[0,-1]
#     #Dphi_an[j] = solution.y[0,-1]
#     #Dphi_ca[j] = solution.y[1,-1]
    
#     print(V_cell[j])
    
#     C_H2_CL[j] = solution.y[4,-1]
#     C_H2O_CL[j] = solution.y[5,-1]
#     #print(solution.y[4,-1],solution.y[5,-1])
    
# plt.figure(0)
# plt.plot(i_array,V_cell,'.')
# plt.xlabel('External Current')
# plt.ylabel('Voltage')
# plt.show()

# plt.figure(1)
# plt.plot(i_array,C_H2_CL,'.')
# plt.plot(i_array,C_H2O_CL,'.')
# plt.xlabel('External Current')
# plt.ylabel('Concentration in Anode')
# plt.show()

########################Pressure Study######################################    
i_array = np.array([1000,5000,10000])
pres_array = np.linspace(100000,500000,50)
V_cell = np.zeros([len(pres_array),len(pres_array)])
C_H2_CL = np.zeros_like(V_cell)
C_H2O_CL = np.zeros_like(V_cell)
ASR = np.zeros_like(V_cell)
OCV = np.zeros_like(pres_array)

for i, current in enumerate(i_array):

    for j, pressure in enumerate(pres_array):
        #print(pressure)
        solution = SOFC_model(1E-10,pressure)
        
        OCV[j] = solution.y[1,-1] - solution.y[0,-1]
        
        solution = SOFC_model(current,pressure)
    
        V_cell[i,j] = solution.y[1,-1] - solution.y[0,-1]
        
        ASR[i,j] = (V_cell[i,j]-OCV[j])/param.i_ext
        #Dphi_an[j] = solution.y[0,-1]
        #Dphi_ca[j] = solution.y[1,-1]
        
        print(ASR[i,j])
        
        # C_H2_CL[j] = solution.y[4,-1]
        # C_H2O_CL[j] = solution.y[5,-1]
        #print(solution.y[4,-1],solution.y[5,-1])

    
plt.figure(0)
plt.plot(pres_array,ASR[0,:])
plt.plot(pres_array,ASR[1,:])
plt.plot(pres_array,ASR[2,:])
plt.xlabel('System Pressure [Pa]')
plt.ylabel('Area Specific Resistance [ohm/m^2]')
plt.title('ASR vs Pressure for Various Current')
plt.legend(('0.1 A/cm^2','0.5 A/cm^2','1 A/m^2'))
plt.show()

# plt.figure(1)
# plt.plot(pres_array,C_H2_CL,'.')
# plt.plot(pres_array,C_H2O_CL,'.')
# plt.xlabel('System Pressure')
# plt.ylabel(' Concentration in Anode')
# plt.show()
