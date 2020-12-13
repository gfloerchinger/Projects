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
    

# for var in solution.y
#     plt.plot(solution.t,var)
    
# plt.legend(['Anode double layer','Cathode double layer'])
    
    
i_array = np.linspace(0.00000001,100000,50)
V_cell = np.zeros_like(i_array)

for j, current in enumerate(i_array):
    print(current)
    solution = SOFC_model(current)
    V_cell[j] = solution.y[1,-1] - solution.y[0,-1]
    print(V_cell[j])

plt.plot(i_array,V_cell,'.')
plt.savefig('results/polarization.png',dpi=350)
plt.show()

plt.plot(i_array,V_cell,'.')

plt.show()