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
from scipy.integrate import solve_ivp #integration function for ODE system.
from matplotlib import pyplot as plt
import numpy as np


from scipy.integrate import solve_ivp #integration function for ODE system.
from SOFC_function_1D import residual # point the model to the residual function
from SOFC_inputs import param, SV_0, ptr
    
    
    
    
solution = solve_ivp(lambda t, y: residual(t, y, param, ptr),param.time_span, SV_0, rtol=1e-9, atol=1e-7, method='BDF')
    
    

for var in solution.y:
    plt.plot(solution.t,var)
    
plt.legend(['Anode double layer','Cathode double layer'])
    
    
# V_elyte = solution.y[0,:]
# V_ca = V_elyte + solution.y[1,:]


# plt.plot(solution.t,V_elyte)
# plt.plot(solution.t,V_ca)

# plt.xlabel('Time (s)',fontsize=14)
# plt.ylabel('Electric Potential (V)',fontsize=14)

# plt.legend([r'$\phi_{\rm elyte}$',r'$\phi_{\rm cathode}$'],fontsize=14,frameon=False)
   




