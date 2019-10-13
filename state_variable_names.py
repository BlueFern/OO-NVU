# -*- coding: utf-8 -*-
"""
Initialise the state variables
"""

# Define a class that will contain attributes = state variables
class state_vars():   
    pass

# Returns a class 'v' containing all state variables
# u: 2D array of state variables vs time
# idx: dictionary containing the state variable names ('keys') and their corresponding index ('values')
# location: for use inside a function ('in function') or an array for use outside of a function ('out function')
def initialise_var_names(u, idx, location):
    
    # Initialise class to contain the variables
    v = state_vars()
    
    if location == 'in function':
        
        # Neuron
        v.E_t =   u[idx['E_t']]         
        v.I_t =   u[idx['I_t']]            
        v.K_e =   u[idx['K_e']]             
        v.Na_sa = u[idx['Na_sa']]          
        v.Na_d =  u[idx['Na_d']]          
        v.O2 =    u[idx['O2']]        
        v.CBV =   u[idx['CBV']]           
        v.HbR =   u[idx['HbR']]            
        v.Ca_n =  u[idx['Ca_n']]                
        v.nNOS =  u[idx['nNOS']]  
        v.NO_n =  u[idx['NO_n']]           
        
        # Astrocyte
        v.v_k = u[idx['v_k']]
        v.K_p = u[idx['K_p']]
        v.Ca_p = u[idx['Ca_p']]
        v.Na_k = u[idx['Na_k']]
        v.K_k = u[idx['K_k']]
        v.Cl_k = u[idx['Cl_k']]
        v.HCO3_k = u[idx['HCO3_k']]
        v.Na_s = u[idx['Na_s']]
        v.K_s = u[idx['K_s']]
        v.HCO3_s = u[idx['HCO3_s']]
        v.w_k = u[idx['w_k']]
        v.I_k = u[idx['I_k']]
        v.Ca_k = u[idx['Ca_k']]
        v.h_k = u[idx['h_k']]
        v.s_k = u[idx['s_k']]
        v.m_k = u[idx['m_k']] 
        v.eet_k = u[idx['eet_k']]
        v.NO_k = u[idx['NO_k']]
        v.AA_k = u[idx['AA_k']]
        
        # SMC
        v.Ca_i = u[idx['Ca_i']]
        v.s_i = u[idx['s_i']]
        v.v_i = u[idx['v_i']]
        v.w_i = u[idx['w_i']]
        v.I_i = u[idx['I_i']]
        v.NO_i = u[idx['NO_i']]
        v.E_b = u[idx['E_b']]
        v.E_6c = u[idx['E_6c']]
        v.cGMP_i = u[idx['cGMP_i']]
        v.H_i = u[idx['H_i']]
        v.AA_i = u[idx['AA_i']]
        
        # EC
        v.Ca_j = u[idx['Ca_j']]
        v.s_j = u[idx['s_j']]
        v.v_j = u[idx['v_j']]
        v.I_j = u[idx['I_j']]
        v.eNOS = u[idx['eNOS']]
        v.NO_j = u[idx['NO_j']]
        
        # Wall Mechanics
        v.Mp = u[idx['Mp']]
        v.AMp = u[idx['AMp']]
        v.AM = u[idx['AM']]
        v.R = u[idx['R']]
        
        # GABA, NPY and Glu
        v.GABA = u[idx['GABA']]
        v.NPY = u[idx['NPY']]
        v.Glu = u[idx['Glu']]
        
    elif location == 'out function':
        
        # Neuron
        v.E_t =   u[:, idx['E_t']]            
        v.I_t =   u[:, idx['I_t']]             
        v.K_e =   u[:, idx['K_e']]            
        v.Na_sa = u[:, idx['Na_sa']]            
        v.Na_d =  u[:, idx['Na_d']]            
        v.O2 =    u[:, idx['O2']]              
        v.CBV =   u[:, idx['CBV']]            
        v.HbR =   u[:, idx['HbR']]      
        v.Ca_n =  u[:, idx['Ca_n']]            
        v.nNOS = u[:, idx['nNOS']] 
        v.NO_n =  u[:, idx['NO_n']]      
        
        # Astrocyte
        v.v_k = u[:, idx['v_k']]
        v.K_p = u[:, idx['K_p']]
        v.Ca_p = u[:, idx['Ca_p']]
        v.Na_k = u[:, idx['Na_k']]
        v.K_k = u[:, idx['K_k']]
        v.Cl_k = u[:, idx['Cl_k']]
        v.HCO3_k = u[:, idx['HCO3_k']]
        v.Na_s = u[:, idx['Na_s']]
        v.K_s = u[:, idx['K_s']]
        v.HCO3_s = u[:, idx['HCO3_s']]
        v.w_k = u[:, idx['w_k']]
        v.I_k = u[:, idx['I_k']]
        v.Ca_k = u[:, idx['Ca_k']]
        v.h_k = u[:, idx['h_k']]
        v.s_k = u[:, idx['s_k']]
        v.m_k = u[:, idx['m_k']] 
        v.eet_k = u[:, idx['eet_k']]
        v.NO_k = u[:, idx['NO_k']]
        v.AA_k = u[:, idx['AA_k']]
        
        # SMC
        v.Ca_i = u[:, idx['Ca_i']]
        v.s_i = u[:, idx['s_i']]
        v.v_i = u[:, idx['v_i']]
        v.w_i = u[:, idx['w_i']]
        v.I_i = u[:, idx['I_i']]
        v.NO_i = u[:, idx['NO_i']]
        v.E_b = u[:, idx['E_b']]
        v.E_6c = u[:, idx['E_6c']]
        v.cGMP_i = u[:, idx['cGMP_i']]
        v.H_i = u[:, idx['H_i']]
        v.AA_i = u[:, idx['AA_i']]
        
        # EC
        v.Ca_j = u[:, idx['Ca_j']]
        v.s_j = u[:, idx['s_j']]
        v.v_j = u[:, idx['v_j']]
        v.I_j = u[:, idx['I_j']]
        v.eNOS = u[:, idx['eNOS']]
        v.NO_j = u[:, idx['NO_j']]
        
        # Wall Mechanics
        v.Mp = u[:, idx['Mp']]
        v.AMp = u[:, idx['AMp']]
        v.AM = u[:, idx['AM']]
        v.R = u[:, idx['R']]
        
        # GABA, NPY and Glu
        v.GABA = u[:, idx['GABA']]
        v.NPY = u[:, idx['NPY']]
        v.Glu = u[:, idx['Glu']]
        
    else:
        print("initialise_var_names needs a location!")
            
    return v