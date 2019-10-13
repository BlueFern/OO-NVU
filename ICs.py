# -*- coding: utf-8 -*-
"""
Initial conditions
"""
# Returns u0, an array containing the initial conditions for the state variables of the NVU model
# idx: dictionary containing the state variable names ('keys') and their corresponding index ('values')
def set_initial_conditions(idx):
    
    u0 = [None]*len(idx) # Initialise with correct size
    
    # Neuron
    u0[idx['E_t']] = 0          # Proportion of excitatory cells firing per unit time [-] 
    u0[idx['I_t']] = 0          # Proportion of inhibitory cells firing per unit time [-]
    u0[idx['K_e']] = 3.5        # Extracellular K+ concentration [mM]
    u0[idx['Na_sa']] = 9.37     # Na+ concentration in the soma/axon [mM]
    u0[idx['Na_d']] = 9.42      # Na+ concentration in the dendrite [mM]
    u0[idx['O2']] = 0.02566     # Tissue oxygen concentration [mM] 
    u0[idx['CBV']] = 0.794      # Nondimensional cerebral blood volume [-] 
    u0[idx['HbR']] = 0.7641     # Normalised deoxyhemoglobin concentration [-] 
    u0[idx['Ca_n']] = 0.1       # Neuronal Ca2+ concentration [uM]
    u0[idx['nNOS']] = 0.01056   # Neuronal nitric oxide synthase concentration [uM]
    u0[idx['NO_n']] = 0.02425   # Neuronal nitric oxide concentration [uM] 
    
    # Astrocyte
    u0[idx['Na_k']] = 18740     # Astrocytic Na+ concentration [uM]
    u0[idx['K_k']] = 92660      # Astrocytic K+ concentration [uM]  
    u0[idx['HCO3_k']] = 9085    # Astrocytic HCO3- concentration [uM]
    u0[idx['Cl_k']] = 8212      # Astrocytic Cl- concentration [uM] 
    u0[idx['Na_s']] = 149200    # Synaptic Na+ concentration [uM] 
    u0[idx['K_s']] = 2932       # Synaptic K+ concentration [uM] 
    u0[idx['HCO3_s']] = 16980   # Synaptic HCO3- concentration [uM] 
    u0[idx['K_p']] = 3039       # Perivascular K+ concentration [uM] 
    u0[idx['w_k']] = 8.26e-5    # Open probability of the astrocytic BK channel [-] 
    u0[idx['Ca_k']] = 0.1435    # Astrocytic Ca2+ concentration [uM] 
    u0[idx['s_k']] = 480.8      # Ca2+ concentration in the astrocytic ER [uM] 
    u0[idx['h_k']] = 0.4107     # Inactivation variable of the astrocytic IP3R channel [-] 
    u0[idx['I_k']] = 0.048299   # Astrocytic IP3 concentration [uM]
    u0[idx['eet_k']] = 0.4350   # Astrocytic EET concentration [uM] 
    u0[idx['m_k']] = 0.513      # Open probability of the astrocytic TRPV4 channel [-] 
    u0[idx['Ca_p']] = 1853      # Perivascular Ca2+ concentration [uM] 
    u0[idx['NO_k']] = 0.02234   # Astrocytic nitric oxide concentration [uM] 
    u0[idx['v_k']] = -88.79     # Astrocytic membrane potential [mV] 
    u0[idx['AA_k']] = 9.3       # Astrocytic arachidonic acid concentration [uM]
    
    # SMC
    u0[idx['Ca_i']] = 0.2641    # SMC Ca2+ concentration [uM]
    u0[idx['s_i']] = 1.1686     # Ca2+ concentration in the SR of the SMC [uM]
    u0[idx['v_i']] = -34.7      # SMC membrane potential [mV]
    u0[idx['w_i']] = 0.2206     # Open probability of the SMC BK channel [-]
    u0[idx['I_i']] = 0.275      # SMC IP3 concentration [uM]
    u0[idx['NO_i']] = 0.02047   # SMC nitric oxide concentration [uM] 
    u0[idx['E_b']] = 0.6372     # Fraction of sGC in the basal state [-]
    u0[idx['E_6c']] = 0.2606    # Fraction of sGC in the intermediate form [-] 
    u0[idx['cGMP_i']] = 6.1     # SMC cGMP concentration [uM]
    u0[idx['H_i']] = 0.069   # SMC 20-HETE concentration [uM] 0.069
    u0[idx['AA_i']] = 9.3       # SMC arachidonic acid concentration [uM]
    
    # EC
    u0[idx['Ca_j']] = 0.8339    # EC Ca2+ concentration [uM]
    u0[idx['s_j']] = 0.6262     # Ca2+ concentration in the ER of the EC [uM]
    u0[idx['v_j']] = -68.39     # EC membrane potential [mV]
    u0[idx['I_j']] = 0.825      # EC IP3 concentration [uM]
    u0[idx['eNOS']] = 0.4451    # Endothelial nitric oxide synthase concentration [uM]
    u0[idx['NO_j']] = 0.02051   # EC nitric oxide concentration [uM] 
    
    # Wall Mechanics
    u0[idx['Mp']] = 0.0842      # Fraction of free phosphorylated cross-bridges [-]
    u0[idx['AMp']] = 0.0622     # Fraction of attached phosphorylated cross-bridges [-]
    u0[idx['AM']] = 0.2746      # Fraction of attached dephosphorylated cross-bridges [-]
    u0[idx['R']] = 22.44        # Arteriolar radius [um]
    
    # GABA, Glu and NPY
    u0[idx['GABA']] = 0         # Normalised GABA concentration [-]
    u0[idx['Glu']] = 0          # Normalised glutamate concentration [-]
    u0[idx['NPY']] = 0          # Normalised NPY concentration [-]
    return u0