# OO-NVU - Version 3
Matlab implementation of the NVU model.

Contains the following pathways:
* Excitatory/inhibitory neuron dynamics based on Wilson & Cowan (1972)
* K+ pathway via K+ release from the neuron into the SC
* Nitric oxide pathway via glutamate release from the neuron into the SC
* Astrocytic calcium pathway via glutamate release from the neuron into the SC
* TRPV4 calcium channel on the astrocytic endfoot
* 20-HETE and Arachidonic Acid dynamics in the astroycte and SMC

How to Run 
----
The code is run in Matlab:
* To solve the system use the script `run_NVU_DE_system.m`
* To change any model parameters see the file `all_parameters.m`
* To adjust any algebraic variables see the file `all_algebraic_variables.m`
* To adjust any differential equations see the file `NVU_DE_system.m`
* To adjust initial conditions see the file `initial_conditions.m`
* When adding any new state variables to `NVU_DE_system.m` these must also be added to `all_indices` and `set_variable_names.m`

Note that time is in milliseconds.

For further information on the NVU model refer to the following papers:
* Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory interactions in localized populations of model neurons. Biophysical Journal, 12, 1–24.
* Mathias, E., Kenny, A., Plank, M. J., & David, T. (2018). Integrated models of neurovascular coupling and BOLD signals: Responses for varying neural activations. NeuroImage, 174(March), 69–86. http://doi.org/10.1016/j.neuroimage.2018.03.010
* Kenny, A., Plank, M. J., & David, T. (2018). The role of astrocytic calcium and TRPV4 channels in neurovascular coupling. Journal of Computational Neuroscience, 44(1), 97–114. http://doi.org/10.1007/s10827-017-0671-7
* Dormanns, K., Brown, R. G. G., & David, T. (2016). The role of nitric oxide in neurovascular coupling. Journal of Theoretical Biology, 394, 1–17. http://doi.org/10.1016/j.jtbi.2016.01.009
* Dormanns, K., van Disseldorp, E. M. J., Brown, R. G., & David, T. (2015). Neurovascular coupling and the influence of luminal agonists via the endothelium. Journal of Theoretical Biology, 364, 49–70. http://doi.org/10.1016/j.jtbi.2014.08.029
* Farr, H., & David, T. (2011). Models of neurovascular coupling via potassium and EET signalling. Journal of Theoretical Biology, 286(1), 13–23. http://doi.org/10.1016/j.jtbi.2011.07.006

