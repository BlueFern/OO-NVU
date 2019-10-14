# OO-NVU - Version 3.2
Python implementation of the NVU model.

Contains the following pathways:
* Excitatory/inhibitory neuron dynamics based on Wilson & Cowan (1972)
* GABA and NPY release from the neuron (inhibitory only) 
* K+ pathway via K+ release from the neuron into the SC
* Nitric oxide pathway via glutamate release from the neuron into the SC
* Astrocytic calcium pathway via glutamate release from the neuron into the SC
* TRPV4 calcium channel on the astrocytic endfoot
* 20-HETE and Arachidonic Acid dynamics in the astroycte and SMC

How to Run 
----
The code is run in Python. All functions are commented. If using Spyder it is recommended that you set the following preferences:
* To have plots open in new windows rather than in the console so they can be edited/saved, go to `Preferences->IPython Console->Graphics` and set `Graphics Backend` to `Automatic`
* To show all variables in Variable explorer (including figures and data objects), go to `Preferences->Variable Explorer` and turn off `Exclude unsupported data types`

The code is laid out as follows:
* To solve the system and plot variables use the script `run.py`
* To change any model parameters see the file `parameters.py` which is loaded as a module in `run.py` (via `import parameters as p`)
* To adjust any algebraic variables see the file `algebraic_variables.py`
* To adjust any differential equations see the file `ODEsystem.py`
* To adjust model functions (e.g. time dependent input functions like P(t) or Q(t)) see the file `model_functions.py`
* To adjust initial conditions see the file `ICs.py`
* When adding any new state variables to `ODEsystem.py` these must also be added to `indices.py` and `state_variable_names.py`
* The normalised hemodynamics variables are defined in `normalised_hemo.py` (solved separately from `ODEsystem.py` as they depend on the steady state conditions) 
* Two custom plotting functions are defined in `plotting_functions.py` for plotting either state or algebraic variables either on the same graph or in subplots; examples of use are shown in `run.py`
* Experimental data stored as mat files in the subfolder `./mat_files` can be extracted and plotted using custom functions defined in `import_mat_files.py`; examples of use are shown in `run.py`

Note that time is in milliseconds.

For further information on the NVU model refer to the following papers:
* Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory interactions in localized populations of model neurons. Biophysical Journal, 12, 1–24. https://doi.org/10.1016/S0006-3495(72)86068-5
* Mathias, E., Kenny, A., Plank, M. J., & David, T. (2018). Integrated models of neurovascular coupling and BOLD signals: Responses for varying neural activations. NeuroImage, 174(March), 69–86. http://doi.org/10.1016/j.neuroimage.2018.03.010
* Kenny, A., Plank, M. J., & David, T. (2018). The role of astrocytic calcium and TRPV4 channels in neurovascular coupling. Journal of Computational Neuroscience, 44(1), 97–114. http://doi.org/10.1007/s10827-017-0671-7
* Dormanns, K., Brown, R. G. G., & David, T. (2016). The role of nitric oxide in neurovascular coupling. Journal of Theoretical Biology, 394, 1–17. http://doi.org/10.1016/j.jtbi.2016.01.009
* Dormanns, K., van Disseldorp, E. M. J., Brown, R. G., & David, T. (2015). Neurovascular coupling and the influence of luminal agonists via the endothelium. Journal of Theoretical Biology, 364, 49–70. http://doi.org/10.1016/j.jtbi.2014.08.029
* Farr, H., & David, T. (2011). Models of neurovascular coupling via potassium and EET signalling. Journal of Theoretical Biology, 286(1), 13–23. http://doi.org/10.1016/j.jtbi.2011.07.006


