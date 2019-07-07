function input_NPY = NPY(t,p) 
    % Normalised NPY concentration, normalised wrt some maximum value NPY_max (not currently known) i.e. NPY/NPY_max
    
    NPY_input = 1; % Max value for normalised NPY, i.e. NPY_max/NPY_max = 1
    pulse_width_P = p.startpulse + p.lengthpulse; 
    input_NPY = NPY_input .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P));

end
