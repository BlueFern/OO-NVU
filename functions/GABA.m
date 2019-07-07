function input_GABA = GABA(t,p) 
    % Normalised GABA concentration, normalised wrt some maximum value GABA_max (not currently known) i.e. GABA/GABA_max
    
    GABA_input = 1; % Max value for normalised GABA, i.e. GABA_max/GABA_max = 1
    pulse_width_P = p.startpulse + p.lengthpulse; 
    input_GABA = GABA_input .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P));

end
