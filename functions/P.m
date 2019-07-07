function input_e = P(t,p,input_data) 
% External input given to excitatory cells, used in arg_e in all_algebraic_variables()

    if p.Whiskerpad == 0 || p.Whiskerpad == 5 % Square pulse
            pulse_width_P = p.startpulse + p.lengthpulse; 
            input_e = p.Pinput .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P));
%             input_e = 0.5*tanh((t - p.startpulse)/100) - 0.5*tanh((t - pulse_width_P)/100);
            
    elseif p.Whiskerpad == 1 || p.Whiskerpad == 2 || p.Whiskerpad == 4 % Use Zheng input data
            index_t = round(t./p.dt + 1);
            input_e = p.Pinput .* input_data(index_t);
            
    elseif p.Whiskerpad == 3 % Two square pulses
            pulse_width_P = p.startpulse + p.lengthpulse;
            pulse_width_P_1 = p.secondpulse + p.secondlength;
            input_e = p.Pinput .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P) + heaviside(t - p.secondpulse) - heaviside(t - pulse_width_P_1));
    end
     
end
