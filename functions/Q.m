function input_i = Q(t,p,input_data) 
% External input given to inhibitory cells, used in arg_i in all_algebraic_variables()

    if p.Whiskerpad == 0 || p.Whiskerpad == 5 % Square pulse
        pulse_width_P = p.startpulse + p.lengthpulse;
        input_i = p.Qinput .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P));

    elseif p.Whiskerpad == 1 || p.Whiskerpad == 2 || p.Whiskerpad == 4
        index_t = round(t./p.dt + 1);
        input_i = p.Qinput .* input_data(index_t);

    elseif p.Whiskerpad == 3
        pulse_width_P = p.startpulse + p.lengthpulse;
        pulse_width_P_1 = p.secondpulse + p.secondlength;
        input_i = p.Qinput .* (heaviside(t - p.startpulse) - heaviside(t - pulse_width_P) + heaviside(t - p.secondpulse) - heaviside(t - pulse_width_P_1));
    end
            
end