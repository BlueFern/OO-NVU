function gamma = Gamma(t, p, E_t, I_t)  
% Component of dKe/dt that is dependent on E_t - I_t. Creates a linear decay from max 
% E-I down to 0 rather than the sharp jump that occurs when Gamma is not used. 
% Set the parameter p.gamma_switch = 0 if not wanted
            
            t_end = p.startpulse + p.lengthpulse + p.gamma_delay; 
            endgamma = t_end + p.gamma_width;
            
            if      t < t_end
                gamma = p.alpha_K_e * p.beta_K_e .* (abs(E_t - I_t) ./ p.EI_relative);
            elseif  t > endgamma
                gamma = p.alpha_K_e * p.beta_K_e .* (abs(E_t - I_t) ./ p.EI_relative);
            else
                gamma = (p.alpha_K_e .* p.beta_K_e) .* (endgamma - t)./(endgamma - t_end); % Decreasing from max (E-I)/EI_relative value (1) to zero at endgamma
%                 gamma = (alpha_K_e .* beta_K_e) .* (0.5 - 0.5*tanh((t - pulse_width_P - 2400)/1400) ); % Sigmoidal
            end
            
end