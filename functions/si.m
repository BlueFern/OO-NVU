function input_si = si(p,x) 
% Sigmoid function
     input_si = 1 ./ (1 + exp(-(p.a_i .* (x - p.theta_i)))) - 1 ./ (1 + exp(p.a_i .* p.theta_i));
end