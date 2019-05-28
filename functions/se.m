function input_se = se(p,x) 
% Sigmoid function
     input_se = 1 ./ (1 + exp(-(p.a_e .* (x - p.theta_e)))) - 1 ./ (1 + exp(p.a_e .* p.theta_e));
end