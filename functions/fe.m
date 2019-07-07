function input_fe = fe(p,x) 
% Sigmoid function
     input_fe = 5.12 ./ ( 1 + exp(-(x-15)./4.16) );
end