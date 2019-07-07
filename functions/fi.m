function input_fi = fi(p,x) 
% Sigmoid function
     input_fi = 11.61 ./ ( 1 + exp(-(x-15)./3.94) );
end