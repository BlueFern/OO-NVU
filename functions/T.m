function input_t = T(t,p,input_data) 
% Thalamic input given to neurons, used in arg_e and arg_i in all_algebraic_variables()
% Taken from input_data which is generated in run_NVU_DE_system.m
        index_t = round(t./p.dt + 1);
        input_t = p.Tinput .* input_data(index_t);     
end