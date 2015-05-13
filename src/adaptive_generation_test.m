
pop_w = 0;
pop_a = 1;
UI_level = 0.05;
% ini_popsize = 20;
for ini_i = 1:ini_NrPar
    ini_randMat(ini_i,:) = vertcat(1,(up_b(ini_i)-low_b(ini_i))*rand(ini_popsize-1,1)+low_b(ini_i));
end

while pop_w == 0;
    
    
    
    for pop_i = 1 : ini_NrPar
        pop_c = ((up_b(pop_i)-low_b(pop_i))/2);
        for pop_j = pop_a : ini_popsize
            pop_d = pop_c  - ini_randMat(pop_i,pop_j);
            pop_vector_1(pop_i,pop_j) = pop_d;
        end
        pop_p(pop_i) = abs(sum(pop_vector_1(pop_i,:))/ini_popsize);
        pop_p_1(pop_i) = (pop_p(pop_i) / abs(pop_c - up_b(pop_i)));
    end
    
    pop_res = (sum(pop_p_1) / ini_NrPar);
    
    if pop_res > UI_level
        pop_a = ini_popsize;
        ini_popsize = ini_popsize +  20;
            
        for ini_i = 1:ini_NrPar
            ini_randMat(ini_i,(pop_a+1):ini_popsize) = vertcat((up_b(ini_i)-low_b(ini_i))*rand((ini_popsize-pop_a),1)+low_b(ini_i));
        end
    else
        pop_w = 1;
    end
end