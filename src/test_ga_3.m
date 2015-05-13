for ini_k = 1:ini_popsize
    
    if ui_choice == 1
        for ini_i = 1:ini_NrPar
            nv.smcec.params.(ui_param{ini_i+1})= ini_value(ini_i) * ini_randMat(ini_i,ini_k);
        end
    end
    
    if ui_choice == 2
        for ini_i = 1:ini_NrPar
            nv.wall.params.(ui_param{ini_i+1})= ini_value(ini_i) * ini_randMat(ini_i,ini_k);
        end
    end
    
    if ui_choice == 3
        for ini_i = 1:ini_NrPar
            nv.astrocyte.params.(ui_param{ini_i+1})= ini_value(ini_i) * ini_randMat(ini_i,ini_k);
        end
    end
    
    nv.simulate()
    
    for ini_l = 1 : ini_NrFitVal
        results_1 (ini_l,:)  = nv.out(ini_FitString{ini_l});
        results(ini_k,ini_l) = results_1 (ini_FitTime(ini_l));              %?? does it work ini_FitTime
    end
end

Max_matrix = [results(:,:),ini_randMat'];

for ini_k = 1 : ini_NrFitVal
    Results_end_max(r_simu_2,ini_k)  = max(results(:,ini_k));                %essential for plotting results
    Results_end_mean(r_simu_2,ini_k) = mean(results(:,ini_k));               %essential for plotting results
end

%% Selection (roulette wheel method)-> sel_...
%probability of each parameterset to be choosen for the offspring
%population -> probability for parameter sets with a high
%calciumconcentration is higher

%probabilities

sel_sum_prob(:,:) = zeros(ini_popsize,ini_NrFitVal);
for sel_j = 1 : ini_NrFitVal                                                %??
    
    for sel_1 = 1:ini_popsize                                               %ini_popsize= size of the population
        sel_prob(sel_1,sel_j)       = abs((abs(ui_fitness{sel_j} - Max_matrix(sel_1,sel_j)))/(ui_fitness{sel_j}));                                  %probability for each individual to be chosen
    end                                                                     %saving the probabilities
    delete_1(:,sel_j) = sel_prob(:,sel_j) ./ sum(sel_prob(:,sel_j));
    sel_sum_prob(:,1) = sel_sum_prob(:,1) + sel_prob(:,sel_j);
    
end

for sel_1 = 1 : ini_popsize 
delete_2 (:,1) = sum(delete_1(sel_1,:)) ./ ini_NrFitVal
end

sel_cum_probvec(:,:) = cumsum((sel_sum_prob(:,1)) ./ (ini_popsize * ini_NrFitVal));

sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)',sel_cum_probvec, Max_matrix]; %numbered max_matrix with probabilities

% decision
sel_output      = zeros(ini_popsize/2 , 2 + ini_NrFitVal + ini_NrPar);   %creates a zero matrix plus 1 row sel_cumprobvec + 1 linspace + # of fitness values

sel_rand_value  = rand(length(sel_output(:,1)),1);                          %random vector with the length of the number of parents and with random numbers between 0-1
for sel_2=1:length(sel_output(:,1))                                         %for the length of the number of the parents, that should be chosen
    
    while sel_output(sel_2,1) == 0                                          %till the row of the sel_output rows does include a individual
        sel_p=1;
        if sel_rand_value(sel_2,1) < sel_max_matrix(sel_p,2);               % if random value at the position n is smaler as the cummultative probability ...
            sel_output(sel_2,:) = sel_max_matrix(sel_p,:);                  % ... save the individum in the sel_output matrix
        end
        
        while sel_rand_value(sel_2,1) > sel_max_matrix(sel_p,2);            % while the random number is bigger ....
            sel_p=sel_p+1;
            if sel_rand_value(sel_2,1) <= sel_max_matrix(sel_p,2);          % test again if the cummultative probability of the next individium (parameter set) is bigger than the random number
                sel_output(sel_2,:) = sel_max_matrix(sel_p,:);              % ... save the individum in the sel_output matrix
            end
        end
        
        
        if sel_2>1                                                                      % prevent, that one individum (parameterset) is taken twice
            sel_duplex = ismember(sel_output(sel_2,1),sel_output(1:(sel_2-1),1));       % proves if the chosen individum is already in the sel_outtput matrix
            if sel_duplex == 1                                                          % if yes ...
                sel_output(sel_2,1)     = 0;                                            % delet it in the sel_otuput matrix
                sel_rand_value(sel_2,1) = rand(1);                                      % create a new random number at this place and repeat all the steps
            end
        end
        
    end
end