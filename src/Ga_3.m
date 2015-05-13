%% Genetic Algorithm
%Initiation of a random population

clear; clc; close all
%% Options
% First the options need to be set for the ODE solver (currently |ode15s|)
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC(), ...
    'odeopts', odeopts);
%% New user input

% Fitness values (Strings) and number of fitness values
ui_promp_FitString = {'How many fitness values should be optimised(min:1 max:20):'
    'Fitness value 1:'
    'Fitness value 2:'
    'Fitness value 3:'
    'Fitness value 4:'
    'Fitness value 5:'
    'Fitness value 6:'
    'Fitness value 7:'
    'Fitness value 8:'
    'Fitness value 9:'
    'Fitness value 10:'
    'Fitness value 11:'
    'Fitness value 12:'
    'Fitness value 13:'
    'Fitness value 14:'
    'Fitness value 15:'
    'Fitness value 16:'
    };

ui_title_FitString = 'What should be optimised?';
ui_lines_FitString = 1;
ui_defaultanswer_FitString  = {'3','c_k','Ca_i','F_r','0','0','0','0','0','0','0','0','0','0','0','0','0'};
ui_FitString= inputdlg(ui_promp_FitString,ui_title_FitString,ui_lines_FitString,ui_defaultanswer_FitString);
ini_NrFitVal  = str2double(ui_FitString{1});
ini_FitString = ui_FitString(2:17);



% Fitness values to the fitness strings
for   ui_i_FitVal = 1 : ini_NrFitVal
    ui_promp_FitValSave = {'Enter value after which you want to optimize for' ini_FitString{ui_i_FitVal}};
    ui_prompt_vec_FitVal{ui_i_FitVal} = strjoin(ui_promp_FitValSave);
end

ui_title_FitVal = 'Add fitness values to the fitness strings';
ui_lines_FitVal = 1;

for ui_i_FitVal = 1 : ini_NrFitVal
    ui_defaultanswer_FitVal(ui_i_FitVal) = {'0.45'};
end
ui_save_FitVal   = ui_prompt_vec_FitVal(1:ini_NrFitVal);
ui_first_FitVal = inputdlg(ui_save_FitVal,ui_title_FitVal,ui_lines_FitVal,ui_defaultanswer_FitVal);

% Timestep to each fitness value
for   ui_i_FitTime = 1 : ini_NrFitVal
    ui_promp_FitTimeSave = {'Enter the timestep for the results for' ini_FitString{ui_i_FitTime}};
    ui_prompt_vec_FitTime{ui_i_FitTime} = strjoin(ui_promp_FitTimeSave);
end

ui_title_FitTime = 'Add timestep for you fitness results';
ui_lines_FitTime = 1;

for ui_i_FitTime = 1 : ini_NrFitVal
    ui_defaultanswer_FitTime(ui_i_FitTime) = {'740'};
end
ui_save_FitTime   = ui_prompt_vec_FitTime(1:ini_NrFitVal);
ui_first_FitTime = inputdlg(ui_save_FitTime,ui_title_FitTime,ui_lines_FitTime,ui_defaultanswer_FitTime);

%FitTime and FitValue saved as double
for ui_j_FitTime =1:length(ui_first_FitTime)
    ini_FitTime(ui_j_FitTime) = str2double(ui_first_FitTime{ui_j_FitTime});
    ini_FitVal(ui_j_FitTime) = str2double(ui_first_FitVal{ui_j_FitTime});
end



%Paramters
%Astrocyte
ui_prompt_param_A = { 'How many parameter do you want to change in the Astrocyte?'
    'Parameter1:'
    'Parameter2:'
    'Parameter3:'
    'Parameter4:'
    'Parameter5:'
    'Parameter6:'
    'Parameter7:'
    'Parameter8:'
    'Parameter9:'
    'Parameter10:'
    'Parameter11:'
    'Parameter12:'
    'Parameter13:'
    'Parameter14:'
    'Parameter15:'
    'Parameter16:'
    };                                                                      % At which timepoint the user wants to know the calcium concentration

ui_title_param_A = 'Which parameters do you want to use?';
ui_lines_param_A = 1;
ui_defaultanswer_param_A = {'16','k_deg','k_on','K_inh','r_h','K_G','K_I','BK_end','K_ex','B_ex','J_max','K_act','P_L','V_max','k_pump','delta','VR_ER_cyt'};
ui_param_A = inputdlg(ui_prompt_param_A,ui_title_param_A,ui_lines_param_A,ui_defaultanswer_param_A);
ini_NrPar(1)= str2double(ui_param_A{1});
ini_ParamString = ui_param_A(2:sum(ini_NrPar)+1);

%SMCEC
ui_prompt_param_S = { 'How many parameter do you want to change in the SMCEC?'
    'Parameter1:'
    'Parameter2:'
    'Parameter3:'
    'Parameter4:'
    'Parameter5:'
    'Parameter6:'
    'Parameter7:'
    'Parameter8:'
    'Parameter9:'
    'Parameter10:'
    'Parameter11:'
    'Parameter12:'
    'Parameter13:'
    'Parameter14:'
    'Parameter15:'
    'Parameter16:'
    
    };                                                                      % At which timepoint the user wants to know the calcium concentration

ui_title_param_S = 'Which parameters do you want to use?';
ui_lines_param_S = 1;
ui_defaultanswer_param_S = {'5','gamma_i','lambda_i','C_m_j','J_PLC','J_0_j','0','0','0','0','0','0','0','0','0','0','0'};
ui_param_S = inputdlg(ui_prompt_param_S,ui_title_param_S,ui_lines_param_S,ui_defaultanswer_param_S);
ini_NrPar(2) = str2double(ui_param_S{1});
ini_ParamString((ini_NrPar(1)+1):sum(ini_NrPar)) = ui_param_S(2:ini_NrPar(2)+1);


%WallMechanics
ui_prompt_param_W = { 'How many parameter do you want to change in the WallMechanics?'
    'Parameter1:'
    'Parameter2:'
    'Parameter3:'
    'Parameter4:'
    'Parameter5:'
    'Parameter6:'
    'Parameter7:'
    'Parameter8:'
    'Parameter9:'
    'Parameter10:'
    'Parameter11:'
    'Parameter12:'
    'Parameter13:'
    'Parameter14:'
    'Parameter15:'
    'Parameter16:'
    };

ui_title_param_W = 'Which parameters do you want to use?';
ui_lines_param_W = 1;
ui_defaultanswer_param_W  = {'7','K_2','K_3','K_4','K_5','K_7','gamma_cross','R_0_passive','0','0','0','0','0','0','0','0','0'};
ui_param_W = inputdlg(ui_prompt_param_W,ui_title_param_W,ui_lines_param_W,ui_defaultanswer_param_W);
ini_NrPar(3) = str2double(ui_param_W{1});
ini_ParamString((ini_NrPar(1)+ini_NrPar(2)+1):sum(ini_NrPar)) = ui_param_W(2:ini_NrPar(3)+1);


%individual boundaries for the different parameters
%Check if user wants individual boundaries

ui_check_1 = {'Do you want individual parameter boundaries (y=yes/n=no)?'};
ui_title_check = 'In which percentage boundaries should the parameter change?';
ui_lines_check = 1;
ui_defaultanswer_check = {'n'};
ui_check_dlg = inputdlg(ui_check_1,ui_title_check,ui_lines_check,ui_defaultanswer_check);
ui_answer            = ui_check_dlg{1};
ui_yes               = 'y';
ui_check = strcmp(ui_answer,ui_yes);

if ui_check == 1
    
    for ui_i_bound = 1 : sum(ini_NrPar)
        ui_A = {'Enter lower boundary of randomisation of', ini_ParamString{ui_i_bound}};
        ui_B = {'Enter upper boundary of randomisation of', ini_ParamString{ui_i_bound}};
        ui_i_bound_vec{ui_i_bound*2-1} = strjoin(ui_A);
        ui_i_bound_vec{ui_i_bound*2}   = strjoin(ui_B);
    end
    
    
    ui_title_bound = 'Individual boundaries for the Parameter';
    ui_lines_bound        = 1;
    for ui_i_bound      = 1 : sum(ini_NrPar)
        ui_defaultanswer_bound(ui_i_bound*2-1) = {'0.8'};
        ui_defaultanswer_bound(ui_i_bound*2)   = {'1.2'};
    end
    
    ui_i_bound_1      = ui_i_bound_vec(1:sum(ini_NrPar));
    ui_bound_1        = inputdlg(ui_i_bound_1,ui_title_bound,ui_lines_bound,ui_defaultanswer_bound,'on');
    
    ui_i_bound_2         = ui_i_bound_vec(sum(ini_NrPar)+1 : sum(ini_NrPar)*2);
    ui_title_bound_2     = 'Individual boundaries for the parameters';
    ui_lines_bound_2     = 1;
    ui_defaultanswer_bound_2 = ui_defaultanswer_bound;
    
    ui_bound_2 = inputdlg(ui_i_bound_2,ui_title_bound_2,ui_lines_bound_2,ui_defaultanswer_bound_2,'on');
    
    for ui_j_bound = 1 : sum(ini_NrPar)/2
        
        ini_low_b(ui_j_bound)                              = str2double(ui_bound_1{ui_j_bound * 2 -1});
        ini_up_b(ui_j_bound)                               = str2double(ui_bound_1{ui_j_bound * 2});
        ini_low_b(sum(ini_NrPar)/2 + ui_j_bound)           = str2double(ui_bound_2{ui_j_bound * 2 -1});
        ini_up_b(sum(ini_NrPar)/2 + ui_j_bound)            = str2double(ui_bound_2{ui_j_bound * 2});
        
    end
    ui_j_bound = 1;
else                                                               % IF user doesn't want to change default boundaries
    
    ui_prompt_samebound = {'Enter lower boundary of randomisation for all parameters:'
        'Enter upper boundary of randomisation for all parameters:'};
    
    ui_title_samebound     = 'Individual boundaries for the parameters';
    
    ui_lines_samebound     = 1;
    ui_defaultanswer_samebound = { '0.8','1.2'};
    
    ui_samebound = inputdlg(ui_prompt_samebound,ui_title_samebound,ui_lines_samebound,ui_defaultanswer_samebound,'on');
    
    for ui_j_bound = 1 : sum(ini_NrPar)
        
        ini_low_b(ui_j_bound)           = str2double(ui_samebound{1});
        ini_up_b (ui_j_bound)           = str2double(ui_samebound{2});
        
    end
    ui_j_bound = 1;
end

%General User Input
ui_prompt_general = { 'Populationsize: With which probability should the generation be ultimate distributed:'                 % If the input is 0.95 the first generation is to 95% ultimate distributed
    'Cross-Over: Enter the initial probability of cross-over:'                                              % This value determines the frequency of co
    'Mutation: Enter the probability of Mutation'                                                           % Mutation probability
    'Elitism: What percentage of members of the first generation should survive without any genetic change' % 0.005 means that a organism must have a bigger c_k than 99.5% of population to be an elite
    'Stopping criterion: Enter the maximum simulation time [s]'                                                                  % Max. number of runs / generation
    'Stopping criterion: Enter the convergenz criteria number (How many generations should not have a change):'                                                                 % How many "close" values have to be in row, that he code stops
    'Stopping criterion: Enter the boundaries for convergenz:'                                                                   % Decision criteria if 2 values are close together or not
    'Stopping criterion: Do you want different accuracies for the different fitness values'};

ui_title_general = 'General User Input';
ui_lines_general = 1;
ui_defaultanswer_general = {'0.9','0.8','0.5','0.1','10','15','0.02','y'};                 % Just a default answer
ui_general = inputdlg(ui_prompt_general,ui_title_general,ui_lines_general,ui_defaultanswer_general,'on');


ini_ulti_dist          = 1-str2double(ui_general{1});
ini_co_p               = str2double(ui_general{2});
ini_mu_p               = str2double(ui_general{3});
ini_el_k               = str2double(ui_general{4});
ini_max_time           = str2double(ui_general{5});
ini_nr_convergence     = str2double(ui_general{6});
ini_bound_convergence  = str2double(ui_general{7});
ui_answer_accuracy     = ui_general{8};
ui_yes_accuracy        = 'y';
ui_check_accuracy = strcmp(ui_answer_accuracy,ui_yes_accuracy);

if ui_check_accuracy == 1
    
    for ui_i_accuracy = 1 : ini_NrFitVal
        ui_promp_accuracy = {'Enter the accuarcy for the fitness value' ini_FitString{ui_i_accuracy}};
        ui_prompt_vec_accuracy{ui_i_accuracy} = strjoin(ui_promp_accuracy);
    end
    
    ui_title_accuracy = 'Add the accuracy for every fitness value to reach the aimed value';
    ui_lines_accuracy = 1;
    
    for ui_i_accuracy = 1 : ini_NrFitVal
        ui_defaultanswer_accuracy(ui_i_accuracy) = {'0.1'};
    end
    ui_save_accuracy   = ui_prompt_vec_accuracy(1:ini_NrFitVal);
    ui_first_accuracy = inputdlg(ui_save_accuracy,ui_title_accuracy,ui_lines_accuracy,ui_defaultanswer_accuracy);
    
    for ui_i_accuracy =1:length(ui_first_FitTime)
        ini_accuracy(ui_i_accuracy) = str2double(ui_first_accuracy{ui_i_accuracy});
    end
    
else                                                               % IF user doesn't want to change default boundaries
    
    ui_prompt_sameaccuracy = {'Enter the accuracy for all parameters:'};
    ui_title_sameaccuracy     = 'Add the accuracy for all fitness values to reach the aimed value' ;
    ui_lines_sameaccuracy     = 1;
    ui_defaultanswer_sameaccuracy = { '0.8'};
    ui_sameaccuracy = inputdlg(ui_prompt_sameaccuracy,ui_title_sameaccuracy,ui_lines_sameaccuracy,ui_defaultanswer_sameaccuracy,'on');
    
    for ui_i_accuracy = 1 : ini_NrFitVal
        ini_accuracy(ui_i_accuracy)           = str2double(ui_sameaccuracy{1});
    end
end



%% initial values

for r_i= 1:ini_NrPar(1)
    ini_a            = nv.astrocyte.params;
    ini_value(r_i) = getfield(ini_a, ini_ParamString{r_i});
end
for r_i= 1:ini_NrPar(2)
    ini_a            = nv.smcec.params;
    ini_value(r_i+ini_NrPar(1)) = getfield(ini_a, ini_ParamString{r_i+ini_NrPar(1)});
end
for r_i= 1:ini_NrPar(3)
    ini_a            = nv.wall.params;
    ini_value(r_i+ini_NrPar(1)+ini_NrPar(2)) = getfield(ini_a, ini_ParamString{r_i+ini_NrPar(1)+ini_NrPar(2)});
end

% Initial population

T           = linspace (0,500,1000);                                                   % opens up the timevector
% # of different parameters
% Ratio between inital population and selected population

%running variables

mpar_1      = 1;
mpar_c      = 1;

% Running variable in the parenting choosing
mu_var      = 1;


count_max   = 1;
conv_crit   = 1;
r_simu_2    = 1;
ini_w       = 1;


%% First Generation -> Initiation of a random population

pop_w = 0;
pop_a = 1;


ini_popsize = 2^4;
for ini_i = 1:sum(ini_NrPar)
    ini_randMat(ini_i,:) = vertcat(1,(ini_up_b(ini_i)-ini_low_b(ini_i))*rand(ini_popsize-1,1)+ini_low_b(ini_i));
end

while pop_w == 0;
    
    for pop_i = 1 : sum(ini_NrPar)
        pop_c       = ((ini_up_b(pop_i)+ini_low_b(pop_i))/2);
        for pop_j   = pop_a : ini_popsize
            pop_d   = pop_c  - ini_randMat(pop_i,pop_j);
            pop_vector_1(pop_i,pop_j) = pop_d;
        end
        pop_p(pop_i)    = abs(sum(pop_vector_1(pop_i,:))/ini_popsize);
        pop_p_1(pop_i)  = (pop_p(pop_i) / abs(pop_c - ini_up_b(pop_i)));
    end
    
    pop_res = (sum(pop_p_1) / sum(ini_NrPar));
    
    if pop_res > ini_ulti_dist
        pop_a       = ini_popsize;
        ini_popsize = ini_popsize*2;
        
        for ini_i = 1:sum(ini_NrPar)
            ini_randMat(ini_i,(pop_a+1):ini_popsize) = vertcat((ini_up_b(ini_i)-ini_low_b(ini_i))*rand((ini_popsize-pop_a),1)+ini_low_b(ini_i));
        end
    else
        pop_w = 1;
    end
end


for r_k = 1:ini_popsize
    
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * ini_randMat(r_i,r_k);
    end
    
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * ini_randMat(r_i,r_k);
    end
    
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * ini_randMat(r_i,r_k);
    end
    
    nv.simulate()
    
    for r_l = 1 : ini_NrFitVal
        results_1 (r_l,:)  = nv.out(ini_FitString{r_l});
        results(r_k,r_l) = results_1(r_l,ini_FitTime(r_l));              %?? does it work ini_FitTime
    end
end

Max_matrix = [results(:,:), ini_randMat'];


%% Selection (roulette wheel method)-> sel_...
%probability of each parameterset to be choosen for the offspring
%population -> probability for parameter sets with a high
%calciumconcentration is higher

%probabilities

sel_sum_prob(:,:) = zeros(ini_popsize,ini_NrFitVal);
for sel_j = 1 : ini_NrFitVal                                               %??
    
    for sel_i = 1:ini_popsize                                               %ini_popsize= size of the population
        sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - Max_matrix(sel_i,sel_j)));                                  %probability for each individual to be chosen
        
        if sel_prob(sel_i,sel_j) == 0;
            sel_prob(sel_i,sel_j) = 0.0000000000000000001;
        end
    end
    %saving the probabilities
    sel_sum(:,sel_j) = sum(sel_prob(:,sel_j));
    
    for sel_i = 1 : ini_popsize
        sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
    end
    
    sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
    
    for sel_i = 1 : ini_popsize
        sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
    end
end

for ini_i = 1 : ini_popsize
    sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
    sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
end


sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec, Max_matrix]; %numbered max_matrix with probabilities


stop_max = max(sel_max_matrix(:,3));
stop_i = 1;
while sel_max_matrix(stop_i,3) < stop_max
    stop_i = stop_i+1;
end

stop_max_line = sel_max_matrix(stop_i,:);
best_results(count_max,:) = stop_max_line(:,:);


% decision
sel_output      = zeros(ini_popsize/2 , 3 + ini_NrFitVal + sum(ini_NrPar));   %creates a zero matrix plus 1 row sel_cumprobvec + 1 linspace + # of fitness values

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


%% Multiplicate Parents -> mpar...

% # of calculations

for mpar_i = 1:sel_2/2                                                          %divided by 2, because there are 2 parents
    
    mpar_P1 = sel_output(mpar_1,:);                                               %choosing parents out of input matrix
    mpar_P2 = sel_output(mpar_1+1,:);
    
    for mpar_j = 1:(ini_popsize/(ini_popsize / 2))                                                         % Repeat until initial population size is reached
        
        parent_matrix(mpar_c,:)   = mpar_P1;
        parent_matrix(mpar_c+1,:) = mpar_P2;
        
        mpar_c = mpar_c + 2;
        
    end
    mpar_1 = mpar_1 + 2;
end


%% Mutation of Parent Matrix -> mu...
% mutation occurs randomly

for mu_1 = 1:length(parent_matrix(:,1))                                       % Running variable for lines
    
    for mu_2 = 1 : sum(ini_NrPar)                                                     % Running variable for number of parameters / columns
        mu_rand = rand(1);                                                     % Random number between 0 and 1
        
        if mu_rand <= ini_mu_p                                                 % decision of mutation occurs or not
            for mu_var_1 = 1 : sum(ini_NrPar)
                if mu_2 == mu_var_1;                                                      % defines parameter specific boundaries
                    mu_r                            = ini_low_b(mu_2) + (ini_up_b(mu_2) - ini_low_b(mu_2))*rand(1,1);
                    parent_matrix(mu_1,mu_2 + 3 + ini_NrFitVal)    = mu_r;
                end
            end
            
        end
        
        mu_2 = mu_2 + 1;
    end
    mu_2 = 1;
    
end

mu_p_total(r_simu_2) = mean(ini_mu_p);

%% Cross-over -> co...

co_s = length(parent_matrix(:,1));                                          %length of the column of the parents_matrix (should be initial population)
for co_var  = 1:(co_s-1)
    co_rand = rand(1);
    
    if co_rand <= ini_co_p                                                          %probability of cross-over
        co_var1 = randi(sum(ini_NrPar))+3+ini_NrFitVal;                                                  %defines the lower boundary of the crossover section
        co_var2 = randi([co_var1,sum(ini_NrPar)+3+ini_NrFitVal]);                                        %defines the upper boundary of the crossover section
        
        co_1 = parent_matrix(co_var,co_var1:co_var2);                               %cutting out the cross-over section
        co_2 = parent_matrix(co_var+1,co_var1:co_var2);
        
        parent_matrix(co_var,co_var1:co_var2)  = co_2;                              %change of the cross-over section between the 2 childrens
        parent_matrix(co_var+1,co_var1:co_var2) = co_1;
        
    end
    co_var = co_var + 2;
end                                                                         %output is a matrix which is 2*mpar_rbigger than matrix_input


%% Elitism
el_max(1)  = max(parent_matrix(:,3)) * (1-ini_el_k);                        %Who is elite? maximal calcium concentration of generation multiplied by a number between 1-0
for el_var = 1:length(sel_output(:,1))  %for the length of sel_output check if the the calcium concentration of the set is bigger than the el_max....
    el_u=0;
    if sel_output(el_var,4) >= el_max(1)%... if it is bigger save it in a matrix and check where the same number of the set is in the parentmatrix, when founded overwirte the line
        el_p(el_var,:) = sel_output(el_var,:);
        for el_i = 1:length(parent_matrix(:,1));
            if el_p(el_var,1) == parent_matrix(el_i,1) && el_u== 0
                parent_matrix(el_i,:) = el_p(el_var,:);
                el_u=1;
            end
        end
    end
end                                                                         %output is a matrix which is 2*mpar_rbigger than matrix_input

%% simulation of the second generation

for r_k = 1:ini_popsize
    
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
    
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
    
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
    
    nv.simulate()
    
    for r_l = 1 : ini_NrFitVal
        results_1 (r_l,:)  = nv.out(ini_FitString{r_l});
        results(r_k,r_l) = results_1(r_l,ini_FitTime(r_l));              %?? does it work ini_FitTime
    end
end

parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(:,1:sum(ini_NrFitVal));


%% Third and following generations


%Maximum running time in ?????
stop_time_cond   = 0;
tic
stop_ac_count = 0;

while  stop_time_cond ==0 && stop_ac_count < ini_NrFitVal
    
    %Reset all running variables so the whole process can run again
    sel_p       = 1;
    
    mpar_1      = 1;
    mpar_c      = 1;
    
    co_var      = 1;
    
    el_i        = 1;
    
    %% Selection (roulette wheel method)-> sel_...
    
    %probability of each parameterset to be choosen for the offspring
    %population -> probability for parametersets with a high
    %calciumconcentration is higher
    
    %probabilities
    
    sel_sum_prob(:,:) = zeros(ini_popsize,ini_NrFitVal);
    for sel_j = 1 : ini_NrFitVal                                               %??
        
        for sel_i = 1:ini_popsize                                               %ini_popsize= size of the population
            sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - parent_matrix(sel_i,sel_j+3)));                                  %probability for each individual to be chosen
            
            if sel_prob(sel_i,sel_j) == 0;
                sel_prob(sel_i,sel_j) = 0.0000000000000000001;
            end
        end
        %saving the probabilities
        sel_sum(:,sel_j) = sum(sel_prob(:,sel_j));
        
        for sel_i = 1 : ini_popsize
            sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
        end
        
        sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
        
        for sel_i = 1 : ini_popsize
            sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
        end
    end
    
    for ini_i = 1 : ini_popsize
        sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
        sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
    end
    
    
    sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
    sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec, parent_matrix(:,4:sum(ini_NrPar)+3+ini_NrFitVal)]; %numbered max_matrix with probabilities
    
    stop_max = max(sel_max_matrix(:,3));
    stop_i = 1;
    while sel_max_matrix(stop_i,3) < stop_max
        stop_i = stop_i+1;
    end
    
    stop_max_line = sel_max_matrix(stop_i,:);
    best_results(count_max+1,:) = stop_max_line(:,:);
    
    
    % decision
    sel_output      = zeros(ini_popsize/2 , 3 + ini_NrFitVal + sum(ini_NrPar));   %creates a zero matrix plus 1 row sel_cumprobvec + 1 linspace + # of fitness values
    
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
    
    
    
    %% Multiplicate Parents
    
    for mpar_i = 1:sel_2/2                                                              %divided by 2, because there are 2 parents
        
        mpar_P1 = sel_output(mpar_1,:);                                               %choosing parents out of input matrix
        mpar_P2 = sel_output(mpar_1+1,:);
        
        for mpar_j = 1:ini_popsize/(ini_popsize / 2);
            
            parent_matrix(mpar_c,:)   = mpar_P1;
            parent_matrix(mpar_c+1,:) = mpar_P2;
            
            mpar_c = mpar_c + 2;
        end
        mpar_1 = mpar_1 + 2;
    end
    
    
    %% Mutation
    
    for mu_1 = 1:length(parent_matrix(:,1))
        
        for mu_2 = 1 : sum(ini_NrPar)
            
            mu_max      = max(parent_matrix(:,3));                          %Row where survival probability is located
            mu_mean     = mean(parent_matrix(:,3));
            mu_current  = parent_matrix(mu_1,3);
            
            if mu_current >= mu_mean                                                     %makes mutation rate adaptive to quality of organism's fitness
                p_mu(mu_var) = ini_mu_p * ((mu_max - mu_current) / (mu_max - mu_mean));     %If organism is fitter than average -> mutation should decrease
            else
                p_mu(mu_var) = ini_mu_p;                                            %if not mutation rate stays the same
            end
            
            mu_rand = rand(1);
            
            if mu_rand <= p_mu(mu_var)                                                 % decision of mutation occurs or not
                for mu_var_1 = 1 : sum(ini_NrPar)
                    if mu_2 == mu_var_1;                                                      % defines parameter specific boundaries
                        mu_r = ini_low_b(mu_2) + (ini_up_b(mu_2) - ini_low_b(mu_2))*rand(1,1);
                        parent_matrix(mu_1,mu_2 + 3 + ini_NrFitVal) = mu_r;
                    end
                end
            end
            mu_var   = mu_var + 1;
        end
        mu_2 = 1;
        
    end
    
    mu_p_total(r_simu_2+1) = mean(p_mu);                                              %in order to plot change of mutation over all generations
    
    %% Cross-over
    co_s = length(parent_matrix(:,1));                                          %length of the column of the parents_matrix (should be initial population)
    
    while co_var      <= (co_s-1)
        
        co_rand        = rand(1);
        cr_max         = max(parent_matrix(:,3));
        cr_mean        = mean(parent_matrix(:,3));
        cr_current     = parent_matrix(co_var,3);
        
        if cr_current >= cr_mean
            co_p = ((cr_max - cr_current) / (cr_max - cr_mean));
        else
            co_p = 1;
        end
        
        if co_rand <= co_p                                                              %probability of cross-over
            
            co_var1 = randi(sum(ini_NrPar))+ 3 + ini_NrFitVal;                                               %defines the lower boundary of the crossover section
            co_var2 = randi([co_var1,sum(ini_NrPar) + 3 + ini_NrFitVal]);                                     %defines the upper boundary of the crossover section
            
            co_1    = parent_matrix(co_var,co_var1:co_var2);                               %cutting out the cross-over section
            co_2    = parent_matrix(co_var+1,co_var1:co_var2);
            
            parent_matrix(co_var,co_var1:co_var2)   = co_2;                             %change of the cross-over section between the 2 childrens
            parent_matrix(co_var+1,co_var1:co_var2) = co_1;
            
        end
        co_var = co_var + 2;
    end                                                                                 %output is a matrix which is 2*mpar_rbigger than matrix_input
    
    %% Elitism
    clear el_p;
    el_max(r_simu_2+1)  = max(parent_matrix(:,3)) * (1-ini_el_k);
    for el_var = 1:length(sel_output(:,1))
        el_u=0;
        if sel_output(el_var,3) >= el_max(1,r_simu_2+1)
            el_p(el_var,:) = sel_output(el_var,:);
            for el_i = 1:length(parent_matrix(:,1));
                if el_p(el_var,1) == parent_matrix(el_i,1) && el_u== 0
                    parent_matrix(el_i,:) = el_p(el_var,:);
                    el_u=1;
                end
            end
        end
    end
    
    %%simulation
    for r_k = 1:ini_popsize
        
        for r_i = 1 : ini_NrPar(1)
            nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
        end
        
        for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
            nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
        end
        
        for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
            nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
        end
        
        nv.simulate()
        
        for r_l = 1 : ini_NrFitVal
            results_1 (r_l,:)  = nv.out(ini_FitString{r_l});
            results(r_k,r_l) = results_1(r_l,ini_FitTime(r_l));              %?? does it work ini_FitTime
        end
    end
    
    parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(:,1:sum(ini_NrFitVal));
    
    
    %time
    
    if ini_max_time <= toc
        stop_time_cond = 1;
    end
    
    
    r_simu_2 = r_simu_2 + 1;
    count_max = count_max+1;
    
    
end

%probablities of the generation for saving the best
sel_sum_prob(:,:) = zeros(ini_popsize,ini_NrFitVal);
for sel_j = 1 : ini_NrFitVal                                               %??
    
    for sel_i = 1:ini_popsize                                               %ini_popsize= size of the population
        sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - parent_matrix(sel_i,sel_j+3)));                                  %probability for each individual to be chosen
        
        if sel_prob(sel_i,sel_j) == 0;
            sel_prob(sel_i,sel_j) = 0.0000000000000000001;
        end
    end
    %saving the probabilities
    sel_sum(:,sel_j) = sum(sel_prob(:,sel_j));
    
    for sel_i = 1 : ini_popsize
        sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
    end
    
    sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
    
    for sel_i = 1 : ini_popsize
        sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
    end
end

for ini_i = 1 : ini_popsize
    sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
    sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
end


sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec, parent_matrix(:,4:sum(ini_NrPar)+3+ini_NrFitVal)]; %numbered max_matrix with probabilities

stop_max = max(sel_max_matrix(:,3));
stop_i = 1;
while sel_max_matrix(stop_i,3) < stop_max
    stop_i = stop_i+1;
end

stop_max_line = sel_max_matrix(stop_i,:);
best_results(count_max+1,:) = stop_max_line(:,:);





%% Plots
plot_w = 1;
plot_j = 1;
plot_k = 1;

while plot_w ==1
    plot_l = 1;
    for plot_i = plot_j : plot_j+3
        plot_inifit_vec = linspace(ini_FitVal(plot_i),ini_FitVal(plot_i),r_simu_2+1); %vector for the initial value of the calcium concentration
        
        figure(plot_k)
        subplot(2,2,plot_l)
        plot(linspace(1,r_simu_2+1,r_simu_2+1)', best_results(:,plot_i+3))
        plot_str_1 = {'Best',ini_FitString{plot_i}, 'value of each generation'};
        plot_strjoin_1 = strjoin(plot_str_1);
        legend(plot_strjoin_1,'User defined value');
        title(plot_strjoin_1);
        xlabel('# of generation');
        ylabel(ini_FitString{plot_i});
        hold on;
        grid on;
        
        subplot(2,2,plot_l)
        plot(linspace(1,r_simu_2+1,r_simu_2+1)',plot_inifit_vec);
        plot_str_2 = {'Best',ini_FitString{plot_i}, 'value of each Generation'};
        plot_strjoin_2 = strjoin(plot_str_2);
        legend(plot_strjoin_2,'User defined value');
        title(plot_strjoin_1);
        xlabel('# of generation');
        ylabel(ini_FitString{plot_i});
        grid on;
        grid minor;
        
        plot_l = plot_l + 1;
        
        if plot_i >= ini_NrFitVal
            plot_w = 0;
            break
        end
        
        
    end
    
    plot_k = plot_k + 1;
    plot_j = plot_j + 4;
end


figure(plot_k + 1)
plot(linspace(1,r_simu_2,r_simu_2)',mu_p_total');
legend('Change of average mutation rate');
title('Change of average mutation rate');
xlabel('# of generation');
ylabel('probability of mutation');
grid on;
grid minor;






