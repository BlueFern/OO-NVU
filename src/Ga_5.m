%% Genetic Algorithm

clear; clc; close all

% The Genetic Algorithm (GA) is a search method based on the natural selection. 
% A fitness value or fitness function, after which a certain problem is
% optimized, is used by the GA. In this code it is also possible to
% optimize after more values.
% The algorithm includes the implementation of a randomised first generation, 
% selection, mutation, cross-over and elitism.
% After these mechanisms a new generation is formed and treated like the first generation. 
% These mechanisms are repeated till one of the stopping critera is fulfilled.

%% General definitions

% Organism/Member   = parameterset consists of the parameters
% Generation        = one population with x members

%% Workspace definitions

%ui_...     -> variables/matricies/vectors of the user input
%ini_...    -> variables/matricies/vectors which are inital settings like the mutation rate
%pop_...    -> variables/matricies/vectors of the inital population 
%sel_...    -> variables/matricies/vectors of the selection
%mpar_...   -> variables/matricies/vectors of the multiplication of the parents
%mu_...     -> variables/matricies/vectors of the mutation 
%co_...     -> variables/matricies/vectors of the crossing over
%stop_...   -> variables/matricies/vectors of the stopping criterion
%plot_...   -> variables/matricies/vectors of the plots
%r_...      -> running variables

%T              -> Timesteps

%parent_matrix  ->  first column: nummeration, 
%                   second column: cumulative probability,
%                   third column: survivial probablility,
%                   4+number fitness value column: results,
%                   behind: all parameters / line = one organism  
 
%best_results   ->  Matrix with the best result of each generation (like
%                   parent_matrix)

%% Options
% First the options need to be set for the ODE solver (currently |ode15s|)
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC(), ...
    'odeopts', odeopts);

%% User input
% number of fitness values and name of fitness values (strings) 
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

 
 
% Input of fitness values to the according "fitness strings"
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
% fitness value are timedepended, so user has to decide at which time he is
% interested in.
for   ui_i_FitTime = 1 : ini_NrFitVal
    ui_promp_FitTimeSave = {'Enter the timestep for the results for' ini_FitString{ui_i_FitTime}};
    ui_prompt_vec_FitTime{ui_i_FitTime} = strjoin(ui_promp_FitTimeSave);
end

ui_title_FitTime = 'Add timesteps for you fitness results';
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

 
 
%Paramters saves the number of paramters in each cell and the strings
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
    };                                                                      
 
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
    };                                                                      
 
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
%in which space the initial parameter value should be changed - 0.8 means
%80% of the initial values 1.2 menas 120% of the initial value 
%so a random number between 0.8 and 1.2 is created

%Checks if user wants individual boundaries for each parameter
ui_check_1 = {'Do you want individual parameter boundaries (y=yes/n=no)?'};
ui_title_check = 'In which percentage boundaries should the parameter change?';
ui_lines_check = 1;
ui_defaultanswer_check = {'n'};
ui_check_dlg = inputdlg(ui_check_1,ui_title_check,ui_lines_check,ui_defaultanswer_check);
ui_answer            = ui_check_dlg{1};
ui_yes               = 'y';
ui_check = strcmp(ui_answer,ui_yes);

%if the user wants individual boundaries:
%depending on the user input (which parameters)
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
    
    %saves inputs as double
    for ui_j_bound = 1 : sum(ini_NrPar)/2
        ini_low_b(ui_j_bound)                              = str2double(ui_bound_1{ui_j_bound * 2 -1});
        ini_up_b(ui_j_bound)                               = str2double(ui_bound_1{ui_j_bound * 2});
        ini_low_b(sum(ini_NrPar)/2 + ui_j_bound)           = str2double(ui_bound_2{ui_j_bound * 2 -1});
        ini_up_b(sum(ini_NrPar)/2 + ui_j_bound)            = str2double(ui_bound_2{ui_j_bound * 2})   
    end
    
    ui_j_bound = 1;
    
% same boundaries for all parameter 
else                                                              
    
    ui_prompt_samebound = {'Enter lower boundary of randomisation for all parameters:'
        'Enter upper boundary of randomisation for all parameters:'};
    ui_title_samebound     = 'Individual boundaries for the parameters'; 
    ui_lines_samebound     = 1;
    ui_defaultanswer_samebound = { '0.8','1.2'};
    ui_samebound = inputdlg(ui_prompt_samebound,ui_title_samebound,ui_lines_samebound,ui_defaultanswer_samebound,'on');
    
    %saves inputs as double
    for ui_j_bound = 1 : sum(ini_NrPar)
        ini_low_b(ui_j_bound)           = str2double(ui_samebound{1});
        ini_up_b (ui_j_bound)           = str2double(ui_samebound{2});
    end
    
    ui_j_bound = 1;
end

%General User Input
ui_prompt_general = { 
    'Populationsize: With which probability should the generation be uniform distributed:' 
    % if the input is 0.95 the first generation is to 95% ultimate distributed                
    'Cross-Over: Enter the initial probability of cross-over:'  
    % this value determines the frequency of cossing over
    % 80% means in 80% of the organisms of one population occurs crossing over
    % size and position is randomly
    'Mutation: Enter the initial probability of Mutation'                                                           
    % Mutation probability
    % 50% means the probablity of beeing randomly changed (in boundaries of
    % the parameters) of each parameter of each member of the FIRST generation
    % is 50%. For the following generations it is calculated (see code)
    'Elitism: What percentage of the members of the first generation should survive without any genetic change' 
    % 90% means that every organism of one gerenation, which survival probability is higher than 90% 
    % of the best of the generation, is allowed to survive
    % without any change in the parameters 
    'Stopping criterion: Enter the maximal simulation time [minutes]'   
    'Stopping criterion: Enter the convergenz criteria number (How many generations should not have a change):'
    % How many "close" values have to be in a row, in order to make the
    % code stop
    'Stopping criterion: Enter the boundaries of convergence:'  
    % Decision criteria if 2 values are close together or not
    % o.1 means if the distance between two following values of the bests of two generations of one fitness
    % value (i.e. c_k) is smaller than 0.1. 
    % this is proved for all entered fitness values seperatly
    % how often this have to occur befor the code stops is determed by the
    % convergenz criteria number
    'Stopping criterion: Do you want different accuracies for the different fitness values'
    % accuarcy means how accurate the aimed value should be reached
    % i.e. a value of 0.1 means the programm should stop if one fitness
    % value is an element the interval of [0.9 1.1] of the aimed value 
    };

ui_title_general = 'General User Input';
ui_lines_general = 1;
ui_defaultanswer_general = {'0.9','0.8','0.5','0.1','10','15','0.02','y'}; 
ui_general = inputdlg(ui_prompt_general,ui_title_general,ui_lines_general,ui_defaultanswer_general,'on');

ini_ulti_dist          = 1-str2double(ui_general{1});
ini_co_p               = str2double(ui_general{2});
ini_mu_p               = str2double(ui_general{3});
ini_el_k               = str2double(ui_general{4});
ini_max_time           = str2double(ui_general{5})*60;
ini_nr_convergence     = str2double(ui_general{6});
ini_bound_convergence  = str2double(ui_general{7});
ui_answer_accuracy     = ui_general{8};
ui_yes_accuracy        = 'y';
ui_check_accuracy = strcmp(ui_answer_accuracy,ui_yes_accuracy);

%different accuracies for the different fitness values
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
%same accuracies for all parameters    
else                                                            
    ui_prompt_sameaccuracy = {'Enter the accuracy for all parameters:'};
    ui_title_sameaccuracy     = 'Add the accuracy for all fitness values to reach the aimed value' ;
    ui_lines_sameaccuracy     = 1;
    ui_defaultanswer_sameaccuracy = { '0.8'};
    ui_sameaccuracy = inputdlg(ui_prompt_sameaccuracy,ui_title_sameaccuracy,ui_lines_sameaccuracy,ui_defaultanswer_sameaccuracy,'on');
    
    for ui_i_accuracy = 1 : ini_NrFitVal
        ini_accuracy(ui_i_accuracy)           = str2double(ui_sameaccuracy{1});
    end
end

%% predefined vectors/matrices/variables apart form the user input
T = linspace (0,500,1000);  
ini_popsize = 2^4;

%% running variables
mpar_1      = 1;
mpar_c      = 1;

mu_var      = 1;

conv_crit   = 1;
r_simu_2    = 1;
ini_w       = 1;

%% Preallocation 1
ini_randMat = zeros(sum(ini_NrPar),ini_popsize); %this preallocation reduces the simulation time

%% initial values of the parameters
% saves the inital parameter values of the NVU model in an vector(ini_value)

for r_i= 1:ini_NrPar(1) % values of the astrocyte
    % ini_NrPar(1) is the number of parameters in the astrocyte which are changed
    ini_a = nv.astrocyte.params; 
    ini_value(r_i) = getfield(ini_a, ini_ParamString{r_i});
end
for r_i= 1:ini_NrPar(2)
    ini_a = nv.smcec.params;
    ini_value(r_i+ini_NrPar(1)) = getfield(ini_a, ini_ParamString{r_i+ini_NrPar(1)}); 
    %saves the values under the values of the astrocyte parameters
end
for r_i= 1:ini_NrPar(3)
    ini_a = nv.wall.params;
    ini_value(r_i+ini_NrPar(1)+ini_NrPar(2)) = getfield(ini_a, ini_ParamString{r_i+ini_NrPar(1)+ini_NrPar(2)});
end

%% inital population 
 
% the first generation is randomly created with a defined population size
for ini_i = 1:sum(ini_NrPar)
    ini_randMat(ini_i,:) = vertcat(1,(ini_up_b(ini_i)-ini_low_b(ini_i))*rand(ini_popsize-1,1)+ini_low_b(ini_i));
   %creats an random number in the user defined boundaries 
end
% ini_randMat matrix [n x m] (n = parameters (e.g frist line k_on) 
% and m = number of organisms

 
% uniform distribution
pop_w = 0;

while pop_w == 0;
    % calculates the distribution of the first generation 
    for pop_i = 1 : sum(ini_NrPar)
        pop_c = ((ini_up_b(pop_i)+ini_low_b(pop_i))/2);
        
        for pop_j   = 1 : ini_popsize
            pop_d   = pop_c  - ini_randMat(pop_i,pop_j);
            pop_vector_1(pop_i,pop_j) = pop_d;
        end
        
        pop_p(pop_i)    = abs(sum(pop_vector_1(pop_i,:))/ini_popsize);
        pop_p_1(pop_i)  = (pop_p(pop_i) / abs(pop_c - ini_up_b(pop_i)));
    end
    
    pop_res = (sum(pop_p_1) / sum(ini_NrPar));
    
    % checks if the population is bigger than the user input 
    if pop_res > ini_ulti_dist
        ini_popsize = ini_popsize*2;
        
        for ini_i = 1:sum(ini_NrPar)
            ini_randMat(ini_i,2:ini_popsize) = vertcat((ini_up_b(ini_i)-ini_low_b(ini_i))*...
            rand((ini_popsize-1),1)+ini_low_b(ini_i)); %creates new random parametersets
        end
        
    else
        pop_w = 1; %jumps out of the loop if the gerneration is uniform distributed
    end
end

%% Preallocation 2
%depending on the population size 
%reduces the simulation time
results             = zeros(ini_popsize,ini_NrFitVal);
sel_sum_prob(:,:)   = zeros(ini_popsize,ini_NrFitVal);
sel_prob            = zeros(ini_popsize,ini_NrFitVal);
sel_sum             = zeros(1,ini_NrFitVal);
sel_prob_1          = zeros(ini_popsize,ini_NrFitVal);
sel_prob_2          = zeros(ini_popsize,ini_NrFitVal);
sel_prob_vec_1      = zeros(ini_popsize,1);
sel_prob_vec        = zeros(ini_popsize,1);
sel_output          = zeros(ini_popsize/2 , 3 + ini_NrFitVal + sum(ini_NrPar));
parent_matrix       = zeros(ini_popsize , 3 + ini_NrFitVal + sum(ini_NrPar));
el_var              = zeros(ini_popsize, 3 + ini_NrFitVal + sum(ini_NrPar));
stop_count          = zeros(1,ini_NrFitVal);

parent_matrix(:,4+ini_NrFitVal:3+ini_NrFitVal+sum(ini_NrPar))= ini_randMat';
% rand matrix is now inclued in the parent_matrix (parent_matrix is now
% constantly used) 
 
%% simulation of the first generation
for r_k = 1:ini_popsize
    
    %changes the parameter values in the astrocyte for one organism
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
   
    %changes the parameter values in the smooth muscle cell for one organism
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the wall mechanics for one organism
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_k,r_i + 3 + ini_NrFitVal);
    end
    
    %simulates the NVU model with the changed parameter values 
    nv.simulate()

    %saves the results
    for r_l = 1 : ini_NrFitVal
        results_1(r_l,:)  = nv.out(ini_FitString{r_l}); %all timesteps of one value for one paraemterset
        results(r_k,r_l) = results_1(r_l,ini_FitTime(r_l));  % one timestep  
    end
end

%saves all results in the parent matrix
parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(:,1:sum(ini_NrFitVal));

%% calculation of the survival probabilities and cumultative probabilities

for sel_j = 1 : ini_NrFitVal                                               
    
    for sel_i = 1:ini_popsize                                               
        sel_prob(sel_i,sel_j) = (abs(ini_FitVal(sel_j) - parent_matrix(sel_i,sel_j+3))); 
      % difference between current value and aimed value
        
        if sel_prob(sel_i,sel_j) == 0;
            sel_prob(sel_i,sel_j) = 0.0000000000000000001;
            % if one fitness value reaches the aimed value the code doesn't
            % work anymore, because it can't divivde trough 0 so the value
            % is nearly settet to 0 
        end
    end
    
    sel_sum(:,sel_j) = sum(sel_prob(:,sel_j)); 
    % sum of all distances for one fitness value 
    
    for sel_i = 1 : ini_popsize
        sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
        %sum of distances divided by every single distance
    end
    
    sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
    % sum of (sum of distances divided by every single distance)
    
    for sel_i = 1 : ini_popsize
        sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
        % (sum of distances divided by every single distance)divided by 
         %(sum of (sum of distances divided by every single distance))
    end
    
    %survival probability for every organism of one generation for ONE
    %Fitness value is calculated
end

%average survival probability 
for ini_i = 1 : ini_popsize
    sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
    sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
end
% average survival probability for every organism of one generation
% calculated by the average of the survival probabilities of each fitness
% value 
 
%cumulative probability 
sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec,...
                   parent_matrix(:,4:sum(ini_NrPar)+3+ini_NrFitVal)]; %shape like parent_matrix

%% save best resut of the generation 
stop_max = max(sel_max_matrix(:,3)); % best survival probability -> best member of the generation 
stop_i = 1;
while sel_max_matrix(stop_i,3) < stop_max
    stop_i = stop_i+1;
end
best_results(r_simu_2,:) = sel_max_matrix(stop_i,:);

 
%% settings and calculations for stopping criterions

% accuracy (stopping criterion)
stop_ac_count = 0;
for stop_i = 1:ini_NrFitVal
    if best_results(r_simu_2,stop_i+3) > (ini_FitVal(stop_i)-ini_FitVal(stop_i)*ini_accuracy(stop_i)) && ...
       best_results(r_simu_2,stop_i+3 ) < (ini_FitVal(stop_i)+ini_FitVal(stop_i)*ini_accuracy(stop_i))
       %checks if the best result is in special boundaries around the aimed
       %value -> depending on the user input for the accuracy 
       stop_ac_count = stop_ac_count+1;
    end
end

%time
stop_time_beginning = clock;
%Start time of simulation
stop_time_elapsed = 0;
%Just a preallocation

 
 
%% Loop for second and following generations

while  stop_time_elapsed < ini_max_time  && stop_ac_count < ini_NrFitVal && stop_count(1)< ini_nr_convergence
    
    %Reset all running variables so the whole process can run again
    sel_p       = 1;
    
    mpar_1      = 1;
    mpar_c      = 1;
    
    co_var      = 1;
    
    el_i        = 1;
    
    
    %% Selection (roulette wheel method)-> sel_..   
    %probability of each parameterset to be choosen for the offspring
    %population -> probability to survice for parametersets with a 
    %smaller distance to the aimed value is higher
    
    % decision
    sel_output =  zeros(ini_popsize/2 , 3 + ini_NrFitVal + sum(ini_NrPar));
    %preallocation is needed again because of the following while loop
    %conditions
    sel_rand_value  = rand(length(sel_output(:,1)),1);
    %random vector with the length of the number of parents and with random numbers between 0-1
    for sel_2=1:length(sel_output(:,1))
        %for the length of the number of the parents, that should be chosen
        
        while sel_output(sel_2,1) == 0
            %till the row of the sel_output rows does include an organism
            sel_p=1;
            if sel_rand_value(sel_2,1) < sel_max_matrix(sel_p,2);
                % if random value at the position n is smaler as the cummultative probability ...
                sel_output(sel_2,:) = sel_max_matrix(sel_p,:);
                % ... save the individum in the sel_output matrix
            end
            
            while sel_rand_value(sel_2,1) > sel_max_matrix(sel_p,2);
                % while the random number is bigger ....
                sel_p=sel_p+1;
                if sel_rand_value(sel_2,1) <= sel_max_matrix(sel_p,2);
                    % test again if the cummultative probability of the next organism
                    % (parameterset) is bigger than the random number
                    sel_output(sel_2,:) = sel_max_matrix(sel_p,:);
                    % ... save the individum in the sel_output matrix
                end
            end
            
            
            if sel_2>1
                % prevent, that one individum (parameterset) is taken twice
                sel_duplex = ismember(sel_output(sel_2,1),sel_output(1:(sel_2-1),1));
                % proves if the chosen individum is already in the sel_outtput matrix
                if sel_duplex == 1
                    % if yes ...
                    sel_output(sel_2,1)     = 0;
                    % delet it in the sel_otuput matrix
                    sel_rand_value(sel_2,1) = rand(1);
                    % create a new random number at this place and repeat all the steps
                end
            end
            
        end
    end
    
    
    %% Multiplicate Parents
    for mpar_i = 1:sel_2/2                                                              
        
        mpar_P1 = sel_output(mpar_1,:);                                               
        mpar_P2 = sel_output(mpar_1+1,:);
        %Taking out two parents (are always located close to each other)
        
        for mpar_j = 1:ini_popsize/(ini_popsize / 2);
        %Repeat in (in this case) two times in order to achieve initial
        %population size.
            
            parent_matrix(mpar_c,:)   = mpar_P1;
            parent_matrix(mpar_c+1,:) = mpar_P2;
            
            mpar_c = mpar_c + 2;
            %Always a pair of parents
        end
        mpar_1 = mpar_1 + 2;
    end
        
    
    %% Mutation
    for mu_1 = 1:length(parent_matrix(:,1))
    %Mutation occurs for ALL organism in ONE generation 
    
        for mu_2 = 1 : sum(ini_NrPar)
        %For every Parameter the GA decides whether Mutation occurs or not   
        
            for  mu_3 = 1 : ini_NrFitVal
            %Mutation rate is depended on how close the fitness value(s) 
            %of the current observed organism are to
            %the user defined threshold(s). The closer they are, the smaller the 
            %mutation probability. The average of all fitness values is used to 
            %determine p_mu
  
                mu_max      = ini_FitVal(mu_3);                          
                mu_mean     = mean(parent_matrix(:,3+mu_3));
                mu_current  = parent_matrix(mu_1,3+mu_3);
                
                if abs(mu_current) <= 2*abs(mu_max)
                %Stops the mutation rate to rise over the inital mutation
                %prob.
                    p_mu_1(mu_3) =  ini_mu_p * abs(((mu_max - mu_current) / (mu_max)));
                else
                   p_mu_1(mu_3) = ini_mu_p;
                end

 
            end
            
            p_mu(mu_1) = mean(p_mu_1);           
            %As already described the mutation rate p_mu is depended on all
            %fitness values. Therefore a mean is taken
            
            mu_rand = rand(1);            
            %Random number between 0 and 1
            
            if mu_rand <= p_mu(mu_1)                 
            %Here is the decision if mutation happens or not. The higher
            %p_mu the higher the prob. that mu_rand is bigger than p_mu.
            
                for mu_var_1 = 1 : sum(ini_NrPar)
                    if mu_2 == mu_var_1;                        
                    %Refers to the position of the current parameter and
                    %changes there the the parameter within in the
                    %predefined boundaries and overwrites it in the
                    %parentmatrix
                    
                        mu_r = ini_low_b(mu_2) + (ini_up_b(mu_2) - ini_low_b(mu_2))*rand(1,1);
                        parent_matrix(mu_1,mu_2 + 3 + ini_NrFitVal) = mu_r;
                    end
                end
            end
           mu_var   = mu_var + 1; 
        end
        
        mu_2 = 1;
        
    end
    
    mu_p_total(r_simu_2) = mean(p_mu);                                             
    %Average mutation rate of current generation. This is used in plots.
    
    %% Cross-over
    co_s = length(parent_matrix(:,1));                                     
    co_i = 1;
    while co_var      <= (co_s-1)
    %Repeat until end of generation
        
        for co_var_1 = 1 : ini_NrFitVal
        %Calculate c.o. prob. for one organism. In order to do so the code 
        %has to calculate for every fitness value the distance between fitness value and user
        %defined threshold. The c.o. prob. is resulting of the average of 
        %these distances.
        
                        
            cr_max           = ini_FitVal(co_var_1);
            cr_mean          = mean(parent_matrix(:,3+co_var_1));
            cr_current_1     = parent_matrix(co_var,3+co_var_1);
            cr_current_2     = parent_matrix(co_var+1,3+co_var_1);
            
            if abs(cr_current_1) <= 2*abs(cr_max) && abs(cr_current_2) <= 2*abs(cr_max)
            %Same principle as in mutation. This should prevent a rise of
            %the c.o. prob. over the initial value
            %in crossing over the "quality" of both parents is important to
            %decide if it occurs or not.
            
                co_p_1(co_var_1) = ini_co_p * abs(((cr_max - cr_current_1) / (cr_max)));
                co_p_2(co_var_1) = ini_co_p * abs(((cr_max - cr_current_2) / (cr_max)));
            else
                co_p_1(co_var_1) = ini_co_p;
                co_p_2(co_var_1) = ini_co_p;
            end

            
        end
        co_mean_1 = mean(co_p_1);
        co_mean_2 = mean(co_p_2);
        co_p_1 = co_mean_1 + co_mean_2;
        co_p(co_i)  = mean(co_p_1);
        co_rand     = rand(1);
        if co_rand <= co_p                                                              
            
            co_var1 = randi(sum(ini_NrPar))+ 3 + ini_NrFitVal;                                               
            co_var2 = randi([co_var1,sum(ini_NrPar) + 3 + ini_NrFitVal]);           
            
            co_1    = parent_matrix(co_var,co_var1:co_var2);                               
            co_2    = parent_matrix(co_var+1,co_var1:co_var2);
            %Defines the size and position of the crossover section
            
            parent_matrix(co_var,co_var1:co_var2)   = co_2;                             
            parent_matrix(co_var+1,co_var1:co_var2) = co_1;
            %Swobs sections between parents
            
        end
        co_i = co_i + 1;
        co_var = co_var + 2;
        
    end
    %output is a matrix which is 2*mpar_rbigger than matrix_input
    
    
    %% Elitism
    clear el_p;
    el_max(r_simu_2)  = max(parent_matrix(:,3)) * (1-ini_el_k);
    %Defines the threshold for being elite. Template is the best organism
    %of current generation
    
    for el_var = 1:length(sel_output(:,1))
        el_u=0;
        if sel_output(el_var,3) >= el_max(1,r_simu_2)
        %Checks each line if it's bigger than threshold for elitism if it
        %is, this line gets safed in order to overwrite it in the parent
        %matrix
        
            el_p(el_var,:) = sel_output(el_var,:);
            for el_i = 1:length(parent_matrix(:,1));
            %Searches line in parent_matrix and overwrite line by elite
            %line
            
                if el_p(el_var,1) == parent_matrix(el_i,1) && el_u== 0
                    parent_matrix(el_i,:) = el_p(el_var,:);
                    el_u=1;
                end
            end
        end
    end

    
    %% simulation
    %same as before -> see comments above
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
            results(r_k,r_l) = results_1(r_l,ini_FitTime(r_l));             
        end
    end
    
    parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(:,1:sum(ini_NrFitVal));
    
    
  %% calculation of the survival probabilities and cumultative probabilities  
  %same as before -> see comments above
    sel_sum_prob(:,:) = zeros(ini_popsize,ini_NrFitVal);
    for sel_j = 1 : ini_NrFitVal                                              
        
        for sel_i = 1:ini_popsize                                               
            sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - parent_matrix(sel_i,sel_j+3)));                                  
            if sel_prob(sel_i,sel_j) == 0;
                sel_prob(sel_i,sel_j) = 0.0000000000000000001;
            end
       end
        
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
    sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec,...
                       parent_matrix(:,4:sum(ini_NrPar)+3+ini_NrFitVal)]; 
 
    
    stop_max = max(sel_max_matrix(:,3));
    stop_i = 1;
    
    while sel_max_matrix(stop_i,3) < stop_max
        stop_i = stop_i+1;
    end
    
    stop_max_line = sel_max_matrix(stop_i,:);
    best_results(r_simu_2+1,:) = stop_max_line(:,:);
    
 %% settings and calculations for stopping criterions  
    
    % accuracy -> see commments above
    
    stop_ac_count = 0;
    for stop_i = 1:ini_NrFitVal
        if best_results(r_simu_2+1,stop_i+3) > (ini_FitVal(stop_i)-ini_FitVal(stop_i)*ini_accuracy(stop_i)) && ...
           best_results(r_simu_2+1,stop_i+3) < (ini_FitVal(stop_i)+ini_FitVal(stop_i)*ini_accuracy(stop_i))
            stop_ac_count = stop_ac_count+1;
        end
    end
    
    
    % convergenz
    
    %if more best results in a row don't have a difference in EVERY fitness
    %value stop the loop 
    for stop_i = 1:ini_NrFitVal
    % check every fitness value if it hasn't changed since the last
    % generation 
        if abs(best_results(r_simu_2+1,stop_i)- best_results(r_simu_2,stop_i))< ini_bound_convergence
            stop_count(stop_i) = stop_count(stop_i)+1;
           %if it hasn't change, count one more
        else
            stop_count(stop_i) = 0;
        end
    end
    
    if 1 == ismember([0],stop_count)
    %if there is a zero in the vector set all to 0 
        stop_count(:)=0;
    end
  
    
    %time
    
    stop_time_end     = clock;
    stop_time_elapsed = etime(stop_time_end,stop_time_beginning);
    %Elapsed time between start of simulation and "now"

    r_simu_2 = r_simu_2 + 1;
end

%% Plots
%This plot adapts to the number of the fitness values. Max. number of
%graphs in one figure is 4. If there are more than 4 fitness values a new
%figure gets created. The input of the plot is completely depended on the
%user input or the results. 
 
plot_w = 1;
plot_j = 1;
plot_k = 1;

while plot_w ==1
    plot_l = 1;
    for plot_i = plot_j : plot_j+3
        plot_inifit_vec = linspace(ini_FitVal(plot_i),ini_FitVal(plot_i),r_simu_2); 
        %Creates the plot line for the user defined theshold
        
        figure(plot_k)
        subplot(2,2,plot_l)
        plot(linspace(1,r_simu_2,r_simu_2)', best_results(:,plot_i+3))
        plot_str_1 = {'Best',ini_FitString{plot_i}, 'value of each generation'};
        plot_strjoin_1 = strjoin(plot_str_1);
        %Makes legend adaptive to different fitness value names
        legend(plot_strjoin_1,'User defined value');
        title(plot_strjoin_1);
        xlabel('# of generation');
        ylabel(ini_FitString{plot_i});
        hold on;
        grid on;
        
        subplot(2,2,plot_l)
        plot(linspace(1,r_simu_2,r_simu_2)',plot_inifit_vec);
        plot_str_2 = {'Best',ini_FitString{plot_i}, 'value of each Generation'};
        plot_strjoin_2 = strjoin(plot_str_2);
        %Makes legend adaptive to different fitness value names
        legend(plot_strjoin_2,'User defined value');
        title(plot_strjoin_1);
        xlabel('# of generation');
        ylabel(ini_FitString{plot_i});
        grid on;
        grid minor;
        
        plot_l = plot_l + 1;
        
        if plot_i >= ini_NrFitVal
        %Check if all fitness values are already plotted or not, if yes
        %loop stops
        
            plot_w = 0;
            break
        end
        
        
    end
    
    plot_k = plot_k + 1;
    plot_j = plot_j + 4;
end

 
%% Sort best results
%The Matrix best_results gets treated like a new generation. The parameters
%stay the same, but the survival probability gets calculated a second time
%in order to find the best of the Matrix. With this result the best_result
%matrix can be sorted again by the survival probability.
clear sel_prob_vec
clear sel_prob
clear sel_sum
clear sel_prob_1
clear sel_prob_2
clear sel_prob_vec_1
clear sel_prob_vec
clear sel_cum_probvec

best_l = length(best_results(:,1));

%% Preallocation 3
sel_sum_prob        = zeros(best_l,ini_NrFitVal);
sel_prob            = zeros(best_l,ini_NrFitVal);
sel_sum             = zeros(1,ini_NrFitVal);
sel_prob_1          = zeros(best_l,ini_NrFitVal);
sel_prob_2          = zeros(best_l,ini_NrFitVal);
sel_prob_vec_1      = zeros(best_l,1);
sel_prob_vec        = zeros(best_l,1);

%% select the best organisms of the best organisms of every generation 
    for sel_j = 1 : ini_NrFitVal                                             
        
        for sel_i = 1:best_l                                
            sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - best_results(sel_i,sel_j+3)));    
            if sel_prob(sel_i,sel_j) == 0;
                sel_prob(sel_i,sel_j) = 0.0000000000000000001;
            end
        end
        
        
        sel_sum(:,sel_j) = sum(sel_prob(:,sel_j));
        
        for sel_i = 1 : best_l
            sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
        end
        
        sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
        
        for sel_i = 1 : best_l
            sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
        end
    end
    
    for ini_i = 1 : best_l
        sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
        sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
    end
    
sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
best_results  = [best_results(:,1), sel_cum_probvec, sel_prob_vec, best_results(:,4:sum(ini_NrPar)+3+ini_NrFitVal)];

%% Plot Mutation rates 
%This is just to visualise the decrease of the mutation rate over all 
%generations

figure(plot_k + 1)
plot(linspace(1,r_simu_2-1,r_simu_2-1)',mu_p_total');
legend('Change of average mutation rate');
title('Change of average mutation rate');
xlabel('# of generation');
ylabel('probability of mutation');
grid on;
grid minor;

%% Plot best results
%Here the Matrix get sorted by the newly calculated survival probs. 
%In this case the best 7 of the sorted Matrix get selected and plotted.

plot_max = sortrows(best_results,3);
result_string = {'Best result','Second best result','Third best result','4th best result','5th best result','6th best result','7th best result'};
color_string = {'b','g','r','y','m','c','k'};
plot_w = 1;
plot_i = 1;
plot_j = 1;
plot_max_unique = zeros(7,length(best_results(1,:)));

while plot_w ==1 && plot_j <=7 && plot_i < best_l
    if plot_max(best_l+1-plot_i,3) == plot_max(best_l-plot_i,3)
        plot_w=1;
    else
       plot_max_unique(plot_j,:) = plot_max(best_l+1-plot_i,:);
       plot_j = plot_j + 1;
    end
    plot_i = plot_i + 1;
end
figure(plot_k + 2)
for plot_i = 1 : 7
plot_max_line_1 = plot_max_unique(plot_i,:);
x = 1:length(ini_ParamString);
y_1 = plot_max_line_1(3+ini_NrFitVal+1 : 3+ini_NrFitVal+sum(ini_NrPar));
plot(x, y_1, color_string{plot_i})
hold on;
end
legend(result_string);
set(gca, 'XTick',1:length(ini_ParamString), 'XTickLabel',ini_ParamString)
grid on;
grid minor;

