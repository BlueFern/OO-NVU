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
    %List of output varaibles in the GA
    ui_prompt_plotfeat1 = {
        'R_k'
        'N_K_s'
        'N_Na_s'
        'N_HCO3_s'
        'N_K_k'
        'N_Na_k'
        'N_HCO3_k'
        'N_Cl_k'
        'w_k'
        'K_p'
        'Ca_i'
        's_i'
        'v_i'
        'w_i'
        'I_i'
        'K_i'
        'Ca_j'
        's_j'
        'v_j'
        'I_j'
        'Mp'
        'AMp'
        'AM'
        'M'
        'F_r'
        'R'};
    
    [out_plotfeat1,out_plotfeat_zero1] = listdlg('PromptString','Do you want to plot some output variables?',...
        'SelectionMode','multiple',...
        'ListString',ui_prompt_plotfeat1,...
        'ListSize',[250 300]);
    
    ini_outputplots1(1)= length(out_plotfeat1);
    
    plot_output1 = 1;
    plot_j = 2;
    
    if out_plotfeat1 ~= 0;
        for r_i = 1 : ini_outputplots1(1)
            plot_output_variables1(r_i) = ui_prompt_plotfeat1(out_plotfeat1(r_i));
        end
    else
        plot_output1 = 0;
    end 
    
    if plot_output1 ==1
        plot_l = 1;
       
        u = 1;
       
nv.simulate        
        
        for r_i = 1 : ini_outputplots1(1)

            figure(r_i)
            plot(nv.T, nv.out(ui_prompt_plotfeat1{r_i}))
            plot_str_1 = {plot_output_variables1{r_i}};
            hold on;
            grid on;
            
           
           
            plot_strjoin_1 = strjoin(plot_str_1);
            %Makes legend adaptive to different fitness value names
            legend('initial value', plot_strjoin_1);
            title(plot_strjoin_1);
            xlabel('Time [s]');
            ylabel(plot_output_variables1{r_i});
            hold on;
            grid on;
        
        end
        
    end

 
 
 
 
 
%Fitness function
ui_j = 1;
ui_prompt_Fit = {'Enter the fitness value:'
    'Do you want to minimise or maximise the fitness value? [min] or [max]'
    'Enter P = timepoint A = average (between 2 timepoints, R = changing rate (between 2 timepoints)'
   'Enter timepoint one'
    'Enter timepoint two (just for average and changing rate)'
    'Do you want to have addtional optimization features (y=yes/n=no)?'};

ui_title_Fit = 'specification of the fitness value';
ui_lines = 1;
ui_defaultanswer_Fit = {'R','max','P','370','350','y'};
ui_Fit = inputdlg(ui_prompt_Fit,ui_title_Fit,ui_lines,ui_defaultanswer_Fit,'on');

 
ini_FitString(1)    = ui_Fit(1);

ini_min_max = strcmp(ui_defaultanswer_Fit{2},ui_Fit{2});
%If ini_mmin_max = 1 means default suggestion is used

ui_What(1) = ui_Fit{3};
ui_P_What     = 'P';
ui_check_What1(1) = strcmp(ui_What(1),ui_P_What);
ui_A_What ='A';
ui_check_What2(1) = strcmp(ui_What(1),ui_A_What);

if ui_check_What1(1) == 1
    ini_FitWhat(1) = 1;
elseif ui_check_What2(1) == 1
    ini_FitWhat(1) = 2;
else
    ini_FitWhat(1) = 3;
end

ini_FitTime(ui_j) = str2double(ui_Fit{4})*2;
ini_FitTime2(ui_j) = str2double(ui_Fit{5})*2;

ui_answer_feature    = ui_Fit{6};
ui_yes_feature       = 'y';
ui_check_feature = strcmp(ui_answer_feature,ui_yes_feature);

if ui_check_feature == 1
    ui_i = 0;
else
    ui_i =1;
    
end
ini_NrFitVal = 1;
% All settings for the values, which should help the optimization
ui_j = ui_j + 1;
while ui_i == 0
    ui_What(length(ui_What)+1) = 0;
    ui_prompt_Feat = {'Enter the additional fitness feature:'
        'Enter P = time point, A = average, R = changing rate, B = Boundaries, I = Interval Boundaries'
        'Enter the aimed value (just for P, A or R'
        'Enter timepoint one'
        'Enter timepoint two (just for A, R, I)'
        'Enter the upper boundary value (just for B, I)'
        'Enter the lower boundary value (just for B, I)'
        'Do you want to have more optimization features (y=yes/n=no)'};
    
    ui_title_Feat = 'What should be added';
    ui_lines = 1;
    ui_defaultanswer_Feat = {'Ca_i','I','0.1','320','370','0.15','0.07','n'};
    ui_Feat = inputdlg(ui_prompt_Feat,ui_title_Feat,ui_lines,ui_defaultanswer_Feat,'on');
    
    
    ini_FitString(ui_j)    = ui_Feat(1);
    
    ui_What(ui_j) = ui_Feat{2};
    ui_P_What     = 'P';
    ui_check_What1(ui_j) = strcmp(ui_What(ui_j),ui_P_What);
    ui_A_What ='A';
    ui_check_What2(ui_j) = strcmp(ui_What(ui_j),ui_A_What);
    ui_R_What ='R';
    ui_check_What3(ui_j) = strcmp(ui_What(ui_j),ui_R_What);
    ui_B_What ='B';
    ui_check_What4(ui_j) = strcmp(ui_What(ui_j),ui_B_What);
    ui_I_What ='I';
    ui_check_What5(ui_j) = strcmp(ui_What(ui_j),ui_I_What);
    
    
    if ui_check_What1(ui_j) == 1
        ini_FitWhat(ui_j) = 1;
    elseif ui_check_What2(ui_j) == 1
        ini_FitWhat(ui_j) = 2;
    elseif ui_check_What3(ui_j) == 1
        ini_FitWhat(ui_j) = 3;
    elseif ui_check_What4(ui_j) == 1
        ini_FitWhat(ui_j) = 4;
    elseif ui_check_What5(ui_j) == 1
        ini_FitWhat(ui_j) = 5;
    else
        break;
    end
    
    ini_FitVal(ui_j) = str2double(ui_Feat{3});
    ini_FitTime(ui_j) = str2double(ui_Feat{4})*2;
    ini_FitTime2(ui_j) = str2double(ui_Feat{5})*2;
    ini_FeatUpBound(ui_j-1) =str2double(ui_Feat{6});
    ini_FeatLowBound(ui_j-1) =str2double(ui_Feat{7});
    ui_answer_more    = ui_Feat{8};
    ui_yes_more       = 'y';
    ui_check_more = strcmp(ui_answer_more,ui_yes_more);
    
    if ui_check_more == 1
        ui_i = 0;
        ui_j = ui_j + 1;
    else
        ui_i =1;
        ini_NrFitVal = ui_j;
    end
end

 
%Paramters saves the number of paramters in each cell and the strings
%Astrocyte
ui_prompt_param_A = {
    'startpulse'
    'lengthpulse'
    'length1'
    %t2
    %t3
    'F_input'
    'alpha'
    'beta'
    'delta_t'
    'L_p'
    'X_k'
    'R_tot'
    
    'k_C'
    
    'VR_pa'
    'VR_ps'
    %'R_decay'
    %'K_p_min'
    
    'F'
    'R_g'
    'T'
    'g_K_k'
    'g_Na_k'
    'g_NBC_k'
    'g_KCC1_k'
    'g_NKCC1_k'
    'J_NaK_max'
    'G_BK_k'
    %Cinput
    'K_Na_k'
    'K_K_s'
    %'C_correction'
    %'A_ef_k'
    
    'g_Cl_k'
    'z_K'
    'z_Na'
    'z_Cl'
    'z_NBC'
    'v_6'
    'v_4'
    'psi_w'
    %'BK_end'
    %'K_ex'
    %'B_ex'
    %'K_G'
    %'v_5'
    };

 
[out_astrocyte,out_astrocyte_zero] = listdlg('PromptString','Astrocyte Parameter',...
    'SelectionMode','multiple',...
    'ListString',ui_prompt_param_A);

ini_NrPar(1)= length(out_astrocyte);
if out_astrocyte ~= 0;
    for r_i = 1 : ini_NrPar(1)
        ini_ParamString(r_i) = ui_prompt_param_A(out_astrocyte(r_i));
    end
else
    ini_NrPar(1) = 0;
end

%SMCEC
ui_prompt_param_S = {
    'gamma_i'
    'lambda_i'
    
    'C_m_j'
    'J_PLC'
    'J_0_j'
    
    'F_i'
    'K_r_i'
    'B_i'
    'c_b_i'
    'C_i'
    's_c_i'
    'c_c_i'
    'D_i'
    'v_d'
    'R_d_i'
    'L_i'
    'G_Ca_i'
    'v_Ca1_i'
    'v_Ca2_i'
    'R_Ca_i'
    'G_NaCa_i'
    'c_NaCa_i'
    'v_NaCa_i'
    'G_stretch'
    'alpha_stretch'
    'delta_p'
    'sigma_0'
    'E_SAC'
    'F_NaK_i'
    'G_Cl_i'
    'v_Cl_i'
    'G_K_i'
    'v_K_i'
    'F_KIR_i'
    'k_d_i'
    
    'F_j'
    'K_r_j'
    'B_j'
    'c_b_j'
    'C_j'
    's_c_j'
    'c_c_j'
    'D_j'
    'L_j'
    'G_cat_j'
    'E_Ca_j'
    'm_3_cat_j'
    'm_4_cat_j'
    'G_tot_j'
    'v_K_j'
    'c'
    'bb_j'
    'a_1_j'
    'a_2_j'
    'm_3b_j'
    'm_4b_j'
    'm_3s_j'
    'm_4s_j'
    'G_R_j'
    'v_rest_j'
    'k_d_j'
    
    'G_coup'
    'P_IP3'
    'P_Ca'
    
    'c_w_i'
    'beta_i'
    'v_Ca3_i'
    'R_K_i'
    'z_1'
    'z_2'
    'z_3'
    'z_4'
    'z_5'
    };

 
[out_smcec,out_smcec_zero] = listdlg('PromptString','SMCEC Parameter?',...
    'SelectionMode','multiple',...
    'ListString',ui_prompt_param_S);

ini_NrPar(2)= length(out_smcec);

if out_smcec ~= 0;
    for r_i = 1 : ini_NrPar(2)
        ini_ParamString(ini_NrPar(1) + r_i) = ui_prompt_param_S(out_smcec(r_i));
    end
else
    ini_NrPar(2) = 0;
end

%WallMechanics
ui_prompt_param_W = {
    'K_2'
    'K_3'
    'K_4'
    'K_5'
    'K_7'
    'gamma_cross'
    'eta'
    'R_0_passive'
    'P_T'
    'E_passive'
    'E_active'
    'alpha'};

[out_wall,out_wall_zero] = listdlg('PromptString','Wallmechanics Parameter?',...
    'SelectionMode','multiple',...
    'ListString',ui_prompt_param_W);

ini_NrPar(3)= length(out_wall);

if out_wall ~= 0;
    for r_i = 1 : ini_NrPar(3)
        ini_ParamString(ini_NrPar(1)+ini_NrPar(2) + r_i) = ui_prompt_param_W(out_wall(r_i));
    end
else
    ini_NrPar(3) = 0;
end

 
%individual boundaries for the different parameters
%in which space the initial parameter value should be changed - 0.8 means
%80% of the initial values 1.2 menas 120% of the initial value
%so a random number between 0.8 and 1.2 is created

%Checks if user wants individual boundaries for each parameter
ui_check_1 = {'Do you want individual parameter boundaries (y=yes/n=no)?'};
ui_title_check = 'In which percentage boundaries should the parameter change?';
ui_lines = 1;
ui_defaultanswer_check = {'n'};
ui_check_dlg = inputdlg(ui_check_1,ui_title_check,ui_lines,ui_defaultanswer_check);
ui_answer            = ui_check_dlg{1};
ui_yes               = 'y';
ui_check = strcmp(ui_answer,ui_yes);

%if the user wants individual boundaries:
%depending on the user input (which parameters)

ini_low_b = zeros(sum(ini_NrPar),1);
ini_up_b = zeros(sum(ini_NrPar),1);

if ui_check == 1
    for r_i = 1 : sum(ini_NrPar)
        ui_A = {'Enter lower boundary of randomisation of', ini_ParamString{r_i}};
        ui_B = {'Enter upper boundary of randomisation of', ini_ParamString{r_i}};
        ui_i_bound_vec{r_i*2-1} = strjoin(ui_A);
        ui_i_bound_vec{r_i*2}   = strjoin(ui_B);
    end
    
    ui_title_bound = 'Individual boundaries for the Parameter';
    ui_lines       = 1;
    
    for r_i      = 1 : sum(ini_NrPar)
        ui_defaultanswer_bound(r_i*2-1) = {'0.8'};
        ui_defaultanswer_bound(r_i*2)   = {'1.2'};
    end
    for r_i      = 1 : sum(ini_NrPar)
        ui_defaultanswer_bound_uneven(r_i*2-1) = {'1.2'};
        ui_defaultanswer_bound_uneven(r_i*2)   = {'0.8'};
    end
    ui_i_bound_1         = ui_i_bound_vec(1:sum(ini_NrPar));
    ui_bound_1           = inputdlg(ui_i_bound_1,ui_title_bound,ui_lines,ui_defaultanswer_bound,'on');
    ui_i_bound_2         = ui_i_bound_vec(sum(ini_NrPar)+1 : sum(ini_NrPar)*2);
    ui_title_bound_2     = 'Individual boundaries for the parameters';
    ui_lines     = 1;
    ui_defaultanswer_bound_2 = ui_defaultanswer_bound_uneven;
    ui_bound_2 = inputdlg(ui_i_bound_2,ui_title_bound_2,ui_lines,ui_defaultanswer_bound_2,'on');
    
    
    
    
    
    if mod(sum(ini_NrPar),2)
        for r_i = 1 : sum(ini_NrPar)/2 + 0.5
            ini_low_b(r_i)                              = str2double(ui_bound_1{r_i * 2 -1});
            ini_up_b(r_i + sum(ini_NrPar)/2 - 0.5)      = str2double(ui_bound_2{r_i * 2 -1});
        end
        for r_j = 1 : sum(ini_NrPar)/2 - 0.5
            ini_low_b(r_j + r_i)                        = str2double(ui_bound_2{r_j * 2});
            ini_up_b(r_j)                               = str2double(ui_bound_1{r_j * 2});
        end
    else
        for r_i = 1 : sum(ini_NrPar)/2
            ini_low_b(r_i)                              = str2double(ui_bound_1{r_i * 2 -1});
            ini_up_b(r_i)                               = str2double(ui_bound_1{r_i * 2});
            ini_low_b(sum(ini_NrPar)/2 + r_i)           = str2double(ui_bound_2{r_i * 2 -1});
            ini_up_b(sum(ini_NrPar)/2 + r_i)            = str2double(ui_bound_2{r_i * 2});
        end
        
    end
    
else
    
    ui_prompt_samebound = {'Enter lower boundary of randomisation for all parameters:'
        'Enter upper boundary of randomisation for all parameters:'};
    ui_title_samebound     = 'Individual boundaries for the parameters';
    ui_lines     = 1;
    ui_defaultanswer_samebound = { '0.8','1.2'};
    ui_samebound = inputdlg(ui_prompt_samebound,ui_title_samebound,ui_lines,ui_defaultanswer_samebound,'on');
    
    %saves inputs as double
    
    for r_i = 1 : sum(ini_NrPar)
        ini_low_b(r_i)           = str2double(ui_samebound{1});
        ini_up_b (r_i)           = str2double(ui_samebound{2});
    end
    
end

%General User Input
ui_prompt_general = {
    'Enter ultimate distribution of first generation'
    'Enter static population size (should be a binary number):'
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
    'Stopping criterion: Enter the convergence criteria number (How many generations should not have a change):'
    % How many "close" values have to be in a row, in order to make the
    % code stop
    'Stopping criterion: Enter the boundaries of convergence:'
    % Decision criteria if 2 values are close together or not
    % o.1 means if the distance between two following values of the bests of two generations of one fitness
    % value (i.e. c_k) is smaller than 0.1.
    % this is proved for all entered fitness values seperatly
    % how often this have to occur befor the code stops is determed by the
    % convergenz criteria number
    };

ui_title_general = 'General User Input';
ui_lines = 1;
ui_defaultanswer_general = {'0.80','16','0.8','0.5','0.1','10','15','0.02'};
ui_general = inputdlg(ui_prompt_general,ui_title_general,ui_lines,ui_defaultanswer_general,'on');
ini_ulti_dist          = 1-str2double(ui_general{1});
ini_static_popsize     = str2double(ui_general{2});
ini_co_p               = str2double(ui_general{3});
ini_mu_p               = str2double(ui_general{4});
ini_el_k               = str2double(ui_general{5});
ini_max_time           = str2double(ui_general{6})*60;
ini_nr_convergence     = str2double(ui_general{7});
ini_bound_convergence  = str2double(ui_general{8});

%% settings and calculations for stopping criteria
%time
stop_time_beginning = clock;
%Start time of simulation
stop_time_elapsed = 0;
%Just a preallocation  
ini_max_time_1           = str2double(ui_general{6});
t=now;
DT = datestr(t);
t1 = t + ini_max_time_1/1440;
DT1 = datestr(t1);
msgstring = {'The code will be certainly finished at', DT1};
msgstring1 = strjoin(msgstring);
h = msgbox(msgstring1,'Finish time','help');

 
%% running variables
mpar_1      = 1;
mpar_c      = 1;
mu_var      = 1;

conv_crit   = 1;
r_simu_2    = 1;
ini_w       = 1;

%% predefined vectors/matrices/variables apart form the user input
T = nv.T;
ini_popsize = 4;

%% Preallocation 1
ini_randMat = zeros(ini_popsize,sum(ini_NrPar)); %this preallocation reduces the simulation time

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
for r_i = 1:sum(ini_NrPar)
    ini_randMat(:,r_i) = vertcat(1,(ini_up_b(r_i)-ini_low_b(r_i))*rand(ini_popsize-1,1)+ini_low_b(r_i));
    %creats an random number in the user defined boundaries
end
% ini_randMat matrix [n x m] (and n = number of organisms
% m = parameters (e.g first column k_on)

%% uniform distribution
pop_w = 0;

while pop_w == 0;
    % calculates the distribution of the first generation
    for r_i = 1 : sum(ini_NrPar)
        pop_c = ((ini_up_b(r_i)+ini_low_b(r_i))/2);
        
        for r_j   = 1 : ini_popsize
            pop_d   = pop_c  - ini_randMat(r_j,r_i);
            pop_vector_1(r_j,r_i) = pop_d;
        end
        
        pop_p(r_i)    = abs(sum(pop_vector_1(:,r_i))/ini_popsize);
        pop_p_1(r_i)  = (pop_p(r_i) / abs(pop_c - ini_up_b(r_i)));
    end
    
    pop_res = (sum(pop_p_1) / sum(ini_NrPar));
    
    % checks if the population is bigger than the user input
    if pop_res > ini_ulti_dist
        ini_popsize = ini_popsize*2;
        ini_randMat = zeros(ini_popsize,sum(ini_NrPar));
        for r_i = 1:sum(ini_NrPar)
            ini_randMat(:,r_i) = vertcat(1,(ini_up_b(r_i)-ini_low_b(r_i))*rand(ini_popsize-1,1)+ini_low_b(r_i)); %creates new random parametersets
        end
        
    else
        pop_w = 1; %jumps out of the loop if the gerneration is uniform distributed
    end
end
vector_ini_popsize(r_simu_2) = ini_popsize;

%% Preallocation 2
%depending on the population size
%reduces the simulation time
results             = zeros(ini_popsize,ini_NrFitVal);
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

parent_matrix(:,4+ini_NrFitVal:3+ini_NrFitVal+sum(ini_NrPar))= ini_randMat;
% rand matrix is now inclued in the parent_matrix (parent_matrix is now
% constantly used)

%% simulation of the first generation
for r_j = 1:ini_popsize
    
    %changes the parameter values in the astrocyte for one organism
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the smooth muscle cell for one organism
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the wall mechanics for one organism
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %simulates the NVU model with the changed parameter values
    nv.simulate()
   r_n = 1;
    %saves the results
    for r_i = 1 : ini_NrFitVal
        results_1(r_i,:)  = nv.out(ini_FitString{r_i}); %all timesteps of one value for one parameterset
        if ini_FitWhat(r_i) == 1
            results(r_j,r_i) = results_1(r_i,ini_FitTime(r_i));% one timestep
            initial_result = results(1,1);
        elseif ini_FitWhat(r_i) == 2
            results(r_j,r_i) = mean(results_1(r_i,ini_FitTime(r_i):ini_FitTime2(r_i)));
            initial_result = results(1,1);
        elseif ini_FitWhat(r_i) ==3
            results(r_j,r_i) = (results_1(r_i,ini_FitTime2(r_i))-results_1(r_i,ini_FitTime(r_i)));
            initial_result = results(1,1);
        elseif ini_FitWhat(r_i) == 4
            
            if results_1(r_i,ini_FitTime(r_i))< ini_FeatUpBound(r_n) && results_1(r_i,ini_FitTime(r_i)) > ini_FeatLowBound(r_n)
                results(r_j,r_i) = 1;
            else
                results(r_j,r_i) = 0;
            end
            r_n = r_n+1;
        elseif ini_FitWhat(r_i) == 5
            results_2 = results_1(r_i,ini_FitTime(r_i):ini_FitTime2(r_i));
            for r_k=1:length(results_2)
                if results_2(r_k)< ini_FeatUpBound(r_n) && results_2(r_k) > ini_FeatLowBound(r_n)
                    results_2(r_k) = 1;
                else
                    results_2(r_k) = 0;
                end                
            end
            r_n = r_n+1;
            if length(results_2) == sum(results_2(:))
                results(r_j,r_i) = 1;
            else
                results(r_j,r_i) = 0;
            end
        end
    end
    if ini_min_max == 0
        results(:,1) = 1./results(:,1);
    end
end

 
%saves all results in the parent matrix
parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(:,1:sum(ini_NrFitVal));

%% calculation of the survival probabilities and cumultative probabilities
sel_sumFit = sum(parent_matrix(:,4));
sel_prob_2 = zeros(length(parent_matrix(:,1)),ini_NrFitVal);

for r_i = 1:length(parent_matrix(:,1))
    sel_prob_2(r_i,1) = parent_matrix(r_i,4)./sel_sumFit;
end

 
for r_i = 2 : ini_NrFitVal
    if ini_FitWhat(r_i) == 1 || ini_FitWhat(r_i)== 2 || ini_FitWhat(r_i) == 3
    for r_j = 1:ini_popsize
        sel_prob(r_j,r_i) = (abs(ini_FitVal(r_i) - parent_matrix(r_j,r_i+3)));
        % difference between current value and aimed value
        
        if sel_prob(r_j,r_i) == 0;
            sel_prob(r_j,r_i) = 0.0000000000000000001;
            % if one fitness value reaches the aimed value the code doesn't
            % work anymore, because it can't divivde trough 0 so the value
            % is nearly settet to 0
        end
    end
    
    sel_sum(:,r_i) = sum(sel_prob(:,r_i));
    % sum of all distances for one fitness value
    
    for r_j = 1 : ini_popsize
        sel_prob_1(r_j,r_i) = sel_sum(r_i) ./ sel_prob(r_j,r_i);
        %sum of distances divided by every single distance
    end
    
    sel_sum_2(:,r_i) = sum(sel_prob_1(:,r_i));
    % sum of (sum of distances divided by every single distance)
    
    for r_j = 1 : ini_popsize
        sel_prob_2(r_j,r_i) = sel_prob_1(r_j,r_i) ./ sel_sum_2(:,r_i);
        % (sum of distances divided by every single distance)divided by
        %(sum of (sum of distances divided by every single distance))
    end
    else
        if sum(parent_matrix(:,3+r_i)) == 0 
            sel_prob_2(:,r_i)= 0;
        else
        sel_prob_2(:,r_i)= parent_matrix(:,r_i+3)./sum(parent_matrix(:,r_i+3));
        end
    end
    %survival probability for every organism of one generation for ONE
    %Fitness value is calculated
end

%average survival probability
for r_i = 1 : ini_popsize
    sel_prob_vec_1(r_i,:) = sum(sel_prob_2(r_i,:));
    sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
end
% average survival probability for every organism of one generation
% calculated by the average of the survival probabilities of each fitness
% value

%cumulative probability
sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
parent_matrix (:,1:3) = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec, sel_prob_vec]; %shape like parent_matrix

%% save best resut of the generation
stop_max = max(parent_matrix(:,3)); % best survival probability -> best member of the generation
stop_i = 1;
while parent_matrix(stop_i,3) < stop_max
    stop_i = stop_i+1;
end
best_results(r_simu_2,:) = parent_matrix(stop_i,:);

 
 
 
 
 
%% Loop for second and following generations

while  stop_time_elapsed < ini_max_time && stop_count(1)< ini_nr_convergence %&& stop_ac_count < ini_NrFitVal
    
    %Reset all running variables so the whole process can run again
    sel_p       = 1;
    
    mpar_1      = 1;
    mpar_c      = 1;
    
    co_var      = 1;
    
    
    
    %% Elitism 1
    clear el_p;
    clear el_p_1;
    el_max(r_simu_2)  = max(parent_matrix(:,3)) * (1-ini_el_k);
    %Defines the threshold for being elite. Template is the best organism
    %of current generation
    r_j = 1;
    for r_i = 1:length(parent_matrix(:,1))
        
        if parent_matrix(r_i,3) >= el_max(r_simu_2)
            %Checks each line if it's bigger than threshold for elitism if it
            %is, this line gets safed in order to overwrite it in the parent
            %matrix
            el_p_1(r_j,:) = parent_matrix(r_i,:);
            r_j = r_j + 1;
        end
    end
    if r_j >=4
        el_nr = 3;
    else
        el_nr = r_j-1;
    end
    %el_p = zeros(length(el_p_1)-el_nr : length(el_p_1),1:length(parent_matrix(1,:)));
    el_p_1 = sortrows(el_p_1,3);
    el_p(:,:) = el_p_1(length(el_p_1(:,1))-el_nr+1 : length(el_p_1(:,1)),:);
    
    
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
            if sel_rand_value(sel_2,1) < parent_matrix(sel_p,2);
                % if random value at the position n is smaler as the cummultative probability ...
                sel_output(sel_2,:) = parent_matrix(sel_p,:);
                % ... save the individum in the sel_output matrix
            end
            
            while sel_rand_value(sel_2,1) > parent_matrix(sel_p,2);
                % while the random number is bigger ....
                sel_p=sel_p+1;
                if sel_rand_value(sel_2,1) <= parent_matrix(sel_p,2);
                    % test again if the cummultative probability of the next organism
                    % (parameterset) is bigger than the random number
                    sel_output(sel_2,:) = parent_matrix(sel_p,:);
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
    
    vector_ini_popsize(r_simu_2+1) = ini_popsize;
    
    
    %% Multiplicate Parents
    if length(sel_output(:,1)) <= (ini_static_popsize*0.5)
        %Waits until population has size: ini_static_popsize/2 and
        %will resize it to ini_static_popsize
        for r_i = 1:sel_2/2
            
            mpar_P1 = sel_output(mpar_1,:);
            mpar_P2 = sel_output(mpar_1+1,:);
            %Taking out two parents (are always located close to each other)
            
            for r_j = 1:ini_popsize/(ini_popsize / 2);
                %Repeat in (in this case) two times in order to achieve initial
                %population size.
                
                parent_matrix(mpar_c,:)   = mpar_P1;
                parent_matrix(mpar_c+1,:) = mpar_P2;
                
                mpar_c = mpar_c + 2;
                %Always a pair of parents
            end
            mpar_1 = mpar_1 + 2;
        end
        
        
    else
        clear parent_matrix
        parent_matrix = sel_output(:,:);
        ini_popsize = ini_popsize / 2;
        sel_prob_2(ini_popsize+1:length(sel_prob_2(:,1)),:) =[];
        sel_prob(ini_popsize+1:length(sel_prob(:,1)),:) =[];
        sel_prob_1(ini_popsize+1:length(sel_prob_1(:,1)),:) =[];
        sel_prob_vec_1(ini_popsize+1:length(sel_prob_vec_1(:,1)),:) =[];
        sel_prob_vec(ini_popsize+1:length(sel_prob_vec(:,1)),:) =[];
        sel_cum_probvec(ini_popsize+1:length(sel_cum_probvec(:,1)),:) =[];
    end
    
    
    %% Mutation
    p_mu_1 = zeros(1,length(parent_matrix(:,1)));
    p_mu = zeros(length(parent_matrix(:,1)),1);
   
    
    for r_i = 1:length(parent_matrix(:,1))
        r_bound = 0; 
        %Mutation occurs for ALL organism in ONE generation
        
        mu_max      = max(parent_matrix(:,4))*1.02;
        mu_mean     = mean(parent_matrix(:,4));
        mu_current  = parent_matrix(r_i,4);
        
        p_mu_1(1,r_i) =  ini_mu_p * abs(((mu_max - mu_current) / (mu_max - mu_mean)));
        
        if p_mu_1(1,r_i) > ini_mu_p
            p_mu_1(1,r_i) =  ini_mu_p;
        end
        
        
        for r_k = 2 : ini_NrFitVal
            if ini_FitWhat(r_k) == 1 || ini_FitWhat(r_k) == 2 || ini_FitWhat(r_k) == 3
            %Mutation rate is depended on how close the fitness value(s)
            %of the current observed organism are to
            %the user defined threshold(s). The closer they are, the smaller the
            %mutation probability. The average of all fitness values is used to
            %determine p_mu
            
            mu_max      = ini_FitVal(r_k);
            mu_current  = parent_matrix(r_i,3+r_k);
            % mu_mean = mean(parent_matrix(:,3+r_k));
            
            p_mu_1(r_k,r_i) =  ini_mu_p * abs(((mu_max - mu_current) / (mu_max)));
            %p_mu_1(r_k,r_i) =  ini_mu_p * abs(((mu_max - mu_current) / (mu_max-mu_mean)));
            if p_mu_1(r_k,r_i) > ini_mu_p
                p_mu_1(r_k,r_i) =  ini_mu_p;
            end
            else           
                r_bound = r_bound + 1;
                p_mu_1(r_k,r_i) = 0;
            end
        end
        
        
        
        for r_j = 1 : sum(ini_NrPar)
            %For every Parameter the GA decides whether Mutation occurs or not
            
            if mod(r_simu_2,30)
                p_mu(r_i) = sum(p_mu_1(:,r_i))/(length(p_mu_1(:,r_i)) - r_bound);
            else
                p_mu(r_i) = ini_mu_p;
            end
            
            %As already described the mutation rate p_mu is depended on all
            %fitness values. Therefore a mean is taken
            
            mu_rand = rand(1);
            %Random number between 0 and 1
            
            if mu_rand <= p_mu(r_i)
                %Here is the decision if mutation happens or not. The higher
                %p_mu the higher the prob. that mu_rand is bigger than p_mu.
                
                %Refers to the position of the current parameter and
                %changes there the the parameter within in the
                %predefined boundaries and overwrites it in the
                %parentmatrix
                
                mu_r = ini_low_b(r_j) + (ini_up_b(r_j) - ini_low_b(r_j))*rand(1,1);
                parent_matrix(r_i,r_j + 3 + ini_NrFitVal) = mu_r;
                
            end
        end
        
    end
    
    
    mu_p_total(r_simu_2) = mean(p_mu);
    %Average mutation rate of current generation. This is used in plots.
    
    %% Cross-over
    co_s = length(parent_matrix(:,1));
    r_j = 1;
    co_p_1 = zeros(1,co_s/2);
    co_p_2 =  zeros(1,co_s/2);
    co_p = zeros(co_s/2,1);
    co_bound = 0;
    
    while co_var      <= (co_s-1)
        %Repeat until end of generation
        
        %Calculate c.o. prob. for one organism. In order to do so the code
        %has to calculate for every fitness value the distance between fitness value and user
        %defined threshold. The c.o. prob. is resulting of the average of
        %these distances.
        
        cr_max           = max(parent_matrix(:,4))*1.02;
        cr_mean          = mean(parent_matrix(:,4));
        cr_current_1     = parent_matrix(co_var,4);
        cr_current_2     = parent_matrix(co_var+1,4);
        
        %in crossing over the "quality" of both parents is important to
        %decide if it occurs or not.
        
        co_p_1(1,r_j) =  abs(((cr_max - cr_current_1) / (cr_max - cr_mean)));
        co_p_2(1,r_j) =  abs(((cr_max - cr_current_2) / (cr_max - cr_mean)));
        
        if co_p_1(1,r_j) > ini_co_p
            co_p_1(1,r_j) = ini_co_p;
        end
        if co_p_2(1,r_j) > ini_co_p
            co_p_2(1,r_j) = ini_co_p;
        end
        
        for r_k = 2 : ini_NrFitVal
            if ini_FitWhat(r_k) == 1 || ini_FitWhat(r_k) == 2 || ini_FitWhat(r_k) == 3
            cr_max           = ini_FitVal(r_k);
            cr_current_1     = parent_matrix(co_var,3+r_k);
            cr_current_2     = parent_matrix(co_var+1,3+r_k);
            
            %in crossing over the "quality" of both parents is important to
            %decide if it occurs or not.
            
            co_p_1(r_k,r_j) = ini_co_p * abs(((cr_max - cr_current_1) / (cr_max)));
            co_p_2(r_k,r_j) = ini_co_p * abs(((cr_max - cr_current_2) / (cr_max)));
            else
                co_bound = co_bound + 1;
                co_p_1(r_k,r_j) = 0;
                co_p_2(r_k,r_j) = 0;
            end
        end
        if co_p_1(r_k,r_j) > ini_co_p
            co_p_1(r_k,r_j) = ini_co_p;
        end
        if co_p_2(r_k,r_j) > ini_co_p
            co_p_2(r_k,r_j) = ini_co_p;
        end
        
        co_mean_1 = sum(co_p_1(:,r_j))/(length(co_p_1(:,r_j))-co_bound);
        co_mean_2 = sum(co_p_2(:,r_j))/(length(co_p_1(:,r_j))-co_bound);
        co_p_1 = co_mean_1 + co_mean_2;
        co_p(r_j)  = mean(co_p_1);
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
        r_j = r_j + 1;
        co_var = co_var + 2;
        
    end
    %output is a matrix which is 2*mpar_rbigger than matrix_input
    
    %% Elitism 2
    
    for r_j = 1:length(el_p(:,1));
        %Searches line in parent_matrix and overwrite line by elite
        %line (if not already existing
        parent_matrix(length(parent_matrix(:,1))-r_j+1,:) = el_p(r_j,:);
    end
    
    
    
    
    %% simulation
    for r_j = 1:ini_popsize
        
        %changes the parameter values in the astrocyte for one organism
        for r_i = 1 : ini_NrPar(1)
           nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
        end
        
        %changes the parameter values in the smooth muscle cell for one organism
        for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
            nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
        end
        
        %changes the parameter values in the wall mechanics for one organism
        for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
            nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * parent_matrix(r_j,r_i + 3 + ini_NrFitVal);
        end
        
        %simulates the NVU model with the changed parameter values
        nv.simulate()
        r_n = 1;
        %saves the results
        for r_i = 1 : ini_NrFitVal
            results_1(r_i,:)  = nv.out(ini_FitString{r_i}); %all timesteps of one value for one paraemterset
            if ini_FitWhat(r_i) == 1
                results(r_j,r_i) = results_1(r_i,ini_FitTime(r_i));% one timestep
            elseif ini_FitWhat(r_i) == 2
                results(r_j,r_i) = mean(results_1(r_i,ini_FitTime(r_i):ini_FitTime2(r_i)));
            elseif ini_FitWhat(r_i) ==3
                results(r_j,r_i) = (results_1(r_i,ini_FitTime2(r_i))-results_1(r_i,ini_FitTime(r_i)));
            elseif ini_FitWhat(r_i) == 4
                
                if results_1(r_i,ini_FitTime(r_i))< ini_FeatUpBound(r_n) && results_1(r_i,ini_FitTime(r_i)) > ini_FeatLowBound(r_n)
                    results(r_j,r_i) = 1;
                else
                    results(r_j,r_i) = 0;
                end
                r_n = r_n+1;
            elseif ini_FitWhat(r_i) == 5
                results_2 = results_1(r_i,ini_FitTime(r_i):ini_FitTime2(r_i));
                
                for r_k = 1:length(results_2)
                    if results_2(r_k)< ini_FeatUpBound(r_n) && results_2(r_k) > ini_FeatLowBound(r_n)
                        results_2(r_k) = 1;
                    else
                        results_2(r_k) = 0;
                    end
                end
                r_n = r_n+1;
                if length(results_2) == sum(results_2(:))
                    results(r_j,r_i) = 1;
                else
                    results(r_j,r_i) = 0;
                end
            end
        end
        if ini_min_max == 0
            results(:,1) = 1./results(:,1);
        end
    end
    if ini_min_max == 0
        results(:,1) = 1./results(:,1);
    end
    %saves all results in the parent matrix
    parent_matrix(:,4:3+sum(ini_NrFitVal)) = results(1:length(parent_matrix(:,1)),1:sum(ini_NrFitVal));
        
    %% calculation of the survival probabilities and cumultative probabilities
   sel_sumFit = sum(parent_matrix(:,4));
    for r_i = 1:length(parent_matrix(:,1))
        sel_prob_2(r_i,1) = parent_matrix(r_i,4)./sel_sumFit;
    end
    for r_i = 2 : ini_NrFitVal
        if ini_FitWhat(r_i) == 1 || ini_FitWhat(r_i)== 2 || ini_FitWhat(r_i) == 3
            for r_j = 1:ini_popsize
              sel_prob(r_j,r_i) = (abs(ini_FitVal(r_i) - parent_matrix(r_j,r_i+3)));
                % difference between current value and aimed value
                
                if sel_prob(r_j,r_i) == 0;
                    sel_prob(r_j,r_i) = 0.0000000000000000001;
                    % if one fitness value reaches the aimed value the code doesn't
                    % work anymore, because it can't divivde trough 0 so the value
                    % is nearly settet to 0
                end
            end
            
            sel_sum(:,r_i) = sum(sel_prob(:,r_i));
            % sum of all distances for one fitness value
            
            for r_j = 1 : ini_popsize
                sel_prob_1(r_j,r_i) = sel_sum(r_i) ./ sel_prob(r_j,r_i);
                %sum of distances divided by every single distance
            end
            
            sel_sum_2(:,r_i) = sum(sel_prob_1(:,r_i));
            % sum of (sum of distances divided by every single distance)
            
            for r_j = 1 : ini_popsize
                sel_prob_2(r_j,r_i) = sel_prob_1(r_j,r_i) ./ sel_sum_2(:,r_i);
                % (sum of distances divided by every single distance)divided by
                %(sum of (sum of distances divided by every single distance))
            end
        else
            sel_prob_2(:,r_i)= parent_matrix(:,r_i+3)./sum(parent_matrix(:,r_i+3));
        end
    %survival probability for every organism of one generation for ONE
    %Fitness value is calculated
    end
    
    %average survival probability
    for r_i = 1 : ini_popsize
        sel_prob_vec_1(r_i,:) = sum(sel_prob_2(r_i,:));
        sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
    end
    % average survival probability for every organism of one generation
    % calculated by the average of the survival probabilities of each fitness
    % value
    
    %cumulative probability
    sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
    parent_matrix (:,1:3) = [linspace(1,ini_popsize,ini_popsize)', sel_cum_probvec(1:length(parent_matrix(:,1))), sel_prob_vec(1:length(parent_matrix(:,1)))]; %shape like parent_matrix
    
    
    stop_max = max(parent_matrix(:,3));
    stop_i = 1;
    
    while parent_matrix(stop_i,3) < stop_max
        stop_i = stop_i+1;
    end
    
    stop_max_line = parent_matrix(stop_i,:);
    best_results(r_simu_2+1,:) = stop_max_line(:,:);
    
    %% settings and calculations for stopping criteria
    
    % convergence
    
    %if more best results in a row don't have a difference in EVERY fitness
    %value stop the loop
  
    
    for r_i = 1:ini_NrFitVal
        
    stop_1 = best_results(r_simu_2-stop_count(r_i),r_i+3)+ best_results(r_simu_2-stop_count(r_i),r_i+3)*ini_bound_convergence;
    stop_2 = best_results(r_simu_2-stop_count(r_i),r_i+3)- best_results(r_simu_2-stop_count(r_i),r_i+3)*ini_bound_convergence;
        % check every fitness value if it hasn't changed since the last
        % generation
        if best_results(r_simu_2+1,3+r_i)< stop_1 && best_results(r_simu_2+1,3+r_i)>stop_2 
            stop_count(r_i) = stop_count(r_i)+1;
            %if it hasn't change, count one more
        else
            stop_count(r_i) = 0;
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

%plot best results over the generation
if ini_min_max == 0
    best_results_plot1 = 1./best_results(:,4);
else
    best_results_plot1 = best_results(:,4);
end
initial_result_plot(1:r_simu_2) = initial_result;

figure(1)
plot(linspace(1,r_simu_2,r_simu_2)', best_results_plot1(:))
plot_str_1 = {'Best',ini_FitString{1}, 'value of each generation'};
plot_strjoin_1 = strjoin(plot_str_1);
title(plot_strjoin_1);
xlabel('# of generation');
ylabel(ini_FitString{1});
hold on;
grid on;

if ini_FitWhat(1) == 1 || ini_FitWhat(1) == 2 || ini_FitWhat(1) == 3
    plot(linspace(1,r_simu_2,r_simu_2)', initial_result_plot(:));
    legend(plot_strjoin_1,'initial value')
end

plot_k = 2;

if ini_NrFitVal > 1
    plot_w = 1;
    plot_j = 2;
    
    while plot_w ==1
        plot_l = 1;
        for r_i = plot_j : plot_j+3
            plot_inifit_vec = linspace(ini_FitVal(r_i),ini_FitVal(r_i),r_simu_2);
            %Creates the plot line for the user defined theshold
            
            if ini_FitWhat(plot_j) == 1 || ini_FitWhat(plot_j) == 2 || ini_FitWhat(plot_j) == 3            
                figure(plot_k)
                subplot(2,2,plot_l)
                plot(linspace(1,r_simu_2,r_simu_2)', best_results(:,r_i+3))
                plot_str_1 = {'Best',ini_FitString{r_i}, 'value of each generation'};
                plot_strjoin_1 = strjoin(plot_str_1);
                %Makes legend adaptive to different fitness value names
                legend(plot_strjoin_1,'User defined value');
                title(plot_strjoin_1);
                xlabel('# of generation');
                ylabel(ini_FitString{r_i});
                hold on;
                grid on;
            
 
                subplot(2,2,plot_l)
                plot(linspace(1,r_simu_2,r_simu_2)',plot_inifit_vec);
                plot_str_2 = {'Best',ini_FitString{r_i}, 'value of each Generation'};
                plot_strjoin_2 = strjoin(plot_str_2);
                %Makes legend adaptive to different fitness value names
                legend(plot_strjoin_2,'User defined value');
                grid on;
                grid minor;
            end
            plot_l = plot_l + 1;
            
            if r_i >= ini_NrFitVal
                %Check if all fitness values are already plotted or not, if yes
                %loop stops
                plot_w = 0;
                break
            end
            
            
        end
        
        plot_k = plot_k + 1;
        plot_j = plot_j + 4;
    end
    
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
% A ne preallocation is needed, because the size has changed (now 1st
% generation is also in the matrix)

sel_sum_prob        = zeros(best_l,ini_NrFitVal);
sel_prob            = zeros(best_l,ini_NrFitVal);
sel_sum             = zeros(1,ini_NrFitVal);
sel_prob_1          = zeros(best_l,ini_NrFitVal);
sel_prob_2          = zeros(best_l,ini_NrFitVal);
sel_prob_vec_1      = zeros(best_l,1);
sel_prob_vec        = zeros(best_l,1);

%% select the best organisms of the best organisms of every generation
sel_sumFit = sum(best_results(:,4));
for r_i = 1:length(best_results(:,1))
    sel_prob_2(r_i,1) = best_results(r_i,4)./sel_sumFit;
end

for r_i = 2 : ini_NrFitVal
    
    for r_j = 1:best_l
        sel_prob(r_j,r_i)       = (abs(ini_FitVal(r_i) - best_results(r_j,r_i+3)));
        if sel_prob(r_j,r_i) == 0;
            sel_prob(r_j,r_i) = 0.0000000000000000001;
        end
    end
    
    
    sel_sum(:,r_i) = sum(sel_prob(:,r_i));
    
    for r_j = 1 : best_l
        sel_prob_1(r_j,r_i) = sel_sum(r_i) ./ sel_prob(r_j,r_i);
    end
    
    sel_sum_2(:,r_i) = sum(sel_prob_1(:,r_i));
    
    for r_j = 1 : best_l
        sel_prob_2(r_j,r_i) = sel_prob_1(r_j,r_i) ./ sel_sum_2(:,r_i);
    end
end

for r_i = 1 : best_l
    sel_prob_vec_1(r_i,:) = sum(sel_prob_2(r_i,:));
    sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
end

sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
best_results  = [best_results(:,1), sel_cum_probvec, sel_prob_vec, best_results(:,4:sum(ini_NrPar)+3+ini_NrFitVal)];

if ini_min_max == 0
    best_results_plot1 = 1./best_results(:,4);
end

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
r_i = 1;
plot_j = 1;
best_of_best(1,:) = plot_max(best_l,:);
plot_max_unique = zeros(7,length(best_results(1,:)));

while plot_w ==1 && plot_j <=7 && r_i < best_l
    if plot_max(best_l+1-r_i,3) == plot_max(best_l-r_i,3)
        plot_w=1;
    else
        plot_max_unique(plot_j,:) = plot_max(best_l+1-r_i,:);
        plot_j = plot_j + 1;
    end
    r_i = r_i + 1;
end
figure(plot_k + 2)
for r_i = 1 : 7
    plot_max_line_1 = (-best_of_best(1,:) + plot_max_unique(r_i,:));
    x = 1:length(ini_ParamString);
    y_1 = plot_max_line_1(4+ini_NrFitVal : 3+ini_NrFitVal+sum(ini_NrPar));
    plot(x, y_1, color_string{r_i})
    hold on;
end
legend(result_string);
title('Distance from best organism [%]');
set(gca, 'XTick',1:length(ini_ParamString), 'XTickLabel',ini_ParamString)
grid on;
grid minor;

 
%% simulation with the inital values an the best one
comp_ini_best =ones(2,3+ini_NrFitVal+sum(ini_NrPar));
comp_ini_best(2,:) = plot_max(r_simu_2,:);
u = 1;
for r_j = 1:2
    
    %changes the parameter values in the astrocyte for one organism
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the smooth muscle cell for one organism
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the wall mechanics for one organism
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %simulates the NVU model with the changed parameter values
    nv.simulate()
   
    %saves the results
    
    for r_i = 1 : ini_NrFitVal
        results_ini(u,:)  = nv.out(ini_FitString{r_i});%all timesteps of one value for one paraemter set
        u=u+1;
    end
end
%% Plot the initial and the best parameterset over the time
figure(plot_k+3)
for r_i = 1:ini_NrFitVal
    subplot(2,2,r_i)
    plot(nv.T,results_ini(r_i,:))
    hold on
    plot(nv.T,results_ini(r_i+ini_NrFitVal,:))
    plot_str_3 = {ini_FitString{r_i},'of the initial parameterset and the best one over the time'};
    plot_strjoin_3 = strjoin(plot_str_3);
    %Makes legend adaptive to different fitness value names
    
    legend('Initial parameter set','Best parameter set');
    title(plot_strjoin_3);
    xlabel('Timesteps');
    plot_str_4 = ini_FitString(r_i);
    plot_strjoin_4 = strjoin(plot_str_4);
    ylabel(plot_strjoin_4 );
    grid on;
    grid minor;
end

%Plot more output values
    
    %List of output varaibles in the GA
    ui_prompt_plotfeat = {
        'R_k'
        'N_K_s'
        'N_Na_s'
        'N_HCO3_s'
        'N_K_k'
        'N_Na_k'
        'N_HCO3_k'
        'N_Cl_k'
        'w_k'
        'K_p'
        'Ca_i'
        's_i'
        'v_i'
        'w_i'
        'I_i'
        'K_i'
        'Ca_j'
        's_j'
        'v_j'
        'I_j'
        'Mp'
        'AMp'
        'AM'
        'M'
        'F_r'
        'R'};
    
    [out_plotfeat,out_plotfeat_zero] = listdlg('PromptString','Plot more output variables?',...
        'SelectionMode','multiple',...
        'ListString',ui_prompt_plotfeat);
    
    ini_outputplots(1)= length(out_plotfeat);
    
    plot_output = 1;
    plot_j = 2;
    
    if out_plotfeat ~= 0;
        for r_i = 1 : ini_outputplots(1)
            plot_output_variables(r_i) = ui_prompt_plotfeat(out_plotfeat(r_i));
        end
    else
        plot_output = 0;
    end 
    
    if plot_output ==1
        plot_l = 1;
       
        u = 1;
for r_j = 1:2
    
    %changes the parameter values in the astrocyte for one organism
    for r_i = 1 : ini_NrPar(1)
        nv.astrocyte.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the smooth muscle cell for one organism
    for r_i = ini_NrPar(1)+1 : ini_NrPar(1)+ini_NrPar(2)
        nv.smcec.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %changes the parameter values in the wall mechanics for one organism
    for r_i = ini_NrPar(1)+ini_NrPar(2) + 1 : sum(ini_NrPar)
        nv.wall.params.(ini_ParamString{r_i})= ini_value(r_i) * comp_ini_best(r_j,r_i + 3 + ini_NrFitVal);
    end
    
    %simulates the NVU model with the changed parameter values
    nv.simulate()
    
    %saves the results
    
    for r_i = 1 : ini_outputplots(1)
        results_plotfeat(u,:)  = nv.out(plot_output_variables{r_i});%all timesteps of one value for one paraemter set
        u=u+1;
    end
end
        
        for r_i = 1 : ini_outputplots(1)

            figure(plot_k+3+r_i)
            plot(nv.T, results_plotfeat(r_i,:))
            plot_str_1 = {'Result of',plot_output_variables{r_i},'due to optimization'};
            hold on;
            grid on;
            
            plot(nv.T, results_plotfeat(r_i+ini_outputplots(1),:))
            plot_str_1 = {'Result of',plot_output_variables{r_i},'due to optimization'};
            plot_strjoin_1 = strjoin(plot_str_1);
            %Makes legend adaptive to different fitness value names
            legend('initial value', plot_strjoin_1);
            title(plot_strjoin_1);
            xlabel('Time [s]');
            ylabel(plot_output_variables{r_i});
            hold on;
            grid on;
        
        end
        
    end

 
 
 
 


