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

%% Determine the initial condition - User Input -> ui...

ui_choice = menu('In which cell do you like to run a parameter sensifity test?','SMCEC','WallMechanisms','Astrocyte');

ui_prompt_param = { 'How many parameter?'
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
                    'Parameter17:'
                    'Parameter18:'
                    'Parameter19:'
                    'Parameter20:'
    };                                                                      % At which timepoint the user wants to know the calcium concentration

ui_title_param = 'Which parameters do you want to use?';

ui_lines_param = 1;

if ui_choice == 1
    ui_defaultanswer_param  = {'5','gamma_i','lambda_i','C_m_j','J_PLC','J_0_j','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0'};
end

if ui_choice == 2
    ui_defaultanswer_param  = {'7','K_2','K_3','K_4','K_5','K_7','gamma_cross','R_0_passive','0','0','0','0','0','0','0','0','0','0','0','0','0'};
end

if ui_choice == 3
    ui_defaultanswer_param  = {'16','k_deg','k_on','K_inh','r_h','K_G','K_I','BK_end','K_ex','B_ex','J_max','K_act','P_L','V_max','k_pump','delta','VR_ER_cyt','0','0','0','0'};                 % Just a default answer
end

 
ui_param                = inputdlg(ui_prompt_param,ui_title_param,ui_lines_param,ui_defaultanswer_param);

ini_NrPar               = str2double(ui_param{1});

for ui_i = 1:ini_NrPar
    if ui_choice == 1
        ui_a            = nv.smcec.params;
        ini_value(ui_i) = getfield(ui_a, ui_param{ui_i+1});
    end
    if ui_choice == 2
        ui_a            = nv.wall.params;
        ini_value(ui_i) = getfield(ui_a, ui_param{ui_i+1});
    end
    if ui_choice == 3
        ui_a            = nv.astrocyte.params;
        ini_value(ui_i) = getfield(ui_a, ui_param{ui_i+1});
    end
end

 
ui_prompt_1 = { 'Enter size of population (if divided by number of organism, result must be out of Z:'      % This population size is constant
    'Enter the initial probability of cross-over:'                                              % This value determines the frequency of co
    'Enter the threshold for elitism:'                                                          % 0.005 means that a organism must have a bigger c_k than 99.5% of population to be an elite
    'Enter the initial probability of mutation:'                                                % Mutation probability
    'Do you want individual parameterboundaries [y] [n]'                                        % If User wants to adjust parameter boundaries
    'Enter the minimal generation number:'                                                      % Number of calculation which the code runs at least
    'Enter the maximal generation number:'                                                      % Max. number of runs / generations
    'Enter the boundaries for convergenz'                                                       % Decision criteria if 2 values are close together or not
    'Enter the convergenz criteria number:'                                                     % How many "close" values have to be in row, that he code stops
    'Enter the a biological treshhold for the optimized value'                                  % At which max. c_k should the code stop
    'Enter the timestep for the results'                                                        % At which timepoint the user wants to know the calcium concentration
    'Enter the percentage level of ultimate distribution:'
    'Enter the fitness value'};

ui_dlg_title_1 = 'User input';

ui_num_lines_1 = 1;

if ui_choice == 1
ui_defaultanswer_1 = {'8','0.8','0.005','0.1','n','5','15','0.02','3','0.55','740','0.5','Ca_i'};                 % Just a default answer
end

if ui_choice == 2
ui_defaultanswer_1 = {'8','0.8','0.005','0.1','n','5','15','0.02','3','0.55','740','0.5','F_r'};                 % Just a default answer
end

if ui_choice == 3
ui_defaultanswer_1 = {'8','0.8','0.005','0.1','n','5','15','0.02','3','0.55','740','0.5','c_k'};                 % Just a default answer
end
ui_1 = inputdlg(ui_prompt_1,ui_dlg_title_1,ui_num_lines_1,ui_defaultanswer_1,'on');

ini_popsize       = str2double(ui_1{1});
co_p              = str2double(ui_1{2});
el_k              = str2double(ui_1{3});
mu_p_ini          = str2double(ui_1{4});
answer            = ui_1{5};
yes               = 'y';
min_nr_generation = str2double(ui_1{6});
max_nr_generation = str2double(ui_1{7});
eps               = str2double(ui_1{8});
nr_convergenz     = str2double(ui_1{9});
conc_nr           = str2double(ui_1{10});
timestep          = str2double(ui_1{11});
ui_nr_org_level   = str2double(ui_1{12});

check = strcmp(answer,yes);                                                                             % Returns answer of user about changing parameter boundaries

 
%individual boundaries for the different parameters

if check == 1
    
    for ini_prompt = 1 : ini_NrPar
        A = {'Enter lower boundary of randomisation of', ui_param{ini_prompt+1}};
        B = {'Enter upper boundary of randomisation of', ui_param{ini_prompt+1}};
        ini_prompt_vec{ini_prompt*2-1} = strjoin(A);
        ini_prompt_vec{ini_prompt*2}   = strjoin(B);
    end
    
    
    dlg_title_2         = 'Ca^{2+} optimization in the astrocyte';
    num_lines_2         = 1;
    for ini_prompt      = 1 : ini_NrPar
        defaultanswer_2(ini_prompt*2-1) = {'0.8'};
        defaultanswer_2(ini_prompt*2)   = {'1.2'};
    end
    ini_prompt_2    = ini_prompt_vec(1:ini_NrPar);
    ui_2            = inputdlg(ini_prompt_2,dlg_title_2,num_lines_2,defaultanswer_2,'on');
    
    ini_prompt_3    = ini_prompt_vec(ini_NrPar+1 : ini_NrPar*2);
    dlg_title_3     = 'Ca^{2+} optimization in the astrocyte';
    num_lines_3     = 1;
    defaultanswer_3 = defaultanswer_2;
    
    ui_3 = inputdlg(ini_prompt_3,dlg_title_3,num_lines_3,defaultanswer_3,'on');
    
    for ini_i = 1 : ini_NrPar/2
        
        low_b(ini_NrPar/2 + ini_i)           = str2double(ui_3{ini_i * 2 -1});
        up_b(ini_NrPar/2 + ini_i)            = str2double(ui_3{ini_i * 2});
        low_b(ini_i)                         = str2double(ui_2{ini_i * 2 -1});
        up_b(ini_i)                          = str2double(ui_2{ini_i * 2});
    end
    ini_i = 1;
else                                                               % IF user doesn't want to change default boundaries
    
    prompt_4 = {'Enter lower boundary of randomisation for all parameters:'
        'Enter upper boundary of randomisation for all parameters:'};
    
    dlg_title_4     = 'Ca^{2+} optimization in the astrocyte';
    
    num_lines_4     = 1;
    defaultanswer_4 = { '0.8','1.2'};
    
    ui_4 = inputdlg(prompt_4,dlg_title_4,num_lines_4,defaultanswer_4,'on');
    
    for ini_i = 1 : ini_NrPar
        
        low_b(ini_i)           = str2double(ui_4{1});
        up_b (ini_i)           = str2double(ui_4{2});
        
    end
    ini_i = 1;
end
% end of code for individual boundaries

%% initial values

% Initial population

T           = linspace (0,500,1000);                                                   % opens up the timevector
% # of different parameters
% Ratio between inital population and selected population

%running variables
mpar_i      = 1;
mpar_1      = 1;
mpar_c      = 1;

% Running variable in the parenting choosing
mu_1        = 1;
mu_2        = 1;
mu_var      = 1;

 
count_max   = 1;
conv_crit   = 1;
r_simu_2    = 1;
ini_w       = 1;

%used variables

%selection:
%sel_1      = variable for the for loop

%% First Generation -> Initiation of a random population

pop_w = 0;
pop_a = 1;

% ini_popsize = 20;
for ini_i = 1:ini_NrPar
    ini_randMat(ini_i,:) = vertcat(1,(up_b(ini_i)-low_b(ini_i))*rand(ini_popsize-1,1)+low_b(ini_i));
end

while pop_w == 0;
    
    for pop_i = 1 : ini_NrPar
        pop_c       = ((up_b(pop_i)+low_b(pop_i))/2);
        for pop_j   = pop_a : ini_popsize
            pop_d   = pop_c  - ini_randMat(pop_i,pop_j);
            pop_vector_1(pop_i,pop_j) = pop_d;
        end
        pop_p(pop_i)    = abs(sum(pop_vector_1(pop_i,:))/ini_popsize);
        pop_p_1(pop_i)  = (pop_p(pop_i) / abs(pop_c - up_b(pop_i)));
    end
    
    pop_res = (sum(pop_p_1) / ini_NrPar);
    
    if pop_res > ui_nr_org_level
        pop_a       = ini_popsize;
        ini_popsize = ini_popsize +  20;
        
        for ini_i = 1:ini_NrPar
            ini_randMat(ini_i,(pop_a+1):ini_popsize) = vertcat((up_b(ini_i)-low_b(ini_i))*rand((ini_popsize-pop_a),1)+low_b(ini_i));
        end
    else
        pop_w = 1;
    end
end

 
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
    results (ini_k,:) = nv.out(ui_1{13});
end

Max_matrix = [results(:,timestep),ini_randMat'];                     %takes c_k value at time 740ms because this is the area of interest and combines it with the responsible param. factors

Results_end_max(r_simu_2)  = max(results(:,timestep));                     %essential for plotting results
Results_end_mean(r_simu_2) = mean(results(:,timestep));                    %essential for plotting results

%% Selection (roulette wheel method)-> sel_...
%probability of each parameterset to be choosen for the offspring
%population -> probability for parameter sets with a high
%calciumconcentration is higher

%probabilities
sel_sum = sum(Max_matrix(1:ini_popsize,1));                                         %sum of all calcium concentrations of the generation

for sel_1 = 1:ini_popsize                                                                %ini_popsize= size of the population
    sel_prob            = Max_matrix(sel_1,1)./sel_sum;                                  %probability for each individual to be chosen
    sel_probvec(sel_1,:)= sel_prob;                                               %saving the probabilities
end

sel_cum_probvec = cumsum(sel_probvec);                                        %cummultative porbability
sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)',sel_cum_probvec,sel_probvec, Max_matrix]; %numberd max_matrix with probabilities

% decision
sel_output      = zeros(ini_popsize / 2,4+ini_NrPar);                                  %creates a zero matrix 4 because of probabilities, numberation

sel_rand_value  = rand(length(sel_output(:,1)),1);                             %random vector with the length of the number of parents and with random numbers between 0-1
for sel_2=1:length(sel_output(:,1))                                            %for the length of the number of the parents, that should be chosen
    
    while sel_output(sel_2,1) == 0                                             %till the row of the sel_output rows does include a individual
        sel_p=1;
        if sel_rand_value(sel_2,1) < sel_max_matrix(sel_p,2);     % if random value at the position n is smaler as the cummultative probability ...
            sel_output(sel_2,:) = sel_max_matrix(sel_p,:);            % ... save the individum in the sel_output matrix
        end
        
        while sel_rand_value(sel_2,1) > sel_max_matrix(sel_p,2);   % while the random number is bigger ....
            sel_p=sel_p+1;
            if sel_rand_value(sel_2,1) <= sel_max_matrix(sel_p,2); % test again if the cummultative probability of the next individium (parameter set) is bigger than the random number
                sel_output(sel_2,:) = sel_max_matrix(sel_p,:);         % ... save the individum in the sel_output matrix
            end
        end
        
        
        if sel_2>1                                                                    % prevent, that one individum (parameterset) is taken twice
            sel_duplex = ismember(sel_output(sel_2,1),sel_output(1:(sel_2-1),1)); % proves if the chosen individum is already in the sel_outtput matrix
            if sel_duplex == 1                                                    % if yes ...
                sel_output(sel_2,1)     = 0;                                          % delet it in the sel_otuput matrix
                sel_rand_value(sel_2,1) = rand(1);                                    % create a new random number at this place and repeat all the steps
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
    
    for mu_2 = 1 : ini_NrPar                                                     % Running variable for number of parameters / columns
        mu_rand = rand(1);                                                     % Random number between 0 and 1
        
        if mu_rand <= mu_p_ini                                                 % decision of mutation occurs or not
            for mu_var_1 = 1 : ini_NrPar
                if mu_2 == mu_var_1;                                                      % defines parameter specific boundaries
                    mu_r                            = low_b(mu_2) + (up_b(mu_2) - low_b(mu_2))*rand(1,1);
                    parent_matrix(mu_1,mu_2 + 4)    = mu_r;
                end
            end
            
        end
        
        mu_2 = mu_2 + 1;
    end
    mu_1 = mu_1 + 1;
    mu_2 = 1;
    
end

mu_p_total(r_simu_2) = mean(mu_p_ini);

%% Cross-over -> co...

co_s = length(parent_matrix(:,1));                                          %length of the column of the parents_matrix (should be initial population)
for co_var  = 1:(co_s-1)
    co_rand = rand(1);
    
    if co_rand <= co_p                                                          %probability of cross-over
        co_var1 = randi(ini_NrPar)+4;                                                  %defines the lower boundary of the crossover section
        co_var2 = randi([co_var1,ini_NrPar+4]);                                        %defines the upper boundary of the crossover section
        
        co_1 = parent_matrix(co_var,co_var1:co_var2);                               %cutting out the cross-over section
        co_2 = parent_matrix(co_var+1,co_var1:co_var2);
        
        parent_matrix(co_var,co_var1:co_var2)  = co_2;                              %change of the cross-over section between the 2 childrens
        parent_matrix(co_var+1,co_var1:co_var2) = co_1;
        
    end
    co_var = co_var + 2;
end                                                                         %output is a matrix which is 2*mpar_rbigger than matrix_input

 
%% Elitism
el_max(1)  = max(parent_matrix(:,4)) * (1-el_k);                        %Who is elite? maximal calcium concentration of generation multiplied by a number between 1-0
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

for ini_k = 1:ini_popsize
    
    if ui_choice == 1
        for ini_i = 1:ini_NrPar
            nv.smcec.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
        end
    end
    
    if ui_choice == 2
        for ini_i = 1:ini_NrPar
            nv.wall.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
        end
    end
    
    if ui_choice == 3
        for ini_i = 1:ini_NrPar
            nv.astrocyte.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
        end
    end
    
    nv.simulate()
    
    results (ini_k,:)= nv.out(ui_1{13});
    parent_matrix(ini_k,4)=  results(ini_k,timestep);
end

Results_end_max(r_simu_2+1)  = max(results(:,timestep));                     %essential for plotting results
Results_end_mean(r_simu_2+1) = mean(results(:,timestep));                    %essential for plotting results

 
%% Third and following generations

%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher
while count_max <= max_nr_generation-2 && Results_end_mean(r_simu_2+1)< conc_nr && conv_crit < nr_convergenz
    
    sel_1       = 1;                                        %Reset all running variables so the whole process can run again
    sel_2       = 1;
    sel_p       = 1;
    
    mpar_1      = 1;
    mpar_j      = 1;
    mpar_i      = 1;
    mpar_c      = 1;
    
    mu_1        = 1;
    mu_2        = 1;
    
    co_var      = 1;
    
    el_var      = 1;
    el_i        = 1;
    
    r_sim_1     = 1;
    
    %% Selection (roulette wheel method)-> sel_...
    
    %probability of each parameterset to be choosen for the offspring
    %population -> probability for parametersets with a high
    %calciumconcentration is higher
    sel_sum = sum(parent_matrix(1:ini_popsize,4));
    
    for sel_1= 1:ini_popsize
        sel_prob            = parent_matrix(sel_1,4)./sel_sum;
        sel_probvec(sel_1,:)= sel_prob;
    end
    
    
    % cummultative probability
    sel_cum_probvec = cumsum(sel_probvec);
    
    sel_max_matrix  = [linspace(1,ini_popsize,ini_popsize)',sel_cum_probvec,sel_probvec, parent_matrix(:,4:ini_NrPar+4)]; %numbered parent_matrix with probabilities
    
    % decision
    sel_output      = zeros((ini_popsize / 2),4+ini_NrPar);                               %4 because of probabilities, numberation
    sel_rand_value  = rand(length(sel_output(:,1)),1);
    
    
    for sel_2=1:length(sel_output(:,1))
        
        while sel_output(sel_2,1) == 0
            sel_p=1;
            if sel_rand_value(sel_2,1) < sel_max_matrix(sel_p,2);
                sel_output(sel_2,:) = sel_max_matrix(sel_p,:);
            end
            
            while sel_rand_value(sel_2,1) > sel_max_matrix(sel_p,2);
                sel_p=sel_p+1;
                if sel_rand_value(sel_2,1) <= sel_max_matrix(sel_p,2);
                    sel_output(sel_2,:) = sel_max_matrix(sel_p,:);
                end
            end
            
            if sel_2>1                                                              %prevent, that one parameterset is taken twice
                sel_duplex = ismember(sel_output(sel_2,1),sel_output(1:(sel_2-1),1));
                if sel_duplex == 1
                    sel_output(sel_2,1)     = 0;
                    sel_rand_value(sel_2,1) = rand(1);
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
            
            mpar_c                    = mpar_c + 2;
            
        end
        mpar_1 = mpar_1 + 2;
    end
    
    
    %% Mutation
    
    for mu_1 = 1:length(parent_matrix(:,1))
        
        for mu_2 = 1 : ini_NrPar
            
            mu_max      = max(parent_matrix(:,4));
            mu_mean     = mean(parent_matrix(:,4));
            mu_current  = parent_matrix(mu_1,4);
            
            if mu_current >= mu_mean                                                     %makes mutation rate adaptive to quality of organism's fitness
                p_mu(mu_var) = 0.1 * ((mu_max - mu_current) / (mu_max - mu_mean));     %If organism is fitter than average -> mutation should decrease
            else
                p_mu(mu_var) = 0.10;                                            %if not mutation rate stays the same
            end
            
            mu_rand = rand(1);
            
            if mu_rand <= p_mu(mu_var)                                                 % decision of mutation occurs or not
                for mu_var_1 = 1 : ini_NrPar
                    if mu_2 == mu_var_1;                                                      % defines parameter specific boundaries
                        mu_r = low_b(mu_2) + (up_b(mu_2) - low_b(mu_2))*rand(1,1);
                        parent_matrix(mu_1,mu_2 + 4) = mu_r;
                    end
                end
                
            end
            
            mu_2     = mu_2 + 1;
            mu_var   = mu_var + 1;
        end
        mu_1 = mu_1 + 1;
        mu_2 = 1;
        
    end
    
    mu_p_total(r_simu_2+1) = mean(p_mu);                                              %in order to plot change of mutation over all generations
    
    %% Cross-over
    co_s = length(parent_matrix(:,1));                                          %length of the column of the parents_matrix (should be initial population)
    
    while co_var      <= (co_s-1)
        
        co_rand        = rand(1);
        cr_max         = max(parent_matrix(:,4));
        cr_mean        = mean(parent_matrix(:,4));
        cr_current     = parent_matrix(co_var,4);
        
        if cr_current >= cr_mean
            co_p = ((cr_max - cr_current) / (cr_max - cr_mean));
        else
            co_p = 1;
        end
        
        if co_rand <= co_p                                                              %probability of cross-over
            
            co_var1 = randi(ini_NrPar)+4;                                               %defines the lower boundary of the crossover section
            co_var2 = randi([co_var1,ini_NrPar+4]);                                     %defines the upper boundary of the crossover section
            
            co_1    = parent_matrix(co_var,co_var1:co_var2);                               %cutting out the cross-over section
            co_2    = parent_matrix(co_var+1,co_var1:co_var2);
            
            parent_matrix(co_var,co_var1:co_var2)   = co_2;                             %change of the cross-over section between the 2 childrens
            parent_matrix(co_var+1,co_var1:co_var2) = co_1;
            
        end
        co_var = co_var + 2;
    end                                                                                 %output is a matrix which is 2*mpar_rbigger than matrix_input
    
    %% Elitism
    clear el_p;
    el_max(r_simu_2+1)  = max(parent_matrix(:,4)) * (1-el_k);
    for el_var = 1:length(sel_output(:,1))
        el_u=0;
        if sel_output(el_var,4) >= el_max(1,r_simu_2+1)
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
    for ini_k = 1:ini_popsize
        
        if ui_choice == 1
            for ini_i = 1:ini_NrPar
                nv.smcec.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
            end
        end
        
        if ui_choice == 2
            for ini_i = 1:ini_NrPar
                nv.wall.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
            end
        end
        
        if ui_choice == 3
            for ini_i = 1:ini_NrPar
                nv.astrocyte.params.(ui_param{ini_i+1}) = ini_value(ini_i) * parent_matrix(ini_k,ini_i+4);
            end
        end
        nv.simulate()
        
        results (ini_k,:)= nv.out(ui_1{13});
        parent_matrix(ini_k,4)=  results(ini_k,timestep);
    end
    
    Results_end_max(r_simu_2+2)  = max(results(:,timestep));                     %essential for plotting results
    Results_end_mean(r_simu_2+2) = mean(results(:,timestep));                    %essential for plotting results
    
    if count_max >= min_nr_generation && Results_end_mean(r_simu_2+2)<= Results_end_mean(r_simu_2+1)*(1+eps) && Results_end_mean(r_simu_2+2) >= Results_end_mean(r_simu_2+1)*(1-eps)
        conv_crit = conv_crit+1;
    else
        conv_crit = 0;
    end
    
    
    r_simu_2 = r_simu_2 + 1;
    count_max = count_max+1;
    
end

%% Plots
v_initial = linspace(Max_matrix(1,1),Max_matrix(1,1),count_max+1); %vector for the initial value of the calcium concentration

figure(1)
plot(linspace(1,count_max+1,count_max+1)', Results_end_max)
legend('Best value of each Generation','Initial value')
title('Best value of every generation');
xlabel('# of generation');
ylabel('');
hold on;
grid on;

plot(linspace(1,count_max+1,count_max+1)',v_initial);
legend('Best value of each Generation','Initial value');
title('Best value of every generation');
xlabel('# of generation');
ylabel('');
grid on;
grid minor;

figure(2)
plot(linspace(1,count_max+1,count_max+1)', Results_end_mean)
legend('Average value of each Generation','Initial value')
title('Average value  of every generation');
xlabel('# of generation');
ylabel('');
hold on;
grid on;

plot(linspace(1,count_max+1,count_max+1)',v_initial);
legend('Average value of each Generation','Initial value');
title('Average value of every generation');
xlabel('# of generation');
ylabel('');
grid on;
grid minor;

figure(3)
plot(linspace(1,r_simu_2,r_simu_2)',mu_p_total');
legend('Change of average mutation rate');
title('Change of average mutation rate');
xlabel('# of generation');
ylabel('probability of mutation');
grid on;
grid minor;


