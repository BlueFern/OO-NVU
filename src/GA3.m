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
%% Run a basic simulation
nv.simulate()
ini_c_k_all = nv.out('c_k'); 
initial_c_k = ini_c_k_all(740,:);

      
T = linspace (0,500,1000); 
 
%% Determine the initial condition - User Input

prompt_1 = {'Enter size of population (if divided by number of organism, result must be out of Z:'      % This population size is constant      
            'Enter number of organism which are allowed to reproduce:'                                  % pop.:24 this can just be: 2,4,6,12      
            'Enter the initial probability of cross-over:'                                              % This value determines the frequency of co
            'Enter the threshold for elitism:'                                                          % 0.005 means that a organism must have a bigger c_k than 99.5% of population to be an elite
            'Enter the initial probability of mutation:'                                                % Mutation probability
            'Do you want to change the parameter boundaries(+/-20%): [y] [n]'                                  % If User wants to adjust parameter boundaries
            'Enter the minimal generation number:'                                                      % Number of calculation which the code runs at least
            'Enter the maximal generation number:'                                                      % Max. number of runs / generations
            'Enter the boundaries for convergenz'                                                       % Decision criteria if 2 values are close together or not
            'Enter the convergenz criteria number:'                                                     % How many "close" values have to be in row, that he code stops                        
            'Enter the maximal CA^(2+) concentration number:'};                                         % At which max. c_k should the code stop

 
        
dlg_title_1 = 'Ca^{2+} optimization in the astrocyte';

num_lines_1 = 1;

defaultanswer_1 = {'32','16','0.8','0.005','0.1','y','5','15','0.02','3','0.55'};                       % Just a default answer

x_1 = inputdlg(prompt_1,dlg_title_1,num_lines_1,defaultanswer_1,'on');

num_pop           = str2double(x_1{1});
nr_parents        = str2double(x_1{2});
p_co              = str2double(x_1{3});
p_el              = str2double(x_1{4}); 
p_mu_ini          = str2double(x_1{5}); 
answer            = x_1{6};
yes               = 'y'
min_nr_generation = str2double(x_1{7}); 
max_nr_generation = str2double(x_1{8});
eps               = str2double(x_1{9});
nr_convergenz     = str2double(x_1{10});
conc_nr           = str2double(x_1{11});

 
check = strcmp(answer,yes);                                                                             % Returns answer of user about changing parameter boundaries

if check == 1
prompt_2 = {'Enter lower boundary of randomisation of k_deg:'
            'Enter upper boundary of randomisation of k_deg:'
           
            'Enter lower boundary of randomisation of k_on:'
            'Enter upper boundary of randomisation of k_on:'
            
            'Enter lower boundary of randomisation of k_inh:'
            'Enter upper boundary of randomisation of k_inh:'
            
            'Enter lower boundary of randomisation of r_h:'
            'Enter upper boundary of randomisation of r_h:'
            
            'Enter lower boundary of randomisation of the K_G:'
            'Enter upper boundary of randomisation of the K_G:'
            
            'Enter lower boundary of randomisation of the K_I:'
            'Enter upper boundary of randomisation of the K_I:'
           
            'Enter lower boundary of randomisation of BK_end:'
            'Enter upper boundary of randomisation of BK_end:'
            
            'Enter lower boundary of randomisation of K_ex:'
            'Enter upper boundary of randomisation of K_ex:'
                        
};

dlg_title_2 = 'Ca^{2+} optimization in the astrocyte';

num_lines_2 = 1;
defaultanswer_2 = { '0.8','1.2','0.8','1.2','0.8','1.2','0.8','1.2','0.8','1.2',...
                    '0.8','1.2','0.8','1.2','0.8','1.2'};

x = inputdlg(prompt_2,dlg_title_2,num_lines_2,defaultanswer_2,'on');

 
 
mu_low_b_k_deg      = str2double(x{1});
mu_up_b_k_deg       = str2double(x{2});

mu_low_b_k_on      = str2double(x{3});
mu_up_b_k_on       = str2double(x{4});

mu_low_b_k_inh      = str2double(x{5});
mu_up_b_k_inh       = str2double(x{6});

mu_low_b_r_h      = str2double(x{7});
mu_up_b_r_h       = str2double(x{8});

mu_low_b_K_G      = str2double(x{9});
mu_up_b_K_G       = str2double(x{10});

mu_low_b_K_I      = str2double(x{11});
mu_up_b_K_I       = str2double(x{12});

mu_low_b_BK_end      = str2double(x{13});
mu_up_b_BK_end       = str2double(x{14});

mu_low_b_K_ex      = str2double(x{15});
mu_up_b_K_ex       = str2double(x{16});

 
prompt_3 = {'Enter lower boundary of randomisation of B_ex:'
            'Enter upper boundary of randomisation of B_ex:'
            
            'Enter lower boundary of randomisation of J_max:'
            'Enter upper boundary of randomisation of J_max:'
            
            'Enter lower boundary of randomisation of K_act:'
            'Enter upper boundary of randomisation of K_act:'
            
            'Enter lower boundary of randomisation of P_L:'
            'Enter upper boundary of randomisation of P_L:'
            
            'Enter lower boundary of randomisation of V_max:'
            'Enter upper boundary of randomisation of V_max:'
            
            'Enter lower boundary of randomisation of k_pump:'
            'Enter upper boundary of randomisation of k_pump:'
            
            'Enter lower boundary of randomisation of delta:'
            'Enter upper boundary of randomisation of delta:'
            
            'Enter lower boundary of randomisation of VR_ER_cyt:'
            'Enter upper boundary of randomisation of VR_ER_cyt:'
            
 
};

dlg_title_3 = 'Ca^{2+} optimization in the astrocyte';

num_lines_3 = 1;
defaultanswer_3 = { '0.8','1.2','0.8','1.2','0.8','1.2','0.8','1.2','0.8','1.2',...
                    '0.8','1.2','0.8','1.2','0.8','1.2'};

x_3 = inputdlg(prompt_3,dlg_title_3,num_lines_3,defaultanswer_3,'on');

mu_low_b_B_ex           = str2double(x_3{1});
mu_up_b_B_ex            = str2double(x_3{2});

mu_low_b_J_max          = str2double(x_3{3});
mu_up_b_J_max           = str2double(x_3{4});

mu_low_b_K_act          = str2double(x_3{5});
mu_up_b_K_act           = str2double(x_3{6});

mu_low_b_P_L            = str2double(x_3{7});
mu_up_b_P_L             = str2double(x_3{8});

mu_low_b_V_max          = str2double(x_3{9});
mu_up_b_V_max           = str2double(x_3{10});

mu_low_b_k_pump         = str2double(x_3{11});
mu_up_b_k_pump          = str2double(x_3{12});

mu_low_b_delta          = str2double(x_3{13});
mu_up_b_delta           = str2double(x_3{14});

mu_low_b_VR_ER_cyt      = str2double(x_3{15});
mu_up_b_VR_ER_cyt       = str2double(x_3{16});

else                                                                                        % IF user doesn't want to change default boundaries
  
mu_low_b_k_deg          = 0.8;
mu_up_b_k_deg           = 1.2;

mu_low_b_k_on           = 0.8;
mu_up_b_k_on            = 1.2;

mu_low_b_k_inh          = 0.8;
mu_up_b_k_inh           = 1.2;

mu_low_b_r_h            = 0.8;
mu_up_b_r_h             = 1.2;

mu_low_b_K_G            = 0.8;
mu_up_b_K_G             = 1.2;

mu_low_b_K_I            = 0.8;
mu_up_b_K_I             = 1.2;

mu_low_b_BK_end     	= 0.8;
mu_up_b_BK_end          = 1.2;

mu_low_b_K_ex           = 0.8;
mu_up_b_K_ex            = 1.2;

mu_low_b_B_ex           = 0.8;
mu_up_b_B_ex            = 1.2;

mu_low_b_J_max          = 0.8;
mu_up_b_J_max           = 1.2;

mu_low_b_K_act          = 0.8;
mu_up_b_K_act           = 1.2;

mu_low_b_P_L            = 0.8;
mu_up_b_P_L             = 1.2;

mu_low_b_V_max          = 0.8;
mu_up_b_V_max           = 1.2;

mu_low_b_k_pump         = 0.8;
mu_up_b_k_pump          = 1.2;

mu_low_b_delta          = 0.8;
mu_up_b_delta           = 1.2;

mu_low_b_VR_ER_cyt      = 0.8;
mu_up_b_VR_ER_cyt       = 1.2;
end

 
%% initial values

 
% First Generation
%num_pop = 30; %number of individuals in the first generation
z = 740; %timestep for the first population 
 
x1  = nv.astrocyte.params.k_deg;                                            %orginal values for the parameters used in the AC
x2  = nv.astrocyte.params.k_on;
x3  = nv.astrocyte.params.K_inh;
x4  = nv.astrocyte.params.r_h;
x5  = nv.astrocyte.params.K_G;
x6  = nv.astrocyte.params.K_I;
x7  = nv.astrocyte.params.BK_end;
x8  = nv.astrocyte.params.K_ex;
x9  = nv.astrocyte.params.B_ex;
x10 = nv.astrocyte.params.J_max;
x11 = nv.astrocyte.params.K_act;
x12 = nv.astrocyte.params.P_L;
x13 = nv.astrocyte.params.V_max;
x14 = nv.astrocyte.params.k_pump;
x15 = nv.astrocyte.params.delta;
x16 = nv.astrocyte.params.VR_ER_cyt;

%j = for loop, size of the first generation -> please change in the code

 
 
%Selection
%nr_parents = 10;
%n = number of parents -> please change in the code

%% Parentchoosing,Crossover

r_co        = num_pop/nr_parents;                                           % Ratio between inital population and selected population
pp_1        = 1;
pp          = 1;
rc          = 1;                                                            % Running variable in the parenting choosing
mu_1        = 1;
mu_2        = 1;
n_el        = 1;
el_parII    = 1;
nr_par      = 16;                                                           % # of different parameters
mu_var      = 1;
count_max   = 1;
conv_crit   = 1;
ga          = 1;                                                           
%p_co = 0; %Probability of cross over
% co = size of parent-matrix

%Mutation Initial values

%p_mu        = 0.2;             %defines how often mutation occur
%mu_low_b    = 0.8;         %defines the range of mutation
%mu_up_b     = 1.2;
%p_el        = 0.005; % Probability that one organism is allowed to skip mutation and live on

 
 
%% First Generation -> Initiation of a random population

%Create r.  Range of r is a to b. Limits for the parameterchange
%low_bound = 0.8; % Lower (min) value
%up_bound = +1.2; % Upper (max) value

for j = 1:num_pop  %size of the first generation 
 
randVec1  = vertcat(1,mu_low_b_k_deg + (mu_up_b_k_deg - mu_low_b_k_deg).*rand(j-1, 1)); 
randVec2  = vertcat(1,mu_low_b_k_on + (mu_up_b_k_on - mu_low_b_k_on).*rand(j-1, 1));
randVec3  = vertcat(1,mu_low_b_k_inh + (mu_up_b_k_inh - mu_low_b_k_inh).*rand(j-1, 1));
randVec4  = vertcat(1,mu_low_b_r_h + (mu_up_b_r_h - mu_low_b_r_h).*rand(j-1, 1));
randVec5  = vertcat(1,mu_low_b_K_G + (mu_up_b_K_G - mu_low_b_K_G).*rand(j-1, 1));
randVec6  = vertcat(1,mu_low_b_K_I + (mu_up_b_K_I - mu_low_b_K_I).*rand(j-1, 1)); 
randVec7  = vertcat(1,mu_low_b_BK_end + (mu_up_b_BK_end - mu_low_b_BK_end).*rand(j-1, 1));
randVec8  = vertcat(1,mu_low_b_K_ex + (mu_up_b_K_ex - mu_low_b_K_ex).*rand(j-1, 1));
randVec9  = vertcat(1,mu_low_b_B_ex + (mu_up_b_B_ex - mu_low_b_B_ex).*rand(j-1, 1));
randVec10 = vertcat(1,mu_low_b_J_max + (mu_up_b_J_max - mu_low_b_J_max).*rand(j-1, 1));
randVec11 = vertcat(1,mu_low_b_K_act + (mu_up_b_K_act - mu_low_b_K_act).*rand(j-1, 1)); 
randVec12 = vertcat(1,mu_low_b_P_L + (mu_up_b_P_L - mu_low_b_P_L).*rand(j-1, 1));
randVec13 = vertcat(1,mu_low_b_V_max + (mu_up_b_V_max - mu_low_b_V_max).*rand(j-1, 1));
randVec14 = vertcat(1,mu_low_b_k_pump + (mu_up_b_k_pump - mu_low_b_k_pump).*rand(j-1, 1));
randVec15 = vertcat(1,mu_low_b_delta + (mu_up_b_delta - mu_low_b_delta).*rand(j-1, 1));
randVec16 = vertcat(1,mu_low_b_VR_ER_cyt + (mu_up_b_VR_ER_cyt - mu_low_b_VR_ER_cyt).*rand(j-1, 1));

 
 
nv.astrocyte.params.k_deg       = randVec1(j,1)*x1;
nv.astrocyte.params.k_on        = randVec2(j,1)*x2;
nv.astrocyte.params.K_inh       = randVec3(j,1)*x3;
nv.astrocyte.params.r_h         = randVec4(j,1)*x4;
nv.astrocyte.params.K_G         = randVec5(j,1)*x5;
nv.astrocyte.params.K_I         = randVec6(j,1)*x6;
nv.astrocyte.params.BK_end      = randVec7(j,1)*x7;
nv.astrocyte.params.K_ex        = randVec8(j,1)*x8;
nv.astrocyte.params.B_ex        = randVec9(j,1)*x9;
nv.astrocyte.params.J_max       = randVec10(j,1)*x10;
nv.astrocyte.params.K_act       = randVec11(j,1)*x11;
nv.astrocyte.params.P_L         = randVec12(j,1)*x12;
nv.astrocyte.params.V_max       = randVec13(j,1)*x13;
nv.astrocyte.params.k_pump      = randVec14(j,1)*x14;
nv.astrocyte.params.delta       = randVec15(j,1)*x15;
nv.astrocyte.params.VR_ER_cyt   = randVec16(j,1)*x16;

nv.simulate()  
 
results_c_k (j,:)   = nv.out('c_k');
time_results_c_k    = [T;results_c_k]';

results_R (j,:)     = nv.out('R');
time_results_R      = [T;results_R]';

end

Big_randMat = [randVec1, randVec2, randVec3, randVec4, randVec5,randVec6, randVec7, randVec8, ...
               randVec9, randVec10,randVec11, randVec12, randVec13, randVec14, randVec15,randVec16];

Max_matrix = [time_results_c_k(z,2:j+1)', Big_randMat];                     %takes c_k value at time 740ms because this is the area of interest and combines it with the responsible param. factors       
 
Results_end_max(ga)  = max(time_results_c_k(z,2:j+1)');                     %essential for plotting results
Results_end_mean(ga) = mean(time_results_c_k(z,2:j+1)');                    %essential for plotting results 
 
%% Selection (roulette wheel method)-> sel_...

                                                                            %probability of each parameterset to be choosen for the offspring
                                                                            %population -> probability for parametersets with a high
                                                                            %calciumconcentration is higher

sel_sum_c_k = sum(Max_matrix(1:j,1));

for i= 1:j
sel_prob        = Max_matrix(i,1)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
end

% cummultative porbability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix  = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, Max_matrix]; %numberd max_matrix with probabilities

% decision
sel_output      = zeros(nr_parents,4+nr_par);                               %4 because of probabilities, numberation
sel_rand_value  = rand(length(sel_output(:,1)),1);

 
for n=1:length(sel_output(:,1))
      
    while sel_output(n,1) == 0 
         p=1;
         if sel_rand_value(n,1) < sel_max_matrix(p,2);
         sel_output(n,:) = sel_max_matrix(p,:);
         end

        while sel_rand_value(n,1) > sel_max_matrix(p,2);
         p=p+1;
            if sel_rand_value(n,1) <= sel_max_matrix(p,2);
            sel_output(n,:) = sel_max_matrix(p,:);
            end
        end
        
        
        if n>1                                                              %prevent, that one parameterset is taken twice
            sel_duplex = ismember(sel_output(n,1),sel_output(1:(n-1),1));
            if sel_duplex == 1
            sel_output(n,1)     = 0;
            sel_rand_value(n,1) = rand(1);
            end
        end
      
    end
end
     
 
%% Selecting Parents

% # of calculations


var_p = 1;
for pp = 1:n/2                                                              %divided by 2, because there are 2 parents   
        
    P1 = sel_output(var_p,:);                                               %choosing parents out of input matrix
    P2 = sel_output(var_p+1,:);
    
    for sp = 1:r_co                                                         % Repeat until initial population size is reached
        
        parent_matrix(rc,:)   = P1;
        parent_matrix(rc+1,:) = P2;
        
        rc = rc + 2;
        
    end 
   var_p = var_p + 2;
end 
 
 
%%Mutation of Parent Matrix

for mu_1 = 1:size(parent_matrix(:,1))                                       % Running variable for lines
    
   for mu_2 = 1 : nr_par                                                    % Running variable for number of parameters / columns
       r_mu_1 = rand(1);                                                    % Random number between 0 and 1
       
       if r_mu_1 <= p_mu_ini                                                % decision of mutation occurs or not 

           if mu_2 == 1;                                                    %defines parameter specific boundaries
                r_mu_2 = mu_low_b_k_deg + (mu_up_b_k_deg - mu_low_b_k_deg)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 2;
                r_mu_2 = mu_low_b_k_on + (mu_up_b_k_on - mu_low_b_k_on)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 3;
                r_mu_2 = mu_low_b_k_inh + (mu_up_b_k_inh - mu_low_b_k_inh)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 4;
                r_mu_2 = mu_low_b_r_h + (mu_up_b_r_h - mu_low_b_r_h)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 5
                r_mu_2 = mu_low_b_K_G + (mu_up_b_K_G - mu_low_b_K_G)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 6;
                r_mu_2 = mu_low_b_K_I + (mu_up_b_K_I - mu_low_b_K_I)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 7;
                r_mu_2 = mu_low_b_BK_end + (mu_up_b_BK_end - mu_low_b_BK_end)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 8;
                r_mu_2 = mu_low_b_K_ex + (mu_up_b_K_ex - mu_low_b_K_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 9;
                r_mu_2 = mu_low_b_B_ex + (mu_up_b_B_ex - mu_low_b_B_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 10;
               r_mu_2 = mu_low_b_J_max + (mu_up_b_J_max - mu_low_b_J_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 11;
               r_mu_2 = mu_low_b_K_act + (mu_up_b_K_act - mu_low_b_K_act)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 12;
               r_mu_2 = mu_low_b_P_L + (mu_up_b_P_L - mu_low_b_P_L)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 13;
               r_mu_2 = mu_low_b_V_max + (mu_up_b_V_max - mu_low_b_V_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 14;
               r_mu_2 = mu_low_b_k_pump + (mu_up_b_k_pump - mu_low_b_k_pump)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 15;
               r_mu_2 = mu_low_b_delta + (mu_up_b_delta - mu_low_b_delta)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 16;
               r_mu_2 = mu_low_b_VR_ER_cyt + (mu_up_b_VR_ER_cyt - mu_low_b_VR_ER_cyt)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                             
           end
           end
       
       mu_2 = mu_2 + 1;
   end
   mu_1 = mu_1 + 1;
   mu_2 = 1;
   
end

p_mu_total(ga) = mean(p_mu_ini);                                            
%% Cross-over

co_s        = size(parent_matrix(:,1));                                     %size of the column of the parents_matrix (should be initial population)
for co_var  = 1:(co_s-1)
cr_rand = rand(1);

if cr_rand <= p_co                                                          %probability of cross-over
cr_var1 = randi(nr_par)+4;                                                  %defines the lower boundary of the crossover section
cr_var2 = randi([cr_var1,nr_par+4]);                                        %defines the upper boundary of the crossover section

CO_1 = parent_matrix(co_var,cr_var1:cr_var2);                               %cutting out the cross-over section
CO_2 = parent_matrix(co_var+1,cr_var1:cr_var2); 
 
parent_matrix(co_var,cr_var1:cr_var2)  = CO_2;                              %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,cr_var1:cr_var2) = CO_1;

end 
co_var = co_var + 2;
end                                                                         %output is a matrix which is 2*r_co bigger than matrix_input 
 
 
%% Elitism 
el_max_c_k(1)  = max(parent_matrix(:,4)) * (1-p_el);                        %Who is elite? 
for var_el = 1:size(sel_output(:,1))
el_u=0;
     if sel_output(var_el,4) >= el_max_c_k(1)
         el_p(var_el,:) = sel_output(var_el,:);
         for el_i = 1:size(parent_matrix(:,1));
             if el_p(var_el,1) == parent_matrix(el_i,1) && el_u== 0
                 parent_matrix(el_i,:) = el_p(var_el,:);
                 el_u=1;
             end
         end
     end
end                                                                         %output is a matrix which is 2*r_co bigger than matrix_input 
 
for g=1:length(parent_matrix(:,1))                                          % simulate for the second generation

nv.astrocyte.params.k_deg       = x1*parent_matrix(g,5);
nv.astrocyte.params.k_on        = x2*parent_matrix(g,6);
nv.astrocyte.params.K_inh       = x3*parent_matrix(g,7);
nv.astrocyte.params.r_h         = x4*parent_matrix(g,8);
nv.astrocyte.params.K_G         = x5*parent_matrix(g,9);
nv.astrocyte.params.K_I         = x6*parent_matrix(g,10);
nv.astrocyte.params.BK_end      = x7*parent_matrix(g,11);
nv.astrocyte.params.K_ex        = x8*parent_matrix(g,12);
nv.astrocyte.params.B_ex        = x9*parent_matrix(g,13);
nv.astrocyte.params.J_max       = x10*parent_matrix(g,14);
nv.astrocyte.params.K_act       = x11*parent_matrix(g,15);
nv.astrocyte.params.P_L         = x12*parent_matrix(g,16);
nv.astrocyte.params.V_max       = x13*parent_matrix(g,17);
nv.astrocyte.params.k_pump      = x14*parent_matrix(g,18);
nv.astrocyte.params.delta       = x15*parent_matrix(g,19);
nv.astrocyte.params.VR_ER_cyt   = x16*parent_matrix(g,20);

nv.simulate()

results_c_k_gen(g,:)        =  nv.out('c_k');
results_R_gen(g,:)          =  nv.out('R');
parent_matrix(g,4)          =  results_c_k_gen(g,z);
end

Results_end_max(ga+1)  = max(results_c_k_gen(:,z)');                        %essential for plotting results
Results_end_mean(ga+1) = mean(results_c_k_gen(:,z)');                       %This is located afeter each generation to save the CA result

 
if count_max >= min_nr_generation && Results_end_mean(ga+1)<= Results_end_mean(ga)*(1+eps) && Results_end_mean(ga+1) >= Results_end_mean(ga)*(1-eps)
       conv_crit = conv_crit+1;
   else 
       conv_crit = 0;
end

 %% Third and following generations                                                     
    
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher
while count_max <= max_nr_generation-2 && Results_end_mean(ga+1)< conc_nr && conv_crit < nr_convergenz
    
    pp_1        = 1;                                                        %Reset all running variables so the whole process can run again
    pp          = 1;
    rc          = 1;                                                        % Running variable in the parenting choosing
    mu_1        = 1;
    mu_2        = 1;
    n_el        = 1;
    el_parII    = 1;
    i           = 1;
    var_p       = 1;
    n           = 1;
    co_var      = 1;
    g           = 1;
    p           = 1;
    var_el      = 1;
    sp          = 1;
    el_i        = 1;
    
 
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher
sel_sum_c_k = sum(parent_matrix(1:j,4));

for i= 1:j
sel_prob        = parent_matrix(i,4)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
end

 
% cummultative probability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix  = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, parent_matrix(:,4:nr_par+4)]; %numbered parent_matrix with probabilities

% decision
sel_output      = zeros(nr_parents,4+nr_par);                               %4 because of probabilities, numberation
sel_rand_value  = rand(length(sel_output(:,1)),1);

 
for n=1:length(sel_output(:,1))
    
    while sel_output(n,1) == 0 
         p=1;  
         if sel_rand_value(n,1) < sel_max_matrix(p,2);
         sel_output(n,:) = sel_max_matrix(p,:);
         end

        while sel_rand_value(n,1) > sel_max_matrix(p,2);
         p=p+1;
            if sel_rand_value(n,1) <= sel_max_matrix(p,2);
            sel_output(n,:) = sel_max_matrix(p,:);
            end
        end
        
        if n>1                                                              %prevent, that one parameterset is taken twice
            sel_duplex = ismember(sel_output(n,1),sel_output(1:(n-1),1));
            if sel_duplex == 1
            sel_output(n,1)=0;
            sel_rand_value(n,1)=rand(1);
            end
        end
      
    end
end
    
 
 
%% Selecting Parents

% # of calculations


var_p = 1;
for pp = 1:n/2                                                              %divided by 2, because there are 2 parents   
        
    P1 = sel_output(var_p,:);                                               %choosing parents out of input matrix
    P2 = sel_output(var_p+1,:);
    
    for sp = 1:r_co 
        
        parent_matrix(rc,:)   = P1;
        parent_matrix(rc+1,:) = P2;
        
        rc = rc + 2;
        
    end 
   var_p = var_p + 2;
end 
 
%% Mutation 
for mu_1 = 1:length(parent_matrix(:,1))
    
   for mu_2 = 1 : nr_par
       
        mu_max      = max(parent_matrix(:,4));
        mu_mean     = mean(parent_matrix(:,4));    
        mu_current  = parent_matrix(mu_1,4);
       
       if mu_current >= mu_mean                                                     %makes mutation rate adaptive to quality of organism's fitness
            p_mu(mu_var) = 0.10 * ((mu_max - mu_current) / (mu_max - mu_mean));     %If organism is fitter than average -> mutation should decrease
       else 
            p_mu(mu_var) = 0.10;                                            %if not mutation rate stays the same
       end
       
%        if mu_current == mu_max                                            %Otherwise the best organism would not mutate
%            p_mu(mu_var) = 0.02;
%        end
        
       r_mu_1 = rand(1);
       
       if r_mu_1 <= p_mu(mu_var)                                            %Decision if Mutation occurs or not
           
%            mu_factor_current    = parent_matrix(mu_1,mu_2 + 4);           %factor which is currently observed
%            mu_low_b             = mu_factor_current - p_mu(mu_var);       %new boundaries around factor
%            mu_up_b              = mu_factor_current + p_mu(mu_var); 
%            
%            if mu_low_b <= mu_low_b_ini                                    %initial boundaries are still important
%                mu_low_b = 0.8;
%            end
%            
%            if mu_up_b >= mu_up_b_ini 
%                mu_up_b = 1.2;
%            end

       if mu_2 == 1;                                                        %defines parameter specific boundaries
                r_mu_2 = mu_low_b_k_deg + (mu_up_b_k_deg - mu_low_b_k_deg)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 2;
                r_mu_2 = mu_low_b_k_on + (mu_up_b_k_on - mu_low_b_k_on)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 3;
                r_mu_2 = mu_low_b_k_inh + (mu_up_b_k_inh - mu_low_b_k_inh)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 4;
                r_mu_2 = mu_low_b_r_h + (mu_up_b_r_h - mu_low_b_r_h)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 5
                r_mu_2 = mu_low_b_K_G + (mu_up_b_K_G - mu_low_b_K_G)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 6;
                r_mu_2 = mu_low_b_K_I + (mu_up_b_K_I - mu_low_b_K_I)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 7;
                r_mu_2 = mu_low_b_BK_end + (mu_up_b_BK_end - mu_low_b_BK_end)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 8;
                r_mu_2 = mu_low_b_K_ex + (mu_up_b_K_ex - mu_low_b_K_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 9;
                r_mu_2 = mu_low_b_B_ex + (mu_up_b_B_ex - mu_low_b_B_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 10;
               r_mu_2 = mu_low_b_J_max + (mu_up_b_J_max - mu_low_b_J_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 11;
               r_mu_2 = mu_low_b_K_act + (mu_up_b_K_act - mu_low_b_K_act)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 12;
               r_mu_2 = mu_low_b_P_L + (mu_up_b_P_L - mu_low_b_P_L)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 13;
               r_mu_2 = mu_low_b_V_max + (mu_up_b_V_max - mu_low_b_V_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 14;
               r_mu_2 = mu_low_b_k_pump + (mu_up_b_k_pump - mu_low_b_k_pump)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 15;
               r_mu_2 = mu_low_b_delta + (mu_up_b_delta - mu_low_b_delta)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                 
           elseif mu_2 == 16;
               r_mu_2 = mu_low_b_VR_ER_cyt + (mu_up_b_VR_ER_cyt - mu_low_b_VR_ER_cyt)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                             
           end    
 
       end
       
       mu_2     = mu_2 + 1;
       mu_var   = mu_var + 1;
   end
   mu_1 = mu_1 + 1;
   mu_2 = 1;
   
end

p_mu_total(ga+1) = mean(p_mu);                                              %in order to plot change of mutation over all generations

%% Cross-over
co_s = length(parent_matrix(:,1));                                          %size of the column of the parents_matrix (should be initial population)

while co_var      <= (co_s-1)
    
    cr_rand        = rand(1);
    cr_max         = max(parent_matrix(:,4));
    cr_mean        = mean(parent_matrix(:,4));    
    cr_current     = parent_matrix(co_var,4);
    
    if cr_current >= cr_mean 
        p_co = ((cr_max - cr_current) / (cr_max - cr_mean));
    else
        p_co = 1;
    end

if cr_rand <= p_co                                                          %probability of cross-over
    
cr_var1 = randi(nr_par)+4;                                                  %defines the lower boundary of the crossover section
cr_var2 = randi([cr_var1,nr_par+4]);                                        %defines the upper boundary of the crossover section

CO_1 = parent_matrix(co_var,cr_var1:cr_var2);                               %cutting out the cross-over section
CO_2 = parent_matrix(co_var+1,cr_var1:cr_var2); 
 
parent_matrix(co_var,cr_var1:cr_var2)   = CO_2;                             %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,cr_var1:cr_var2) = CO_1;

end 
co_var = co_var + 2;
end                                                                         %output is a matrix which is 2*r_co bigger than matrix_input 
 
%% Elitism 
clear el_p;
el_max_c_k(ga+1)  = max(parent_matrix(:,4)) * (1-p_el); 
for var_el = 1:size(sel_output(:,1))
el_u=0;
     if sel_output(var_el,4) >= el_max_c_k(1,ga+1)
         el_p(var_el,:) = sel_output(var_el,:);
         for el_i = 1:size(parent_matrix(:,1));
             if el_p(var_el,1) == parent_matrix(el_i,1) && el_u== 0
                 parent_matrix(el_i,:) = el_p(var_el,:);
                 el_u=1;
             end
         end       
     end
end

for g=1:length(parent_matrix(:,1))

% nv.simulate()
    
nv.astrocyte.params.k_deg       = x1*parent_matrix(g,5);
nv.astrocyte.params.k_on        = x2*parent_matrix(g,6);
nv.astrocyte.params.K_inh       = x3*parent_matrix(g,7);
nv.astrocyte.params.r_h         = x4*parent_matrix(g,8);
nv.astrocyte.params.K_G         = x5*parent_matrix(g,9);
nv.astrocyte.params.K_I         = x6*parent_matrix(g,10);
nv.astrocyte.params.BK_end      = x7*parent_matrix(g,11);
nv.astrocyte.params.K_ex        = x8*parent_matrix(g,12);
nv.astrocyte.params.B_ex        = x9*parent_matrix(g,13);
nv.astrocyte.params.J_max       = x10*parent_matrix(g,14);
nv.astrocyte.params.K_act       = x11*parent_matrix(g,15);
nv.astrocyte.params.P_L         = x12*parent_matrix(g,16);
nv.astrocyte.params.V_max       = x13*parent_matrix(g,17);
nv.astrocyte.params.k_pump      = x14*parent_matrix(g,18);
nv.astrocyte.params.delta       = x15*parent_matrix(g,19);
nv.astrocyte.params.VR_ER_cyt   = x16*parent_matrix(g,20);

nv.simulate()

results_c_k_gen(g,:)= nv.out('c_k');
parent_matrix(g,4)  = results_c_k_gen(g,z);

Results_end_max(ga+2)  = max  (parent_matrix(:,4));                         %essential for plotting results
Results_end_mean(ga+2) = mean (parent_matrix(:,4));                         %essential for plotting results

end

 
   if count_max >= min_nr_generation && Results_end_mean(ga+2)<= Results_end_mean(ga+1)*(1+eps) && Results_end_mean(ga+2) >= Results_end_mean(ga+1)*(1-eps)
       conv_crit = conv_crit+1;
   else 
       conv_crit = 0;
   end

 
ga = ga + 1;
count_max = count_max+1;

end

 
 
figure(1)
v_initial_c_k = linspace(initial_c_k,initial_c_k,count_max+1);
plot(linspace(1,count_max+1,count_max+1)', Results_end_max)
    legend('Best [Ca2+] of each Generation','Initial [Ca2+]')
    title('Best CA^2^+ concentration in AC of every generation'); 
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    hold on; 
    grid on;

    plot(linspace(1,count_max+1,count_max+1)',v_initial_c_k);
    legend('Best [Ca2+] of each Generation','Initial [Ca2+]');
    title('Best CA^2^+ concentration in AC of every generation');
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    grid on;
    grid minor;
    
figure(2)
plot(linspace(1,count_max+1,count_max+1)', Results_end_mean)
    legend('Average [Ca2+] of each Generation','Initial [Ca2+]')
    title('Average CA^2^+ concentration in AC of every generation'); 
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    hold on; 
    grid on;

    plot(linspace(1,count_max+1,count_max+1)',v_initial_c_k);
    legend('Average [Ca2+] of each Generation','Initial [Ca2+]');
    title('Average CA^2^+ concentration in AC of every generation');
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    grid on;
    grid minor;

figure(3)
plot(linspace(1,ga,ga)',p_mu_total');
    legend('Change of average mutation rate');
    title('Change of average mutation rate');
    xlabel('# of generation'); 
    ylabel('probability of mutation'); 
    grid on;
    grid minor;

