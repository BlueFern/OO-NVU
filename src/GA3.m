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

%% Determine the initial condition - User Input -> us...

us_prompt_1 = {'Enter size of population (if divided by number of organism, result must be out of Z:'      % This population size is constant      
            'Enter number of organism which are allowed to reproduce:'                                  % pop.:24 this can just be: 2,4,6,12      
            'Enter the initial probability of cross-over:'                                              % This value determines the frequency of co
            'Enter the threshold for elitism:'                                                          % 0.005 means that a organism must have a bigger c_k than 99.5% of population to be an elite
            'Enter the initial probability of mutation:'                                                % Mutation probability
            'Do you want individual parameterboundaries [y] [n]'                           % If User wants to adjust parameter boundaries
            'Enter the minimal generation number:'                                                      % Number of calculation which the code runs at least
            'Enter the maximal generation number:'                                                      % Max. number of runs / generations
            'Enter the boundaries for convergenz'                                                       % Decision criteria if 2 values are close together or not
            'Enter the convergenz criteria number:'                                                     % How many "close" values have to be in row, that he code stops                        
            'Enter the maximal CA^(2+) concentration number:'                                           % At which max. c_k should the code stop
            'Enter the timestep for the CA^(2+)results'};                                               % At which timepoint the user wants to know the calcium concentration 
 
us_dlg_title_1 = 'Ca^{2+} optimization in the astrocyte';

us_num_lines_1 = 1;

us_defaultanswer_1 = {'8','4','0.8','0.005','0.1','n','5','15','0.02','3','0.55','740'};                 % Just a default answer

us_x1 = inputdlg(us_prompt_1,us_dlg_title_1,us_num_lines_1,us_defaultanswer_1,'on');

ini_popsize       = str2double(us_x1{1});
nr_parents        = str2double(us_x1{2});
co_p              = str2double(us_x1{3});
el_k              = str2double(us_x1{4}); 
mu_p_ini          = str2double(us_x1{5}); 
answer            = us_x1{6};
yes               = 'y'
min_nr_generation = str2double(us_x1{7}); 
max_nr_generation = str2double(us_x1{8});
eps               = str2double(us_x1{9});
nr_convergenz     = str2double(us_x1{10});
conc_nr           = str2double(us_x1{11});
timestep          = str2double(us_x1{12});

check = strcmp(answer,yes);                                                                             % Returns answer of user about changing parameter boundaries

 
%individual boundaries for the different parameters

if check == 1
us_prompt_2 = {'Enter lower boundary of randomisation of k_deg:'
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

x = inputdlg(us_prompt_2,dlg_title_2,num_lines_2,defaultanswer_2,'on');



low_b_k_deg      = str2double(x{1});
up_b_k_deg       = str2double(x{2});

low_b_k_on       = str2double(x{3});
up_b_k_on        = str2double(x{4});

low_b_k_inh      = str2double(x{5});
up_b_k_inh       = str2double(x{6});

low_b_r_h        = str2double(x{7});
up_b_r_h         = str2double(x{8});

low_b_K_G        = str2double(x{9});
up_b_K_G         = str2double(x{10});

low_b_K_I        = str2double(x{11});
up_b_K_I         = str2double(x{12});

low_b_BK_end     = str2double(x{13});
up_b_BK_end      = str2double(x{14});

low_b_K_ex       = str2double(x{15});
up_b_K_ex        = str2double(x{16});

 
us_us_prompt_3 = {'Enter lower boundary of randomisation of B_ex:'
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

x_3 = inputdlg(us_prompt_3,dlg_title_3,num_lines_3,defaultanswer_3,'on');

low_b_B_ex           = str2double(x_3{1});
up_b_B_ex            = str2double(x_3{2});

low_b_J_max          = str2double(x_3{3});
up_b_J_max           = str2double(x_3{4});

low_b_K_act          = str2double(x_3{5});
up_b_K_act           = str2double(x_3{6});

low_b_P_L            = str2double(x_3{7});
up_b_P_L             = str2double(x_3{8});

low_b_V_max          = str2double(x_3{9});
up_b_V_max           = str2double(x_3{10});

low_b_k_pump         = str2double(x_3{11});
up_b_k_pump          = str2double(x_3{12});

low_b_delta          = str2double(x_3{13});
up_b_delta           = str2double(x_3{14});

low_b_VR_ER_cyt      = str2double(x_3{15});
up_b_VR_ER_cyt       = str2double(x_3{16});

else                                                               % IF user doesn't want to change default boundaries

prompt_4 = {'Enter lower boundary of randomisation for all parameters:'
            'Enter upper boundary of randomisation for all parameters:'};
        
dlg_title_4 = 'Ca^{2+} optimization in the astrocyte';

num_lines_4 = 1;
defaultanswer_4 = { '0.8','1.2'};

x_4 = inputdlg(prompt_4,dlg_title_4,num_lines_4,defaultanswer_4,'on');

low_bound           = str2double(x_4{1});
up_bound            = str2double(x_4{2}); 
 
low_b_k_deg          = low_bound;
up_b_k_deg           = up_bound;

low_b_k_on           = low_bound;
up_b_k_on            = up_bound;

low_b_k_inh          = low_bound;
up_b_k_inh           = up_bound;

low_b_r_h            = low_bound;
up_b_r_h             = up_bound;

low_b_K_G            = low_bound;
up_b_K_G             = up_bound;

low_b_K_I            = low_bound;
up_b_K_I             = up_bound;

low_b_BK_end         = low_bound;
up_b_BK_end          = up_bound;

low_b_K_ex           = low_bound;
up_b_K_ex            = up_bound;

low_b_B_ex           = low_bound;
up_b_B_ex            = up_bound;

low_b_J_max          = low_bound;
up_b_J_max           = up_bound;

low_b_K_act          = low_bound;
up_b_K_act           = up_bound;

low_b_P_L            = low_bound;
up_b_P_L             = up_bound;

low_b_V_max          = low_bound;
up_b_V_max           = up_bound;

low_b_k_pump         = low_bound;
up_b_k_pump          = up_bound;

low_b_delta          = low_bound;
up_b_delta           = up_bound;

low_b_VR_ER_cyt      = low_bound;
up_b_VR_ER_cyt       = up_bound;
end
% end of code for individual boundaries

%% initial values

% Initial population

x1          = nv.astrocyte.params.k_deg;                                            %orginal values for the parameters used in the AC
x2          = nv.astrocyte.params.k_on;
x3          = nv.astrocyte.params.K_inh;
x4          = nv.astrocyte.params.r_h;
x5          = nv.astrocyte.params.K_G;
x6          = nv.astrocyte.params.K_I;
x7          = nv.astrocyte.params.BK_end;
x8          = nv.astrocyte.params.K_ex;
x9          = nv.astrocyte.params.B_ex;
x10         = nv.astrocyte.params.J_max;
x11         = nv.astrocyte.params.K_act;
x12         = nv.astrocyte.params.P_L;
x13         = nv.astrocyte.params.V_max;
x14         = nv.astrocyte.params.k_pump;
x15         = nv.astrocyte.params.delta;
x16         = nv.astrocyte.params.VR_ER_cyt;

T           = linspace (0,500,1000);                                                   % opens up the timevector
ini_NrPar     = 16;                                                           % # of different parameters
mpar_r      = ini_popsize/nr_parents;                                           % Ratio between inital population and selected population

%running variables
mpar_i          = 1;
mpar_1           = 1;
mpar_c          = 1;  
 
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

%Create r.  Range of r is a to b. Limits for the parameterchange
%low_bound = 0.8; % Lower (min) value
%up_bound = +1.2; % Upper (max) value

for j = 1:ini_popsize  %size of the first generation 

%     while ini_w == 1;
    
randVec1  = vertcat(1,low_b_k_deg + (up_b_k_deg - low_b_k_deg).*rand(j-1, 1)); 
randVec2  = vertcat(1,low_b_k_on + (up_b_k_on - low_b_k_on).*rand(j-1, 1));
randVec3  = vertcat(1,low_b_k_inh + (up_b_k_inh - low_b_k_inh).*rand(j-1, 1));
randVec4  = vertcat(1,low_b_r_h + (up_b_r_h - low_b_r_h).*rand(j-1, 1));
randVec5  = vertcat(1,low_b_K_G + (up_b_K_G - low_b_K_G).*rand(j-1, 1));
randVec6  = vertcat(1,low_b_K_I + (up_b_K_I - low_b_K_I).*rand(j-1, 1)); 
randVec7  = vertcat(1,low_b_BK_end + (up_b_BK_end - low_b_BK_end).*rand(j-1, 1));
randVec8  = vertcat(1,low_b_K_ex + (up_b_K_ex - low_b_K_ex).*rand(j-1, 1));
randVec9  = vertcat(1,low_b_B_ex + (up_b_B_ex - low_b_B_ex).*rand(j-1, 1));
randVec10 = vertcat(1,low_b_J_max + (up_b_J_max - low_b_J_max).*rand(j-1, 1));
randVec11 = vertcat(1,low_b_K_act + (up_b_K_act - low_b_K_act).*rand(j-1, 1)); 
randVec12 = vertcat(1,low_b_P_L + (up_b_P_L - low_b_P_L).*rand(j-1, 1));
randVec13 = vertcat(1,low_b_V_max + (up_b_V_max - low_b_V_max).*rand(j-1, 1));
randVec14 = vertcat(1,low_b_k_pump + (up_b_k_pump - low_b_k_pump).*rand(j-1, 1));
randVec15 = vertcat(1,low_b_delta + (up_b_delta - low_b_delta).*rand(j-1, 1));
randVec16 = vertcat(1,low_b_VR_ER_cyt + (up_b_VR_ER_cyt - low_b_VR_ER_cyt).*rand(j-1, 1));

% for ini_i = 1 : ini_nr_param
%     
%     for ini_j = a : ini_nr_organ
%    ini_d = (
%     end
%     
% end
%  
%     end
 
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

end

Max_matrix = [time_results_c_k(timestep,2:j+1)',randVec1, randVec2, randVec3, randVec4, randVec5,randVec6, randVec7,... 
              randVec8,randVec9, randVec10,randVec11, randVec12, randVec13, randVec14, randVec15,randVec16 ];                     %takes c_k value at time 740ms because this is the area of interest and combines it with the responsible param. factors       
                                             
Results_end_max(r_simu_2)  = max(time_results_c_k(timestep,2:j+1)');                     %essential for plotting results
Results_end_mean(r_simu_2) = mean(time_results_c_k(timestep,2:j+1)');                    %essential for plotting results 
 
%% Selection (roulette wheel method)-> sel_...
%probability of each parameterset to be choosen for the offspring                                                              
%population -> probability for parameter sets with a high
%calciumconcentration is higher

%probabilities
sel_sum_c_k = sum(Max_matrix(1:j,1));                                         %sum of all calcium concentrations of the generation

for sel_1= 1:j                                                                %j= size of the population 
sel_prob = Max_matrix(sel_1,1)./sel_sum_c_k;                                  %probability for each individual to be chosen
sel_probvec(sel_1,:)= sel_prob;                                               %saving the probabilities
end

sel_cum_probvec = cumsum(sel_probvec);                                        %cummultative porbability
sel_max_matrix  = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, Max_matrix]; %numberd max_matrix with probabilities

% decision
sel_output      = zeros(nr_parents,4+ini_NrPar);                                  %creates a zero matrix 4 because of probabilities, numberation
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
   
    for mpar_j = 1:mpar_r                                                        % Repeat until initial population size is reached
        
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
 
           if mu_2 == 1;                                                      % defines parameter specific boundaries
                mu_r = low_b_k_deg + (up_b_k_deg - low_b_k_deg)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 2;
                mu_r = low_b_k_on + (up_b_k_on - low_b_k_on)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 3;
                mu_r = low_b_k_inh + (up_b_k_inh - low_b_k_inh)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 4;
                mu_r = low_b_r_h + (up_b_r_h - low_b_r_h)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 5
                mu_r = low_b_K_G + (up_b_K_G - low_b_K_G)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 6;
                mu_r = low_b_K_I + (up_b_K_I - low_b_K_I)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 7;
                mu_r = low_b_BK_end + (up_b_BK_end - low_b_BK_end)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 8;
                mu_r = low_b_K_ex + (up_b_K_ex - low_b_K_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 9;
                mu_r = low_b_B_ex + (up_b_B_ex - low_b_B_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 10;
               mu_r = low_b_J_max + (up_b_J_max - low_b_J_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 11;
               mu_r = low_b_K_act + (up_b_K_act - low_b_K_act)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 12;
               mu_r = low_b_P_L + (up_b_P_L - low_b_P_L)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 13;
               mu_r = low_b_V_max + (up_b_V_max - low_b_V_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 14;
               mu_r = low_b_k_pump + (up_b_k_pump - low_b_k_pump)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 15;
               mu_r = low_b_delta + (up_b_delta - low_b_delta)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 16;
               mu_r = low_b_VR_ER_cyt + (up_b_VR_ER_cyt - low_b_VR_ER_cyt)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                             
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
el_max_c_k(1)  = max(parent_matrix(:,4)) * (1-el_k);                        %Who is elite? maximal calcium concentration of generation multiplied by a number between 1-0
for el_var = 1:length(sel_output(:,1))  %for the length of sel_output check if the the calcium concentration of the set is bigger than the el_max....
el_u=0;
     if sel_output(el_var,4) >= el_max_c_k(1)%... if it is bigger save it in a matrix and check where the same number of the set is in the parentmatrix, when founded overwirte the line
         el_p(el_var,:) = sel_output(el_var,:);
         for el_i = 1:length(parent_matrix(:,1));
             if el_p(el_var,1) == parent_matrix(el_i,1) && el_u== 0
                 parent_matrix(el_i,:) = el_p(el_var,:);
                 el_u=1;
             end
         end
     end
end                                                                         %output is a matrix which is 2*mpar_rbigger than matrix_input 
 
for r_sim_1=1:length(parent_matrix(:,1))                                          % simulate for the second generation

nv.astrocyte.params.k_deg       = x1*parent_matrix(r_sim_1,5);
nv.astrocyte.params.k_on        = x2*parent_matrix(r_sim_1,6);
nv.astrocyte.params.K_inh       = x3*parent_matrix(r_sim_1,7);
nv.astrocyte.params.r_h         = x4*parent_matrix(r_sim_1,8);
nv.astrocyte.params.K_G         = x5*parent_matrix(r_sim_1,9);
nv.astrocyte.params.K_I         = x6*parent_matrix(r_sim_1,10);
nv.astrocyte.params.BK_end      = x7*parent_matrix(r_sim_1,11);
nv.astrocyte.params.K_ex        = x8*parent_matrix(r_sim_1,12);
nv.astrocyte.params.B_ex        = x9*parent_matrix(r_sim_1,13);
nv.astrocyte.params.J_max       = x10*parent_matrix(r_sim_1,14);
nv.astrocyte.params.K_act       = x11*parent_matrix(r_sim_1,15);
nv.astrocyte.params.P_L         = x12*parent_matrix(r_sim_1,16);
nv.astrocyte.params.V_max       = x13*parent_matrix(r_sim_1,17);
nv.astrocyte.params.k_pump      = x14*parent_matrix(r_sim_1,18);
nv.astrocyte.params.delta       = x15*parent_matrix(r_sim_1,19);
nv.astrocyte.params.VR_ER_cyt   = x16*parent_matrix(r_sim_1,20);

nv.simulate()

results_c_k_gen(r_sim_1,:)  =  nv.out('c_k');
parent_matrix(r_sim_1,4)    =  results_c_k_gen(r_sim_1,timestep);
end

Results_end_max(r_simu_2+1)  = max(results_c_k_gen(:,timestep)');                        %essential for plotting results
Results_end_mean(r_simu_2+1) = mean(results_c_k_gen(:,timestep)');                       %This is located afeter each generation to save the CA result

 
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
sel_sum_c_k = sum(parent_matrix(1:j,4));

for sel_1= 1:j
sel_prob = parent_matrix(sel_1,4)./sel_sum_c_k;
sel_probvec(sel_1,:)= sel_prob;
end

 
% cummultative probability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix  = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, parent_matrix(:,4:ini_NrPar+4)]; %numbered parent_matrix with probabilities

% decision
sel_output      = zeros(nr_parents,4+ini_NrPar);                               %4 because of probabilities, numberation
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
            sel_output(sel_2,1)=0;
            sel_rand_value(sel_2,1)=rand(1);
            end
        end
      
    end
end
    
 
 
%% Multiplicate Parents

 
for mpar_i = 1:sel_2/2                                                              %divided by 2, because there are 2 parents   
        
    mpar_P1 = sel_output(mpar_1,:);                                               %choosing parents out of input matrix
    mpar_P2 = sel_output(mpar_1+1,:);
    
    for mpar_j = 1:mpar_r
        
        parent_matrix(mpar_c,:)   = mpar_P1;
        parent_matrix(mpar_c+1,:) = mpar_P2;
        
        mpar_c = mpar_c + 2;
        
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
       
%        if mu_current == mu_max                                            %Otherwise the best organism would not mutate
%            p_mu(mu_var) = 0.02;
%        end
        
       mu_rand = rand(1);
       
       if mu_rand <= p_mu(mu_var)                                            %Decision if Mutation occurs or not
           
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
                mu_r = low_b_k_deg + (up_b_k_deg - low_b_k_deg)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 2;
                mu_r = low_b_k_on + (up_b_k_on - low_b_k_on)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 3;
                mu_r = low_b_k_inh + (up_b_k_inh - low_b_k_inh)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 4;
                mu_r = low_b_r_h + (up_b_r_h - low_b_r_h)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 5
                mu_r = low_b_K_G + (up_b_K_G - low_b_K_G)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 6;
                mu_r = low_b_K_I + (up_b_K_I - low_b_K_I)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 7;
                mu_r = low_b_BK_end + (up_b_BK_end - low_b_BK_end)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 8;
                mu_r = low_b_K_ex + (up_b_K_ex - low_b_K_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 9;
                mu_r = low_b_B_ex + (up_b_B_ex - low_b_B_ex)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 10;
               mu_r = low_b_J_max + (up_b_J_max - low_b_J_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 11;
               mu_r = low_b_K_act + (up_b_K_act - low_b_K_act)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 12;
               mu_r = low_b_P_L + (up_b_P_L - low_b_P_L)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 13;
               mu_r = low_b_V_max + (up_b_V_max - low_b_V_max)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 14;
               mu_r = low_b_k_pump + (up_b_k_pump - low_b_k_pump)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 15;
               mu_r = low_b_delta + (up_b_delta - low_b_delta)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                 
           elseif mu_2 == 16;
               mu_r = low_b_VR_ER_cyt + (up_b_VR_ER_cyt - low_b_VR_ER_cyt)*rand(1,1);      
                parent_matrix(mu_1,mu_2 + 4) = mu_r;                                             
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

if co_rand <= co_p                                                          %probability of cross-over
    
co_var1 = randi(ini_NrPar)+4;                                                  %defines the lower boundary of the crossover section
co_var2 = randi([co_var1,ini_NrPar+4]);                                        %defines the upper boundary of the crossover section

co_1 = parent_matrix(co_var,co_var1:co_var2);                               %cutting out the cross-over section
co_2 = parent_matrix(co_var+1,co_var1:co_var2); 
 
parent_matrix(co_var,co_var1:co_var2)   = co_2;                             %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,co_var1:co_var2) = co_1;

end 
co_var = co_var + 2;
end                                                                         %output is a matrix which is 2*mpar_rbigger than matrix_input 
 
%% Elitism 
clear el_p;
el_max_c_k(r_simu_2+1)  = max(parent_matrix(:,4)) * (1-el_k); 
for el_var = 1:length(sel_output(:,1))
el_u=0;
     if sel_output(el_var,4) >= el_max_c_k(1,r_simu_2+1)
         el_p(el_var,:) = sel_output(el_var,:);
         for el_i = 1:length(parent_matrix(:,1));
             if el_p(el_var,1) == parent_matrix(el_i,1) && el_u== 0
                 parent_matrix(el_i,:) = el_p(el_var,:);
                 el_u=1;
             end
         end       
     end
end

for r_sim_1=1:length(parent_matrix(:,1))

% nv.simulate()
    
nv.astrocyte.params.k_deg       = x1*parent_matrix(r_sim_1,5);
nv.astrocyte.params.k_on        = x2*parent_matrix(r_sim_1,6);
nv.astrocyte.params.K_inh       = x3*parent_matrix(r_sim_1,7);
nv.astrocyte.params.r_h         = x4*parent_matrix(r_sim_1,8);
nv.astrocyte.params.K_G         = x5*parent_matrix(r_sim_1,9);
nv.astrocyte.params.K_I         = x6*parent_matrix(r_sim_1,10);
nv.astrocyte.params.BK_end      = x7*parent_matrix(r_sim_1,11);
nv.astrocyte.params.K_ex        = x8*parent_matrix(r_sim_1,12);
nv.astrocyte.params.B_ex        = x9*parent_matrix(r_sim_1,13);
nv.astrocyte.params.J_max       = x10*parent_matrix(r_sim_1,14);
nv.astrocyte.params.K_act       = x11*parent_matrix(r_sim_1,15);
nv.astrocyte.params.P_L         = x12*parent_matrix(r_sim_1,16);
nv.astrocyte.params.V_max       = x13*parent_matrix(r_sim_1,17);
nv.astrocyte.params.k_pump      = x14*parent_matrix(r_sim_1,18);
nv.astrocyte.params.delta       = x15*parent_matrix(r_sim_1,19);
nv.astrocyte.params.VR_ER_cyt   = x16*parent_matrix(r_sim_1,20);

nv.simulate()

results_c_k_gen(r_sim_1,:)= nv.out('c_k');
parent_matrix(r_sim_1,4)  = results_c_k_gen(r_sim_1,timestep);

Results_end_max(r_simu_2+2)  = max  (parent_matrix(:,4));                         %essential for plotting results
Results_end_mean(r_simu_2+2) = mean (parent_matrix(:,4));                         %essential for plotting results

end

 
   if count_max >= min_nr_generation && Results_end_mean(r_simu_2+2)<= Results_end_mean(r_simu_2+1)*(1+eps) && Results_end_mean(r_simu_2+2) >= Results_end_mean(r_simu_2+1)*(1-eps)
       conv_crit = conv_crit+1;
   else 
       conv_crit = 0;
   end

 
r_simu_2 = r_simu_2 + 1;
count_max = count_max+1;

end

%% Plots
v_initial_c_k = linspace(Max_matrix(1,1),Max_matrix(1,1),count_max+1); %vector for the initial value of the calcium concentration

figure(1)
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
    plot(linspace(1,r_simu_2,r_simu_2)',mu_p_total');
    legend('Change of average mutation rate');
    title('Change of average mutation rate');
    xlabel('# of generation'); 
    ylabel('probability of mutation'); 
    grid on;
    grid minor;
