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

      
T = linspace (0,500,1000); 
 
 
%Position in the beginning

prompt = {'Enter size of population:',
            'Enter number of organism which are allowed to reproduce:',
            'Enter lower boundary of randomisation of the initial population:',
            'Enter upper boundary of randomisation of the initial population:',
            'Enter the probability of mutation:',
            'Enter lower boundary of mutation:',
            'Enter upper boundary of mutation:',
            'Enter the probability of cross-over:',
            'Enter the threshold for elitism:'};
dlg_title = 'Ca^{2+} optimization in the astrocyte';
num_lines = 1;
defaultanswer = {'30','10','0.8','1.2','0.1','0.8','1.2','0.8','0.005'};

x = inputdlg(prompt,dlg_title,num_lines,defaultanswer);

num_pop     = str2num(x{1});
nr_parents  = str2num(x{2});
low_bound 	= str2num(x{3});
up_bound    = str2num(x{4});
p_mu        = str2num(x{5});
mu_low_b    = str2num(x{6});
mu_up_b     = str2num(x{7});
p_co        = str2num(x{8});
p_el        = str2num(x{9});

%% initial values

 
% First Generation
%num_pop = 30; %number of individuals in the first generation

x1  = nv.astrocyte.params.k_deg; %orginal values for the parameters used in the AC
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

rc   = 1; % Running variable in the parenting choosing
r_co = num_pop/nr_parents; %Ratio between inital population and selected population
pp_1 = 1;
pp   = 1;
nr_par = 16; % # of different parameters
%p_co = 0; %Probability of cross over
% co = size of parent-matrix

%Mutation Initial values
mu_1        = 1;
mu_2        = 1;
%p_mu        = 0.2;             %defines how often mutation occur
%mu_low_b    = 0.8;         %defines the range of mutation
%mu_up_b     = 1.2;
%p_el        = 0.005; % Probability that one organism is allowed to skip mutation and live on
n_el        = 1;
el_parII    = 1;

 
%% First Generation -> Initiation of a random population

%Create r.  Range of r is a to b. Limits for the parameterchange
%low_bound = 0.8; % Lower (min) value
%up_bound = +1.2; % Upper (max) value

for j = 1:num_pop
    %size of the first generation 
 
randVec1    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec2    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec3    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec4    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec5    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec6    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec7    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec8    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec9    = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec10   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec11   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec12   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec13   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec14   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec15   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec16   = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));

   
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

Big_randMat = [randVec1, randVec2, randVec3, randVec4,...
randVec5,randVec6, randVec7, randVec8, randVec9, randVec10,randVec11, randVec12, randVec13, randVec14, randVec15,randVec16];

 
 
Max_matrix = [time_results_c_k(740,2:j+1)', Big_randMat];

Sort_max_matrix = sortrows(Max_matrix);
 
 
 %% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher

sel_sum_c_k = sum(Max_matrix(1:j,1));

for i= 1:j
sel_prob = Max_matrix(i,1)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
end

% cummultative porbability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, Max_matrix]; %numberd max_matrix with probabilities

% decision
sel_output = zeros(nr_parents,4+nr_par); %4 because of probabilities, numberation
sel_rand_value = rand(length(sel_output(:,1)),1);

 
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
        
        
        if n>1                                                          %prevent, that one parameterset is taken twice
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

% r = randperm(n);%creates a shuffled linspace vector from 1 to n, but without duplicates -> so a parent can just be parent once
var_p = 1;
for pp = 1:n/2    %divided by 2, because there are 2 parents   
        
    P1 = sel_output(var_p,:);  %choosing parents out of input matrix
    P2 = sel_output(var_p+1,:);
    
    for sp = 1:r_co 
        
        parent_matrix(rc,:)   = P1;
        parent_matrix(rc+1,:) = P2;
        
        rc = rc + 2;
        
    end 
   var_p = var_p + 2;
end 
 
 
%%Mutation of Parent Matrix
%Cross-over
co_s = size(parent_matrix(:,1));                                %size of the column of the parents_matrix (should be initial population)
for co_var = 1:(co_s-1)
cr_rand = rand(1);

if cr_rand <= p_co                                               %probability of cross-over
cr_var1 = randi(nr_par)+4;                                         %defines the lower boundary of the crossover section
cr_var2 = randi([cr_var1,nr_par+4]);                               %defines the upper boundary of the crossover section

CO_1 = parent_matrix(co_var,cr_var1:cr_var2);                    %cutting out the cross-over section
CO_2 = parent_matrix(co_var+1,cr_var1:cr_var2); 
 
parent_matrix(co_var,cr_var1:cr_var2)  = CO_2;                   %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,cr_var1:cr_var2) = CO_1;

end 
co_var = co_var + 2;
end                                                              %output is a matrix which is 2*r_co bigger than matrix_input 
 
%%Mutation
%Initial values

 
for mu_1 = 1:size(parent_matrix(:,1))
    
   for mu_2 = 1 : nr_par
       r_mu_1 = rand(1);
       
       if r_mu_1 <= p_mu
                r_mu_2 = mu_low_b + (mu_up_b-mu_low_b)*rand(1,1);                       %randomised mutation factor with boundaries
                parent_matrix(mu_1,mu_2 + 4) = parent_matrix(mu_1,mu_2 + 4) * r_mu_2;   %Mutation changes parameter at current position
       end
       
       mu_2 = mu_2 + 1;
   end
   mu_1 = mu_1 + 1;
   mu_2 = 1;
   
end

%Elitism Part I - Who is elite? ------- position after Parent choosing
el_max_c_k  = max(parent_matrix(:,4) * (1-p_el));

for var_el = 1:size(parent_matrix(:,1))
   if parent_matrix(var_el,4) >= el_max_c_k         %Decision if organism is good enough to be elite - [Ca^2+] is 4 column
       v_el(n_el,:) = var_el;                       %Saving the position of the elite organism in order to find the line after mutation
       M_el(n_el,:) = parent_matrix(var_el,:);      %Save elite Set of parameter to save this information
       n_el = n_el + 1;
   end 
end

 
%Elitism Part II - Delete mutation changes of elite sets and reset the mutated line
%by its initial line
for el_parII = 1:length(v_el)
    parent_matrix(v_el(el_parII),:) = M_el(el_parII,:);
end

%%Output 
                                                      %output is a matrix which is 2*r_co bigger than matrix_input 
 
 
    
    
