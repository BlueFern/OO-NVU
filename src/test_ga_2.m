%% Genetic Algorithm
%Initiation of a random population

clear; clc; close all
%% Load NVU model
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
 
%Position in the beginning

prompt = {  'Enter size of population (if divided by number of organism, result must be out of Z:',
            'Enter number of organism which are allowed to reproduce:',
            'Enter lower boundary of randomisation of the initial population:',
            'Enter upper boundary of randomisation of the initial population:',
            'Enter the initial probability of mutation:',
            'Enter lower boundary of mutation:',
            'Enter upper boundary of mutation:',
            'Enter the initial probability of cross-over:',
            'Enter the threshold for elitism:'
            'Enter the number of generations'};
        
dlg_title = 'Ca^{2+} optimization in the astrocyte';

num_lines = 1;

defaultanswer = {'20','10','0.8','1.2','0.1','0.8','1.2','0.8','0.005','5'};

x = inputdlg(prompt,dlg_title,num_lines,defaultanswer);

num_pop           = str2num(x{1});
nr_parents        = str2num(x{2});
low_bound         = str2num(x{3});
up_bound          = str2num(x{4});
p_mu_ini          = str2num(x{5});
mu_low_b_ini      = str2num(x{6});
mu_up_b_ini       = str2num(x{7});
p_co              = str2num(x{8});
p_el              = str2num(x{9}); 
nr_generation     = str2num(x{10}); 
 
%% initial values

 
% First Generation
%num_pop = 30; %number of individuals in the first generation
z = 740; %timestep for the first population 
 
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

r_co        = num_pop/nr_parents; %Ratio between inital population and selected population
pp_1        = 1;
pp          = 1;
rc          = 1; % Running variable in the parenting choosing
mu_1        = 1;
mu_2        = 1;
n_el        = 1;
el_parII    = 1;
nr_par      = 16; % # of different parameters
mu_var      = 1;

 
%% First Generation -> Initiation of a random population

%Create r.  Range of r is a to b. Limits for the parameterchange
%low_bound = 0.8; % Lower (min) value
%up_bound = +1.2; % Upper (max) value

for j = 1:num_pop  %size of the first generation 
 
randVec1  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec2  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec3  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec4  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec5  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec6  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec7  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec8  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec9  = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec10 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec11 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1)); 
randVec12 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec13 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec14 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec15 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));
randVec16 = vertcat(1,low_bound + (up_bound-low_bound).*rand(j-1, 1));

 
nv.astrocyte.params.k_deg = randVec1(j,1)*x1;
nv.astrocyte.params.k_on = randVec2(j,1)*x2;
nv.astrocyte.params.K_inh = randVec3(j,1)*x3;
nv.astrocyte.params.r_h = randVec4(j,1)*x4;
nv.astrocyte.params.K_G = randVec5(j,1)*x5;
nv.astrocyte.params.K_I = randVec6(j,1)*x6;
nv.astrocyte.params.BK_end = randVec7(j,1)*x7;
nv.astrocyte.params.K_ex = randVec8(j,1)*x8;
nv.astrocyte.params.B_ex = randVec9(j,1)*x9;
nv.astrocyte.params.J_max = randVec10(j,1)*x10;
nv.astrocyte.params.K_act = randVec11(j,1)*x11;
nv.astrocyte.params.P_L = randVec12(j,1)*x12;
nv.astrocyte.params.V_max = randVec13(j,1)*x13;
nv.astrocyte.params.k_pump = randVec14(j,1)*x14;
nv.astrocyte.params.delta = randVec15(j,1)*x15;
nv.astrocyte.params.VR_ER_cyt = randVec16(j,1)*x16;

nv.simulate()  
 
results_c_k (j,:) = nv.out('c_k');
time_results_c_k = [T;results_c_k]';

results_R (j,:) = nv.out('R');
time_results_R = [T;results_R]';

end

Big_randMat = [randVec1, randVec2, randVec3, randVec4,...
randVec5,randVec6, randVec7, randVec8, randVec9, randVec10,randVec11, randVec12, randVec13, randVec14, randVec15,randVec16];
 
Max_matrix = [time_results_c_k(z,2:j+1)', Big_randMat];   
 
Results_end_max(1)  = max(time_results_c_k(z,2:j+1)');        %essential for plotting results
Results_end_mean(1) = mean(time_results_c_k(z,2:j+1)');       %essential for plotting results 
 
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher

sel_sum_c_k = sum(Max_matrix(1:j,1));

for i= 1:j
sel_prob = Max_matrix(i,1)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
end

% cummultative probability
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


%% Mutation of Parent Matrix

for mu_1 = 1:size(parent_matrix(:,1))
    
   for mu_2 = 1 : nr_par
       r_mu_1 = rand(1);
       
       if r_mu_1 <= p_mu_ini
                r_mu_2 = mu_low_b_ini + (mu_up_b_ini - mu_low_b_ini)*rand(1,1);                       %randomised mutation factor with boundaries
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                  %Mutation changes parameter at current position
       end
       
       mu_2 = mu_2 + 1;
   end
   mu_1 = mu_1 + 1;
   mu_2 = 1;
   
end

%% Cross-over

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
 

%% Elitism 
el_max_c_k(1)  = max(parent_matrix(:,4)) * (1-p_el); %Who is elite? ------- position after Parent choosing
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
end                                                      %output is a matrix which is 2*r_co bigger than matrix_input 
 
 
 %% Third and following generations                                                     
    
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher
for ga = 1 : nr_generation
    pp_1        = 1;           %Reset all running variables
    pp          = 1;
    rc          = 1;           % Running variable in the parenting choosing
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

end

 
 
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher
sel_sum_c_k = sum(parent_matrix(1:j,4));

for i= 1:j
sel_prob = parent_matrix(i,4)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
end

 
% cummultative probability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, parent_matrix(:,4:nr_par+4)]; %numbered parent_matrix with probabilities

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
 
%% Mutation 
for mu_1 = 1:length(parent_matrix(:,1))
    
   for mu_2 = 1 : nr_par
       
        mu_max      = max(parent_matrix(:,4));
        mu_mean     = mean(parent_matrix(:,4));    
        mu_current  = parent_matrix(mu_1,4);
       
       if mu_current >= mu_mean                                                     %makes mutation rate adaptive to quality of organism's fitness
            p_mu(mu_var) = 0.10 * ((mu_max - mu_current) / (mu_max - mu_mean));      %If organism is fitter than average -> mutation should decrease
       else 
            p_mu(mu_var) = 0.10;                                                     %if not mutation rate stays the same
       end
       
%        if mu_current == mu_max                                                      %Otherwise the best organism would not be mutated
%            p_mu(mu_var) = 0.02;
%        end
        
       r_mu_1 = rand(1);
       
       if r_mu_1 <= p_mu(mu_var)                                      %Decision if Mutation occurs or not
           
%            mu_factor_current    = parent_matrix(mu_1,mu_2 + 4);                %factor which is currently observed
%            mu_low_b             = mu_factor_current - p_mu(mu_var);                     %new boundaries are around factor
%            mu_up_b              = mu_factor_current + p_mu(mu_var); 
%            
%            if mu_low_b <= mu_low_b_ini                                      %initial boundaries are still important
%                mu_low_b = 0.8;
%            end
%            
%            if mu_up_b >= mu_up_b_ini 
%                mu_up_b = 1.2;
%            end
           
                r_mu_2 = mu_low_b_ini + (mu_up_b_ini - mu_low_b_ini)*rand(1,1);                       %randomised mutation factor with boundaries
                parent_matrix(mu_1,mu_2 + 4) = r_mu_2;                                  %Mutation changes parameter at current position
       end
       
       mu_2 = mu_2 + 1;
       mu_var = mu_var + 1;
   end
   mu_1 = mu_1 + 1;
   mu_2 = 1;
   
end
p_mu_total(ga) = mean(p_mu);                                                %in order to plot change of mutation over all generations
 
%% Cross-over
co_s = length(parent_matrix(:,1));                                %size of the column of the parents_matrix (should be initial population)

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

if cr_rand <= p_co                                         %probability of cross-over
    
cr_var1 = randi(nr_par)+4;                                         %defines the lower boundary of the crossover section
cr_var2 = randi([cr_var1,nr_par+4]);                               %defines the upper boundary of the crossover section

CO_1 = parent_matrix(co_var,cr_var1:cr_var2);                    %cutting out the cross-over section
CO_2 = parent_matrix(co_var+1,cr_var1:cr_var2); 
 
parent_matrix(co_var,cr_var1:cr_var2)   = CO_2;                   %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,cr_var1:cr_var2) = CO_1;

end 
co_var = co_var + 2;
end                                                              %output is a matrix which is 2*r_co bigger than matrix_input 

%% Elitism 
clear el_p;
el_max_c_k(ga+1)  = max(parent_matrix(:,4)) * (1-p_el); 
for var_el = 1:size(sel_output(:,1))
el_u=0;
     if sel_output(var_el,4) >= el_max_c_k(ga+1)
         el_p(var_el,:) = sel_output(var_el,:);
         for el_i = 1:size(parent_matrix(:,1));
             if el_p(var_el,1) == parent_matrix(el_i,1) && el_u== 0
                 parent_matrix(el_i,:) = el_p(var_el,:);
                 el_u=1;
             end
         end       
     end
end

Results_end_max(ga+1)  = max  (parent_matrix(:,4));        %essential for plotting results
Results_end_mean(ga+1) = mean (parent_matrix(:,4));        %essential for plotting results

end

%% Getting Ca2+ 

g = 1;

for g=1:length(parent_matrix(:,1))

 
    
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
parent_matrix(g,nr_par+5)   =  results_R_gen(g,z);
end

%% Plotting results

Results_end_max(ga+2)  = max(results_c_k_gen(:,z)');        %essential for plotting results
Results_end_mean(ga+2) = mean(results_c_k_gen(:,z)');       %This is located afeter each generation to save the CA result


figure(1)
v_initial_c_k = linspace(initial_c_k,initial_c_k,nr_generation+2);
plot(linspace(1,nr_generation+2,nr_generation+2)', Results_end_max)
    legend('Best [Ca2+] of each Generation','Initial [Ca2+]')
    title('Best CA^2^+ concentration in AC of every generation'); 
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    hold on; 
    grid on;

    plot(linspace(1,nr_generation+2,nr_generation+2)',v_initial_c_k);
    legend('Best [Ca2+] of each Generation','Initial [Ca2+]');
    title('Best CA^2^+ concentration in AC of every generation');
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    grid on;
    grid minor;
    
figure(2)
plot(linspace(1,nr_generation+2,nr_generation+2)', Results_end_mean)
    legend('Average [Ca2+] of each Generation','Initial [Ca2+]')
    title('Average CA^2^+ concentration in AC of every generation'); 
    xlabel('# of generation'); 
    ylabel('[Ca^2^+][\muM]'); 
    hold on; 
    grid on;

    plot(linspace(1,nr_generation+2,nr_generation+2)',v_initial_c_k);
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
