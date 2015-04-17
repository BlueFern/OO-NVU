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
 
 
 
%% initial values

 
% First Generation

x1 = nv.astrocyte.params.k_deg; %orginal values for the parameters used in the AC
x2 = nv.astrocyte.params.k_on;
x3 = nv.astrocyte.params.K_inh;
x4 = nv.astrocyte.params.r_h;
x5 = nv.astrocyte.params.K_G;
%j = for loop, size of the first generation

 
%Selection

%n = number of parents

%% Parentchoosing,Crossover

rc   = 1; % Running variable in the parenting choosing
r_co = 10; %Ratio between inital population and selected population
pp_1 = 1;
pp   = 1;
nr_par = 16; % # of different parameters
p_co = 0.7; %Probability of cross over
% co = size of parent-matrix

 
 
%% First Generation -> Initiation of a random population

%Create r.  Range of r is a to b. Limits for the parameterchange
low_bound = 0.8; % Lower (min) value
b = +1.2; % Upper (max) value

for j = 1:30 %size of the first generation 
 
randVec1 = vertcat(1,low_bound + (b-low_bound).*rand(j-1, 1)); 
randVec2 = vertcat(1,low_bound + (b-low_bound).*rand(j-1, 1));
randVec3 = vertcat(1,low_bound + (b-low_bound).*rand(j-1, 1));
randVec4 = vertcat(1,low_bound + (b-low_bound).*rand(j-1, 1));
randVec5 = vertcat(1,low_bound + (b-low_bound).*rand(j-1, 1));
    
nv.simulate()  
 
nv.astrocyte.params.k_deg = randVec1(j,1)*x1;
r1(j,:) = nv.astrocyte.params.k_deg;

nv.astrocyte.params.k_on = randVec2(j,1)*x2;
r2(j,:) = nv.astrocyte.params.k_on;

nv.astrocyte.params.K_inh = randVec3(j,1)*x3;
r3(j,:) = nv.astrocyte.params.K_inh;

nv.astrocyte.params.r_h = randVec4(j,1)*x4;
r4(j,:) = nv.astrocyte.params.r_h;

nv.astrocyte.params.K_G = randVec5(j,1)*x5;
r5(j,:) = nv.astrocyte.params.K_G;

results_c_k (j,:) = nv.out('c_k');
time_results_c_k = [T;results_c_k]';

results_R (j,:) = nv.out('R');
time_results_R = [T;results_R]';

% figure (1)
% plot(nv.T,nv.out('c_k'))
%     %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
%     title('Different Ca^{2+} concentrations in AC due to Variations of Parameters'); 
%     xlabel('time [s]'); 
%     ylabel('[Ca^2^+][\muM]'); 
%     hold on; 
%     grid on;

 
end
resultsparams = [r1,r2,r3,r4,r5];
Big_randMat = [randVec1, randVec2, randVec3, randVec4, randVec5];

G= 1:j;

% figure (2)
% subplot(3,1,1)
% plot(G,Big_randMat(:,:),'-x')
%     
%     legend('k_{deg}','k_{on}','K_{inh}','r_h','k_{G}')
%     title('percentual change of each parameter');
%     xlabel('# of calculation [-]'); 
%     ylabel('Percentage [-]'); 
%     grid on;
%     grid minor;
% 
% subplot(3,1,2)
% plot(G, time_results_c_k(740,2:j+1),'-x')
%     %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
%     title('Different C^{2+} concentrations in AC due to variations of the parameters'); 
%     xlabel('# of calculation [-]'); 
%     ylabel('[Ca^2^+] [\muM]'); 
%     hold on; 
%     grid on;
%     grid minor;
% 
% subplot(3,1,3)
% plot(G, time_results_R(740,2:j+1),'-x')
%     %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
%     title('Different diameters R of the vessel due to variations of the parameters'); 
%     xlabel('# of calculation [-]'); 
%     ylabel('Diameter R in [m]'); 
%     hold on; 
%     grid on;
%     grid minor;

Max_matrix = [time_results_c_k(740,2:j+1)', randVec1, randVec2, randVec3, randVec4, randVec5 ];

Sort_max_matrix = sortrows(Max_matrix);
    
% q = 3;  %number of highest concentrations you want
% figure (3)
% pp = j-q+1;
% z=1;
% for z = 1:5 %number of parameters
%     plot(Sort_max_matrix(pp:j,1),Sort_max_matrix(pp:j,z+1),'x') 
%     hold on
%     z = z+1;
% end

[sorted_c_k,] = sort(time_results_c_k(740,2:j+1),'descend');  %/sortingIndices
Max_c_k = sorted_c_k(:,1:q)';

y1 = Sort_max_matrix((j-q):j,2);
y2 = Sort_max_matrix((j-q):j,2);
y3 = Sort_max_matrix((j-q):j,2);
y4 = Sort_max_matrix((j-q):j,2);
y5 = Sort_max_matrix((j-q):j,2);

% legend('k_{deg}','k_{on}','K_{inh}','r_h','k_{deg}')
% title('percentual change of each parameter');
% xlabel('CA^{2+} []'); 
% ylabel('Percentage of parameter change [-]'); 
% grid on;
% grid minor;

 
%% Selection (roulette wheel method)-> sel_...

%probability of each parameterset to be choosen for the offspring
%population -> probability for parametersets with a high
%calciumconcentration is higher

sel_sum_c_k = sum(Max_matrix(1:j,1));

for i= 1:j
sel_prob = Max_matrix(i,1)./sel_sum_c_k;
sel_probvec(i,:)= sel_prob;
i = i+1;
end

% cummultative porbability
sel_cum_probvec = cumsum(sel_probvec);

sel_max_matrix = [linspace(1,j,j)',sel_cum_probvec,sel_probvec, Max_matrix]; %numberd max_matrix with probabilities

% decision

for n = 1:10 %number of the potential parents
p=1  
sel_rand_value = rand(1,1)%random vector between 0 and 1 
sel_rand_value_save(n,:)= sel_rand_value %save the vector to compare

if sel_rand_value < sel_max_matrix(p,2);
        sel_output(n,:) = sel_max_matrix(p,:);
end
while sel_rand_value > sel_max_matrix(p,2);
    p=p+1;
    if sel_rand_value < sel_max_matrix(p,2);
        sel_output(n,:) = sel_max_matrix(p,:);
    end
end
end

%% Selecting Parents
% # of calculations
r = randperm(n);%creates a shuffled linspace vector from 1 to n, but without duplicates -> so a parent can just be parent once
for pp = 1:n/2    %divided by 2, because there are 2 parents   
    pp_2 = pp_1 + 1;
    
    P1 = sel_output(r(pp_1),:);  %choosing parents out of input matrix
    P2 = sel_output(r(pp_2),:);
    
    for sp = 1:r_co 
        
        parent_matrix(rc,:)   = P1;
        parent_matrix(rc+1,:) = P2;
        
        rc = rc + 2;
        
    end 
   pp_1 = pp_1 + 2;
end 
 
%%Mutation of Parent Matrix

%Cross-over
co_s = size(parent_matrix(:,1));                                %size of the column of the parents_matrix (should be initial population)
for co_var = 1:co_s
cr_rand = rand(1);

if cr_rand <= p_co                                               %probability of cross-over
cr_var1 = randi(nr_par);                                         %defines the lower boundary of the crossover section
cr_var2 = randi(numel(linspace(cr_var1,nr_par,nr_par-cr_var1))); %defines the upper boundary of the crossover section

CO_1 = parent_matrix(co_var,cr_var1:cr_var2);                    %cutting out the cross-over section
CO_2 = parent_matrix(co_var+1,cr_var1:cr_var2); 
 
parent_matrix(co_var,cr_var1:cr_var2)  = CO_2;                   %change of the cross-over section between the 2 childrens
parent_matrix(co_var+1,cr_var1:cr_var2) = CO_1;

end 
co_var = co_var + 2;
end                                                              %output is a matrix which is 2*r_co bigger than matrix_input 




    
   

