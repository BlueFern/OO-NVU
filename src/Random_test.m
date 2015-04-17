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
 
% Create r.  Range of r is a to b.
a = 0.8; % Lower (min) value.
b = +1.2; % Upper (max) value

x1 = nv.astrocyte.params.k_deg;
x2 = nv.astrocyte.params.k_on;
x3 = nv.astrocyte.params.K_inh;
x4 = nv.astrocyte.params.r_h;
x5 = nv.astrocyte.params.K_I;

for j = 1:50

randVec1 = vertcat(1,a + (b-a).*rand(j-1, 1)); 
randVec2 = vertcat(1,a + (b-a).*rand(j-1, 1));
randVec3 = vertcat(1,a + (b-a).*rand(j-1, 1));
randVec4 = vertcat(1,a + (b-a).*rand(j-1, 1));
randVec5 = vertcat(1,a + (b-a).*rand(j-1, 1));
    
nv.simulate()  
 
nv.astrocyte.params.k_deg = randVec1(j,1)*x1;
r1(j,:) = nv.astrocyte.params.k_deg;

nv.astrocyte.params.k_on = randVec2(j,1)*x2;
r2(j,:) = nv.astrocyte.params.k_on;

nv.astrocyte.params.K_inh = randVec3(j,1)*x3;
r3(j,:) = nv.astrocyte.params.K_inh;

nv.astrocyte.params.r_h = randVec4(j,1)*x4;
r4(j,:) = nv.astrocyte.params.r_h;

nv.astrocyte.params.k_deg = randVec5(j,1)*x5;
r5(j,:) = nv.astrocyte.params.K_I;

results_c_k (j,:) = nv.out('c_k');
time_results_c_k = [T;results_c_k]';

results_R (j,:) = nv.out('R');
time_results_R = [T;results_R]';

figure (1)
plot(nv.T,nv.out('c_k'))
    %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different Ca^{2+} concentrations in AC due to Variations of Parameters'); 
    xlabel('time [s]'); 
    ylabel('[Ca^2^+][\muM]'); 
    hold on; 
    grid on;


end
resultsparams = [r1,r2,r3,r4,r5];

G= 1:j;

figure (2)
subplot(3,1,1)
plot(G,randVec1(:),'-x')
    hold on
plot(G,randVec2(:),'-x')
    hold on
plot(G,randVec3(:),'-x')
    hold on
plot(G,randVec4(:),'-x')
    hold on
plot(G,randVec5(:),'-x')
    hold on
    
    legend('k_{deg}','k_{on}','K_{inh}','r_h','K_I')
    title('percentual change of each parameter');
    xlabel('# of calculation [-]'); 
    ylabel('Percentage [-]'); 
    grid on;
    grid minor;

subplot(3,1,2)
plot(G, time_results_c_k(740,2:j+1),'-x')
    %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different C^{2+} concentrations in AC due to variations of the parameters'); 
    xlabel('# of calculation [-]'); 
    ylabel('[Ca^2^+] [\muM]'); 
    hold on; 
    grid on;
    grid minor;

subplot(3,1,3)
plot(G, time_results_R(740,2:j+1),'-x')
    %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different diameters R of the vessel due to variations of the parameters'); 
    xlabel('# of calculation [-]'); 
    ylabel('Diameter R in [m]'); 
    hold on; 
    grid on;
    grid minor;

Max_matrix = [time_results_c_k(740,2:j+1)', randVec1, randVec2, randVec3, randVec4, randVec5 ];

Sort_max_matrix = sortrows(Max_matrix);
    
q = 10;  %number of highest concentrations you want
figure (3)
pp = j-q+1;
z=1;
for z = 1:5
    plot(Sort_max_matrix(pp:j,1),Sort_max_matrix(pp:j,z+1),'-x') 
    hold on
    z = z+1;
end


    
legend('k_{deg}','k_{on}','K_{inh}','r_h','K_I')
title('percentual change of each parameter');
xlabel('Increasing order of [Ca^{2+}] [\muM]'); 
ylabel('Percentage [-]'); 
grid on;
grid minor;

