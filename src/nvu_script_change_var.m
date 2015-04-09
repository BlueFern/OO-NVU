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
a = -1; % Lower (min) value.
b = 1; % Upper (max) value
c = 1;
d =2;
x1 = nv.astrocyte.params.g_BK_k;
x2 = nv.astrocyte.params.g_NBC_k;
x3 = nv.astrocyte.params.g_K_k;
x4 = nv.astrocyte.params.g_Cl_k;
x5 = nv.astrocyte.params.g_TRPV_k;
for j = 1:30
randVec1 = vertcat(1,1*10.^(a + (b-a).*rand(j-1, 1)));
randVec2 = vertcat(1,1*10.^(a + (b-a).*rand(j-1, 1)));
randVec5 = vertcat(1,1*10.^(a + (b-a).*rand(j-1, 1)));
randVec4 = vertcat(1,1*10.^(a + (b-a).*rand(j-1, 1)));
randVec3 = vertcat(1,1*10.^(a + (b-a).*rand(j-1, 1)));

   
nv.simulate() 
 
nv.astrocyte.params.g_BK_k = randVec1(j,1)*x1;
r1(j,:) = nv.astrocyte.params.g_BK_k;
nv.astrocyte.params.g_NBC_k = randVec2(j,1)*x2;
r2(j,:) = nv.astrocyte.params.g_NBC_k;
nv.astrocyte.params.g_K_k = randVec3(j,1)*x3;
r3(j,:) = nv.astrocyte.params.g_K_k;
nv.astrocyte.params.g_Cl_k = randVec4(j,1)*x4;
r4(j,:) = nv.astrocyte.params.g_Cl_k;
nv.astrocyte.params.g_TRPV_k = randVec5(j,1)*x5;
r5(j,:) = nv.astrocyte.params.g_TRPV_k;
results_c_k (j,:) = nv.out('v_k');
time_results_c_k = [T;results_c_k]';
results_R (j,:) = nv.out('R');
time_results_R = [T;results_R]';
figure (1)
plot(nv.T,nv.out('v_k'))
    legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different membrane voltage of the Astrocyte due to Variations of Parameters');
    xlabel('time [s]');
    ylabel('[v_K][V]');
    hold all;
    grid on;
figure(2)
        %suptitle('The Input Signal')
        subplot(3,2,1)
        hold all;
        plot(nv.T, nv.astrocyte.params.g_TRPV_k*nv.out('z_k')'.*nv.out('E_TRPV_k'))
        title('E_TRPV_k')
        xlabel('Time [s]'); ylabel('[V]')

        subplot(3,2,2)
        hold all;
       plot(nv.T, nv.astrocyte.params.g_Na_k*nv.out('E_Na_k'))
        title('E_Na_k')
        xlabel('Time [s]'); ylabel('[V]')

        subplot(3,2,3)
        hold all;
       plot(nv.T, nv.astrocyte.params.g_BK_k*nv.out('w_k')'.*nv.out('E_BK_k'))
        title('E_BK_k')
        xlabel('Time [s]'); ylabel('[V]')
        
        subplot(3,2,4)
        hold all;
       plot(nv.T, nv.astrocyte.params.g_NBC_k*nv.out('E_NBC_k'))
        title('E_NBC_k')
        xlabel('Time [s]'); ylabel('[V]')

        
        subplot(3,2,5)
        hold all;
     plot(nv.T, nv.astrocyte.params.g_Cl_k*nv.out('E_Cl_k'))
        title('E_Cl_k')
        xlabel('Time [s]'); ylabel('[V]')
        subplot(3,2,6)
       plot(nv.T, nv.astrocyte.params.g_K_k*nv.out('E_K_k'))
       hold all;
        title('E_K_k')
        xlabel('Time [s]'); ylabel('[V]')
       
end
resultsparams = [r1,r2,r3,r4,r5];
G= 1:j;
figure (3)
subplot(3,1,1)
plot(G,randVec1(:),'-x')
    hold all
plot(G,randVec2(:),'-x')
    hold all
plot(G,randVec3(:),'-x')
    hold all
plot(G,randVec4(:),'-x')
    hold all
plot(G,randVec5(:),'-x')
    hold all
   
    legend('g_BK','g_NBC','g_K','g_Cl','g_Na')
    title('percentual change of each parameter');
    xlabel('# of calculation [-]');
    ylabel('Percentage [-]');
    grid on;
    grid minor;
subplot(3,1,2)
plot(G, time_results_c_k(740,2:j+1),'-x')
    %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different v_k AC due to variations of the parameters');
    xlabel('# of calculation [-]');
    ylabel('[v_k] [V]');
    hold all;
    grid on;
    grid minor;
subplot(3,1,3)
plot(G, time_results_R(740,2:j+1),'-x')
    %legend('variation 1','variation 2','variation 3','variation 4','variation 5','variation 6','variation 7','variation 8','variation 9')
    title('Different diameters R of the vessel due to variations of the parameters');
    xlabel('# of calculation [-]');
    ylabel('Diameter R in [m]');
    hold all;
    grid on;
    grid minor;
Max_matrix = [time_results_c_k(740,2:j+1)', randVec1, randVec2, randVec3, randVec4, randVec5 ];
Sort_max_matrix = sortrows(Max_matrix);
   
q = 10;  %number of highest concentrations you want
figure (4)
pp = j-q+1;
z=1;
for z = 1:5
    plot(Sort_max_matrix(pp:j,1),Sort_max_matrix(pp:j,z+1),'-x')
    hold all
    z = z+1;
end
 
    
legend('g_BK','g_NBC','g_K','g_Cl','g_Na')
title('percentual change of each parameter');
xlabel('Increasing order of [v_k] [V]');
ylabel('Percentage [-]');
grid on;
grid minor;