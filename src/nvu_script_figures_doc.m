%% Demonstration script for 00-NVU model + Calcium Model
% This script runs the NVU Model under pre-defined base conditions. The
% purpose of this script is to give a quick overview of the outcomes of the
% code for first time viewers. 

% Please note this is not an exhaustive instruction list. Please refer to
% nvu_script for more detail.

%% NVU Model Overview
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are specified when the 
% modules are constructed. Note: if theses are to be changed please refer 
% to the nvu_script for more detail.
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


%% Plots
% Lists of quantities that can be retrieved from the NVU model after
% simulation are given in the documentation pages of the individual model
% components.

%%

%Moritz' Part of plotting
figure(1) 
set(gcf, 'Position', [70 70 1800 900]);
set(gcf, 'Name', 'Neurovascular Coupling Overview');
subplot(3,3,1)
plot(nv.T, 1e3 * nv.out('K_s')); xlim([0 500])
xlabel('Time [s]'); ylabel('K_s [mM]'); grid on
subplot(3,3,2)
plot(nv.T, 1e3 * nv.out('v_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('v_k [mV]'); grid on
subplot(3,3,3)
plot(nv.T, nv.out('J_BK_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{BK_k} [\muM/s]'); grid on
subplot(3,3,4)
plot(nv.T, nv.out('K_p')); xlim([0 500])
xlabel('Time [s]'); ylabel('[K_p] [\muM]'); grid on
subplot(3,3,5)
plot(nv.T, nv.out('J_KIR_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{KIR_i} [\muM/s]'); grid on
subplot(3,3,6)
plot(nv.T, nv.out('v_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('v_i [mV]'); grid on
subplot(3,3,7)
plot(nv.T, nv.out('J_VOCC_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{VOCC_i} [\muM/s]'); grid on
subplot(3,3,8)
plot(nv.T, nv.out('Ca_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('[Ca_i] [\muM]'); grid on
subplot(3,3,9)
plot(nv.T, nv.out('R')*1e6); xlim([0 500])
xlabel('Time [s]'); ylabel('R [\mum]'); grid on

%%


figure(2) 
set(gcf, 'Position', [70 40 1800 900]);
set(gcf, 'Name', 'AC State Variables')
subplot(4,2,1)
plot(nv.T, 0.001*nv.out('K_s')); xlim([0 500])
xlabel('Time [s]'); ylabel('[K^+]_s [mM]'); grid on
subplot(4,2,2)
plot(nv.T, 0.001*nv.out('N_K_k')./(nv.out('R_k'))); xlim([0 500])
xlabel('Time [s]'); ylabel('[K^+]_k (mM)') ; grid on
subplot(4,2,3)
plot(nv.T, nv.out('v_k')*1e3); xlim([0 500])
xlabel('Time [s]');ylabel('v_k [mV]'); grid on
subplot(4,2,4)
plot(nv.T, 0.001*nv.out('N_Na_k')./(nv.out('R_k'))); xlim([0 500])
xlabel('Time [s]'); ylabel('[Na^+]_k (mM)'); grid on
subplot(4,2,5)
plot(nv.T, 0.001*nv.out('N_HCO3_k')./(nv.out('R_k'))); xlim([0 500])
xlabel('Time [s]'); ylabel('[HCO_3^-]_k (mM)'); grid on
subplot(4,2,6)
plot(nv.T, 0.001*nv.out('N_Cl_k')./(nv.out('R_k'))); xlim([0 500])
xlabel('Time [s]'); ylabel('[Cl^-]_k (mM)'); grid on
subplot(4,2,7)
plot(nv.T, 0.001*nv.out('K_p')); xlim([0 500])
xlabel('Time [s]'); ylabel('[K^+]_p (mM)'); grid on
%%

figure(3)
set(gcf, 'Position', [70 10 1800 900]);
set(gcf, 'Name', 'AC Fluxes');
subplot(2,4,1)
plot(nv.T, nv.out('J_KCC1_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{KCC1_k} [muM/s]'); grid on
subplot(2,4,2)
plot(nv.T, nv.out('J_NKCC1_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{NKCC1_k} [muM/s]'); grid on
subplot(2,4,3)
plot(nv.T, nv.out('J_K_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{K_k} [muM/s]'); grid on
subplot(2,4,4)
plot(nv.T, nv.out('J_NaK_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{NaK_k} [muM/s]'); grid on
subplot(2,4,5)
plot(nv.T, nv.out('J_Na_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{Na_k} [muM/s]'); grid on
subplot(2,4,6)
plot(nv.T, nv.out('J_NBC_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{NBC_k} [muM/s]'); grid on
subplot(2,4,7)
plot(nv.T, nv.out('J_BK_k')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{BK_k} [muM/s]'); grid on

%%
figure(4) 
set(gcf, 'Position', [70 -20 1800 900]);
set(gcf, 'Name', 'SMC/EC State Variables and Coupling');
subplot(4,3,1)
plot(nv.T, nv.out('Ca_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('Ca_i [muM]'); grid on
subplot(4,3,2)
plot(nv.T, nv.out('s_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('s_{i} [\muM]'); grid on
subplot(4,3,3)
plot(nv.T, nv.out('v_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('v_i [mV]'); grid on
subplot(4,3,4)
plot(nv.T, nv.out('I_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('I_i [\muM/s]'); grid on
subplot(4,3,5)
plot(nv.T, nv.out('K_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('K_i [muM]'); grid on
subplot(4,3,6)
plot(nv.T, nv.out('Ca_j')); xlim([0 500])
xlabel('Time [s]'); ylabel('[Ca^{2+}]_j [\muM]'); grid on
subplot(4,3,7)
plot(nv.T, nv.out('s_j')); xlim([0 500])
xlabel('Time [s]'); ylabel('s_j [\muM]'); grid on
subplot(4,3,8)
plot(nv.T, nv.out('v_j')); xlim([0 500])
xlabel('Time [s]'); ylabel('v_j [mV]'); grid on
subplot(4,3,9)
plot(nv.T, nv.out('I_j')); xlim([0 500])
xlabel('Time [s]'); ylabel('I_j [\muM]'); grid on
subplot(4,3,10)
plot(nv.T, nv.out('V_coup_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('V_{coup_i} [mV/s]'); grid on
subplot(4,3,11)
plot(nv.T, nv.out('J_IP3_coup_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{IP3 coup_i} [\muM/s]'); grid on
subplot(4,3,12)
plot(nv.T, nv.out('J_Ca_coup_i')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{Ca coup_i} [\muM/s]'); grid on



figure(5)
set(gcf, 'Position', [70 -50 1800 900]);
set(gcf, 'Name', 'SMC Fluxes');
subplot(4,3,1)
plot(nv.T, nv.out('J_KIR_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{KIR_i} [\muM/s]'); grid on
subplot(4,3,2)
plot(nv.T, nv.out('J_NaK_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{Na/K_i} [\muM/s]'); grid on
subplot(4,3,3)
plot(nv.T, nv.out('J_K_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{K_i} [\muM/s]'); grid on
subplot(4,3,4)
plot(nv.T, nv.out('J_extrusion_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{extrusion_i} [\muM/s]'); grid on
subplot(4,3,5)
plot(nv.T, nv.out('J_VOCC_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{VOCC_i} [\muM/s]'); grid on
subplot(4,3,6)
plot(nv.T, nv.out('J_stretch_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{stretch_i} [\muM/s]'); grid on
subplot(4,3,7)
plot(nv.T, nv.out('J_NaCa_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{Na/Ca_i} [\muM/s]'); grid on
subplot(4,3,8)
plot(nv.T, nv.out('J_Cl_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{Cl_i} [\muM/s]'); grid on
subplot(4,3,9)
plot(nv.T, nv.out('J_CICR_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{CICR_i} [\muM/s]'); grid on
subplot(4,3,10)
plot(nv.T, nv.out('J_IP3_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{IP3_i} [\muM/s]'); grid on
subplot(4,3,11)
plot(nv.T, nv.out('J_SR_leak_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{SR\_leak_i} [\muM/s]'); grid on
subplot(4,3,12)
plot(nv.T, nv.out('J_SR_uptake_i')); xlim([0 500])
xlabel('Time [s]'); ylabel('J_{SR\_uptake_i} [\muM/s]'); grid on

%%
 
figure(6)
set(gcf, 'Position', [70 -80 1800 900]);
set(gcf, 'Name', 'EC Fluxes')
subplot(4,3,1)
plot(nv.T, nv.out('J_cation_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{cation_j} [\muM/s]'); grid on
subplot(4,3,2)
plot(nv.T, nv.out('J_extrusion_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{extrusion_j} [\muM/s]'); grid on
subplot(4,3,3)
plot(nv.T, nv.out('J_stretch_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{stretch_j} [\muM/s]'); grid on
subplot(4,3,4)
plot(nv.T, nv.out('J_K_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{K_j} [\muM/s]'); grid on
subplot(4,3,5)
plot(nv.T, nv.out('J_R_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{R_j} [\muM/s]'); grid on 
subplot(4,3,6)
plot(nv.T, nv.out('J_CICR_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{CICR_j} [\muM/s]'); grid on
subplot(4,3,7)
plot(nv.T, nv.out('J_IP3_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{_IP3_j} [\muM/s]');grid on
subplot(4,3,8)
plot(nv.T, nv.out('J_ER_leak_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{ER leak_j} [\muM/s]'); grid on
subplot(4,3,9)
plot(nv.T, nv.out('J_ER_uptake_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{ER uptake_j} [\muM/s]'); grid on
subplot(4,3,10)
plot(nv.T, nv.out('J_degrad_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' J_{degrad_j} [\muM/s]'); grid on

% NEED TO CHECK THESE 2 IN DOCUMENTATION! --> Koenigsberger
subplot(4,3,12)
plot(nv.T, nv.out('J_BK_Ca_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' I_{BK Ca_j} [\muM/s]'); grid on
subplot(4,3,11)
plot(nv.T, nv.out('J_SK_Ca_j')); xlim([0 500])
xlabel('Time [s]'); ylabel(' I_{SK Ca_j} [\muM/s]'); grid on
%%

figure(7)
set(gcf, 'Position', [70 -110 1800 900]);
set(gcf, 'Name', 'Contraction Model and Radius')
subplot(3, 2, 1)
plot(nv.T, nv.out('M')); xlim([0 500])
xlabel('Time [s]'); ylabel('M [-]'); grid on
subplot(3, 2, 2)
plot(nv.T, nv.out('Mp')); xlim([0 500])
xlabel('Time [s]'); ylabel('Mp [-]'); grid on
subplot(3, 2, 3)
plot(nv.T, nv.out('AMp')); xlim([0 500])
xlabel('Time [s]'); ylabel('AMp [-]'); grid on
subplot(3, 2, 4)
plot(nv.T, nv.out('AM')); xlim([0 500])
xlabel('Time [s]'); ylabel('AM [-]'); grid on
subplot(3, 2, 5)
plot(nv.T, nv.out('F_r')); xlim([0 500])
xlabel('Time [s]'); ylabel('Fraction [-]'); grid on
subplot(3, 2, 6)
plot(nv.T, 1e6*nv.out('R')); xlim([0 500])
xlabel('Time [s]'); ylabel('R [\mum]'); grid on
%%


save1 = input('Would you like to save all figures and parameters in a separate folder? (yes [y] / no [n] / png [p]): ', 's');
if strcmp(save1,'y') || strcmp(save1,'yes')
        fprintf('Please wait!');
        cl = clock;
        filename = ['Results-',date,'-',num2str(cl(4),'%.2d'),':',num2str(cl(5),'%.2d'),':',num2str(cl(6),'%.0f')];
        mkdir(filename);
        cd(filename);
        saveas(figure(1),'1 Neurovascular Coupling Overview.fig');
        saveas(figure(2),'2 AC State Variables.fig');
        saveas(figure(3),'3 AC Fluxes.fig');
        saveas(figure(4),'4 SMC EC State Variables and Coupling.fig');
        saveas(figure(5),'5 SMC Fluxes.fig');
        saveas(figure(6),'6 EC Fluxes.fig');
        saveas(figure(7),'7 Contraction Model and Radius.fig');
        
        file1=fopen('Parameters.txt','w');
        fprintf(file1, ['Parameters of "' filename '":\r\n\r\n']);
        fprintf(file1, 'J_PLC      = %.3f      [uM/s]\r\n', nv.smcec.params.J_PLC);
        fclose(file1);
        clc; fprintf('Figures and parameters have been saved in the folder "%s".\n', filename);
        value=1;
        cd ..
elseif strcmp(save1,'n') || strcmp(save1,'no')
        clc; fprintf('Figures and parameters have not been saved.\n');
        value=1;
elseif strcmp(save1,'p') || strcmp(save1,'p')
        fprintf('Please wait!');
        cl = clock;
        filename = ['Results-',date,'-',num2str(cl(4),'%.2d'),':',num2str(cl(5),'%.2d'),':',num2str(cl(6),'%.0f')];
        mkdir(filename);
        cd(filename);
        saveas(figure(1),'1 Neurovascular Coupling Overview.png');
        saveas(figure(2),'2 AC State Variables.png');
        saveas(figure(3),'3 AC Fluxes.png');
        saveas(figure(4),'4 SMC EC State Variables and Coupling.png');
        saveas(figure(5),'5 SMC Fluxes.png');
        saveas(figure(6),'6 EC Fluxes.png');
        saveas(figure(7),'7 Contraction Model and Radius.png');
        file1=fopen('Parameters.txt','w');
        fprintf(file1, ['Parameters of "' filename '":\r\n\r\n']);
        fprintf(file1, 'J_PLC      = %.3f      [uM/s]\r\n', nv.smcec.params.J_PLC);
        fclose(file1);
        clc; fprintf('Figures and parameters have been saved in the folder "%s".\n', filename);
        value=1;
        cd ..
else
        clc; fprintf('Warning: "%s" is not a valid input!\n', save1);
end

