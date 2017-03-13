%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model:
% K+, NO, Astrocytic Ca2+, TRPV4, ECS

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all
clc

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 0; XLIM2 = 500;
FIG_NUM = 10;


NEURONAL_START  = 100;      % Start of neuronal stimulation
NEURONAL_END    = 300;      % End of neuronal stimulation 
ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 200000000;      % End of ECS K+ input
J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

ECS             = 1;        % Include ECS compartment or not
CA_SWITCH       = 1;        % Turn on Ca pathway
NO_INPUT_SWITCH = 1;        % Turn on NO stimulation
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
K_SWITCH        = 0;        % Turn on K+ input into SC
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
J_STRETCH_FIX   = -1;       % 1: old (wrong) , -1: new (correct), J_stretch is - instead

nv = NVU(Neuron('F_input', 2.67, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'GluSwitch', NO_INPUT_SWITCH, 'NOswitch', NO_PROD_SWITCH, 'KSwitch', K_SWITCH), ...
    Astrocyte('J_max', 2880*1.1, 'r_h', 4.8, 'ECS_input', 9, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END, 'GluSwitch', CA_SWITCH, 'ECSswitch', ECS, 'PVStoECS', 0, 'SCtoECS', 1), ...
    WallMechanics(), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH, 'stretchFix', J_STRETCH_FIX), 'odeopts', odeopts);

nv.T = linspace(0, XLIM2*2, 5000);

nv.simulate()

% Run this to take the ICs as the end of the last simulation run
ICs = (nv.U(end, :))';
nv.u0 = ICs;
nv.simulateManualICs()

%legend('Ca^{2+}', 'TRPV4', 'Ca^{2+} and TRPV4')
% legend('K^+', 'K^+ and Ca^{2+}', 'K^+ and Ca^{2+} and TRPV4')
% legend('NO', 'NO and Ca^{2+}', 'NO and Ca^{2+} and TRPV4')
% legend('K^+', 'NO', 'K^+ and NO', 'K^+ and NO and Ca^{2+} and TRPV4')
% 
% 
% figure(FIG_NUM);
% 
% subplot(2,2,1)
% hold all;
% % plot(nv.T, nv.out('v_k')*1e3, 'k', 'LineWidth', 1);
% % ylabel('v_k [mV]');
% plot(nv.T, nv.out('I_k'), 'k', 'LineWidth', 1);
% ylabel('IP3_k [mM]');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,2,2)
% hold all;
% plot(nv.T, nv.out('Ca_k'), 'k', 'LineWidth', 1);
% ylabel('Ca_k [mM]');
% % plot(nv.T, nv.out('J_VOCC_i'));
% % ylabel('VOCC flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,2,3)
% hold all;
% plot(nv.T, nv.out('eet_k'), 'k', 'LineWidth', 1);
% ylabel('eet_k [mM]');
% % plot(nv.T, nv.out('J_KIR_i'));
% % ylabel('KIR flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,2,4)
% hold all;
% plot(nv.T, nv.out('w_k'), 'k', 'LineWidth', 1);
% ylabel('w_k [-]');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 


% 
% 
figure(FIG_NUM+1);
% 
% 'k:' 'k--' 'k'

subplot(2,3,1)
hold all;
% plot(nv.T, nv.out('Cak'), 'k--', 'LineWidth', 1);
 plot(nv.T, nv.out('Ca_k'), 'k-', 'LineWidth', 1);
ylabel('Ca_k [\muM]');
xlim([XLIM1 XLIM2])
xlabel('Time [s]'); 
% % 
% subplot(2,3,2)
% hold all;
% plot(nv.T, nv.out('eet_k'), 'k', 'LineWidth', 1);
% ylabel('eet_k [\muM]');
% % plot(nv.T, nv.out('J_KIR_i'));
% % ylabel('KIR flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,3,3)
% hold all;
% %plot(nv.T, nv.out('v_k')*1e3, 'k', 'LineWidth', 1);
% %ylabel('v_k [mV]');
%  plot(nv.T, nv.out('I_k'), 'k', 'LineWidth', 1);
%  ylabel('IP3_k [\muM]');
%  xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,3,4)
% hold all;
% plot(nv.T, nv.out('w_k'), 'k', 'LineWidth', 1);
% ylabel('w_k [-]');
% % plot(nv.T, nv.out('J_BK_k'));
% % ylabel('BK flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,3,5)
% hold all;
% plot(nv.T, nv.out('K_p')/1e3, 'k', 'LineWidth', 1);
% xlabel('Time [s]'); 
% ylabel('K_p [mM]');
% xlim([XLIM1 XLIM2])
% 
% subplot(2,3,6)
% hold all;
% plot(nv.T, nv.out('R')*1e6, 'k', 'LineWidth', 1)
% xlabel('Time [s]'); 
% ylabel('Radius [\mum]');
% xlim([XLIM1 XLIM2])

% 
% 
% 
% 

% figure(FIG_NUM);
% hold all;
% plot(nv.T, nv.out('J_BK_k'))
% ylabel('BK flux')

figure(FIG_NUM+100);

subplot(2,2,1)
hold all;
plot(nv.T, nv.out('Ca_k'), 'k-', 'LineWidth', 1);
ylabel('Ca_k [\muM]');
% plot(nv.T, nv.out('J_VOCC_i'));
% ylabel('VOCC flux');
xlim([XLIM1 XLIM2])
xlabel('Time [s]'); 
% % 
% subplot(2,2,3)
% hold all;
% plot(nv.T, nv.out('K_p')/1e3, 'k-', 'LineWidth', 1);
% ylabel('K_p [mM]');
% % plot(nv.T, nv.out('J_KIR_i'));
% % ylabel('KIR flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% ylim([0 15])
% 
% subplot(2,2,2)
% hold all;
% plot(nv.T, nv.out('w_k'), 'k-', 'LineWidth', 1);
% ylabel('w_k [-]');
% % plot(nv.T, nv.out('J_BK_k'));
% % ylabel('BK flux');
% xlim([XLIM1 XLIM2])
% xlabel('Time [s]'); 
% 
% subplot(2,2,4)
% hold all;
% plot(nv.T, nv.out('R')*1e6, 'k-', 'LineWidth', 1)
% xlabel('Time [s]'); 
% ylabel('Radius [\mum]');
% xlim([XLIM1 XLIM2])
% ylim([18 30])

% figure(FIG_NUM+2);
% 
% subplot(2,2,4)
% hold all;
% plot(nv.T, nv.out('w_i'));
% ylabel('w_i');
% xlim([XLIM1 XLIM2])
% 
% subplot(2,2,1)
% hold all;
% plot(nv.T, nv.out('J_KIR_i'));
% ylabel('J_{KIR}');
% xlim([XLIM1 XLIM2])
% 
% subplot(2,2,3)
% hold all;
% plot(nv.T, nv.out('eet_k'));
% ylabel('eet_k');
% xlim([XLIM1 XLIM2])

% subplot(2,2,2)
% hold all;
% plot(nv.T, nv.out('E_BK_k')*1e3);
% ylabel('E_{BK}');
% xlim([XLIM1 XLIM2])

% figure(1);
% hold all;
% plot(nv.T, nv.out('R')*1e6)
% xlabel('Time [s]'); ylabel('Radius [\mum]');
% xlim([XLIM1 XLIM2])

% 

% figure(FIG_NUM+10)
% hold all;
% plot(nv.T, nv.out('m_c'))
% xlabel('Time [s]'); ylabel('m_c')
% xlim([XLIM1 XLIM2])
% 
% figure(FIG_NUM+11)
% hold all;
% plot(nv.T, nv.out('CaM'))
% xlabel('Time [s]'); ylabel('CaM')
% xlim([XLIM1 XLIM2])
% 
% 
% figure(FIG_NUM+1)
% hold all;
% plot(nv.T, nv.out('R')*1e6)
% xlabel('Time [s]'); ylabel('Radius [\mum]')
% xlim([XLIM1 XLIM2])
% 
% % 
% figure(FIG_NUM+2)
% 
% subplot(1,3,3)
% hold all;
% plot(nv.T, nv.out('R')*1e6)
% xlabel('Time [s]'); ylabel('Radius [\mum]')
% xlim([XLIM1 XLIM2])
% ylim([16 30])
% 
% subplot(1,3,2)
% hold all;
% plot(nv.T, nv.out('K_p')/1e3)
% xlabel('Time [s]'); ylabel('K_p (mM)')
% xlim([XLIM1 XLIM2])
% ylim([0 35])
% 
% subplot(1,3,1)
% hold all;
% plot(nv.T, nv.out('K_e')/1e3)
% xlabel('Time [s]'); ylabel('K_e (mM)')
% xlim([XLIM1 XLIM2])
%ylim([0 55])

% subplot(2,2,4)
% hold all;
% plot(nv.T, nv.out('Ca_n'))
% xlabel('Time [s]'); title('Ca_n')
% xlim([XLIM1 XLIM2])
% 
% subplot(2,3,6)
% hold all;
% plot(nv.T, nv.out('J_BK_k'))
% xlabel('Time [s]'); title('BK flux')
% xlim([XLIM1 XLIM2])
% 
% subplot(1,2,3)
% hold all;
% plot(nv.T, nv.out('K_s')/1e3);
% xlabel('Time [s]'); title('K+ in SC [mM]')
% xlim([XLIM1 XLIM2])
% 
