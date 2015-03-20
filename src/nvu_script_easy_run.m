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

%%
for foo = [0,1]
    nv.astrocyte.params.glu_switch = foo;
    
    for pie = [0,1]
        nv.astrocyte.params.trpv_switch = pie;
        %% Run a basic simulation
        nv.simulate()

        %% Plots
        % Lists of quantities that can be retrieved from the NVU model after
        % simulation are given in the documentation pages of the individual model
        % components.

        figure(1)
        %suptitle('The Input Signal')
        subplot(4,1,1)
        hold all;
        plot(nv.T, nv.out('ft'))
        title('Input signal from the neuron into the synaptic cleft')
        xlabel('Time [s]'); ylabel('Input Signal f(t) [-]')

        subplot(4,1,2)
        hold all;
        plot(nv.T, nv.out('v_k')*1e3)
        title('Membrane Potential of the astrocyte')
        xlabel('Time [s]');ylabel('v_k [mV]')

        subplot(4,1,3)
        hold all;
        [ax,p1,p2] = plotyy(nv.T, nv.out('J_BK_k'), nv.T, nv.out('J_KIR_i'));
        title('The contribution of the BK- and KIR-channel to K_p')
        xlabel('Time [s]');
        ylabel(ax(1), 'J_{BK_k} [\muM ms^{-1}]'); ylabel(ax(2), 'J_{KIR_i} [\muM s^{-1}]') 

        subplot(4,1,4)
        hold all;
        plot(nv.T, nv.out('K_p'))
        title('Potassium concentration in perivascular space')
        xlabel('Time [s]');ylabel('K_p [\muM]')
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
        %%
        figure(2) % Plot the Radius
        plot(nv.T, 1e6 * nv.out('R'))
        hold all;
        title('Radius (R)');   xlabel('time (s)'); ylabel('(\mum)'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
        
        figure(3) % Plot variables from the Astrocyte 
        %suptitle('Astrocyte')
        subplot(3,2,1)
        hold all;
        plot(nv.T, nv.out('N_K_k')./(nv.out('R_k')))

        title('Potassium in AC (K_k)');   xlabel('Time [s]'); ylabel('[K^+]   (mM)') ; grid on
        subplot(3,2,2)
        hold all;
        plot(nv.T, 0.001*nv.out('N_Na_k')./(nv.out('R_k')))
        title('Sodium in AC (Na_k)');   xlabel('Time [s]'); ylabel('[Na^+]   (mM)'); grid on
        subplot(3,2,3)
        hold all;
        % plot(nv.T, 0.001*nv.out('N_HCO3_k')./(nv.out('R_k')))
        % title('HCO_3 in AC (HCO3_k)');   xlabel('Time [s]'); ylabel('[HCO_3]   (mM)'); grid on
        % subplot(3,2,4)
        % hold all;
        plot(nv.T, nv.out('v_k'))
        title('Membrane Voltage Astrocyte');   xlabel('Time [s]'); ylabel('[v_k]   (mV)'); grid on
        subplot(3,2,4)
        hold all;
        plot(nv.T, 0.001*nv.out('N_Cl_k')./(nv.out('R_k')))
        title('Chlorine in AC (Cl_k)');   xlabel('Time [s]'); ylabel('[Cl]   (mM)'); grid on
        subplot(3,2,5)
        hold all;
        plot(nv.T, nv.out('w_k'))
        title('Open probability of the BK Channel in AC (N_K_k)');   xlabel('Time [s]'); ylabel('(-)'); grid on
        subplot(3,2,6)
        hold all;
        plot(nv.T, (nv.out('J_BK_k'))'./(nv.out('R_k')))
        title('K^+ flux through the BK channel in AC (J\_BK\_k/R\_k)');   xlabel('Time [s]'); ylabel('K^+ flux [\muM/s]'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')

        figure(4) % Plot variables from the Perivascular Space and Synaptic Cleft
        subplot(1,2,1)
        hold all;
        plot(nv.T, 0.001*nv.out('K_p'))
        title('Potassium in the Perivascular Space (K_p)');   xlabel('Time [s]'); ylabel('[K^+]  (mM)'); grid on
        subplot(1,2,2)
        hold all;
        plot(nv.T, 0.001*nv.out('K_s'))
        title('[K^+] in synaptic cleft: K_s'); xlabel('Time [s]'); ylabel('[K^+]_s [mM]'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')

        figure(5)
        %suptitle('Smooth Muscle Cell')
        subplot(2,2,1)
        hold all;
        plot(nv.T, nv.out('Ca_i'))
        title('[Ca^{2+}] in smooth muscle cell: Ca_i'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_i [\muM]'); grid on
        subplot(2,2,2)
        hold all;
        plot(nv.T, nv.out('v_i'))
        xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage smooth muscle cell: v_i'); grid on
        subplot(2,2,3)
        hold all;
        plot(nv.T, nv.out('J_VOCC_i'))
        title('Ca^{2+} flux through the VOCC channel: J\_VOCC\_i'); xlabel('Time [s]'); ylabel('Ca^{2+} flux [\muM/s]'); grid on
        subplot(2,2,4)
        hold all;
        plot(nv.T,  nv.out('J_KIR_i'))
        title('K^+ flux through the KIR channel: J\_KIR\_i'); xlabel('Time [s]'); ylabel('K^+ flux [\muM m/s]'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')

        figure(6)
        %suptitle('Endothlial Cell')
        subplot(2,2,1)
        hold all;
        plot(nv.T, nv.out('Ca_j'))
        title('[Ca^{2+}] in Endothlial cell: Ca_j'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_i [\muM]'); grid on
        subplot(2,2,2)
        hold all;
        plot(nv.T, nv.out('v_j'))
        xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage Endothlial cell: v_j'); grid on
        subplot(2,2,[3:4])
        hold all;
        plot(nv.T, nv.out('J_IP3_j'))
        title('IP3 release in EC: J\_IP_3\_j'); xlabel('Time [s]'); ylabel(' [\muM/s]'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')

        figure(7) % Plot the new Calcium Equation Variables
        %suptitle('Calcium Equation Variables')
        subplot(1,2,1)
        hold all;
        plot(nv.T, nv.out('s_k'))
        title('[Ca^{2+}] in astrocytic ER: s_k'); xlabel('Time [s]'); ylabel('[Ca^{2+}] [\muM]'); grid on
        subplot(1,2,2)
        hold all;
        plot(nv.T, nv.out('c_k'))
        title('[Ca^{2+}] in astrocytic Cytosole: c_k'); xlabel('Time [s]'); ylabel('[Ca^{2+}] [\muM]'); grid on
        legend('NO GLU TRPV4 deactivated','NO GLU TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
    end
end