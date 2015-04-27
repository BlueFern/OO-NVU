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
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-3, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC(), ...
    'odeopts', odeopts);
 
%%
 for clu =[2880]
        nv.astrocyte.params.vk_switch =1;
        nv.astrocyte.params.J_max = clu;
        
        for foo =[8e-3] % [0.152,0.3,0.5,0.7,1,1.5]
            nv.astrocyte.params.CAP_switch = 1;
            nv.astrocyte.params.v_5 = foo;
            for pie = [0,1]
                nv.astrocyte.params.trpv_switch = pie;
                
                %% Run a basic simulation
                %nv.astrocyte.u0(nv.astrocyte.index.w_k) = clu;
                nv.simulate()
                figure(10)
                %suptitle('The Input Signal')
                
                hold all;
                subplot(3,1,1)
                plot(nv.T, nv.out('J_BK_k')'./nv.out('R_k'))
                title('J_(BK)/R_k')
                xlabel('Time [s]'); ylabel('K^+ flux [\muM/s]'); grid on
                
                hold all;
                subplot(3,1,2)
                plot(nv.T, nv.out('K_p'))
                title('Potassium concentration in PVS')
                xlabel('Time [s]'); ylabel('K^+ concentration in PVS [\muM]'); grid on
                
                hold all;
                subplot(3,1,3)
                plot(nv.T,1e6 * nv.out('R'))
                title('Radius')
                xlabel('Time [s]'); ylabel('(\mum)'); grid on
                legend('w_k=0.1e-3','w_k=0.1e-2','w_k=0.1e-1','w_k=0.1','w_k=1')
                %% Plots
                % Lists of quantities that can be retrieved from the NVU model after
                % simulation are given in the documentation pages of the individual model
                % components.
                
                figure(1)
                %suptitle('The Input Signal')
                subplot(4,1,1)
                hold all;
                plot(nv.T, nv.out('Ca_p'))
                title('--')
                xlabel('Time [s]'); ylabel('[\muM')
                
                subplot(4,1,2)
                hold all;
                plot(nv.T, nv.out('v_k')*1e3)
                title('Membrane Potential of the astrocyte')
                xlabel('Time [s]');ylabel('v_k [mV]')
                
                subplot(4,1,3)
                hold all;
                [ax,p1,p2] = plotyy(nv.T, nv.out('J_BK_k'),nv.T, nv.out('J_KIR_i'));
                title('The contribution of the BK- and KIR-channel to K_p')
                xlabel('Time [s]');
                ylabel(ax(1), 'J_{BK_k} [\muM ms^{-1}]'); ylabel(ax(2), 'J_{KIR_i} [\muM s^{-1}]')
                
                subplot(4,1,4)
                hold all;
                plot(nv.T, nv.out('K_p'))
                title('Potassium concentration in perivascular space')
                xlabel('Time [s]');ylabel('K_p [\muM]')
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                %%
                figure(2) % Plot the Radius
                plot(nv.T, 1e6 * nv.out('R'))
                hold all;
                title('Radius (R)');   xlabel('time (s)'); ylabel('(\mum)'); grid on
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                
                figure(3) % Plot variables from the Astrocyte
                
                subplot(3,2,1)
                hold all;
                plot(nv.T, nv.out('E_BK_k'),nv.T, nv.out('v_k'),'g')
                
                title('E_BK_k and V-k');   xlabel('Time [s]'); ylabel('[E_BK_k] and vk  (V') ; grid on
                subplot(3,2,2)
                hold all;
                plot(nv.T,nv.out('v_k')- nv.out('E_BK_k'))
                title('v_k-E_BK_k');   xlabel('Time [s]'); ylabel('[v_k-E_BK_k]   (V)'); grid on
                subplot(3,2,3)
                hold all;
                plot(nv.T, 0.001*nv.out('N_HCO3_k')./(nv.out('R_k')))
                title('HCO_3 in AC (HCO3_k)');   xlabel('Time [s]'); ylabel('[HCO_3]   (mM)'); grid on
                subplot(3,2,4)
                hold all;
                plot(nv.T, nv.out('K_k'))
                title('K+ in AC');   xlabel('Time [s]'); ylabel('[K+]   (\muM)'); grid on
                subplot(3,2,4)
                hold all;
                plot(nv.T, nv.out('K_p'))
                title('K+ in PVS');   xlabel('Time [s]'); ylabel('[K]   (\muM)'); grid on
                subplot(3,2,5)
                hold all;
                plot(nv.T, nv.out('J_KIR_i'))
                title('J_KIR');   xlabel('Time [s]'); ylabel('(-)'); grid on
                subplot(3,2,6)
                hold all;
                plot(nv.T, (nv.out('J_BK_k')))
                title('K^+ flux through the BK channel in AC J\_BK\_k');   xlabel('Time [s]'); ylabel('K^+ flux [\muM/s]'); grid on
                
                
                figure(4) % Plot variables from the Perivascular Space and Synaptic Cleft
                subplot(1,2,1)
                hold all;
                plot(nv.T, 0.001*nv.out('K_p'))
                title('Potassium in the Perivascular Space (K_p)');   xlabel('Time [s]'); ylabel('[K^+]  (mM)'); grid on
                subplot(1,2,2)
                hold all;
                plot(nv.T, 0.001*nv.out('K_s'))
                title('[K^+] in synaptic cleft: K_s'); xlabel('Time [s]'); ylabel('[K^+]_s [mM]'); grid on
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                
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
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                
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
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                
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
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                
                figure(8) % Plot variables from the Astrocyte
                
                positionVector1=[ 0.05, 0.8-0.03, 0.4 ,0.2];
                subplot('Position',positionVector1)
                hold all;
                plot(nv.T, nv.out('K_s'))
                title('Potassium in SC(K_s)');
                xlabel('Time [s]'); ylabel('[K+]   (\muM)') ; grid on
                
                positionVector2=[ 0.5, 0.8-0.03, 0.4, 0.2];
                subplot('Position',positionVector2)
                hold all;
                plot(nv.T, 1000*nv.out('v_k'))
                title('Membrane Voltage Astrocyte');
                xlabel('Time [s]'); ylabel('[v_k]   (mV)'); grid on
                
                positionVector3=[ 0.05, 0.55-0.03,0.4,0.2];
                subplot('Position',positionVector3)
                hold all;
                plot(nv.T, nv.out('i_k'))
                title('IP3');
                xlabel('Time [s]'); ylabel('[R_k]   (m)'); grid on
                
                positionVector4=[ 0.5, 0.55-0.03,0.4,0.2];
                subplot('Position',positionVector4)
                hold all;
                plot(nv.T, 0.001*nv.out('K_p'))
                title('Potassium in the Perivascular Space (K_p)');
                xlabel('Time [s]'); ylabel('[K+]   (mM)'); grid on
                
                positionVector5=[ 0.05, 0.30-0.03,0.4,0.2];
                subplot('Position',positionVector5)
                hold all;
                plot(nv.T, nv.out('c_k'))
                title('[Ca^{2+}] in astrocytic Cytosole: c_k');
                xlabel('Time [s]'); ylabel('[Ca^{2+}] [\muM]'); grid on
                
                positionVector6=[ 0.5, 0.30-0.03,0.4,0.2];
                subplot('Position',positionVector6)
                hold all;
                plot(nv.T, nv.out('Ca_i'))
                title('[Ca^{2+}] in smooth muscle cell: Ca_i');
                xlabel('Time [s]'); ylabel('[Ca^{2+}]_i [\muM]'); grid on
                
                positionVector7=[ 0.05, 0.05-0.03,0.4,0.2];
                subplot('Position',positionVector7)
                hold all;
                plot(nv.T, nv.out('eet_k'))
                title('EET');
                xlabel('Time [s]'); ylabel('[EET]   (microM)'); grid on
                
                positionVector8=[ 0.5, 0.05-0.03,0.4,0.2];
                subplot('Position',positionVector8)
                plot(nv.T, 1e6 * nv.out('R'))
                hold all;
                title('Radius (R)');
                xlabel('time (s)'); ylabel('(\mum) Radius'); grid on
                legend('deactivated TRPV4','activated TRPV4')
                
                figure(11)
                
                subplot(4,2,1)
                hold all;
                plot(nv.T, nv.out('J_TRPV_k'))
                title('TRPV flux')
                xlabel('Time [s]'); ylabel('[ Ca2+ \muM s-1]')
                
                subplot(4,2,2)
                hold all;
                plot(nv.T, nv.out('J_IP3'))
                title(' J_IP3 ')
                xlabel('Time [s]'); ylabel('[Ca2+ \muM s-1]')
                
                subplot(4,2,3)
                hold all;
                plot(nv.T, nv.out('J_pump'));
                title('J_pump')
                xlabel('Time [s]');
                ylabel('[Ca2+ \muM s-1]')
                
                subplot(4,2,4)
                hold all;
                plot(nv.T, nv.out('J_ER_leak'));
                title('J_ER_leak')
                xlabel('Time [s]');
                ylabel('[Ca2+ \muM s-1]')
                
                
                subplot(4,2,5)
                hold all;
                plot(nv.T,(nv.out('J_K_k')'./nv.out('R_k')))
                title('J_K_k/R_k')
                xlabel('Time [s]');ylabel('[K+ \muM s-1]')
                
                subplot(4,2,6)
                hold all;
                plot(nv.T,(nv.out('J_BK_k')'./nv.out('R_k')))
                title('J_BK_k/R_k')
                xlabel('Time [s]');ylabel('[K+ \muM s-1]')
                
                subplot(4,2,7)
                hold all;
                plot(nv.T,(nv.out('J_Na_k')'./nv.out('R_k')))
                title('J_Na_k/R_k')
                xlabel('Time [s]');ylabel('[Na+ \muM s-1]')
                
                subplot(4,2,8)
                hold all;
                plot(nv.T,(nv.out('J_NBC_k')'./nv.out('R_k')))
                title('J_NBC_k/R_k')
                xlabel('Time [s]');ylabel('[ \muM s-1]')
                legend('TRPV4 deactivated','TRPV4 activated')
                figure(9)
                %suptitle('The Input Signal')
                subplot(3,2,1)
                hold all;
                plot(nv.T,nv.out('E_TRPV_k'))
                title('E_TRPV_k')
                xlabel('Time [s]'); ylabel('[V]')
                
                subplot(3,2,2)
                hold all;
                plot(nv.T,nv.out('E_Na_k'))
                title('E_Na_k')
                xlabel('Time [s]'); ylabel('[V]')
                
                subplot(3,2,3)
                hold all;
                plot(nv.T, nv.out('E_BK_k'))
                title('E_BK_k')
                xlabel('Time [s]'); ylabel('[V]')
                
                subplot(3,2,4)
                hold all;
                plot(nv.T, nv.out('E_NBC_k'))
                title('E_NBC_k')
                xlabel('Time [s]'); ylabel('[V]')
                
                
                subplot(3,2,5)
                hold all;
                plot(nv.T, nv.out('E_Cl_k'))
                title('E_Cl_k')
                xlabel('Time [s]'); ylabel('[V]')
                subplot(3,2,6)
                plot(nv.T, nv.out('E_K_k'))
                hold all;
                title('E_K_k')
                xlabel('Time [s]'); ylabel('[V]')
                legend('constant Ca_p TRPV4 deactivated','constant Ca_p TRPV4 activated','TRPV4 deactivated','TRPV4 activated')
                %
                
                
            end
        
    end
    
end
