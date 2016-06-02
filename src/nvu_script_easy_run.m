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
clear all 
close all

%% Options 
% First the options need to be set for the ODE solver (currently |ode15s|) 
odeopts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC(), ...
    'odeopts', odeopts);
 for BK = [1,0] ;  %BK = 1 --> original BK channel
    nv.astrocyte.params.switchBK = BK;
    nv.astrocyte.params.trpv_switch=0;
    if BK == 0;
        TRPV12 = 0;
             for TRPV = [1,1,0];
           
                    nv.astrocyte.params.C_astr_k = 40;%pF
                    nv.astrocyte.params.gamma_k =1970;%mV/uM
                    nv.astrocyte.params.gam_cae_k= 200; %uM
                    nv.astrocyte.params.gam_cai_k=0.2; %uM
                    nv.astrocyte.params.epshalf_k= 0.16;
                    nv.astrocyte.params.kappa_k= 0.04;
                    nv.astrocyte.params.v1_TRPV_k= 0.120; %V
                    nv.astrocyte.params.v2_TRPV_k= 0.013; %V
                    nv.astrocyte.params.t_TRPV_k= 0.9; %mV
                    nv.astrocyte.params.R_0_passive_k= 20e-6;
                           
                    nv.astrocyte.params.G_TRPV_k= 200;
           
               if TRPV12 == 0;
                   
                    nv.astrocyte.params.C_astr_k = 40;%pF
                    nv.astrocyte.params.gamma_k =834.3;%mV/uM
                    nv.astrocyte.params.gam_cae_k= 200; %uM
                    nv.astrocyte.params.gam_cai_k=0.01; %uM
                    nv.astrocyte.params.epshalf_k= 0.1;
                    nv.astrocyte.params.kappa_k= 0.1;
                    nv.astrocyte.params.v1_TRPV_k= 0.120; %V
                    nv.astrocyte.params.v2_TRPV_k= 0.013; %V
                    nv.astrocyte.params.t_TRPV_k= 0.9; %mV
                    nv.astrocyte.params.R_0_passive_k= 20e-6;
                               
                    nv.astrocyte.params.G_TRPV_k= 50;
               end
          
           
           nv.astrocyte.params.G_BK_k =225; 
           nv.astrocyte.params.reverseBK =-0.08135;
           nv.astrocyte.params.switchBK = BK;
           nv.astrocyte.params.trpv_switch=TRPV;
             nv.astrocyte.params.Ca_4 = 0.35;
            nv.astrocyte.params.v_7 = -13.57e-3;
          TRPV12 = TRPV12 +1;
          %% Run a basic simulation
                %nv.astrocyte.u0(nv.astrocyte.index.w_k) = clu;
                nv.simulate()
                figure(10)
                %suptitle('The Input Signal')
                
                hold all;
              subplot(3,1,1)
                plot(nv.T, 0.001*nv.out('N_K_k')./(nv.out('R_k')))
            title('Potassium in AC (K_k)');   xlabel('Time [s]'); ylabel('[K^+]   (mM)') ; grid on
             hold all;
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
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
               plot(nv.T, nv.out('K_p')'./nv.out('K_k'))
                title('Potassium concentration in perivascular space / potassium concentration in AC')
                xlabel('Time [s]');ylabel('K_p [\muM]')
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                %%
                figure(2) % Plot the Radius
                plot(nv.T, 1e6 * nv.out('R'))
                hold all;
                title('Radius (R)');   xlabel('time (s)'); ylabel('(\mum)'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(3) % Plot variables from the Astrocyte
                
                subplot(3,2,1)
                hold all;
                plotyy(nv.T, nv.out('E_BK_k'),nv.T, nv.out('v_k'))
                
                title('E_{BK_k} and V-k');   xlabel('Time [s]') ; grid on
                subplot(3,2,2)
                hold all;
                plot(nv.T,nv.out('v_k')- nv.out('E_BK_k'))
                title('v_k-E_{BK_k}');   xlabel('Time [s]'); ylabel('[v_k-E_{BK_k}]   (V)'); grid on
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
                plot(nv.T, 0.001*nv.out('K_k'))
                title('[K^+] in astrocyte [K_k'); xlabel('Time [s]'); ylabel('[K^+]_k [mM]'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(5)
                %suptitle('Smooth Muscle Cell')
                subplot(2,2,1)
                hold all;
                plot(nv.T, nv.out('Ca_p'))
                title('[Ca^{2+}] PVS: Ca_p'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_p [\muM]'); grid on
                subplot(2,2,2)
                hold all;
                plot(nv.T, nv.out('v_i'))
                xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage smooth muscle cell: v_i'); grid on
                subplot(2,2,3)
                hold all;
                plot(nv.T, nv.out('J_VOCC_k'))
                title('Ca^{2+} flux through the VOCC channel: J\_VOCC\_i'); xlabel('Time [s]'); ylabel('Ca^{2+} flux [\muM/s]'); grid on
                subplot(2,2,4)
                hold all;
                plot(nv.T,  nv.out('J_TRPV_k'))
                title('TRPV4'); xlabel('Time [s]'); ylabel('Ca^2+ flux [\muM m/s]'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(8) % Plot variables from the Astrocyte
                 
                positionVector1=[ 0.05, 0.8-0.03, 0.35 ,0.17];
                subplot('Position',positionVector1)
                hold all;
                plot(nv.T, nv.out('K_s'),'LineWidth',1)
                title('(A) [K^+] in SC (K_s)', 'FontSize', 12);
                ylabel('[K+]   (\muM)','FontSize', 12) ; grid on
                set(gca,'fontsize',12,'xticklabel',[])
                xlim([50 450])
                
                positionVector2=[ 0.5, 0.8-0.03, 0.35, 0.17];
                subplot('Position',positionVector2)
                hold all;
                plot(nv.T, 1000*nv.out('v_k'),'LineWidth',1)
                title(' (B) Membrane voltage Astrocyte','FontSize', 12);
                ylabel('[v_k]   (mV)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
               positionVector3=[ 0.05, 0.55-0.02,0.35,0.17];
                subplot('Position',positionVector3)
                hold all;
                plot(nv.T, nv.out('w_k'),'LineWidth',1)
                title('(C) w_k Open probability BK channel','FontSize', 12);
                 ylabel('w_k (-)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector4=[ 0.5, 0.55-0.02,0.35,0.17];
                subplot('Position',positionVector4)
                hold all;
                plot(nv.T, 0.001*nv.out('K_p'),'LineWidth',1)
                title('(D) [K^+] in the Perivascular Space (K_p)','FontSize', 12);
                 ylabel('[K+]   (mM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector5=[ 0.05, 0.30-0.02,0.35,0.17];
                subplot('Position',positionVector5)
                hold all;
                plot(nv.T, nv.out('c_k'),'LineWidth',1)
                title('(E) [Ca^{2+}] in astrocytic Cytosole (c_k)','FontSize', 12);
                ylabel('[Ca^{2+}] (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                 
                positionVector6=[ 0.5, 0.30-0.02,0.35,0.17];
                subplot('Position',positionVector6)
                hold all;
                plot(nv.T, nv.out('Ca_i'),'LineWidth',1)
                title('(F) [Ca^{2+}] in Smooth muscle cell (Ca_i)','FontSize', 12);
                ylabel('[Ca^{2+}]_i (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector7=[ 0.05, 0.05-0.00,0.35,0.17];
                subplot('Position',positionVector7)
                hold all;
                plot(nv.T, nv.out('eet_k'),'LineWidth',1)
                title(' (G) [EET] in Astrocyte','FontSize', 12);
                xlabel('Time [s]','FontSize', 12); ylabel('[EET] (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12)
                 xlim([50 450])
                 
                positionVector8=[ 0.5, 0.05-0.00,0.35,0.17];
                subplot('Position',positionVector8)
                plot(nv.T, 1e6 * nv.out('R'),'LineWidth',1)
                hold all;
                title(' Radius (R)','FontSize', 12);
                xlabel('(H) Time (s)','FontSize', 12); ylabel('Radius (\mum)','FontSize', 12); grid on
               
                set(gca,'fontsize',12)
                 xlim([50 450])
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
           end  
    else
    
                %% Run a basic simulation
                %nv.astrocyte.u0(nv.astrocyte.index.w_k) = clu;
                nv.simulate()
                figure(10)
                %suptitle('The Input Signal')
                
                hold all;
              subplot(3,1,1)
                plot(nv.T, 0.001*nv.out('N_K_k')./(nv.out('R_k')))
            title('Potassium in AC (K_k)');   xlabel('Time [s]'); ylabel('[K^+]   (mM)') ; grid on
             hold all;
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
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
               plot(nv.T, nv.out('K_p')'./nv.out('K_k'))
                title('Potassium concentration in perivascular space / potassium concentration in AC')
                xlabel('Time [s]');ylabel('K_p [\muM]')
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                %%
                figure(2) % Plot the Radius
                plot(nv.T, 1e6 * nv.out('R'))
                hold all;
                title('Radius (R)');   xlabel('time (s)'); ylabel('(\mum)'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(3) % Plot variables from the Astrocyte
                
                subplot(3,2,1)
                hold all;
                plotyy(nv.T, nv.out('E_BK_k'),nv.T, nv.out('v_k'))
                
                title('E_{BK_k} and V-k');   xlabel('Time [s]') ; grid on
                subplot(3,2,2)
                hold all;
                plot(nv.T,nv.out('v_k')- nv.out('E_BK_k'))
                title('v_k-E_{BK_k}');   xlabel('Time [s]'); ylabel('[v_k-E_{BK_k}]   (V)'); grid on
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
                plot(nv.T, 0.001*nv.out('K_k'))
                title('[K^+] in astrocyte [K_k'); xlabel('Time [s]'); ylabel('[K^+]_k [mM]'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(5)
                %suptitle('Smooth Muscle Cell')
                subplot(2,2,1)
                hold all;
                plot(nv.T, nv.out('Ca_p'))
                title('[Ca^{2+}] PVS: Ca_p'); xlabel('Time [s]'); ylabel('[Ca^{2+}]_p [\muM]'); grid on
                subplot(2,2,2)
                hold all;
                plot(nv.T, nv.out('v_i'))
                xlabel('Time [s]'); ylabel('v_i [mV]'); title('Membrane voltage smooth muscle cell: v_i'); grid on
                subplot(2,2,3)
                hold all;
                plot(nv.T, nv.out('J_VOCC_k'))
                title('Ca^{2+} flux through the VOCC channel: J\_VOCC\_i'); xlabel('Time [s]'); ylabel('Ca^{2+} flux [\muM/s]'); grid on
                subplot(2,2,4)
                hold all;
                plot(nv.T,  nv.out('J_TRPV_k'))
                title('TRPV4'); xlabel('Time [s]'); ylabel('Ca^2+ flux [\muM m/s]'); grid on
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
                
                figure(8) % Plot variables from the Astrocyte
                 
                positionVector1=[ 0.05, 0.8-0.03, 0.35 ,0.17];
                subplot('Position',positionVector1)
                hold all;
                plot(nv.T, nv.out('K_s'),'LineWidth',1)
                title('(A) [K^+] in SC (K_s)', 'FontSize', 12);
                ylabel('[K+]   (\muM)','FontSize', 12) ; grid on
                set(gca,'fontsize',12,'xticklabel',[])
                xlim([50 450])
                
                positionVector2=[ 0.5, 0.8-0.03, 0.35, 0.17];
                subplot('Position',positionVector2)
                hold all;
                plot(nv.T, 1000*nv.out('v_k'),'LineWidth',1)
                title(' (B) Membrane voltage Astrocyte','FontSize', 12);
                ylabel('[v_k]   (mV)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
               positionVector3=[ 0.05, 0.55-0.02,0.35,0.17];
                subplot('Position',positionVector3)
                hold all;
                plot(nv.T, nv.out('w_k'),'LineWidth',1)
                title('(C) w_k Open probability BK channel','FontSize', 12);
                 ylabel('w_k (-)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector4=[ 0.5, 0.55-0.02,0.35,0.17];
                subplot('Position',positionVector4)
                hold all;
                plot(nv.T, 0.001*nv.out('K_p'),'LineWidth',1)
                title('(D) [K^+] in the Perivascular Space (K_p)','FontSize', 12);
                 ylabel('[K+]   (mM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector5=[ 0.05, 0.30-0.02,0.35,0.17];
                subplot('Position',positionVector5)
                hold all;
                plot(nv.T, nv.out('c_k'),'LineWidth',1)
                title('(E) [Ca^{2+}] in astrocytic Cytosole (c_k)','FontSize', 12);
                ylabel('[Ca^{2+}] (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                 
                positionVector6=[ 0.5, 0.30-0.02,0.35,0.17];
                subplot('Position',positionVector6)
                hold all;
                plot(nv.T, nv.out('Ca_i'),'LineWidth',1)
                title('(F) [Ca^{2+}] in Smooth muscle cell (Ca_i)','FontSize', 12);
                ylabel('[Ca^{2+}]_i (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12,'xticklabel',[])
                 xlim([50 450])
                
                positionVector7=[ 0.05, 0.05-0.00,0.35,0.17];
                subplot('Position',positionVector7)
                hold all;
                plot(nv.T, nv.out('eet_k'),'LineWidth',1)
                title(' (G) [EET] in Astrocyte','FontSize', 12);
                xlabel('Time [s]','FontSize', 12); ylabel('[EET] (\muM)','FontSize', 12); grid on
                set(gca,'fontsize',12)
                 xlim([50 450])
                 
                positionVector8=[ 0.5, 0.05-0.00,0.35,0.17];
                subplot('Position',positionVector8)
                plot(nv.T, 1e6 * nv.out('R'),'LineWidth',1)
                hold all;
                title(' Radius (R)','FontSize', 12);
                xlabel('(H) Time (s)','FontSize', 12); ylabel('Radius (\mum)','FontSize', 12); grid on
               
                set(gca,'fontsize',12)
                 xlim([50 450])
                
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
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
                legend('NVU', 'TRPV4 set 1', 'TRPV4 set 2', 'constant E_{BK}')
    end       %
                
   
        
end 
