%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 2.1
% K+, NO, Astrocytic Ca2+, TRPV4, ECS, Neuron

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

timeStart = datetime('now');
fprintf('Start time is %s\n', char(timeStart));



odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);
FIG_NUM = 1;
XLIM1 = 50;
XLIM2 = 5000; % End of simulation

% For current type 1 or 2 use max current strength 0.022
% For current type 3 use max current strength 0.042
% For current type 4 use max current strength 0.035

CURRENT_STRENGTH    = 0.022;    % Max strength of current input in mA/cm2     %%%%%%%%%%%%%%%%% 0.006 or 0.022 sometimes 0042?
NEURONAL_START      = 100000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ;      % Start of neuronal stimulation
CURRENT_TYPE        = 1;        % Types of current input. 1: normal, 2: two stimulations (second stimulation is 8 sec after and 1 sec long), 3: obtained from experimental input data, 4: whisker pad (from experiment) + locus coeruleus (pain pathway)

% Used if CURRENT_STRENGTH = 1 or 2
NEURONAL_END        = 100010;      % End of neuronal stimulation 

% Used if CURRENT_STRENGTH = 3 or 4
ISI = 7;                        % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
stim = 2;                       % INDEX for length of initial stimulation [2,8,16]

% Used if CURRENT_STRENGTH = 4: scaling for the two stimulation components,
% alpha for whisker pad and beta for locus coeruleus/pain. Default for both
% is 1. 
% I_total = alpha * I_Wh + beta * I_LC
alpha = 1;
beta = 1;

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations
GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)
R_DECAY         = 0.15;
    
% Load initial NVU ksyn11.5 when replicating experiment use 2.4 10.0... 
nv = NVU(Neuron('k_syn', 2.4, 'CurrentType', CURRENT_TYPE, 'O2switch', O2SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH), ...
    Astrocyte('trpv_switch', TRPV_SWITCH, 'R_decay', R_DECAY), ...
    WallMechanics('wallMech', 3.0), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), ...
    ANLS('startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START), 'odeopts', odeopts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 11.5 goes to 1 (k_syn)

% Adjust time vector
nv.neuron.params.dt = 0.01; dt = nv.neuron.params.dt;
nv.T = 0:dt:XLIM2;
numTimeSteps = length(nv.T);

%% Load whisker stimulation input data from file, save as I_Wh and put into Neuron.m
if nv.neuron.params.CurrentType ~= 1
    
    load Zheng2010_data.mat
    actual_ISI = info.isi_duration(ISI);
    actual_stim = info.condition_stim_duration(stim);
    sum_neural_wh = zeros(size(neural_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_neural_wh = sum_neural_wh+neural_data(:,ISI,stim,experiment,animal)'; % Sum all data 
        end
    end
    mean_neural_wh = sum_neural_wh./110;  % Average the neural data over all animals and experiments, animals*experiments=110
    neural_tim_vector_shifted = neural_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
    interp_neural_wh = interp1(neural_tim_vector_shifted, mean_neural_wh, nv.T); % Interpolate so there is data for all timesteps for NVU
    interp_neural_wh(isnan(interp_neural_wh))=0.02;   % Remove NaNs     
    
    nv.neuron.input_data = interp_neural_wh;   % Replace dummy in Neuron with input data, leave it at that for CURRENT_TYPE = 3
    I_Wh = interp_neural_wh;                   % Save as whisker pad current I_Wh
end

%% Construct additional pain pathway stimulation I_LC and input total current I_total to Neuron.m
if nv.neuron.params.CurrentType == 4
    
    % Construct first stimulation, set to work for any NEURONAL_START or any
    % initial stimulus duration (2,8,16 s)
    time_1 = linspace(0, actual_stim, 10000); 
    I_1 = 0.0006 * time_1.^2 + 0.1;             % Chosen for the shape! Can be modified if needed
    time_1 = time_1 + NEURONAL_START;           % Shift so starts at NEURONAL_START
    I_1 = interp1(time_1, I_1, nv.T);           % Make same size as nv.T
    I_1(isnan(I_1)) = 0;                        % Replace NaNs with zeros
    
    % Construct second stimulation with duration 2 sec
    time_2 = linspace(0, 1, 10000);
    I_2 = 0.1*ones(size(time_2));
    time_2 = time_2 + NEURONAL_START + actual_stim + actual_ISI;
    I_2 = interp1(time_2, I_2, nv.T);
    I_2(isnan(I_2)) = 0;
    
    % Add together
    I_LC = I_1 + I_2;
    
    % Total current (whisker pad plus LC)
    I_total = alpha*I_Wh + beta*I_LC;
    
    % Input to Neuron.m
    nv.neuron.input_data = I_total;
end

%% Run the simulation
nv.simulate() 




figure(2);
hold all
set(gcf,'Name', 'ANLS State Variables')
anls_vars = fieldnames(nv.anls.index);
i_anls = size(anls_vars, 1);
for i = 1:1:i_anls
    subplot(5,6,i)
    hold all
    plot(nv.T, nv.out(char(anls_vars(i))));
    xlabel('Time [s]'); ylabel(strcat(anls_vars(i), ' [mM]'));
    xlim([XLIM1 XLIM2])
    
end
% 
% 
% 
figure(33667);
set(gcf,'Name', 'ATPase pumps')
subplot(3,1,1)
hold all
plot(nv.T, nv.out('Vn_pump')); 
xlabel('Time [s]'); ylabel({'Neuron ATPase'; 'pump (mM/s)'});
xlim([XLIM1 XLIM2])
subplot(3,1,2)
hold all
plot(nv.T, nv.out('Vk_pump')); 
xlabel('Time [s]'); ylabel({'Astrocyte ATPase';' pump (mM/s)'});
xlim([XLIM1 XLIM2])
subplot(3,1,3)
hold all
plot(nv.T, nv.out('CBF')); 
xlabel('Time [s]'); ylabel('CBF');
xlim([XLIM1 XLIM2])
% figure(55);
% set(gcf,'Name', 'stuff')
% subplot(2,3,1)
% hold all
% plot(nv.T, nv.out('K_p')); 
% xlabel('Time [s]'); ylabel('K_p');
% xlim([95 150])
% subplot(2,3,2)
% hold all
% plot(nv.T, nv.out('J_BK_k')); 
% xlabel('Time [s]'); ylabel('J_BK_k');
% xlim([95 150])
% subplot(2,3,3)
% hold all
% plot(nv.T, nv.out('J_KIR_i')); 
% xlabel('Time [s]'); ylabel('J_KIR_i');
% xlim([95 150])
% subplot(2,3,4)
% hold all
% plot(nv.T, nv.out('v_k')); 
% xlabel('Time [s]'); ylabel('v_k');
% xlim([95 150])
% subplot(2,3,5)
% hold all
% plot(nv.T, nv.out('E_BK_k')); 
% xlabel('Time [s]'); ylabel('E_BK_k');
% xlim([95 150])
% subplot(2,3,6)
% hold all
% plot(nv.T, nv.out('v_i')); 
% xlabel('Time [s]'); ylabel('v_i');
% xlim([95 150])
% 

% figure(556);
% set(gcf,'Name', 'Astrocytic K+')
% subplot(2,3,1)
% hold all
% plot(nv.T, nv.out('K_k')); 
% xlabel('Time [s]'); ylabel('K_k');
% xlim([XLIM1 XLIM2])
% subplot(2,3,2)
% hold all
% plot(nv.T, nv.out('K_s')); 
% xlabel('Time [s]'); ylabel('K_s');
% xlim([XLIM1 XLIM2])
% subplot(2,3,3)
% hold all
% plot(nv.T, nv.out('J_BK_k')); 
% xlabel('Time [s]'); ylabel('J_BK_k');
% xlim([XLIM1 XLIM2])
% subplot(2,3,4)
% hold all
% plot(nv.T, nv.out('J_K_k')); 
% xlabel('Time [s]'); ylabel('J_K_k');
% xlim([XLIM1 XLIM2])
% subplot(2,3,5)
% hold all
% plot(nv.T, nv.out('E_K_k')); 
% xlabel('Time [s]'); ylabel('E_K_k');
% xlim([XLIM1 XLIM2])
% subplot(2,3,6)
% hold all
% plot(nv.T, nv.out('v_k')); 
% xlabel('Time [s]'); ylabel('v_k');
% xlim([XLIM1 XLIM2])


XLIM1 = 20;
XLIM2 = 500;
figure(337);
set(gcf,'Name', 'input stuff cloutier')
subplot(3,1,1)
hold all
plot(nv.T, nv.out('Vc_GLC')); 
xlabel('Time [s]'); ylabel('Vc_GLC (mM/s)');
xlim([XLIM1 XLIM2])
subplot(3,1,2)
hold all
plot(nv.T, nv.out('Vc_O2')); 
xlabel('Time [s]'); ylabel('Vc_O2 (mM/s)');
xlim([XLIM1 XLIM2])
subplot(3,1,3)
hold all
plot(nv.T, nv.out('Vc_LAC')); 
xlabel('Time [s]'); ylabel('Vc_LAC (mM/s)');
xlim([XLIM1 XLIM2])


% figure(20000);
% set(gcf,'Name', 'Dropping arteriole GLC')
% subplot(2,2,1)
% hold all
% plot(nv.T, nv.out('Vn_pump'));
% xlabel('Time [s]'); ylabel('ATPase pump (mM/s)');
% subplot(2,2,2)
% hold all
% plot(nv.T, nv.out('ATPn'));
% xlabel('Time [s]'); ylabel('ATPn (mM)');
% subplot(2,2,3)
% hold all
% plot(nv.T, nv.out('GLCc'));
% xlabel('Time [s]'); ylabel('GLCc (mM)');
% subplot(2,2,4)
% hold all
% plot(nv.T, nv.out('GLCn'));
% xlabel('Time [s]'); ylabel('GLCn (mM)');



% 
% 
% figure(2342343);
% set(gcf,'Name', 'ATPn fluxes')
% subplot(2,2,1)
% hold all
% plot(nv.T, nv.out('Vn_ldh'));
% xlabel('Time [s]'); ylabel('Vn_ldh');
% subplot(2,2,2)
% hold all
% plot(nv.T, nv.out('Vn_ldh'));
% xlabel('Time [s]'); ylabel('Vn_ldh');



% anls_vars = fieldnames(nv.anls.index);
% i_anls = size(anls_vars, 1);
% for i = 1:1:i_anls
%     h = figure(100+i);
%     hold all
%     plot(nv.T, nv.out(char(anls_vars(i))));
%     xlabel('Time [s]'); ylabel(strcat(anls_vars(i), ' [mM]'));
%     xlim([XLIM1 XLIM2])
%     
%     set(h,'Units','Inches');
%     pos = get(h,'Position');
%     set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% 
%     print(h, '-dpdf', strcat('/home/s/SD_PhD/thesis/images/anls/test_results/ours/',char(anls_vars(i)),'.pdf'), '-r0')
% end
% hold off


% figure(123414);
% subplot(4,2,1);
%     hold all;
%     plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
%     ylabel('v_{sa} [mV]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,2);
%     hold all;
%     plot(nv.T, nv.out('GLCn'), 'LineWidth', 1);
%     ylabel('GLC_n [mM]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,3);
%     hold all;
%     plot(nv.T, nv.out('R'), 'LineWidth', 1);
%     ylabel('Radius [\mum]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,4);
%     hold all;
%     plot(nv.T, nv.out('K_s')/1e3, 'LineWidth', 1);
%     ylabel('K_s [mM]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,5);
%     hold all;
%     plot(nv.T, nv.out('ATPn'), 'LineWidth', 1);
%     ylabel('ATP_n [mM]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,6);
%     hold all;
%     plot(nv.T, nv.out('Vk_pump'), 'LineWidth', 1);
%     ylabel('ATpase Pump (astrocyte) [mMs^{-1}]');
%     xlim([0 1000])
% subplot(4,2,7);
%     hold all;
%     plot(nv.T, nv.out('Vn_pump'), 'LineWidth', 1);
%     ylabel('ATpase Pump (neuron) [mMs^{-1}]');
%     xlim([XLIM1 XLIM2])
% subplot(4,2,8);
%     hold all;
%     plot(nv.T, nv.out('ATPg'), 'LineWidth', 1);
%     ylabel('ATP_g [mM]');
%     xlim([0 1000])



% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Find a point in time 20 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and the point in time when stimulation begins (preNeuronalStimTime2) and get the values for 
% CBF, CBV, HBR, CMRO2 at that time, then find the midpoint between the min and max and normalise with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
np = nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = (nv.out('CBV'))'; HBR = (nv.out('HbR'))'; CMRO2 = nv.out('CMRO2');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
HBR_0 = 0.5*( max(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + min(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
CMRO2_0 = 0.5*( max(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

%Normalised variables
CBF_N = CBF./CBF_0;
CBV_N = CBV./CBV_0;
HBR_N = HBR./HBR_0;
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N .* HBR_N ./ CMRO2_N;                                          % Total hemoglobin (normalised)
HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1;                                      % Oxyhemoglobin (normalised)
BOLD_N = 100 * np.V_0 * ( np.a_1 * (1 - HBR_N) - np.a_2 * (1 - CBV_N) );    % BOLD (percentage increase from 0)

%% Plot experimental and model CBF from data file
if nv.neuron.params.CurrentType ~= 1
    sum_cbf = zeros(size(cbf_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_cbf = sum_cbf+cbf_data(:,ISI,stim,experiment,animal)';
        end
    end
    mean_cbf = (sum_cbf./110) - 1;
    cbf_tim_vector_shifted = cbf_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
    % Plot
    figure;
    plot(cbf_tim_vector_shifted, mean_cbf, ':k', nv.T, (nv.out('CBF')-CBF_0)./CBF_0, '-k', 'LineWidth', 1);
    ylabel('\Delta CBF')
    xlabel('Time [s]')
    xlim([90 150])
    %ylim([-0.05 0.3])
    legend('experiment','model')
    title(['CBF with initial duration ' num2str(actual_stim) ', ISI ' num2str(actual_ISI)] );
    p1=patch([100 100+actual_stim 100+actual_stim 100],[-0.05 -0.05 0.3 0.3],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[-0.05 -0.05 0.3 0.3],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
end











%% Plot BOLD
% figure;
%     plot(nv.T, BOLD_N, 'k', 'LineWidth', 1)
%     ylabel('BOLD (%)')
%     xlabel('Time [s]')
%     xlim([90 150])
%     p1=patch([100 100+actual_stim 100+actual_stim 100],[-0.6 -0.6 1.2 1.2],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
%     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[-0.6 -0.6 1.2 1.2],'k');
%     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
%     
%% Plot Hemoglobin
% figure;
%     plot(nv.T, HBO_N, ':k', nv.T, HBR_N, '--k', nv.T, HBT_N, '-k', 'LineWidth', 1)
%     xlabel ('Time [s]')
%     xlim([XLIM1 XLIM2])
%     legend('HbO','HbR','HbT')
%     p1=patch([100 100+actual_stim 100+actual_stim 100],[0.9 0.9 1.2 1.2],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
%     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0.9 0.9 1.2 1.2],'k');
%     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');

%% Plot current input
%     figure;
%     plot(nv.T, I_total, 'k', 'LineWidth', 1);
%     ylabel('Neural input')
%     xlabel('Time [s]')
%     xlim([90 150])
    %ylim([0 1])
%     p1=patch([100 100+actual_stim 100+actual_stim 100],[0 0 1 1],'k');
%     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
%     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0 0 1 1],'k');
%     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');


%% Plot CBF
% figure(FIG_NUM);
% hold all;
% plot(nv.T, nv.out('CBF')./CBF_0, 'k', 'LineWidth', 1);
% ylabel('CBF [-]');
% xlim([XLIM1 XLIM2])
% ylim([0.97 1.4])
% title('Normalised CBF');
% xlabel('time (s)');
% p1=patch([100 116 116 100],[0.97 0.97 1.4 1.4],'k');
% set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% p2=patch([124 125 125 124],[0.97 0.97 1.4 1.4],'k');
% set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% 
% %% Plot current
% figure(200);
% hold all;
% plot(nv.T, nv.out('current'));
% ylabel('current');

%% Plot hemoglobin
% 
% figure(FIG_NUM);
% plot(nv.T, HBO_N-1, 'r', nv.T, HBR_N-1, 'b', nv.T, HBT_N-1, 'g', 'LineWidth', 2);
% p=patch([100 104 104 100],[-0.15 -0.15 0.25 0.25],'k');
% set(p,'FaceAlpha',0.1,'EdgeColor', 'none');
% xlim([95 130]);
% ylabel('\Delta Concentration');
% xlabel('Time [s]');
% legend('HbO','HbR','HbT');
% 
% figure(FIG_NUM+1);
% plot(nv.T, BOLD_N, 'k', 'LineWidth', 2);
% p=patch([100 104 104 100],[-0.6 -0.6 1.7 1.7],'k');
% set(p,'FaceAlpha',0.1,'EdgeColor', 'none');
% xlim([95 130]);
% ylim([-0.6 1.7]);
% ylabel('\Delta BOLD');
% xlabel('Time [s]');

% figure(FIG_NUM+2);
% plot( nv.T, nv.out('J_KDR_sa'),  nv.T, nv.out('J_KA_sa'), nv.T, nv.out('J_Kleak_sa'), nv.T, nv.out('J_Kpump_sa'), nv.T, nv.out('J_KDR_d'), nv.T, nv.out('J_KA_d'), nv.T, nv.out('J_Kleak_d'), nv.T, nv.out('J_Kpump_d'), nv.T, nv.out('J_NMDA_K_d'));
% xlim([XLIM1 XLIM2]);
% legend('KDR_s','KA_s','Kleak_s','Kpump_s','KDR_d','KA_d','Kleak_d','Kpump_d','NMDA_d')

% figure(FIG_NUM+1);
% subplot(2,2,1);
%     hold all;
%     plot(nv.T, nv.out('nNOS_act_n'), 'LineWidth', 1);
%     ylabel('nNOS');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]')
% %     p1=patch([100 100+actual_stim 100+actual_stim 100],[0.3 0.3 0.9 0.9],'k');
% %     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% %     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0.3 0.3 0.9 0.9],'k');
% %     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% subplot(2,2,2);
%     hold all;
%     plot(nv.T, nv.out('NO_n'), 'LineWidth', 1);
%     ylabel('NO_n');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]')
% %     p1=patch([100 100+actual_stim 100+actual_stim 100],[0.15 0.15 0.45 0.45],'k');
% %     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% %     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0.15 0.15 0.45 0.45],'k');
% %     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% subplot(2,2,3);
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
%     ylabel('Radius');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]')
% %     p1=patch([100 100+actual_stim 100+actual_stim 100],[22.5 22.5 25.5 25.5],'k');
% %     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% %     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[22.5 22.5 25.5 25.5],'k');
% %     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% subplot(2,2,4);
%     hold all;
%     plot(nv.T, CBF_N, 'LineWidth', 1);
%     ylabel('CBF');
%     xlim([XLIM1 XLIM2])
%     xlabel('Time [s]')
% %     p1=patch([100 100+actual_stim 100+actual_stim 100],[1 1 1.4 1.4],'k');
% %     set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% %     p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[1 1 1.4 1.4],'k');
% %     set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
    
%% Plot some variables
% % 

figure(12341234);
subplot(3,3,1);
    hold all;
    plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
    ylabel('v_{sa} [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,2);
    hold all;
    plot(nv.T, nv.out('v_k'), 'LineWidth', 1);
    ylabel('v_k [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,3);
    hold all;
    plot(nv.T, nv.out('O2')*1e3, 'LineWidth', 1);
    ylabel('O2 [\muM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,4);
    hold all;
    plot(nv.T, nv.out('Nag'), 'LineWidth', 1);
    ylabel('Nag');
    xlim([XLIM1 XLIM2])
subplot(3,3,5);
    hold all;
    plot(nv.T, nv.out('R'), 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(3,3,6);
    hold all;
    plot(nv.T, nv.out('K_s')/1e3, 'LineWidth', 1);
    ylabel('K_s [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, nv.out('K_e'), 'LineWidth', 1);
    ylabel('K_e [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, nv.out('K_p')/1e3, 'LineWidth', 1);
    ylabel('K_p [mM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    % BOLD using normalised variables
    plot(nv.T, nv.out('Na_sa'), 'LineWidth', 1);  
    ylabel('Na_sa');
    xlim([XLIM1 XLIM2])
% %  % Plot figures - whatever you want

%% Plot for paper
% % % % % figure(1);
% % % % % subplot(4,2,1);
% % % % %     hold all;
% % % % %     plot(nv.T, nv.out('K_s')/1e3, 'k', 'LineWidth', 1);
% % % % %     ylabel('K_s [mM]');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([-100 20])
% % % % %     title('A. SC potassium');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, -100, 16, 5],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, -100, 2, 5],'FaceColor',[0 0 0])
% % % % % subplot(4,2,2);
% % % % %     hold all;
% % % % %     plot(nv.T, nv.out('K_e'), 'k', 'LineWidth', 1);
% % % % %     ylabel('K_e [mM]');
% % % % %     xlim([XLIM1 XLIM2])
% % % % %     title('B. Extracellular potassium');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, 2, 16, 0.2],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, 2, 2, 0.2],'FaceColor',[0 0 0])
% % % % % subplot(4,2,3);
% % % % %     hold all;
% % % % %     plot(nv.T, nv.out('R')*1e6, 'k', 'LineWidth', 1);
% % % % %     ylabel('Radius [\mum]');
% % % % % %     ylim([22.4 25])
% % % % %     xlim([XLIM1 XLIM2])
% % % % %     title('C. Radius');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, 22.4, 16, 0.1],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, 22.4, 2, 0.1],'FaceColor',[0 0 0])
% % % % % subplot(4,2,4);
% % % % %     hold all;
% % % % %     plot(nv.T, nv.out('CMRO2'), 'k', 'LineWidth', 1);
% % % % %     ylabel('CMRO2');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([-0.1 0.4])
% % % % %     title('D. Normalised CMRO_2');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, -0.1, 16, 0.02],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, -0.1, 2, 0.02],'FaceColor',[0 0 0])
% % % % % subplot(4,2,5);
% % % % %     hold all;
% % % % %     plot(nv.T, nv.out('O2'), 'k', 'LineWidth', 1);
% % % % %     ylabel('O_2 [mM]');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([0.0265 0.03])
% % % % %     title('E. Oxygen concentration');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, 0.0265, 16, 0.00015],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, 0.0265, 2, 0.00015],'FaceColor',[0 0 0])
% % % % % subplot(4,2,6);
% % % % %     hold all;
% % % % %     plot(nv.T, CBV_N, 'k', 'LineWidth', 1);
% % % % %     ylabel('CBV [-]');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([0.97 1.1])
% % % % %     title('F. Normalised CBV');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, 0.97, 16, 0.005],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, 0.97, 2, 0.005],'FaceColor',[0 0 0])
% % % % % subplot(4,2,7);
% % % % %     hold all;
% % % % %     plot(nv.T, HBR_N, 'k', 'LineWidth', 1);
% % % % %     ylabel('HBR [-]');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([0.86 1.05])
% % % % %     title('G. Normalised deoxyhemoglobin');
% % % % %     xlabel('time (s)');
% % % % % %     rectangle('Position',[300, 0.86, 16, 0.008],'FaceColor',[0 0 0])
% % % % % %     rectangle('Position',[320, 0.86, 2, 0.008],'FaceColor',[0 0 0])
% % % % % subplot(4,2,8);
% % % % %     hold all;
% % % % %     plot(nv.T, BOLD_N, 'k', 'LineWidth', 1);
% % % % %     ylabel('\Delta BOLD (%)');
% % % % %     xlim([XLIM1 XLIM2])
% % % % % %     ylim([-0.9 1.5])
% % % % %     title('H. BOLD response');
% % % % %     xlabel('time (s)');
%     rectangle('Position',[300, -0.9, 16, 0.09],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -0.9, 2, 0.09],'FaceColor',[0 0 0])
 
% % %%
% figure(FIG_NUM+10);
% subplot(2,3,1);
%     hold all;
%     plot(nv.T, nv.out('Ca_i'), 'LineWidth', 1);
%     ylabel('Ca_i [\muM]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,2);
%     hold all;
%     plot(nv.T, nv.out('cGMP_i'), 'LineWidth', 1);
%     ylabel('cGMP [\muM]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,3);
%     hold all;
%     plot(nv.T, nv.out('AMp'), 'LineWidth', 1);
%     ylabel('AMp [-]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,4);
%     hold all;
%     plot(nv.T, nv.out('AM'), 'LineWidth', 1);
%     ylabel('AM [-]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,5);
%     hold all;
%     plot(nv.T, nv.out('F_r'), 'LineWidth', 1);
%     ylabel('F_r [-]');
%     xlim([XLIM1 XLIM2])
% subplot(2,3,6);
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
%     ylabel('R [\mum]');
%     xlim([XLIM1 XLIM2])
% 
% figure(FIG_NUM+11);
%     hold all;
%     plotyy(nv.T, -nv.out('Ca_i'), nv.T, nv.out('R')*1e6)
%     %ylabel('Ca_i [\muM]');
%     %xlim([XLIM1 XLIM2])
%     
    
timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
