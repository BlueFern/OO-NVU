clear;

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

% Time
T = linspace(0, 600, 20000);

% Important constants:
J_PLC_1 = 0.14;
J_PLC_2 = 0.17;

% SMC coupling
D_Ca_i = 2;
D_IP3_i = 0;        % Optional IP3 coupling
D_v_i = 0;          % Optional membrane potential coupling

% EC coupling - optional!!
D_Ca_j = 0;
D_IP3_j = D_Ca_j;   % Optional IP3 coupling
D_v_j = D_Ca_j;     % Optional membrane potential coupling

%K_p = 9300;     % Signal on
%K_p = 3400;     % Signal off

K_p = 8800;

nv = NVU_coupled_red( WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1_red('J_PLC_1', J_PLC_1, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_1', K_p), ...
    SMCEC_2_red('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_2', K_p), ...
    'odeopts', odeopts, 'T', T);

nv.simulate();

% Plot calcium of each SMC
%Calcium
figure(100);
plot(nv.T, nv.out('Ca_i_1'));
xlabel('t'); ylabel('Ca_i');
title('Cell 1 Ca_i');
ylim([0 1]);

figure(200);
plot(nv.T, nv.out('Ca_i_2'));
xlabel('t'); ylabel('Ca_i');
title('Cell 2 Ca_i');
ylim([0 1]);

% Plot another case on the same graph %
% figure(1);
% hold on
% plot(nv.T, nv.out('Ca_i_1'),'r');
% title('Cell 1 Ca_i');
% hold off
% 
% figure(2);
% hold on
% plot(nv.T, nv.out('Ca_i_2'),'r');
% title('Cell 2 Ca_i');
% hold off

% Trajectories: cytosolic calcium vs SR calcium

% For plotting trajectory only when signal on
% starttime = floor((nv.astrocyte_1.params.startpulse+80)/nv.T(length(nv.T))*length(nv.T));
% endtime = floor((nv.astrocyte_1.params.startpulse+nv.astrocyte_1.params.lengthpulse)/nv.T(length(nv.T))*length(nv.T));

starttime = floor(300/nv.T(length(nv.T))*length(nv.T));
endtime = length(nv.T);

Ca_i_1 = nv.out('Ca_i_1'); s_i_1 = nv.out('s_i_1');
figure(5);
plot(Ca_i_1(starttime:endtime,:), s_i_1(starttime:endtime,:));
title('Cell 1 trajectory');
xlabel('Calcium in cytosol'); ylabel('Calcium in store');

Ca_i_2 = nv.out('Ca_i_2'); s_i_2 = nv.out('s_i_2');
figure(6);
plot(Ca_i_2(starttime:endtime,:), s_i_2(starttime:endtime,:));
title('Cell 2 trajectory');
xlabel('Calcium in cytosol'); ylabel('Calcium in store');

%% IP3
figure(3);
plot(nv.T, nv.out('I_i_1'));
title('Cell 1 IP3');

figure(4);
plot(nv.T, nv.out('I_i_2'));
title('Cell 2 IP3');

%% Fluxes
figure(7);
plot(nv.T, nv.out('J_Ca_SMCcoup_i_2'));
title('Calcium coupling from cell 1 to 2');
xlabel('t'); ylabel('J_Ca_SMCcoup_i_2');

figure(8);
plot(nv.T, nv.out('J_Ca_SMCcoup_i_1'));
title('Calcium coupling from cell 2 to 1');
xlabel('t'); ylabel('J_Ca_SMCcoup_i_1');

%% IP3 in EC/SMC
figure(9);
plot(nv.T, nv.out('I_j_1'), nv.T, nv.out('I_j_2'), nv.T, nv.out('I_i_1'), nv.T, nv.out('I_i_2'));
title('IP3');
legend('EC 1', 'EC 2', 'SMC 1', 'SMC 2', 'Location','NorthEastOutside');

%% Membrane potential of SMCs
figure(10);
plot(nv.T, nv.out('v_i_1'), nv.T, nv.out('v_i_2'));

%% Other
figure(11);
plot(nv.T, nv.out('s_i_1'));
title('Cell 1 calcium in store');

figure(12);
plot(nv.T, nv.out('s_i_2'));
title('Cell 2 calcium in store');

%% fgdfgg
figure(13);
plot(nv.T, nv.out('s_i_1')+nv.out('Ca_i_1'));
