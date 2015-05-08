clear;

odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

% Important constants:
J_PLC_1 = 0.2;
J_PLC_2 = 0.3;
D_Ca = 0.1;
D_IP3 = D_Ca; % Optional IP3 coupling
D_v = D_Ca; % Optional membrane potential coupling

nv = NVU_coupled(Astrocyte_1(), Astrocyte_2(), ...
    WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1('J_PLC_1', J_PLC_1, 'D_Ca', D_Ca, 'D_IP3', D_IP3, 'D_v', D_v), ...
    SMCEC_2('J_PLC_2', J_PLC_2, 'D_Ca', D_Ca, 'D_IP3', D_IP3, 'D_v', D_v), ...
    'odeopts', odeopts);

nv.simulate();

% Calcium
figure(1);
plot(nv.T, nv.out('Ca_i_1'));
title('Cell 1 Ca_i');

figure(2);
plot(nv.T, nv.out('Ca_i_2'));
title('Cell 2 Ca_i');

%% IP3
figure(3);
plot(nv.T, nv.out('I_i_1'));
title('Cell 1 IP3');

figure(4);
plot(nv.T, nv.out('I_i_2'));
title('Cell 2 IP3');

%% Trajectories

% For plotting trajectory only when signal on
starttime = floor((nv.astrocyte_1.params.startpulse+80)/nv.T(length(nv.T))*length(nv.T));
endtime = floor((nv.astrocyte_1.params.startpulse+nv.astrocyte_1.params.lengthpulse)/nv.T(length(nv.T))*length(nv.T));

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

%% Fluxes
figure(7);
plot(nv.T, nv.out('J_Ca_SMCcoup_i_2'));
title('Calcium coupling from cell 1 to 2');
xlabel('t'); ylabel('J_Ca_SMCcoup_i_2');

figure(8);
plot(nv.T, nv.out('J_IP3_SMCcoup_i_2'));
title('IP3 coupling from cell 1 to 2');
xlabel('t'); ylabel('J_IP3_SMCcoup_i_2');