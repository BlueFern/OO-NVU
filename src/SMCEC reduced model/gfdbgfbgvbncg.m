odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 0.5, 'Vectorized', 1);

nv = NVU_red(WallMechanics(), ...
    SMCEC_red('J_PLC', 0.4, 'K_p', 9300), ...
    'odeopts', odeopts);

nv.simulate;

figure(400);
plot(nv.T, nv.out('Ca_i'));
title('Ca_i - reduced model');

figure(401);
plot(nv.out('Ca_i'), nv.out('s_i'));
xlabel('Ca_i');
ylabel('s_i');

% figure(41);
% plot(nv.T, nv.out('J_KIR_i'));
% title('J_KIR_i - reduced model');