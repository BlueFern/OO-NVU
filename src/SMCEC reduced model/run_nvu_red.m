clear;

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.05, 'Vectorized', 1);

%K_p = 9300;     % Signal on
%K_p = 3400;     % Signal off

K_p =9200;

JPLC = 0.05;

% Time
T = 0:0.1:5000;

nv = NVU_red(WallMechanics(), ...
    SMCEC_red('J_PLC', JPLC, 'K_p', K_p), ...
    'odeopts', odeopts, 'T', T);

nv.simulate;

figure(400);
%hold on
plot(nv.T, nv.out('Ca_i'));
title('Ca_i - reduced model');
%xlim([200 300]);
%hold off
% 
% figure(100);
% hold on
% plot(nv.T, nv.out('I_i'));
% title('IP3');
% xlim([200 300]);
% hold off

% figure(200);
% hold on
% plot(nv.T, nv.out('J_stretch_i'));
% title('J_stretch_i');
% xlim([900 1000]);
% hold off




% 
% figure(401);
% plot(nv.out('Ca_i'), nv.out('s_i'));
% xlabel('Ca_i');
% ylabel('s_i');
% ylim([1 1.4]);

% figure(41);
% plot(nv.T, nv.out('J_KIR_i'));
% title('J_KIR_i - reduced model');