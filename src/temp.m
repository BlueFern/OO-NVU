odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Astrocyte('startpulse', 10, 'lengthpulse', 600), ...
    WallMechanics(), ...
    SMCEC('J_PLC',0.4), ...
    'odeopts', odeopts);

nv.simulate()

K_p = nv.out('K_p');
J_KIR_i = nv.out('J_KIR_i')';
J_BK_k = nv.out('J_BK_k')';
TT = nv.T;

starttime = floor(1/nv.T(length(nv.T))*length(nv.T));
endtime = length(nv.T);

figure(102);
hold on
plot(TT(starttime:endtime,:), J_KIR_i(starttime:endtime,:));
title('J_KIR_i');
xlim([1 600]);
hold off

figure(300);
plot(nv.T, nv.out('Ca_i'));
title('Ca_i - full model');

% figure(101);
% %hold on
% plot(TT(starttime:endtime,:), K_p(starttime:endtime,:));
% title('K_p');
% xlim([1 600]);
% %hold off

% figure(3);
% hold on
% plot(TT(starttime:endtime,:), J_BK_k(starttime:endtime,:));
% title('J_BK_k');
% xlim([1 600]);
% hold off