% Bifurcation diagram for single NVU, with J_PLC as parameter
% REDUCED MODEL

clear;   

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.1, 'Vectorized', 1);

coords = zeros(1,3);
counter = 0;

% Time
T = 0:0.1:5000;
startnum = 1000;

%K_p = 9300;     % Signal On
%K_p = 3400;   % Signal Off

K_p = 9100;

jplc_it = 0.09:0.005:0.125;

num_iterations = length(jplc_it);
time_per_iteration = zeros(1,length(jplc_it));
i = 1;

fprintf('%d total iterations \n', num_iterations);

for JPLC=jplc_it
    
    nv = NVU_red(WallMechanics(), ...
    SMCEC_red('J_PLC', JPLC, 'K_p', K_p), ...
    'odeopts', odeopts, 'T', T);

    nv.simulate();
    
    starttime = floor(startnum/nv.T(length(nv.T))*length(nv.T));
    endtime = length(nv.T);
    
    Ca_i = nv.out('Ca_i'); s_i = nv.out('s_i');
    
    JJ = JPLC*ones(length(Ca_i),1); 
    JJ = JJ(starttime:endtime,:);
    
    m = horzcat(Ca_i(starttime:endtime,:), ...
    s_i(starttime:endtime,:), JJ);
    coords = vertcat(coords,m);

    counter = counter + 1;
    amount_done = counter/num_iterations;
    percentage_done = amount_done*100;
    time_per_iteration(i) = nv.elapsedtime;
    time_left = mean(nonzeros(time_per_iteration))*(num_iterations-counter);
    min = floor(time_left/60);
    sec = rem(time_left, 60);
    i = i + 1;
    
    fprintf('%.2f%% done, %d iterations left, ETA: %.f min %.f sec \n', percentage_done, num_iterations-counter, min, sec);

end

figure(1111);
hold on
for i=1:counter
    plot3( coords((i-1)*length(m)+i+1 : i*length(m) , 1) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 3) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 2) , 'b', 'LineWidth', 5);
 end
hold off
xlabel('Calcium in cytosol Ca_i'); ylabel('J_{PLC}'); zlabel('Calcium in store s_i');
%grid on
title('Bifurcation diagram');
view([90 -90])
% xlim([0 1]); 
%zlim([0.1 0.24]); 
% ylim([0 1]);

%xlim([0 1]); ylim([0.54 0.66]);

save jplcbifurc.mat

