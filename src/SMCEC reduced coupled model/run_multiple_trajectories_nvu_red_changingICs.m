% Plot trajectories for multiple D for reduced model
% ICs taken from previous trajectory

% Initial coordinates - comment out if continuing a plot
clear;  
coordinates = [0.25,0.25,0.25,0.000015,0.1,0.1,-60,0.1,0.1,100e3,0.1,0.1, ...
             -75,0.1,0.25,0.25,0.25,0.000015,0.1,0.1,-60,0.1,0.1,100e3,0.1,0.1,-75,0.1];

odeopts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.1, 'Vectorized', 1);

% Important constants:
J_PLC_1 = 0.125;
J_PLC_2 = 0.125;
K_p_1 = 9100;
K_p_2 = 9000;

t = 0:0.1:1000;
startnum = 800;     % Starting time to plot trajectories (so ignore transient behaviour)

coords = zeros(1,5);
counter = 0;


startofD = 0.6;
endofD = 1;
space = 0.05;
DD = startofD:space:endofD;
num_iterations = length(DD);
time_per_iteration = zeros(1,length(DD));
i = 1;

fprintf('%d total iterations \n', num_iterations);


for D=DD   
    
    % SMC coupling
    D_Ca_i = D;
    D_IP3_i = 0;   % Optional IP3 coupling
    D_v_i = 0;     % Optional membrane potential coupling
    
    % EC coupling - optional!!
    D_Ca_j = 0;
    D_IP3_j = 0;   % Optional IP3 coupling
    D_v_j = 0;     % Optional membrane potential coupling

    nv = NVU_coupled_red( WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1_red('J_PLC_1', J_PLC_1, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_1', K_p_1), ...
    SMCEC_2_red('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_2', K_p_2), ...
    'odeopts', odeopts, 'T', t);

    % ICs - last point of previous trajectory
    for i = 1:4
        nv.wall_1.u0(i) = coordinates(i);
        nv.wall_2.u0(i) = coordinates(i+14);
    end
    for i = 1:10
        nv.smcec_1.u0(i) = coordinates(i+4);
        nv.smcec_2.u0(i) = coordinates(i+18);
    end

    nv.simulate();
    
    starttime = floor(startnum/nv.T(length(nv.T))*length(nv.T));
    endtime = length(nv.T);
    
    Ca_i_1 = nv.out('Ca_i_1'); s_i_1 = nv.out('s_i_1');
    Ca_i_2 = nv.out('Ca_i_2'); s_i_2 = nv.out('s_i_2');
    
    DD = D*ones(length(Ca_i_1),1); 
    DD = DD(starttime:endtime,:);
    
    m = horzcat(Ca_i_1(starttime:endtime,:), ...
    s_i_1(starttime:endtime,:), ...
    Ca_i_2(starttime:endtime,:), ...
    s_i_2(starttime:endtime,:), DD);
    coords = vertcat(coords,m);

    get_coordinates();
    
    counter = counter + 1;
    amount_done = counter/num_iterations;
    percentage_done = amount_done*100;
    time_per_iteration(i) = nv.elapsedtime;
    time_left = mean(nonzeros(time_per_iteration))*(num_iterations-counter);
    min = floor(time_left/60);
    sec = rem(time_left, 60);
    i = i + 1;
    
    fprintf('%.1f%% done, %d iterations left, ETA: %.f min %.f sec \n', percentage_done, num_iterations-counter, min, sec);
end

figure(3);
%clf
hold on
for i=1:counter
    plot3( coords((i-1)*length(m)+i+1 : i*length(m) , 1) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 5) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 2) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 1');
view([70 25])
%xlim([0 1]); zlim([0 1.5]);

figure(28);
%clf
hold on
for i=1:counter
    plot3( coords((i-1)*length(m)+i+1 : i*length(m) , 3) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 5) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 4) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 2');
view([70 25])

