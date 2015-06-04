% Plot trajectories for multiple D for reduced model

clear;   

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.1, 'Vectorized', 1);

% Important constants:
J_PLC_1 = 0.5;
J_PLC_2 = 0.6;

t = 0:0.1:1000;

coords = zeros(1,5);
counter = 0;

K_p = 9300;

startofD = 0;
endofD = 0.15;
space = 0.025;
DD = startofD:space:endofD;
num_iterations = length(DD);

for D=startofD:space:endofD   % Choose values of D to plot trajectories for
    
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
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_1', K_p), ...
    SMCEC_2_red('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j, 'K_p_2', K_p), ...
    'odeopts', odeopts, 'T', t);

    nv.simulate();
    
    starttime = floor(300/nv.T(length(nv.T))*length(nv.T));
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

    counter = counter + 1;
    amount_done = counter/num_iterations;
    percentage_done = amount_done*100;
    fprintf('%.2f%% done, %d iterations left \n', percentage_done, num_iterations-counter);
end

figure(11);
clf
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
xlim([0 1]); zlim([0 1.5]);

figure(21);
clf
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
xlim([0 1]); zlim([0 1.5]);

