% Plot trajectories for multiple D

clear;

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.1, 'Vectorized', 1);

% Important constants:
J_PLC_1 = 0.35;
J_PLC_2 = 0.4;

coords1 = zeros(1,5);
coords2 = zeros(1,5);
counter = 0;

for D=0:0.04:0.4     % Choose values of D to plot trajectories for
    
    % SMC coupling
    D_Ca_i = D;
    D_IP3_i = 0;   % Optional IP3 coupling
    D_v_i = 0;     % Optional membrane potential coupling
    
    % EC coupling - optional!!
    D_Ca_j = 0;
    D_IP3_j = D_Ca_j;   % Optional IP3 coupling
    D_v_j = D_Ca_j;     % Optional membrane potential coupling

    nv = NVU_coupled(Astrocyte_1(), Astrocyte_2(), ...
    WallMechanics_1(), WallMechanics_2(), ...
    SMCEC_1('J_PLC_1', J_PLC_1, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j), ...
    SMCEC_2('J_PLC_2', J_PLC_2, 'D_Ca_i', D_Ca_i, 'D_IP3_i', D_IP3_i, ...
    'D_v_i', D_v_i, 'D_Ca_j', D_Ca_j, 'D_IP3_j', D_IP3_j, 'D_v_j', D_v_j), ...
    'odeopts', odeopts);
    nv.simulate();

    starttime = floor((nv.astrocyte_1.params.startpulse+50)/nv.T(length(nv.T))*length(nv.T));
    endtime = floor((nv.astrocyte_1.params.startpulse+nv.astrocyte_1.params.lengthpulse)/nv.T(length(nv.T))*length(nv.T));
    starttime2 = floor((nv.astrocyte_1.params.startpulse+nv.astrocyte_1.params.lengthpulse+20)/nv.T(length(nv.T))*length(nv.T));
    
    Ca_i_1 = nv.out('Ca_i_1'); s_i_1 = nv.out('s_i_1');
    Ca_i_2 = nv.out('Ca_i_2'); s_i_2 = nv.out('s_i_2');
    
    DD = D*ones(length(Ca_i_1),1); 
    DD1 = DD(starttime:endtime,:);
    DD2 = DD(starttime2:length(Ca_i_1),:);
    
    m1 = horzcat(Ca_i_1(starttime:endtime,:), ...
    s_i_1(starttime:endtime,:), ...
    Ca_i_2(starttime:endtime,:), ...
    s_i_2(starttime:endtime,:), DD1);
    coords1 = vertcat(coords1,m1);
    
    m2 = horzcat(Ca_i_1(starttime2:length(Ca_i_1),:), ...
    s_i_1(starttime2:length(s_i_1),:), ...
    Ca_i_2(starttime2:length(Ca_i_2),:), ...
    s_i_2(starttime2:length(s_i_2),:), DD2);
    coords2 = vertcat(coords2, m2);
    
    counter = counter+1
end

figure(110);
clf
hold on
for i=1:counter
    plot3( coords1((i-1)*length(m1)+i+1 : i*length(m1) , 1) , ...
        coords1((i-1)*length(m1)+i+1 : i*length(m1) , 5) , ...
        coords1((i-1)*length(m1)+i+1 : i*length(m1) , 2) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 1 signal on');
view([70 25])
xlim([0 1]); zlim([0 1.5]);

figure(151);
clf
hold on
for i=1:counter
    plot3( coords1((i-1)*length(m1)+i+1 : i*length(m1) , 3) , ...
        coords1((i-1)*length(m1)+i+1 : i*length(m1) , 5) , ...
        coords1((i-1)*length(m1)+i+1 : i*length(m1) , 4) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 2 signal on');
view([70 25])
xlim([0 1]); zlim([0 1.5]);

figure(120);
clf
hold on
for i=1:counter
    plot3( coords2((i-1)*length(m2)+i+1 : i*length(m2) , 1) , ...
        coords2((i-1)*length(m2)+i+1 : i*length(m2) , 5) , ...
        coords2((i-1)*length(m2)+i+1 : i*length(m2) , 2) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 1 signal off');
view([70 25])
xlim([0 1]); zlim([0 1.5]);

figure(161);
clf
hold on
for i=1:counter
    plot3( coords2((i-1)*length(m2)+i+1 : i*length(m2) , 3) , ...
        coords2((i-1)*length(m2)+i+1 : i*length(m2) , 5) , ...
        coords2((i-1)*length(m2)+i+1 : i*length(m2) , 4) );
end
hold off
xlabel('Calcium in cytosol'); ylabel('D'); zlabel('Calcium in store');
grid on
title('Cell 2 signal off');
view([70 25])
xlim([0 1]); zlim([0 1.5]);
