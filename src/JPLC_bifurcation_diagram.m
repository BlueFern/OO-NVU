% Bifurcation diagram for single NVU, with J_PLC as parameter

clear;   

odeopts = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.1, 'Vectorized', 1);

coords = zeros(1,3);
counter = 0;

for JPLC=0.55:0.005:0.6
    
    % Signal off
    %nv = NVU(Astrocyte('startpulse', 600), ...
    %WallMechanics(), ...
    %SMCEC('J_PLC', JPLC), ...
    %'odeopts', odeopts);

    % Signal on
    nv = NVU(Astrocyte('startpulse', 10, 'lengthpulse', 590), ...
    WallMechanics(), ...
    SMCEC('J_PLC', JPLC), ...
    'odeopts', odeopts);

    nv.simulate();
    
    starttime = floor(300/nv.T(length(nv.T))*length(nv.T));
    endtime = length(nv.T);
    
    Ca_i = nv.out('Ca_i'); s_i = nv.out('s_i');
    
    JJ = JPLC*ones(length(Ca_i),1); 
    JJ = JJ(starttime:endtime,:);
    
    m = horzcat(Ca_i(starttime:endtime,:), ...
    s_i(starttime:endtime,:), JJ);
    coords = vertcat(coords,m);

    counter = counter+1
end

figure(1);
hold on
for i=1:counter
    plot3( coords((i-1)*length(m)+i+1 : i*length(m) , 1) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 3) , ...
        coords((i-1)*length(m)+i+1 : i*length(m) , 2) , 'r', 'LineWidth', 1, 'Marker', ...
        's', 'MarkerSize', 1);
 end
hold off
xlabel('Calcium in cytosol Ca_i'); ylabel('J_{PLC}'); zlabel('Calcium in store s_i');
%grid on
title('Bifurcation diagram');
view([90 -90])
% xlim([0 1]); zlim([0 1.5]); 
% ylim([0 1]);

%xlim([0 1]); ylim([0.54 0.66]);



