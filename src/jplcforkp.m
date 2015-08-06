% Find variation in Kp with JPLC changing

clear;
odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 1, 'Vectorized', 1);

values = zeros(1,1);
counter = 0;

range = 0:0.01:0.8;

for JPLC=range
    
    % Signal off
%     nv = NVU(Astrocyte('startpulse', 700), ...
%     WallMechanics(), ...
%     SMCEC('J_PLC', JPLC), ...
%     'odeopts', odeopts);

    % Signal on
    nv = NVU(Astrocyte('startpulse', 10, 'lengthpulse', 690), ...
    WallMechanics(), ...
    SMCEC('J_PLC', JPLC), ...
    'odeopts', odeopts);

    nv.simulate();
    
    Kp = nv.out('K_p');
    values = vertcat(values, Kp(length(Kp)));
    counter = counter + 1
end

values = values(2:length(values));

figure(124);
plot(range, values);
maxi = max(values)
mini = min(values)
average = mean(values)