% Find variation in Kp with JPLC changing

odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

values = zeros(1,1);

for JPLC=0:0.05:0.8
    
    % Signal off
    %nv = NVU(Astrocyte('startpulse', 700), ...
    %WallMechanics(), ...
    %SMCEC('J_PLC', JPLC), ...
    %'odeopts', odeopts);

    % Signal on
    nv = NVU(Astrocyte('startpulse', 10, 'lengthpulse', 690), ...
    WallMechanics(), ...
    SMCEC('J_PLC', JPLC), ...
    'odeopts', odeopts);

    nv.simulate();
    
    Kp = nv.out('K_p');
    values = vertcat(values, Kp(length(Kp)));
end

