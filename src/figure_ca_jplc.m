%% Loop made for creating a figure that shows calcium concentration over 
%% increasing J_PLC for new NVU model version 1.2

simulatie.nummer=[25 26];%
simulatie.description   = ['old - normal    ';     %1
                           'new - normal    ';     %2
                           'new - VOCC*0.5  ';     %3
                           'new - VOCC*0.60 ';     %4
                           'new - VOCC*0.75 ';     %5
                           'new - VOCC*0.90 ';     %6
                           'new - VOCC*0.10 ';     %7
                           'new - VOCC*0.01 ';     %8
                           'new - PMCA*0.50 ';     %9
                           'new - PMCA*0.10 ';     %10
                           'new - PMCA*0.05 ';     %11
                           'new - PMCA*0.001';     %12
                           'new - PMCA*0.005';     %13
                           'old - extr*0.50 ';     %14
                           'old - extr*0.10 ';     %15
                           'old - extr*0.05 ';     %16
                           'old - extr*0.001';     %17
                           'old - extr*0.005';     %18
                           'old - VOCC*0.5  ';     %19
                           'old - VOCC*0.60 ';     %20
                           'old - VOCC*0.75 ';     %21
                           'old - VOCC*0.90 ';     %22
                           'old - VOCC*0.10 ';     %23
                           'old - VOCC*0.01 ';     %24
                           'old - no neural ';     %25
                           'new - no neural '];    %26   
                       
simulatie.settings      = ['nv.SMCEC.params.newextrusionswitch';
                           'nv.SMCEC.params.G_Ca_i            ';
                           'nv.SMCEC.params.I_PMCA            ';
                           'nv.SMCEC.params.D_i               ';
                           'startpulse                        '];
 VOCC=1.29e-3;
 IPMCA=5.37;
 D_i=0.24;
 normalstart=100;
simulatie.settingsvalue.newextrusion = [0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1];
simulatie.settingsvalue.G_Ca_i       = [VOCC VOCC VOCC*0.5 VOCC*0.6 VOCC*0.75 VOCC*0.90 VOCC*0.10 VOCC*0.01 VOCC VOCC VOCC VOCC VOCC ...
                                        VOCC VOCC VOCC VOCC VOCC VOCC*0.5 VOCC*0.6 VOCC*0.75 VOCC*0.90 VOCC*0.10 VOCC*0.01 VOCC VOCC];
simulatie.settingsvalue.I_PMCA       = [IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA*0.5 IPMCA*0.1 IPMCA*0.05 IPMCA*0.001 IPMCA*0.005 ...
                                        IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA IPMCA];
simulatie.settingsvalue.D_i          = [D_i D_i D_i D_i D_i D_i D_i D_i D_i D_i D_i D_i D_i ...
                                        D_i*0.50 D_i*0.10 D_i*0.05 D_i*0.001 D_i*0.005 D_i D_i D_i D_i D_i D_i D_i D_i];
simulatie.settingsvalue.startpulse   = [normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart ...
                                        normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart normalstart+400 normalstart+400];

for xxx=simulatie.nummer
    eenmaximumlist=[];
    eenminimumlist=[];
    J_PLC_value=[0.1:0.025:0.8];
    %simulatie.settings(1,:)=simulatie.settingsvalue.newextrusion(xxx);
    %simulatie.settings(2,:)=simulatie.settingsvalue.G_Ca_i(xxx);
    Calciumconcentration=[];
    for xx=J_PLC_value
    
        % Things to change:
        %XLIM1 = 50; XLIM2 = 500;            % x limits of any plots
        PRESSURE = 4000;                     % Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********, 4000 Pa is default
        PRESSURE_CHANGE=2000;                % Pressure at 100 sec till 300 seconds when switch is on, 2000 Pa is default
        START_PULSE  = simulatie.settingsvalue.startpulse(xxx);                  % Start time of neuronal stimulation
        LENGTH_PULSE = 300;       	         % Length of stimulation 
        GLU_SWITCH   = 1;                    % Turn on glutamate input to SC, default 1
        NO_SWITCH    = 1;                    % Turn on Nitric Oxide production, default 1
        TRPV4_SWITCH = 1;                    % Turn on TRPV4 channel flux in AC, default 1
        PLC_SWITCH   = 0;                    % Turn on PLC change in time, default 1
        PressureSwitch =0;                   % Turn on pressure change in time, default 1
        J_PLC_steadystate = xx;     % Jplc value in EC: 0.18 for steady state (default)
        J_PLC_vasomotion  = xx;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds
        NEWEXTRUSIONSWITCH = simulatie.settingsvalue.newextrusion(xxx);      %Turn on new extrusion based on Kapela, default 1
        
        odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);


        % Initiate the NVU
        nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
            Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
            WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
            SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion',J_PLC_vasomotion, 'newextrusionswitch', NEWEXTRUSIONSWITCH,  ...
            'G_Ca_i',simulatie.settingsvalue.G_Ca_i(xxx),  'I_PMCA', simulatie.settingsvalue.I_PMCA(xxx), 'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, ...
            'Pressure_change', PRESSURE_CHANGE, 'D_i',simulatie.settingsvalue.D_i(xxx) ), 'odeopts', odeopts);


        nv.T = linspace(0, 500, 5000);     % Initiate time vector
        nv.simulate(); % Run simulation


        Calciumconcentration=nv.out('Ca_i');
        eenmaximumlist=[eenmaximumlist max(Calciumconcentration(3000:4000))];
        eenminimumlist=[eenminimumlist min(Calciumconcentration(3000:4000))];

    end
    simulatie.maximumlist(xxx,:)=eenmaximumlist;
    simulatie.minimumlist(xxx,:)=eenminimumlist;
    simulatie.J_PLC_value(xxx,:)=J_PLC_value;
    
end


for nummer=[6]
    figure(1)
    hold on
    plot(simulatie.J_PLC_value(nummer,:),simulatie.maximumlist(nummer,:), 'color','k')
    plot(simulatie.J_PLC_value(nummer,:),simulatie.minimumlist(nummer,:), 'color','k')
end


