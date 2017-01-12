%% Loop made for creating a figure that shows calcium concentration over 
%% increasing J_PLC for new NVU model version 1.2
% use this program to follow the movement of the Hopf bifurcation with
% changing pressure
% can be run multiple times with different settings

simulatiehopf.description=['old - pressure            ';     %1
                           'old - pressure no neural  ';     %2
                           'old - pressure no NO      ';     %3
                           'old - extrusion           ';     %4
                           'new - CICR                ';     %5
                           'new - stretch activated   ';     %6
                           'new - IP3 channel EC      '];    %7    
                       
simulatiehopf.settings  = ['pressure          ';
                           'VOCC              ';
                           'PMCA              ';
                           'extrusion         ';
                           'CICR              ';
                           'stretch activated ';
                           'IP3 mediated EC   '];


simulatiehopf.nummer=[1 2];
range=[500:500:8000];
for xxx=simulatiehopf.nummer;
    xxrange=1;
    for xxxx=range
        
        simpres=        [xxxx xxxx xxxx 1 1 1 1]; %P_T
        gluswitch=      [1 1 1 1 1 1 1];    %G_Ca_i   
        NOswitch=       [1 1 0 1 1 1 1]; %I_PMCA
        extrusion=0.24* [1 1 1 xxxx 1 1 1]; %D_i
        CICR= 55*       [1 1 1 1 xxxx 1 1]; %C_i
        stretch=6.1e-3* [1 1 1 1 1 xxxx 1]; %G_stretch
        IPendo=0.23*    [1 1 1 1 1 1 xxxx]; %F_j
        simulatiehopf.settingsvalue.newextrusion = [0 0 0 0 1 1 1]; 
        simulatiehopf.settingsvalue.startpulse= [100 2000 100 100 100 100];
        J_PLC_value=[0.0:0.0025:0.7];
        eenmaximumlist=zeros(1,length(J_PLC_value));
        eenminimumlist=zeros(1,length(J_PLC_value));
        meanlist=zeros(1,length(J_PLC_value));
        
        xxplc=1;
        for xx=J_PLC_value

            % Things to change:
            %XLIM1 = 50; XLIM2 = 500;            % x limits of any plots
            PRESSURE = simpres(xxx);             % Pressure in Pa, change here!! Because the parameter is used in both SMCEC and WallMechanics ***********, 4000 Pa is default
            PRESSURE_CHANGE=2000;                % Pressure at 100 sec till 300 seconds when switch is on, 2000 Pa is default
            START_PULSE  = simulatiehopf.settingsvalue.startpulse(xxx);% Start time of neuronal stimulation
            LENGTH_PULSE = 1900;       	         % Length of stimulation 
            GLU_SWITCH   = gluswitch(xxx);       % Turn on glutamate input to SC, default 1
            NO_SWITCH    = NOswitch(xxx);        % Turn on Nitric Oxide production, default 1
            TRPV4_SWITCH = 1;                    % Turn on TRPV4 channel flux in AC, default 1
            PLC_SWITCH   = 0;                    % Turn on PLC change in time, default 1
            PressureSwitch =0;                   % Turn on pressure change in time, default 1
            J_PLC_steadystate = xx;     % Jplc value in EC: 0.18 for steady state (default)
            J_PLC_vasomotion  = xx;     % Jplc value in EC:  0.4 for vasomotion   (default) from 100 to 300 seconds
            NEWEXTRUSIONSWITCH = simulatiehopf.settingsvalue.newextrusion(xxx);      %Turn on new extrusion based on Kapela, default 1

            odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);


            % Initiate the NVU
            nv = NVU(Neuron('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_SWITCH), ...
                Astrocyte('startpulse', START_PULSE, 'lengthpulse', LENGTH_PULSE, 'rhoSwitch', GLU_SWITCH, 'TRPV4switch', TRPV4_SWITCH), ...
                WallMechanics('PressureSwitch', PressureSwitch,'P_T', PRESSURE,'Pressure_change', PRESSURE_CHANGE), ...
                SMCEC('PressureSwitch', PressureSwitch, 'J_PLC_steadystate', J_PLC_steadystate, 'J_PLC_vasomotion',J_PLC_vasomotion, 'newextrusionswitch', NEWEXTRUSIONSWITCH,  ...
                'NOswitch', NO_SWITCH, 'P_T', PRESSURE, 'PLCSwitch', PLC_SWITCH, ...
                'Pressure_change', PRESSURE_CHANGE, 'D_i',extrusion(xxx), 'C_i',CICR(xxx),'G_stretch', stretch(xxx),'F_j',IPendo(xxx) ), 'odeopts', odeopts);


            nv.T = linspace(0, 1000, 10000);     % Initiate time vector
            nv.simulate(); % Run simulation


            Calciumconcentration=nv.out('Ca_i');
            eenmaximumlist(xxplc)=max(Calciumconcentration(9000:9950));
            eenminimumlist(xxplc)=min(Calciumconcentration(9000:9950));
            meanlist(xxplc)=mean(Calciumconcentration(9000:9950));
            xxplc=xxplc+1;
        end
        absolutlist=eenmaximumlist-eenminimumlist;
        first=find(absolutlist>0.001,1);
        
%         if isempty(first)==0
%             if first==1
%                 simulatiehopf.firsthopf(xxx,xxrange)=J_PLC_value(first);
%             else
%                 simulatiehopf.firsthopf(xxx,xxrange)=J_PLC_value(first-1);
%             end
%             simulatiehopf.factorfirsthopf(xxx,xxrange)=xxxx;
%             xxrange=xxrange+1;
%             
%         end
        
        if isempty(first)==0
            second=find(absolutlist(first+1:end)<0.001,1);
            if isempty(second)==0
                if first==1
                    simulatiehopf.firsthopf(xxx,xxrange)=J_PLC_value(first);
                else
                    simulatiehopf.firsthopf(xxx,xxrange)=J_PLC_value(first-1);
                end
                simulatiehopf.secondhopf(xxx,xxrange)=J_PLC_value(first+second);
                simulatiehopf.factorfirsthopf(xxx,xxrange)=xxxx;
                simulatiehopf.factorsecondhopf(xxx,xxrange)=xxxx;
                xxrange=xxrange+1;
            end
        end
        
    end

    
    
end

% save('simulatiehopf_metPa_old_stretchandersom_noneural_01stretchEC','simulatiehopf')
% clear all

lijst=[1];
legendlist=[];
for i=[1:length(lijst)]
    nummer=lijst(i);
    figure(7)
    kleur=['k','b','r','g'];
    hold on
    plot(simulatiehopf.firsthopf(nummer,:),simulatiehopf.factorfirsthopf(nummer,:), 'color',kleur(i),'LineWidth',2.5)
    plot(simulatiehopf.secondhopf(nummer,:),simulatiehopf.factorsecondhopf(nummer,:),'color',kleur(i),'LineWidth',2.5)
   
    simulatiehopf.description(nummer,:);
    if length(legendlist)==0;
        legendlist=[simulatiehopf.description(nummer,:)];
    else
        legendlist=[legendlist ',' simulatiehopf.description(nummer,:)];
    end
    ylim([500 8000])
    xlim([0.1 0.8])
    xlabel('J PLC [\muM/s]', 'FontSize', 22); 
    ylabel('Pressure [Pa]', 'FontSize', 22)
    box on
    title('The position of the Hopf bifurcations under changing pressure', 'FontSize', 22)
    set(gca,'FontSize',22,'LineWidth',2)
    
    plot([0.1 0.8],[4000 4000],'--', 'color','r','LineWidth',2.5)

end
legenda = strsplit(legendlist,',');
legend(legenda);