%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model. It includes:
% - basic K+ / Ca2+ 

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.
clear all
close all
clc
%odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1, 'Vectorized', 1);

nv = NVU(Neuron(), ...
    Astrocyte(), ...
    WallMechanics(), ...
    SMCEC('J_PLC', 0.18));
    %'odeopts', odeopts);

nv.T = linspace(0, 1000, 2000);

%nv.wall.params.P_T = 4000;
nv.astrocyte.params.switchBK = 0;

nv.astrocyte.params.trpv_switch=0;
nv.simulate()


%% plot all state variables:
for i = 1:2
    figure(1)
    set(gcf,'name','Neuron and Astrocyte')
    subplot(5,4,i)
    plot(nv.T, nv.out(nv.neuron.varnames{i}))
    h1(i) = gca();
    ylabel(nv.neuron.varnames{i})
    hold on
end
    linkaxes(h1, 'x');

for i = 1:17
    figure(1)
    subplot(5,4,2+i)
    plot(nv.T, nv.out(nv.astrocyte.varnames{i}))
    h1(i) = gca();
    ylabel(nv.astrocyte.varnames{i})
    hold on
end
    linkaxes(h1, 'x');


for i = 1:10
    figure(2)
    set(gcf,'name','ECSMC & WallMechanics')
    subplot(3,5,i)
    plot(nv.T, nv.out(nv.smcec.varnames{i}))
    h2(i) = gca();
    ylabel(nv.smcec.varnames{i})
    hold on
end

    
for i = 1:4
    figure(2)
    subplot(3,5,10+i)
    plot(nv.T, nv.out(nv.wall.varnames{i}))
    h2(10+i) = gca();
    ylabel(nv.wall.varnames{i})
    hold on
end
    linkaxes(h2, 'x');
    hold on
    
    

    
    
    figure(3)
    set(gcf,'name','Stretch Figures')
    subplot(3,4,1)
    plot(nv.T, nv.outputs{3}(nv.smcec.idx_out.J_stretch_i,:))
    h3(1) = gca();
    ylabel('J_stretch_i')
    hold on
    
    figure(3)
    subplot(3,4,2)
    plot(nv.T, nv.outputs{4}(nv.wall.idx_out.F_r,:))
    h3(2) = gca();
    ylabel('F_r')
    hold on   

    linkaxes(h3, 'x');