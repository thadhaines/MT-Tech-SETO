%% test to use ode solvers to step PST-esq model
close all;clear;format compact;clc
% %% Octave specific solvers
% pkg load control
% pkg load signal
% solverSelection = {'ode45', 'ode23', 'ode15s'}; % OCTAVE
%%    %
printFigs = 1;
solverSelection = {'ode45', 'ode23', 'ode113','ode15s','ode23s','ode23t','ode23tb'}; % MATLAB

% 'ode15i' requires derivative at t=0... more thought required - availabe in octave aswell

%% pss model definition (miniWECC)
%           1   2   3   4       5   6   7       8       9       10
pss_con = [ 1  	1   20  2     0.25 0.04 0.2   0.03      1.0     -1.0];

%% MATLAB model
tend = 60;
block1 = tf([pss_con(3)*pss_con(4), 0],[pss_con(4), 1]);
block2 = tf([pss_con(5), 1],[pss_con(6), 1]);
block3 = tf([pss_con(7), 1],[pss_con(8), 1]);

G= block1*block2*block3;
tL = 0:1/60/4:tend; % quarter cycle steps
modSig = zeros(size(tL,1),1);
modSig(tL>=1) = .001; % very small input to avoid limiter
yL = lsim(G,modSig,tL);



%% stock solver attempt with statespace
for slnNum = 1:length(solverSelection)
    clear t1 t2 y1 y2
    odeName =solverSelection{slnNum};
    
    % Configure ODE settings
    %options = odeset('RelTol',1e-3,'AbsTol',1e-6); % default settings
    options = odeset('RelTol',1e-5,'AbsTol',1e-8,'InitialStep', 1/60/4, 'MaxStep',20);
    
    % manipulate test sytem to statespace
    [num,den] = tfdata(G);
    global A B U
    [A,B,C,D] = tf2ss(num{1},den{1});
    % initial conditions
    x = zeros(size(A,1),1);
    y0 = x;
    U = 0;
    
    % Pre-perturbance
    [t1,y1] = feval(odeName, @getXdot, [0,1-1/60/4],y0, options);
    yOut1 = C*y1'+D*U;
    
    % Step input
    U = modSig(end);
    [t2,y2] = feval(odeName, @getXdot, [1,tend],y1(end,:)', options);
    yOut2 = C*y2'+D*U;
    
    % combining output
    tCombined = [t1;t2];
    yCombined = [yOut1, yOut2];
    
    %% plotting fixed vs variable results
    figure
    % detail of event
    subplot(3,1,1)
    plot(tL,yL, 'r')
    hold on
    plot(tCombined,yCombined,'--','color','black')
    nPfix = length(tL);
    nPvar = length(tCombined);
    legend( {['Fixed time step (',int2str(nPfix),' points)'];...
        ['Variable time step (',int2str(nPvar),' points)']},'location','east')
    title(['System Output  - Event- ',odeName])
    xlim([0.95,1.25])
    set(gca,'fontsize',13); % font size
    
    % full time
    subplot(3,1,2)
    plot(tL,yL, 'r')
    hold on
    plot(tCombined,yCombined,'--','color','black')
    nPfix = length(tL);
    nPvar = length(tCombined);
%     legend( {['Fixed time step (',int2str(nPfix),' points)'];...
%         ['Variable time step (',int2str(nPvar),' points)']})
    title(['System Output - ',odeName])
    set(gca,'fontsize',13); % font size
    
    %% Plot Time step size
    tStep = zeros(length(tCombined),1);
    for tNdx = 2:length(tCombined)
        tStep(tNdx-1) = tCombined(tNdx)-tCombined(tNdx-1);
    end
    
    subplot(3,1,3)
    stairs(tCombined,tStep)
    title(['Variable Time Step Size- ',odeName])
    set(gca,'fontsize',13); % font size
    
    temp = gcf;
    newPos = temp.Position + [0, 0, temp.Position(3), 0]; % double width
    set(gcf,'Position',newPos)
    
    % pdf output code
    if printFigs
        set(gcf,'color','w'); % to remove border of figure
        export_fig(['step',odeName],'-pdf'); % to print fig
    end
    %% plot detail of figure....
    detailTime = [tend-1, tend];
    figure
    subplot(2,1,1)
    plot(tL,yL, 'r')
    hold on
    plot(tCombined,yCombined,'--','color','black')
    legend( {['Fixed time step'];...
        ['Variable time step']})
    title(['System Output - ',odeName,' - Detail'])
    xlim([detailTime])
    
    set(gca,'fontsize',13); % font size
    
    subplot(2,1,2)
    stairs(tCombined,tStep)
    xlim([detailTime])
    title(['Variable Time Step Size - ',odeName,' - Detail'])
    ymax = max( max(tStep(tCombined>tend-2)), max(tStep(end-2:end))  )*1.25;
    ylim([0,ymax])
    set(gca,'fontsize',13); % font size
    
    temp = gcf;
    newPos = temp.Position + [0, 0, temp.Position(3), 0]; % double width
    set(gcf,'Position',newPos)
    
    
    % pdf output code
    if printFigs
        set(gcf,'color','w'); % to remove border of figure
        export_fig(['stepDetail',odeName],'-pdf'); % to print fig
    end
    
end
%% plot absolute dif
% figure
% [tdif, ydif] = calcVarDif(tL, yL, tCombined, yCombined);
% plot(tdif,abs(ydif))
%
% title('Absolute difference')


%% Results
%{

Manipulation of solver output is required for correct state operation and
model output handling.

Solvers typically have a 'OutputFcn' option that may be useful in indexing,
time advancement, required state/output handling, and/or
network solution calls.

Other options exist that may also be useful in future developement.
%}