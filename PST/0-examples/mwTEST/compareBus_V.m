function  compareBus_V( caseA, caseB, printFigs )
%COMPAREBUS_V plotbus voltage comparison from PST data
%   Plot bus voltage from two case of PST data. Assumes time vector is the 
%   same, but accomodates for global variable.

% load data from case A
load(caseA)
% if exist('g', 'var')
%     % handle global variable
%     tA = g.sys.t;
%     if isfield(g,'mac')
%         varA = g.mac.mac_spd;
%     else
%         varA = mac_spd;
%     end
% else
    tA = t;
    varA = bus_v;
% end

% clear unrequried vars
clearvars -except tA varA caseA caseB printFigs

% load data from case B
load(caseB)
% if exist('g', 'var')
%     % handle global variable
%     tB = g.sys.t;
%     if isfield(g,'mac')
%         varB = g.mac.mac_spd;
%     else
%         varB = mac_spd;
%     end
% else
    tB = t;
    varB = bus_v;
% end

% manipulate case names for labels
nameA = strsplit(caseA,'.');
nameA = nameA{1};

nameB = strsplit(caseB,'.');
nameB = nameB{1};

figure
subplot(2,1,1)
plot(tA, abs(varA), 'linewidth', 1)
hold on
plot(tB, abs(varB), '--', 'linewidth',1.5)

title({'Bus Voltage Comparison'; [nameA,' (solid) - ', nameB, ' (dashed)']})

ylabel('Voltage [PU]')
xlabel('Time [sec]')

% plot absolute difference
subplot(2,1,2)
plot(tA, abs(varA-varB) )
title({ 'Absolute Voltage Difference'; [nameA,' - ', nameB] } )
ylabel('Voltage Difference [PU]')
xlabel('Time [sec]')

% pdf output code
if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([nameA,nameB],'-pdf'); % to print fig
end

end

