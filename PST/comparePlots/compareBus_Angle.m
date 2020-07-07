function  compareBus_Angle( caseA, caseB, printFigs )
%COMPAREBUS_Angle plot bus angle comparison from PST data
%   Plot bus voltage from two case of PST data. Assumes time vector is the 
%   same, but accomodates for global variable.

% load data from case A
load(caseA)
if exist('g', 'var')
    % handle global variable
    tA = g.sys.t;
    if isfield(g,'mac')
        varA = g.sys.theta;
    else
        varA = theta;
    end
else
    tA = t;
    varA = theta;
end

% clear unrequried vars
clearvars -except tA varA caseA caseB printFigs

% load data from case B
load(caseB)
if exist('g', 'var')
    % handle global variable
    tB = g.sys.t;
    if isfield(g,'mac')
        varB = g.sys.theta;
    else
        varB = theta;
    end
else
    tB = t;
    varB = theta;
end

% manipulate case names for labels
nameA = strsplit(caseA,'.');
nameA = nameA{1};

nameB = strsplit(caseB,'.');
nameB = nameB{1};

figure
subplot(2,1,1)
plot(tA, varA, 'linewidth', 1)
hold on
plot(tB, varB, '--', 'linewidth',1.5)

title({'Bus Angle Comparison'; [nameA,' (solid) - ', nameB, ' (dashed)']})

ylabel('Angle [PU]')
xlabel('Time [sec]')

% plot absolute difference
subplot(2,1,2)
plot(tA, abs(varA-varB) )
title({ 'Absolute Angle Difference'; [nameA,' - ', nameB] } )
ylabel('Angle Difference [PU]')
xlabel('Time [sec]')

% pdf output code
if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([nameA,nameB,'Angle'],'-pdf'); % to print fig
end

end

