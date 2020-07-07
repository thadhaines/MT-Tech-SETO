function  compareMac_Spd( caseA, caseB, printFigs )
%COMPAREMAC_SPD plot machine speed comparison from PST data
%   Plot machine speed from two case of PST data. Assumes time vector is the 
%   same, but accomodates for global variable.

% load data from case A
load(caseA)
if exist('g', 'var')
    % handle global variable
    tA = g.sys.t;
    if isfield(g,'mac')
        spdA = g.mac.mac_spd;
    else
        spdA = mac_spd;
    end
else
    tA = t;
    spdA = mac_spd;
end

% clear unrequried vars
clearvars -except tA spdA caseA caseB printFigs

% load data from case B
load(caseB)
if exist('g', 'var')
    % handle global variable
    tB = g.sys.t;
    if isfield(g,'mac')
        spdB = g.mac.mac_spd;
    else
        spdB = mac_spd;
    end
else
    tB = t;
    spdB = mac_spd;
end

% manipulate case names for labels
nameA = strsplit(caseA,'.');
nameA = nameA{1};

nameB = strsplit(caseB,'.');
nameB = nameB{1};

figure
subplot(2,1,1)
plot(tA, spdA, 'linewidth', 1)
hold on
plot(tB, spdB, '--', 'linewidth',1.5)

title({'Machine Speed Comparison'; [nameA,' (solid) - ', nameB, ' (dashed)']})

ylabel('Speed [PU]')
xlabel('Time [sec]')

% plot absolute difference
subplot(2,1,2)
plot(tA, abs(spdA-spdB) )
title({ 'Absolute Speed Difference'; [nameA,' - ', nameB] } )
ylabel('Speed Difference [PU]')
xlabel('Time [sec]')

% pdf output code
if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([nameA,nameB,'MacSpd'],'-pdf'); % to print fig
end

end

