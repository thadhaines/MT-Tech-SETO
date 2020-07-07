function  comparePSS_Out( caseA, caseB, printFigs )
%COMPAREPSS_OUT plot pss_out comparison from PST data
%   Plot pss_out from two case of PST data. Assumes time vector is the 
%   same, but accomodates for global variable.

% load data from case A
load(caseA)
if exist('g', 'var')
    % handle global variable
    if isfield(g,'sys')
        tA = g.sys.t;
    else
        tA = t;
    end
    if isfield(g,'pss')
        varA = g.pss.pss_out;
    else
        varA = pss_out;
    end
else
    tA = t;
    varA = pss_out;
end

% clear unrequried vars
clearvars -except tA varA caseA caseB printFigs

% load data from case B
load(caseB)
if exist('g', 'var')
    % handle global variable
    if isfield(g,'sys')
        tB = g.sys.t;
    else
        tB = t;
    end
    if isfield(g,'pss')
        varB = g.pss.pss_out;
    else
        varB = pss_out;
    end
else
    tB = t;
    varB = pss_out;
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

title({'PSS Out Comparison'; [nameA,' (solid) - ', nameB, ' (dashed)']})

ylabel('Output')
xlabel('Time [sec]')

% plot absolute difference
subplot(2,1,2)
plot(tA, abs(varA-varB) )
title({ 'Absolute PSS Out Difference'; [nameA,' - ', nameB] } )
ylabel('Difference [PU]')
xlabel('Time [sec]')

% pdf output code
if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([nameA,nameB,'PSSout'],'-png'); % to print fig
end

end

