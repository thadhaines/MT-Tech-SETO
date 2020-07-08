function  compareB_cv( caseA, caseB, printFigs )
%COMPAREB_CV plot SVC B_cv comparison from PST data
%   Plot B_cv from two case of PST data. Assumes time vector is the 
%   same, but accomodates for global variable.

% load data from case A
load(caseA)
if exist('g', 'var')
    % handle global variable
    tA = g.sys.t;
    if isfield(g,'svc')
        varA = g.svc.B_cv;
    else
        varA = B_cv;
    end
else
    tA = t;
    varA = B_cv;
end

% clear unrequried vars
clearvars -except tA varA caseA caseB printFigs

% load data from case B
load(caseB)
if exist('g', 'var')
    % handle global variable
    tB = g.sys.t;
    if isfield(g,'svc')
        varB = g.svc.B_cv;
    else
        varB = B_cv;
    end
else
    tB = t;
    varB = B_cv;
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

title({'SVC Susceptance Comparison'; [nameA,' (solid) - ', nameB, ' (dashed)']})

ylabel('Susceptance [PU]')
xlabel('Time [sec]')

% plot absolute difference
subplot(2,1,2)
plot(tA, abs(varA-varB) )
title({ 'Absolute Susceptance Difference'; [nameA,' - ', nameB] } )
ylabel('B\_cv Difference [PU]')
xlabel('Time [sec]')

% pdf output code
if printFigs
    set(gcf,'color','w'); % to remove border of figure
    export_fig([nameA,nameB,'B-cv'],'-pdf'); % to print fig
end

end

