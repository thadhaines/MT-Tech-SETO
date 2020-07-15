%%  Simple script to plot voltage and angle from linear tests to see if
%   linear results are the same between versions of PST.
%   All lines should plot directly on top of eachother
%   Not the most accurate, but definitely the simplest.

%% load data
clear; close all; clc

dataName = 'linResults';
a = ['pstSETO', dataName];
b = ['pstV2P3', dataName];
c = ['pstV3p1', dataName];

nameCell = {a,b,c};

for n=1:length(nameCell)
    load(nameCell{1})
    data(n).v = linV;
    data(n).a = linAng;
    data(n).t = tL;
    clear linAng linV modSig tL
end

dataSetL = length(data);
sytles = {'-','--',':'};
figure

for n=1:dataSetL
    
    subplot(1,2,1)
    hold on
    plot(data(n).t, data(n).v, sytles{n},'linewidth',3/n)
    title('Bus Voltage')
    ylabel('Voltage [PU]')
    xlabel('Time [sec]')
    xlim([min(data(n).t),max(data(n).t)])
    
    subplot(1,2,2)
    hold on
    plot(data(n).t, data(n).a, sytles{n},'linewidth',3/n)
    title('Voltage Angle')
    ylabel('Angle [rad?]')
    xlabel('Time [sec]')
    xlim([min(data(n).t),max(data(n).t)])
end