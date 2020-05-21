function [ ] = stability_plot( data_A, A, data_B, B, data_C, C, title_str, plot_end )
%stability_plot Make plot given PST data of generator speed and power
%   data_X = string to .mat file with t, pelect and mac_spd data - plots
%   data in #2.  Assumes A stable, B borderline, and C unstable
%   X = label for data loaded
%   plot_end = end of x-axis in plots

bfz = 13;   % base font size
lt = 1.5;   % thickness of plotted lines
loc = 'southeast'; % location of all legends

cm = [ ...
    .75 .75 .75; % stable
    0 0 0;  % borderline
    1 0 1]; % unstable

load(data_A)
A_p = pelect(2,:);
A_s = mac_spd(2,:);
A_t = t;

load(data_B)
B_p = pelect(2,:);
B_s = mac_spd(2,:);
B_t = t;

load(data_C)
C_p = pelect(2,:);
C_s = mac_spd(2,:);
C_t = t;

% Plots
figure
subplot(2,1,1)
grid on; hold on
plot(A_t,A_p,'LineWidth',lt, 'color',cm(1,:))
plot(B_t,B_p,'LineWidth',lt, 'color',cm(2,:))
plot(C_t,C_p,'LineWidth',lt, 'color',cm(3,:))
legend({A,B,C},'location',loc)
ylabel('P_E [pu]','fontsize',bfz);
xlabel('Time [sec]','fontsize',bfz);
xlim([0, plot_end])
ylim([-.32,.32])
set(gca,'fontsize',bfz)
title(title_str,'fontsize',bfz)

% Generator Speed
subplot(2,1,2)
grid on; hold on
plot(A_t,A_s,'LineWidth',lt, 'color',cm(1,:))
plot(B_t,B_s,'LineWidth',lt, 'color',cm(2,:))
plot(C_t,C_s,'LineWidth',lt, 'color',cm(3,:))
xlim([0 plot_end])
legend({A,B,C},'location','northeast')
ylabel('Speed [pu]','fontsize',bfz);
xlabel('Time [sec]','fontsize',bfz);
set(gca,'fontsize',bfz)
ylim([.99,1.02])
%set(gcf,'Position',[680 250 750 720])

end

