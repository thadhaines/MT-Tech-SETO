function y_switch

% y_switch
% 9:45 am July 20, 1998
% script file for calculating the reduced Y matrices
% for the various switching options
% part of s_simu
% calls red_ybus

% Author: Graham Rogers
% (c) copyright Cherry Tree Scientific Software/ J.H. Chow 1997 - 1998
% All rights reserved

% modified to add option 7, clear fault with out loss of line
% suggested by Trong Nguyen and M.A. Pai University of Illinois
% errors in line to line, line to ground and line to line to ground corrected

%{
functionalized
added returns from red_ybus to g.int as they are used in the i_simu (network interface code)
added appropriate lines and buses to g.line and g.bus respectively
%}

global g

%% create pre-fault admittance matrix
if ~isempty(g.ncl.load_con)
    g.ncl.nload = length(g.ncl.load_con(:,1));
end

[g.int.Y_gprf,g.int.Y_gncprf,g.int.Y_ncgprf,g.int.Y_ncprf,g.int.V_rgprf,g.int.V_rncprf,g.int.boprf] = red_ybus(g.bus.bus,g.line.line);

g.bus.bus_intprf = g.bus.bus_int;% store the internal bus numbers for the pre_fault system

%nbus = length(bus(:,1)); unused? -thad 07/17/20

%% create fault on matrices according to type
g.bus.bus_f = g.bus.bus;
g.line.line_f = g.line.line;

f_type = g.sys.sw_con(2,6); 

if f_type < 4 || f_type==7
    f_nearbus = g.sys.sw_con(2,2);
    
    g.bus.bus_idx = find(g.bus.bus_f(:,1)==f_nearbus); % global - used in livePlot
    
    if isempty(g.bus.bus_idx)
        error('faulted bus not specified correctly')
    end
    
    if f_type == 0 || f_type==7
        %three phase fault zero impedance to ground
        bf = 1.0/1e-7;
    elseif f_type ==1
        % line to ground fault zf = zn*z0/(zn+z0)
        xf = g.sys.sw_con(2,4)*g.sys.sw_con(2,5)/(g.sys.sw_con(2,4) + g.sys.sw_con(2,5));
        xf = max(xf,1e-7);
        bf = 1.0/xf;
    elseif f_type == 2
        % line to line to ground  zf = zn + z0
        xf = g.sys.sw_con(2,4) + g.sys.sw_con(2,5);
        xf = max(xf,1e-7);
        bf = 1.0/xf;
    elseif f_type==3
        %line to line  zf = zn
        xf = g.sys.sw_con(2,4);
        xf = max(xf,1e-7);
        bf = 1.0/xf;
    end
    g.bus.bus_f(g.bus.bus_idx(1),9) = -bf; % set fault bus B 
end

if f_type == 4
    % remove line with no fault
    f_nearbus = g.sys.sw_con(2,2);
    g.bus.bus_idx = find(g.bus.bus_f(:,1)==f_nearbus);
    f_farbus =  g.sys.sw_con(2,3);
    line_idx = find((g.line.line_f(:,1)==f_nearbus & g.line.line_f(:,2)==f_farbus)...
        | (g.line.line_f(:,2)==f_nearbus & g.line.line_f(:,1)==f_farbus));
    % choose first instance of line
    if ~isempty(line_idx)
        g.line.line_f(line_idx(1),4) = 1.0e7; %make line reactance 'infinite'
    else
        error('can not find faulted line in line data')
    end
end

if f_type == 5
    %loss of load
    f_nearbus = g.sys.sw_con(2,2);
    g.bus.bus_idx = find(g.bus.bus_f(:,1)==f_nearbus);
    if isempty(g.bus.bus_idx)
        error('can not find faulted bus in bus data')
    end
    % set p and q load to zero
    g.bus.bus_f(g.bus.bus_idx(1),6) = 0.0;
    g.bus.bus_f(g.bus.bus_idx(1),7) = 0.0;
end

if f_type == 6;
    %no fault
    f_nearbus = g.sys.sw_con(2,2);
    g.bus.bus_idx = find(g.bus.bus_f(:,1)==f_nearbus);
    if isempty(g.bus.bus_idx)
        error('can not find bus in bus data to do nothing to...')
    end
end

% form fault on reduced matrices
[g.int.Y_gf,g.int.Y_gncf,g.int.Y_ncgf,g.int.Y_ncf,g.int.V_rgf,g.int.V_rncf,g.int.bof] = red_ybus(g.bus.bus_f,g.line.line_f);   % fault-on admittance matrix 
g.bus.bus_intf = g.bus.bus_int;

%% second switching point, clear fault at near end/add new line
if f_type < 4
    f_farbus = g.sys.sw_con(2,3);
    g.line.line_pf1 = g.line.line;
    g.bus.bus_pf1 = g.bus.bus;
    line_idx = find( (g.line.line_pf1(:,1)==f_nearbus & g.line.line_pf1(:,2)==f_farbus)...
        | (g.line.line_pf1(:,2)==f_nearbus & g.line.line_pf1(:,1)==f_farbus));
    if isempty(line_idx)
        fb_str = num2str(f_farbus);
        fn_str = num2str(f_nearbus);
        disp(['can not find line between ',fn_str,' & ',fb_str,' in line data'])
        error('faulted line not specified correctly')
    end
    
    g.line.line_pf1(line_idx(1),4) = 1.0e7; %make faulted line reactance 'infinite'
    new_bus = max(g.bus.bus(:,1))+10; % add new faulted bus
    max_pf1b = length(g.bus.bus(:,1))+1;
    
    g.bus.bus_pf1(max_pf1b,1) = new_bus; % fault bus number
    g.bus.bus_pf1(max_pf1b,2) = 1.0;
    g.bus.bus_pf1(max_pf1b,3:7)=zeros(1,5);
    g.bus.bus_pf1(max_pf1b,9) = g.bus.bus_f(g.bus.bus_idx(1),9); % B
    g.bus.bus_pf1(max_pf1b,10) = 3; % bus type
    
    dlpf1 = length(g.line.line_pf1(:,1))+1; % add new line
    g.line.line_pf1(dlpf1,1)=new_bus;
    g.line.line_pf1(dlpf1,2)=f_farbus;
    g.line.line_pf1(dlpf1,3:6) = g.line.line(line_idx(1),3:6);
    
    [g.int.Y_gpf1,g.int.Y_gncpf1,g.int.Y_ncgpf1,g.int.Y_ncpf1,g.int.V_rgpf1,g.int.V_rncpf1,g.int.bopf1]...
        = red_ybus(g.bus.bus_pf1,g.line.line_pf1);  % post-fault
    % admittance matrix
    g.bus.bus_intpf1 = g.bus.bus_int;
    
elseif f_type==4 || f_type==5 || f_type==6
    % fault type is 4 or 5, 6 no change in system structure
    % set post fault data to fault data (bus_pf1 = bus_f)
    g.bus.bus_pf1 = g.bus.bus_f;
    g.line.line_pf1 = g.line.line_f;
    
    g.int.Y_gpf1 = g.int.Y_gf;
    g.int.Y_gncpf1 = g.int.Y_gncf;
    g.int.Y_ncgpf1 = g.int.Y_ncgf;
    g.int.Y_ncpf1 = g.int.Y_ncf;
    g.int.V_rgpf1 = g.int.V_rgf;
    g.int.V_rncpf1 = g.int.V_rncf;
    g.int.bopf1 = g.int.bof;
    
    g.bus.bus_intpf1 = g.bus.bus_intf;
    
elseif f_type == 7
    % clear fault
    g.bus.bus_pf1 = g.bus.bus;
    g.line.line_pf1 = g.line.line;
    
    g.int.Y_gpf1 = g.int.Y_gprf;
    g.int.Y_gncpf1 = g.int.Y_gncprf;
    g.int.Y_ncgpf1 = g.int.Y_ncgprf;
    g.int.Y_ncpf1 = g.int.Y_ncprf;
    g.int.V_rgpf1 = g.int.V_rgprf;
    g.int.V_rncpf1 = g.int.V_rncprf;
    g.int.bopf1 = g.int.boprf;
    
    g.bus.bus_intpf1 = g.bus.bus_intprf;
end

%% third switching point, clear fault at remote end
if f_type < 4
    g.line.line_pf2 = g.line.line_pf1;
    g.line.line_pf2(dlpf1,4) = 1.0e7;%open line
    
    g.bus.bus_pf2 = g.bus.bus_pf1;
    g.bus.bus_pf2(max_pf1b,9)=0.0;%remove short
    
    [g.int.Y_gpf2,g.int.Y_gncpf2,g.int.Y_ncgpf2,g.int.Y_ncpf2,g.int.V_rgpf2,g.int.V_rncpf2,g.int.bopf2]...
        = red_ybus(g.bus.bus_pf2,g.line.line_pf2);  % post-fault
    % admittance matrix
    g.bus.bus_intpf2 = g.bus.bus_int;
else
    % load type = 4 or 5, 6 or 7
    % no change in system structure
    g.bus.bus_pf2 = g.bus.bus_pf1;
    g.line.line_pf2 = g.line.line_pf1;
    
    g.int.Y_gpf2 = g.int.Y_gpf1;
    g.int.Y_gncpf2 = g.int.Y_gncpf1;
    g.int.Y_ncgpf2 = g.int.Y_ncgpf1;
    g.int.Y_ncpf2 = g.int.Y_ncpf1;
    g.int.V_rgpf2 = g.int.V_rgpf1;
    g.int.V_rncpf2 = g.int.V_rncpf1;
    g.int.bopf2 = g.int.bopf1;
    
    g.bus.bus_intpf2 = g.bus.bus_intpf1;
end
end