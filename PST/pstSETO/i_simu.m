function h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
%I_SIMU Purpose forms the network interface variables
% I_SIMU  Purpose forms the network interface variables
%
% Syntax: h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
%
%   NOTES:  Called by s_simu
%           Calls nc_load
% 
%   Input: 
%   k - the current time step
%  	ks - indicates the switching times
%  	k_inc - the number of time seps between switching points
%  	h - vector of time steps
%  	bus_sim - value of bus matrix at this switching time
%  	Y_g - reduced Y matrix for generators
%  	Y_gnc - mutual reduced Y generators-nc loads
%  	Y_ncg - mutual reduced Y nc loads generators
%   Y_nc - reduced Y matrix nc loads
% 	rec_V1 - voltage recovery matrix generators
% 	rec_V2 - voltage recovery matrix nc loads
%  	bo - bus order for this switching time
%
%   Output: 
%   h_sol - the time step at this value of ks
%
%   History:
%   Date        Time    Engineer        Description
%   03/xx/97    XX:XX   Graham Rogers  	Version 1.0
%   08/xx/97    XX:XX   xxx           	Version 1.1 - add induction generator
%   Copyright (c) Joe Chow All Rights Reserved
%   07/15/20    11:53   Thad Haines     Revised format of globals and internal function documentation
%   07/29/20    15:20   Thad Haines     jay -> 1j

global g

flag =1;
kdc=10*(k-1)+1;

if isempty(g.ind.n_mot)
    g.ind.n_mot = 0;
end

if isempty(g.igen.n_ig)
    g.igen.n_ig =0;
end


psi = g.mac.psi_re(:,k) + 1j*g.mac.psi_im(:,k);
vmp = g.ind.vdp(:,k) + 1j*g.ind.vqp(:,k);
vmpig = g.igen.vdpig(:,k) + 1j*g.igen.vqpig(:,k);

if (g.ind.n_mot~=0 && g.igen.n_ig==0) 
    ntot = g.mac.n_mac + g.ind.n_mot;
    ngm = g.mac.n_mac + g.ind.n_mot;
    int_volt=[psi; vmp]; % internal voltages of generators and motors
elseif (g.ind.n_mot==0 && g.igen.n_ig~=0)
    ntot = g.mac.n_mac + g.igen.n_ig;
    ngm = g.mac.n_mac;
    int_volt=[psi; vmpig]; % internal voltages of generators and ind. generators
elseif  (g.ind.n_mot~=0 && g.igen.n_ig~=0)
    ntot = g.mac.n_mac + g.ind.n_mot + g.igen.n_ig;
    ngm = g.mac.n_mac + g.ind.n_mot;
    int_volt=[psi; vmp; vmpig]; % internal voltages of generators, motors & ind. generators
else
    int_volt = psi;
end
h_sol = h(ks);
nbus = length(bus_sim(:,1));
cur = Y_g*int_volt; % network solution currents into generators
b_v(bo(g.ncl.nload+1:nbus),1) = rec_V1*int_volt;   % bus voltage reconstruction

if g.ncl.nload~=0
    if k~=1
        kk = k-1;
    else
        kk=k;
    end
    vnc = g.bus.bus_v(g.bus.bus_int(g.ncl.load_con(:,1)),kk);% initial value
    vnc = nc_load(g.bus.bus_sim,flag,Y_nc,Y_ncg,int_volt,vnc,1e-5,k,kdc);
    
    % set nc load voltages
    b_v(bo(1:g.ncl.nload),1)=vnc;
    b_v(bo(g.ncl.nload+1:nbus),1) = b_v(bo(g.ncl.nload+1:nbus),1)+rec_V2*vnc;
    cur = cur + Y_gnc*vnc;% modify generator currents for nc loads
end

% note: the dc bus voltages are the equivalent HT bus voltages
%       and not the LT bus voltages
g.bus.bus_v(g.bus.bus_int(bus_sim(:,1)),k) = b_v;
g.bus.theta(g.bus.bus_int(bus_sim(:,1)),k) = angle(b_v);
g.mac.cur_re(:,k) = real(cur(1:g.mac.n_mac));
g.mac.cur_im(:,k) = imag(cur(1:g.mac.n_mac)); % generator currents

if g.ind.n_mot~=0
    g.ind.idmot(:,k) = -real(cur(g.mac.n_mac+1:ngm));%induction motor currents
    g.ind.iqmot(:,k) = -imag(cur(g.mac.n_mac+1:ngm));%current out of network
    g.ind.s_mot(:,k) = g.bus.bus_v(g.bus.bus_int(g.ind.ind_con(:,2)),k).*(g.ind.idmot(:,k)-1j*g.ind.iqmot(:,k));
    g.ind.p_mot(:,k) = real(g.ind.s_mot(:,k));
    g.ind.q_mot(:,k) = imag(g.ind.s_mot(:,k));
end
if g.igen.n_ig~=0
    g.igen.idig(:,k) = -real(cur(ngm+1:ntot));%induction generator currents
    g.igen.iqig(:,k) = -imag(cur(ngm+1:ntot));%current out of network
    g.igen.s_igen(:,k) = g.bus.bus_v(g.bus.bus_int(g.igen.igen_con(:,2)),k)...
        .*(g.igen.idig(:,k)-1j*g.igen.iqig(:,k));
    g.igen.pig(:,k) = real(g.igen.s_igen(:,k));
    g.igen.qig(:,k) = imag(g.igen.s_igen(:,k));
end

if g.dc.n_conv ~=0
    % calculate dc voltage and current
    V0(g.dc.r_idx,1) = abs(g.bus.bus_v(g.dc.rec_ac_bus,k)).*g.dc.dcc_pot(:,7);
    V0(g.dc.i_idx,1) = abs(g.bus.bus_v(g.dc.inv_ac_bus,k)).*g.dc.dcc_pot(:,8);
    g.dc.Vdc(g.dc.r_idx,kdc) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,kdc)) - g.dc.i_dcr(:,kdc).*g.dc.dcc_pot(:,3);
    g.dc.Vdc(g.dc.i_idx,kdc) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,kdc)) - g.dc.i_dci(:,kdc).*g.dc.dcc_pot(:,5);
    g.dc.i_dc(g.dc.r_idx,kdc) = g.dc.i_dcr(:,kdc);
    g.dc.i_dc(g.dc.i_idx,kdc) = g.dc.i_dci(:,kdc);
end
