function h_sol = i_simu(k,ks,k_inc,h,bus_sim,Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
%Syntax: h_sol = i_simu(k,ks,k_inc,h,bus_sim,...
%                 Y_g,Y_gnc,Y_ncg,Y_nc,rec_V1,rec_V2,bo)
% 11:13 AM 18/08/97
%Purpose forms the network interface variables
% Inputs: k - the current time step
%         ks - indicates the switching times
%         k_inc - the number of time seps between switching points
%         h vector of time steps
%         bus_sim value of bus matrix at this switching time
%         Y_g - reduced Y matrix for generators
%         Y_gnc - mutual reduced Y generators-nc loads
%         Y_ncg - mutual reduced Y nc loads generators
%         Y_nc - reduced Y matrix nc loads
%         rec_V1 - voltage recovery matrix generators
%         rec_V2 - voltage recovery matrix nc loads
%         bo bus order for this switching time
% Output: h_sol - the time step at this value of ks
% Called by: s_simu
% Version 1.1
% Date:   August 1997
% Modification: add induction generator
% Version: 1.0
% Date:    March 1997
% Author:   Graham Rogers
% Copyright (c) Joe Chow All Rights Reserved


%global load_con nload

global vdp vqp n_mot idmot iqmot p_mot q_mot ind_con

global vdpig vqpig n_ig idig iqig  pig qig igen_con

global i_dc Vdc alpha gamma dcc_pot i_dcr  i_dci Pdc
global r_idx i_idx ac_bus rec_ac_bus inv_ac_bus n_conv

    
global g

flag =1;
kdc=10*(k-1)+1;
if isempty(n_mot); n_mot = 0;end
if isempty(n_ig); n_ig =0;end
jay = sqrt(-1);
psi = g.sys.psi_re(:,k) + jay*g.sys.psi_im(:,k);
vmp = vdp(:,k) + jay*vqp(:,k);
vmpig = vdpig(:,k) + jay*vqpig(:,k);
if (n_mot~=0&n_ig==0)
    ntot = g.mac.n_mac + n_mot;
    ngm = g.mac.n_mac + n_mot;
    int_volt=[psi; vmp]; % internal voltages of generators and motors
elseif (n_mot==0 & n_ig~=0)
    ntot = g.mac.n_mac + n_ig;
    ngm = g.mac.n_mac;
    int_volt=[psi; vmpig]; % internal voltages of generators and ind. generators
elseif  (n_mot~=0&n_ig~=0)
    ntot = g.mac.n_mac + n_mot + n_ig;
    ngm = g.mac.n_mac + n_mot;
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
    vnc = g.sys.bus_v(g.sys.bus_int(g.ncl.load_con(:,1)),kk);% initial value
    vnc = nc_load(bus_sim,flag,Y_nc,Y_ncg,int_volt,vnc,1e-5,k,kdc);
    
    % set nc load voltages
    b_v(bo(1:g.ncl.nload),1)=vnc;
    b_v(bo(g.ncl.nload+1:nbus),1) = b_v(bo(g.ncl.nload+1:nbus),1)+rec_V2*vnc;
    cur = cur + Y_gnc*vnc;% modify generator currents for nc loads
end
% note: the dc bus voltages are the equivalent HT bus voltages
%       and not the LT bus voltages
g.sys.bus_v(g.sys.bus_int(bus_sim(:,1)),k) = b_v;
g.sys.theta(g.sys.bus_int(bus_sim(:,1)),k) = angle(b_v);
g.sys.cur_re(:,k) = real(cur(1:g.mac.n_mac));
g.sys.cur_im(:,k) = imag(cur(1:g.mac.n_mac)); % generator currents
if n_mot~=0
    idmot(:,k) = -real(cur(g.mac.n_mac+1:ngm));%induction motor currents
    iqmot(:,k) = -imag(cur(g.mac.n_mac+1:ngm));%current out of network
    s_mot(:,k) = g.sys.bus_v(g.sys.bus_int(ind_con(:,2)),k).*(idmot(:,k)-jay*iqmot(:,k));
    p_mot(:,k) = real(s_mot(:,k));
    q_mot(:,k) = imag(s_mot(:,k));
end
if n_ig~=0
    idig(:,k) = -real(cur(ngm+1:ntot));%induction generator currents
    iqig(:,k) = -imag(cur(ngm+1:ntot));%current out of network
    s_igen(:,k) = g.sys.bus_v(g.sys.bus_int(igen_con(:,2)),k).*(idig(:,k)-jay*iqig(:,k));
    pig(:,k) = real(s_igen(:,k));
    qig(:,k) = imag(s_igen(:,k));
end

if n_conv ~=0
    % calculate dc voltage and current
    V0(r_idx,1) = abs(g.sys.bus_v(rec_ac_bus,k)).*dcc_pot(:,7);
    V0(i_idx,1) = abs(g.sys.bus_v(inv_ac_bus,k)).*dcc_pot(:,8);
    Vdc(r_idx,kdc) = V0(r_idx,1).*cos(alpha(:,kdc)) - i_dcr(:,kdc).*dcc_pot(:,3);
    Vdc(i_idx,kdc) = V0(i_idx,1).*cos(gamma(:,kdc)) - i_dci(:,kdc).*dcc_pot(:,5);
    i_dc(r_idx,kdc) = i_dcr(:,kdc);
    i_dc(i_idx,kdc) = i_dci(:,kdc);
end
