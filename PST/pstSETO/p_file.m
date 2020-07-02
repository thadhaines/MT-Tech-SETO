% m.file  for  computing perturbations  
% 5:53 pm 29/12/98
% for svm_mgen.m  
% forms state space model of system
% Author Graham Rogers
% (c) Copyright Joe Chow/ Cherry Tree Scientific Software  1991-1997
% Added code for pwrmod, D. Trudnowski, 2015
% All Rights Reserved
% step 3a: network solution

%% Original globals block condensed into g...

flag = 1;
%generators
mac_ib(0,2,bus,flag);
mac_em(0,2,bus,flag);
mac_tra(0,2,bus,flag);
mac_sub(0,2,bus,flag);
mac_ind(0,2,bus,flag); 
mac_igen(0,2,bus,flag);
psi = g.sys.psi_re(:,2) + jay*g.sys.psi_im(:,2);

if n_mot~=0&&n_ig==0
   vmp = vdp(:,2) + jay*vqp(:,2);
   int_volt=[psi; vmp]; % internal voltages of generators and motors 
elseif n_mot==0&&n_ig~=0
   vmpig = vdpig(:,2) + jay*vqpig(:,2);
   int_volt=[psi; vmpig]; % internal voltages of sync and ind. generators  
elseif n_mot~=0&&n_ig~=0
   vmp = vdp(:,2) + jay*vqp(:,2);
   vmpig = vdpig(:,2) + jay*vqpig(:,2);
   int_volt = [psi;vmp;vmpig];
else
   int_volt = psi;
end
if n_conv~=0
   dc_cont(0,2,2,bus,flag);
end

cur(:,2) = Y_gprf*int_volt; % network solution currents into generators       
b_v(boprf(g.ncl.nload+1:nbus),1) = V_rgprf*int_volt;   % bus voltage reconstruction
if g.ncl.nload~=0
   vnc = v(boprf(1:g.ncl.nload),1);
   vnc = nc_load(bus,flag,Y_ncprf,Y_ncgprf,int_volt,vnc,1e-6,2,2);
   bvnc = full(V_rncprf*vnc);
   b_v(boprf(1:g.ncl.nload),1) = vnc;
   cur(:,2) = cur(:,2) + Y_gncprf*vnc;% modify currents for nc loads
   b_v(boprf(g.ncl.nload+1:nbus),1) =  b_v(boprf(g.ncl.nload+1:nbus),1) + bvnc; % adjust voltages for nc loads
end
v(g.sys.bus_int(bus(:,1)),2) = b_v;
g.sys.bus_v(g.sys.bus_int(bus(:,1)),2) = b_v;
g.sys.theta(g.sys.bus_int(bus(:,1)),2) = angle(b_v); 
g.sys.cur_re(1:g.mac.n_mac,2) = real(cur(1:g.mac.n_mac,2)); 
g.sys.cur_im(1:g.mac.n_mac,2) = imag(cur(1:g.mac.n_mac,2));
cur_mag(1:g.mac.n_mac,2) = abs(cur(1:g.mac.n_mac,2)).*g.mac.mac_pot(:,1);

if n_mot~=0
   idmot = -real(cur(g.mac.n_mac+1:ngm,:));%induction motor currents
   iqmot = -imag(cur(g.mac.n_mac+1:ngm,:));%current out of network
end
if n_ig~=0
   idig = -real(cur(ngm+1:ntot,:));%induction generator currents
   iqig = -imag(cur(ngm+1:ntot,:));%current out of network
end

if n_conv ~=0
   % calculate dc voltage and current
   V0(r_idx,1) = abs(v(rec_ac_bus,2)).*dcc_pot(:,7);
   V0(i_idx,1) = abs(v(inv_ac_bus,2)).*dcc_pot(:,8);
   Vdc(r_idx,2) = V0(r_idx,1).*cos(alpha(:,2)) - i_dcr(:,2).*dcc_pot(:,3);
   Vdc(i_idx,2) = V0(i_idx,1).*cos(gamma(:,2)) - i_dci(:,2).*dcc_pot(:,5);
   i_dc(r_idx,2) = i_dcr(:,2);
   i_dc(i_idx,2) = i_dci(:,2);
end
% DeltaP/omega filter
dpwf(0,2,bus,flag);
% pss
pss(0,2,flag);
% exciters
smpexc(0,2,flag);
smppi(0,2,flag);
exc_dc12(0,2,flag);
exc_st3(0,2,flag);
% turbine/governor
tg(0,2,flag);
tg_hydro(0,2,bus,flag);
% calculate rates of change
flag = 2;
mac_em(0,2,bus,flag);
mac_tra(0,2,bus,flag);
mac_sub(0,2,bus,flag); 
mac_ind(0,2,bus,flag); 
mac_igen(0,2,bus,flag);
dpwf(0,2,bus,flag);
pss(0,2,flag);

smpexc(0,2,flag);
smppi(0,2,flag);
exc_dc12(0,2,flag);
exc_st3(0,2,flag);

tg(0,2,flag);
tg_hydro(0,2,bus,flag);

if n_svc~=0 
   v_svc = abs(v(g.sys.bus_int(svc_con(:,2)),2));
   svc(0,2,bus,flag,v_svc);
end
if n_tcsc~=0
   tcsc(0,2,bus,flag);
end
if g.lmod.n_lmod~=0 
   lmod(0,2,flag);
end
if g.rlmod.n_rlmod~=0 
   rlmod(0,2,flag);
end
if g.pwr.n_pwrmod~=0 
   pwrmod_p(0,2,bus,flag);
   pwrmod_q(0,2,bus,flag);
end

if n_conv ~=0
   dc_cont(0,2,2,bus,flag);
   dc_line(0,2,2,bus,flag);
end
telect(:,2) = g.mac.pelect(:,2).*g.mac.mac_pot(:,1) ...
    +g.mac.mac_con(:,5).*cur_mag(:,2).*cur_mag(:,2);
% form state matrix
% form vector of d states
d_vector = zeros(max_state,1);
mac_state = 6*g.mac.n_mac;
exc_state = mac_state+5*g.exc.n_exc;
pss_state = exc_state + 3*n_pss;
dpw_state = pss_state +6*n_dpw;
d_vector(1:g.mac.n_mac) = g.mac.dmac_ang(:,2);
d_vector(g.mac.n_mac+1:2*g.mac.n_mac) = g.mac.dmac_spd(:,2);
d_vector(2*g.mac.n_mac+1:3*g.mac.n_mac) = g.mac.deqprime(:,2);
d_vector(3*g.mac.n_mac+1:4*g.mac.n_mac) = g.mac.dpsikd(:,2);
d_vector(4*g.mac.n_mac+1:5*g.mac.n_mac) = g.mac.dedprime(:,2);
d_vector(5*g.mac.n_mac+1:6*g.mac.n_mac) = g.mac.dpsikq(:,2);

if g.exc.n_exc~=0
    nEXC = g.exc.n_exc;
   d_vector(mac_state+ 1 : mac_state+ nEXC ) = g.exc.dV_TR(:,2);
   d_vector(mac_state+ 1 + nEXC : mac_state+ 2*nEXC) = g.exc.dV_As(:,2);
   d_vector(mac_state+ 1 + 2*nEXC :mac_state+ 3*nEXC) = g.exc.dV_R(:,2);
   d_vector(mac_state+ 1 + 3*nEXC :mac_state+ 4*nEXC) = g.exc.dEfd(:,2);
   d_vector(mac_state+ 1 + 4*nEXC :mac_state+ 5*nEXC) = g.exc.dR_f(:,2);
end

if n_pss~=0
   d_vector(exc_state+1:exc_state+n_pss) = dpss1(:,2);
   d_vector(exc_state+n_pss+1:exc_state+2*n_pss) = dpss2(:,2);
   d_vector(exc_state+2*n_pss+1:exc_state+3*n_pss) = dpss3(:,2);
end
if n_dpw~=0
   d_vector(pss_state+1:pss_state+n_dpw) = dsdpw1(:,2);
   d_vector(pss_state+n_dpw+1:pss_state+2*n_dpw) = dsdpw2(:,2);
   d_vector(pss_state+2*n_dpw+1:pss_state+3*n_dpw) = dsdpw3(:,2);
   d_vector(pss_state+3*n_dpw+1:pss_state+4*n_dpw) = dsdpw4(:,2);
   d_vector(pss_state+4*n_dpw+1:pss_state+5*n_dpw) = dsdpw5(:,2);
   d_vector(pss_state+5*n_dpw+1:pss_state+6*n_dpw) = dsdpw6(:,2);
end

if g.tg.n_tg~=0 || g.tg.n_tgh~=0
   ngt = g.tg.n_tg + g.tg.n_tgh;
   d_vector(dpw_state+1:dpw_state+ngt) = g.tg.dtg1(:,2);
   d_vector(dpw_state+ngt+1:dpw_state+2*ngt) = g.tg.dtg2(:,2);
   d_vector(dpw_state+2*ngt+1:dpw_state+3*ngt) = g.tg.dtg3(:,2);
   d_vector(dpw_state+3*ngt+1:dpw_state+4*ngt) = g.tg.dtg4(:,2);
   d_vector(dpw_state+4*ngt+1:dpw_state+5*ngt) = g.tg.dtg5(:,2);
end
if n_mot~=0
   mot_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh); % times turbine number? 06/11/20 -thad
   d_vector(mot_start+1:mot_start+n_mot) = dvdp(:,2);
   d_vector(mot_start+n_mot+1:mot_start+2*n_mot) = dvqp(:,2);
   d_vector(mot_start+2*n_mot+1:mot_start+3*n_mot) = dslip(:,2);
end
if n_ig~=0
   ig_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot;
   d_vector(ig_start+1:ig_start+n_ig) = dvdpig(:,2);
   d_vector(ig_start+n_ig+1:ig_start+2*n_ig) = dvqpig(:,2);
   d_vector(ig_start+2*n_ig+1:ig_start+3*n_ig) = dslig(:,2);
end

if n_svc ~= 0
   svc_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig;
   d_vector(svc_start+1:svc_start+n_svc) = dB_cv(:,2);
   d_vector(svc_start+n_svc+1:svc_start+2*n_svc) = dB_con(:,2);
end

if n_tcsc~=0
   tcsc_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig+2*n_svc;
   d_vector(tcsc_start+1:tcsc_start+n_tcsc)=dB_tcsc(:,2);
end
if g.lmod.n_lmod ~= 0
   lmod_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig+2*n_svc+n_tcsc;
   d_vector(lmod_start+1:lmod_start+g.lmod.n_lmod) = g.lmod.dlmod_st(:,2);
end
if g.rlmod.n_rlmod ~= 0
   rlmod_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig+2*n_svc+n_tcsc+g.lmod.n_lmod;
   d_vector(rlmod_start+1:rlmod_start+g.rlmod.n_rlmod) = g.rlmod.drlmod_st(:,2);
end
if g.pwr.n_pwrmod ~= 0
   pwrmod_p_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig+2*n_svc+n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod;
   d_vector(pwrmod_p_start+1:pwrmod_p_start+g.pwr.n_pwrmod) = g.pwr.dpwrmod_p_st(:,2);
   pwrmod_q_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig+2*n_svc+n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod+g.pwr.n_pwrmod;
   d_vector(pwrmod_q_start+1:pwrmod_q_start+g.pwr.n_pwrmod) = g.pwr.dpwrmod_q_st(:,2);
end

if n_conv~=0
   dc_start = dpw_state+5*(g.tg.n_tg + g.tg.n_tgh)+3*n_mot+3*n_ig + 2*n_svc +n_tcsc+ g.lmod.n_lmod+g.rlmod.n_rlmod+2*g.pwr.n_pwrmod;
   d_vector(dc_start+1: dc_start+n_dcl) = dv_conr(:,2);
   d_vector(dc_start+n_dcl+1: dc_start+2*n_dcl) = dv_coni(:,2);
   d_vector(dc_start+2*n_dcl+1: dc_start+3*n_dcl) = di_dcr(:,2);
   d_vector(dc_start+3*n_dcl+1: dc_start+4*n_dcl) = di_dci(:,2);
   d_vector(dc_start+4*n_dcl+1: dc_start+5*n_dcl) = dv_dcc(:,2);
end 

% form state matrix
if c_state == 0
   if k==1
      j_state = j;
   else
      j_state = j + sum(state(1:k-1));
   end
   
   if g.mac.n_ib~=0
      k_nib_idx = find(g.mac.not_ib_idx==k);
   else
      k_nib_idx = k;
   end
   
   if j == 2;  
      if ~isempty(k_nib_idx)
         c_spd(k_nib_idx,j_state) = 1;
      end
   end
   
   a_mat(:,j_state) = p_mat*d_vector/pert;
   
   % form output matrices 
   c_p(g.mac.not_ib_idx,j_state) = (g.mac.pelect(g.mac.not_ib_idx,2)-g.mac.pelect(g.mac.not_ib_idx,1))...
      .*g.mac.mac_pot(g.mac.not_ib_idx,1)/pert;
   c_t(g.mac.not_ib_idx,j_state) = (telect(g.mac.not_ib_idx,2)-telect(g.mac.not_ib_idx,1))/pert;
   c_pm(g.mac.not_ib_idx,j_state) = (g.mac.pmech(g.mac.not_ib_idx,2)-g.mac.pmech(g.mac.not_ib_idx,1))/pert;
   c_v(:,j_state) = (abs(v(:,2)) -abs(v(:,1)))/pert;
   c_ang(:,j_state) = (g.sys.theta(:,2) - g.sys.theta(:,1))/pert;
   c_curd(:,j_state) = (g.mac.curd(:,2) - g.mac.curd(:,1))/pert;  % JHC 12/17/2015
   c_curq(:,j_state) = (g.mac.curq(:,2) - g.mac.curq(:,1))/pert;  % JHC 12/17/2015
   if g.exc.n_exc~=0
      c_Efd(:,j_state) = (g.exc.Efd(:,2)-g.exc.Efd(:,1))/pert;
   end
   if ~isempty(g.sys.lmon_con) 
      from_idx = g.sys.bus_int(line(g.sys.lmon_con,1));
      to_idx = g.sys.bus_int(line(g.sys.lmon_con,2));
      V1 = v(from_idx,1);
      V2 = v(to_idx,1);
      [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
      [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
      V1 = v(from_idx,2);
      V2 = v(to_idx,2);
      [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
      [l_if2,l_it2]=line_cur(V1,V2,R,X,B,tap,phi);
      c_pf1(:,j_state) = (real(s12-s11))/pert; 
      c_qf1(:,j_state) = (imag(s12-s11))/pert;
      c_pf2(:,j_state) = (real(s22-s21))/pert;
      c_qf2(:,j_state) = (imag(s22-s21))/pert;
      c_ilmf(:,j_state) = (abs(l_if2)-abs(l_if1))/pert;
      c_ilmt(:,j_state) = (abs(l_it2)-abs(l_it1))/pert;
      c_ilrf(:,j_state) = real(l_if2-l_if1)/pert;
      c_ilif(:,j_state) = imag(l_if2-l_if1)/pert;
      c_ilrt(:,j_state) = real(l_it2-l_it1)/pert;
      c_ilit(:,j_state) = imag(l_it2-l_it1)/pert;
   end
   if n_conv~=0
      c_dcir(:,j_state) = (i_dcr(:,2)-i_dcr(:,1))/pert;
      c_dcii(:,j_state) = (i_dci(:,2)-i_dci(:,1))/pert;
      c_Vdcr(:,j_state) = (Vdc(r_idx,2)-Vdc(r_idx,1))/pert;
      c_Vdci(:,j_state) = (Vdc(i_idx,2)-Vdc(i_idx,1))/pert;
   end
else
   % form b and d matrices
   if c_state == 1
      b_vr(:,vr_input) = p_mat*d_vector/pert;
      d_pvr(:,vr_input) = (g.mac.pelect(:,2)-g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_vvr(:,vr_input) = abs(v(:,2) - v(:,1))/pert;
      d_angvr(:,vr_input) = (g.sys.theta(:,2)-g.sys.theta(:,1))/pert;
   elseif c_state==2
      b_pr(:,pr_input) = p_mat*d_vector/pert;
      d_ppr(:,pr_input) = (g.mac.pelect(:,2) - g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_vpr(:,pr_input) = abs(v(:,2) - v(:,1))/pert;
      d_angpr(:,pr_input) = (g.sys.theta(:,2)-g.sys.theta(:,1))/pert; 
   elseif c_state==3
      b_svc(:,svc_input) = p_mat*d_vector/pert;
      % note: d_svc is zero because of the time constant
   elseif c_state==4
      b_tcsc(:,tcsc_input)=p_mat*d_vector/pert;
   elseif c_state == 5
      b_lmod(:,lmod_input) = p_mat*d_vector/pert;
      % note: d_lmod is zero because of the time constant
   elseif c_state == 6
      b_rlmod(:,rlmod_input) = p_mat*d_vector/pert;
      % note: d_lmod is zero because of the time constant
   elseif c_state == 7
      b_pwrmod_p(:,pwrmod_p_input) = p_mat*d_vector/pert;
      % note: d_pwrmod is zero because of the time constant
   elseif c_state == 8
      b_pwrmod_q(:,pwrmod_q_input) = p_mat*d_vector/pert;
      % note: d_pwrmod is zero because of the time constant
   elseif c_state == 9
      b_dcr(:,dcmod_input) = p_mat*d_vector/pert;
      d_pdcr(:,dcmod_input) = (g.mac.pelect(:,2)-g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_vdcr(:,dcmod_input) = abs(v(:,2) - v(:,1))/pert;
      d_angdcr(:,dcmod_input) = (g.sys.theta(:,2)-g.sys.theta(:,1))/pert;
      d_pdcr(:,dcmod_input)=(g.mac.pelect(:,2) - g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_idcdcr(:,dcmod_input) = (i_dcr(:,2)-i_dcr(:,1))/pert;
      d_Vdcrdcr(:,dcmod_input) = (Vdc(r_idx,2)-Vdc(r_idx,1))/pert;
      d_Vdcidcr(:,dcmod_input) = (Vdc(i_idx,2)-Vdc(i_idx,1))/pert;
      if ~isempty(g.sys.lmon_con) 
         from_idx = g.sys.bus_int(line(g.sys.lmon_con,1));
         to_idx = g.sys.bus_int(line(g.sys.lmon_con,2));
         V1 = v(from_idx,1);
         V2 = v(to_idx,1);
         [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
         [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
         V1 = v(from_idx,2);
         V2 = v(to_idx,2);
         [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
         [l_if2,l_it2]=line_cur(V1,V2,R,X,B,tap,phi);
         d_pf1cdr(:,dcmod_input) = (real(s12-s11))/pert; 
         d_qf1dcr(:,dcmod_input) = (imag(s12-s11))/pert;
         d_pf2dcr(:,dcmod_input) = (real(s22-s21))/pert;
         d_qf2dcr(:,dcmod_input) = (imag(s22-s21))/pert;
         d_ilmfdcr(:,dcmod_input) = (abs(l_if2)-abs(l_if1))/pert;
         d_ilmtdcr(:,dcmod_input) = (abs(l_it2)-abs(l_it1))/pert;
         d_ilrfdcr(:,dcmod_input) = real(l_if2-l_if1)/pert;
         d_ilifdcr(:,dcmod_input) = imag(l_if2-l_if1)/pert;
         d_ilrtdcr(:,dcmod_input) = real(l_it2-l_it1)/pert;
         d_ilitdcr(:,dcmod_input) = imag(l_it2-l_it1)/pert;
      end      
   elseif c_state == 10
      b_dci(:,dcmod_input) = p_mat*d_vector/pert;
      d_pdci(:,dcmod_input) = (g.mac.pelect(:,2)-g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_vdci(:,dcmod_input) = abs(v(:,2) - v(:,1))/pert;
      d_angdci(:,dcmod_input) = (g.sys.theta(:,2)-g.sys.theta(:,1))/pert;
      d_pdci(:,dcmod_input)=(g.mac.pelect(:,2) - g.mac.pelect(:,1)).*g.mac.mac_pot(:,1)/pert;
      d_idcdci(:,dcmod_input) = (i_dci(:,2)-i_dci(:,1))/pert;
      d_Vdcrdci(:,dcmod_input) = (Vdc(r_idx,2)-Vdc(r_idx,1))/pert;
      d_Vdcdci(:,dcmod_input) = (Vdc(i_idx,2)-Vdc(i_idx,1))/pert;
      if ~isempty(g.sys.lmon_con) 
         from_idx = g.sys.bus_int(line(g.sys.lmon_con,1));
         to_idx = g.sys.bus_int(line(g.sys.lmon_con,2));
         V1 = v(from_idx,1);
         V2 = v(to_idx,1);
         [s11,s21] = line_pq(V1,V2,R,X,B,tap,phi);
         [l_if1,l_it1] = line_cur(V1,V2,R,X,B,tap,phi);
         V1 = v(from_idx,2);
         V2 = v(to_idx,2);
         [s12,s22] = line_pq(V1,V2,R,X,B,tap,phi);
         [l_if2,l_it2]=line_cur(V1,V2,R,X,B,tap,phi);
         d_pf1cdi(:,dcmod_input) = (real(s12-s11))/pert; 
         d_qf1dci(:,dcmod_input) = (imag(s12-s11))/pert;
         d_pf2dci(:,dcmod_input) = (real(s22-s21))/pert;
         d_qf2dci(:,dcmod_input) = (imag(s22-s21))/pert;
         d_ilmfdci(:,dcmod_input) = (abs(l_if2)-abs(l_if1))/pert;
         d_ilmtdci(:,dcmod_input) = (abs(l_it2)-abs(l_it1))/pert;
         d_ilrfdci(:,dcmod_input) = real(l_if2-l_if1)/pert;
         d_ilifdci(:,dcmod_input) = imag(l_if2-l_if1)/pert;
         d_ilrtdci(:,dcmod_input) = real(l_it2-l_it1)/pert;
         d_ilitdci(:,dcmod_input) = imag(l_it2-l_it1)/pert;
      end  
   end
end
%reset states to initial values
g.sys.psi_re(:,2) = g.sys.psi_re(:,1);
g.sys.psi_im(:,2) = g.sys.psi_im(:,1);
v(:,2) = v(:,1);
g.sys.bus_v(:,2)=g.sys.bus_v(:,1);
telect(:,2) = telect(:,1); % what is this? -thad 06/18/20

g.mac.eterm(:,2) = g.mac.eterm(:,1);
g.mac.pelect(:,2) = g.mac.pelect(:,1);
g.mac.qelect(:,2) = g.mac.qelect(:,1);

g.sys.theta(:,2) = g.sys.theta(:,1);

g.mac.pmech(:,2) = g.mac.pmech(:,1);
g.mac.mac_ang(:,2) = g.mac.mac_ang(:,1);
g.mac.dmac_ang(:,2) = g.mac.dmac_ang(:,1);
g.mac.mac_spd(:,2) = g.mac.mac_spd(:,1);
g.mac.dmac_spd(:,2) = g.mac.dmac_spd(:,1);
g.mac.eqprime(:,2) = g.mac.eqprime(:,1);
g.mac.deqprime(:,2) = g.mac.deqprime(:,1);
g.mac.psikd(:,2) = g.mac.psikd(:,1);
g.mac.dpsikd(:,2) = g.mac.dpsikd(:,1);
g.mac.edprime(:,2) = g.mac.edprime(:,1);
g.mac.dedprime(:,2) = g.mac.dedprime(:,1);
g.mac.psikq(:,2)= g.mac.psikq(:,1);

if g.exc.n_exc ~= 0
   g.exc.V_TR(:,2)=g.exc.V_TR(:,1);
   g.exc.dV_TR(:,2)=g.exc.dV_TR(:,1);
   g.exc.V_As(:,2) = g.exc.V_As(:,1);
   g.exc.dV_As(:,2) = g.exc.dV_As(:,1);
   g.exc.V_A(:,2) = g.exc.V_A(:,1);
   g.exc.dV_R(:,2) = g.exc.dV_R(:,1);
   g.exc.V_R(:,2)=g.exc.V_R(:,1);
   g.exc.Efd(:,2)=g.exc.Efd(:,1);
   g.exc.dEfd(:,2) = g.exc.dEfd(:,1);
   g.exc.R_f(:,2)=g.exc.R_f(:,1);
   g.exc.dR_f(:,2) = g.exc.dR_f(:,1);
end
if n_pss~=0
   pss1(:,2)=pss1(:,1);
   pss2(:,2)=pss2(:,1);
   pss3(:,2)=pss3(:,1);
   dpss1(:,2)=dpss1(:,1);
   dpss2(:,2)=dpss2(:,1);
   dpss3(:,2)=dpss3(:,1);
end
if n_dpw~=0
   sdpw1(:,2)=sdpw1(:,1);
   sdpw2(:,2)=sdpw2(:,1);
   sdpw3(:,2)=sdpw3(:,1);
   sdpw4(:,2)=sdpw4(:,1);
   sdpw5(:,2)=sdpw5(:,1);
   sdpw6(:,2)=sdpw6(:,1);   
   dsdpw1(:,2)=dsdpw1(:,1);
   dsdpw2(:,2)=dsdpw2(:,1);
   dsdpw3(:,2)=dsdpw3(:,1);
   dsdpw4(:,2)=dsdpw4(:,1);
   dsdpw5(:,2)=dsdpw5(:,1);
   dsdpw6(:,2)=dsdpw6(:,1);

end

if g.tg.n_tg ~=0 || g.tg.n_tgh ~= 0
   g.tg.tg1(:,2)= g.tg.tg1(:,1);
   g.tg.tg2(:,2)= g.tg.tg2(:,1);
   g.tg.tg3(:,2)= g.tg.tg3(:,1);
   g.tg.tg4(:,2)= g.tg.tg4(:,1);
   g.tg.tg5(:,2)= g.tg.tg5(:,1);
   g.tg.dtg1(:,2)= g.tg.dtg1(:,1);
   g.tg.dtg2(:,2)= g.tg.dtg2(:,1);
   g.tg.dtg3(:,2)= g.tg.dtg3(:,1);
   g.tg.dtg4(:,2)= g.tg.dtg4(:,1);
   g.tg.dtg5(:,2)= g.tg.dtg5(:,1);
end
if n_mot~=0
   vdp(:,2) = vdp(:,1);
   vqp(:,2) = vqp(:,1);
   slip(:,2) = slip(:,1);
   dvdp(:,2) = dvdp(:,1);
   dvqp(:,2) = dvqp(:,1);
   dslip(:,2) = dslip(:,1);
end
if n_ig~=0
   vdpig(:,2) = vdpig(:,1);
   vqpig(:,2) = vqpig(:,1);
   slig(:,2) = slig(:,1);
   dvdpig(:,2) = dvdpig(:,1);
   dvqpig(:,2) = dvqpig(:,1);
   dslig(:,2) = dslig(:,1);
end
if n_svc ~=0
   B_cv(:,2) = B_cv(:,1);
   dB_cv(:,2) = dB_cv(:,1);
   B_con(:,2) = B_con(:,1);
   dB_con(:,2) = dB_con(:,1);
   svc_sig(:,2) = svc_sig(:,1);
end
if n_tcsc~=0
   B_tcsc(:,2)=B_tcsc(:,1);
   dB_tcsc(:,2)=dB_tcsc(:,1);
   tcsc_sig(:,2)=tcsc_sig(:,1);
end
if g.lmod.n_lmod ~=0
   g.lmod.lmod_st(:,2) = g.lmod.lmod_st(:,1);
   g.lmod.dlmod_st(:,2) = g.lmod.dlmod_st(:,1);
   g.lmod.lmod_sig(:,2) = g.lmod.lmod_sig(:,1);
end
if g.rlmod.n_rlmod ~=0
   g.rlmod.rlmod_st(:,2) = g.rlmod.rlmod_st(:,1);
   g.rlmod.drlmod_st(:,2) = g.rlmod.drlmod_st(:,1);
   g.rlmod.rlmod_sig(:,2) = g.rlmod.rlmod_sig(:,1);
end
if g.pwr.n_pwrmod ~=0
   g.pwr.pwrmod_p_st(:,2) = g.pwr.pwrmod_p_st(:,1);
   g.pwr.dpwrmod_p_st(:,2) = g.pwr.dpwrmod_p_st(:,1);
   g.pwr.pwrmod_p_sig(:,2) = g.pwr.pwrmod_p_sig(:,1);
   g.pwr.pwrmod_q_st(:,2) = g.pwr.pwrmod_q_st(:,1);
   g.pwr.dpwrmod_q_st(:,2) = g.pwr.dpwrmod_q_st(:,1);
   g.pwr.pwrmod_q_sig(:,2) = g.pwr.pwrmod_q_sig(:,1);
end

if n_conv~=0
   v_conr(:,2) = v_conr(:,1);
   v_coni(:,2) = v_coni(:,1);
   i_dcr(:,2) = i_dcr(:,1);
   i_dci(:,2) = i_dci(:,1);
   v_dcc(:,2) = v_dcc(:,1);
   dv_conr(:,2) = dv_conr(:,1);
   dv_coni(:,2) = dv_coni(:,1);
   di_dcr(:,2) = di_dcr(:,1);
   di_dci(:,2) = di_dci(:,1);
   dv_dcc(:,2) = dv_dcc(:,1);
   Vdc(:,2) = Vdc(:,1);
   i_dc(:,2) = i_dc(:,1);
   alpha(:,2) = alpha(:,1);
   gamma(:,2) = gamma(:,1); 
   dc_sig(:,2)=dc_sig(:,1);
end

