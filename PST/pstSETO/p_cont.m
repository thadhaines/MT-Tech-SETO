%perturbation control file  
% 5:40 pm 29/12/1996
% modified to include pwrmod, D. Trudnowski, 2015
% modified to include induction generators
% modified to include load modulation
% input disturbance modulation added to svc and lmod
% Author: Graham Rogers
% (c) Copyright Joe Chow/Cherry Tree Scientific Software 1991-1998
% All right reserved
% step 3: perform 0.1% perturbation on each state in turn


vr_input = 0;
pr_input = 0;
c_state = 0;
p_ratio = 1e-5;
for k = 1:g.mac.n_mac
   not_ib = 1;
   if ~isempty(g.mac.mac_ib_idx)
      ib_chk = find(g.mac.mac_ib_idx==k);
      if ~isempty(ib_chk);not_ib = 0;end
   end
   if not_ib==1
      j = 1;
      disp('disturb generator')
      pert = p_ratio*abs(g.mac.mac_ang(k,1));   
      pert = max(pert,p_ratio);
      g.mac.mac_ang(k,2) = g.mac.mac_ang(k,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 1;
      
      j = j+1;
      pert = p_ratio*abs(g.mac.mac_spd(k,1));   
      pert = max(pert,p_ratio);
      g.mac.mac_spd(k,2) = g.mac.mac_spd(k,1) + pert;
      p_file  % m file of perturbations
      st_name(k,j) = 2;
      k_tra = 0;
      k_sub = 0;
      if ~isempty(g.mac.mac_tra_idx)
         k_idx = find(g.mac.mac_tra_idx==k);
         if ~isempty(k_idx);k_tra = 1;end
      end
      if ~isempty(g.mac.mac_sub_idx)
         k_idx = find(g.mac.mac_sub_idx==k);
         if ~isempty(k_idx);k_sub=1;end
      end
      if k_tra==1 || k_sub==1
         j=j+1;
         pert = p_ratio*abs(g.mac.eqprime(k,1));   
         pert = max(pert,p_ratio);
         g.mac.eqprime(k,2) = g.mac.eqprime(k,1) + pert;
         p_file   % m file of perturbations
         st_name(k,j) = 3; 
      end
      if k_sub==1
         j=j+1;
         pert = p_ratio*abs(g.mac.psikd(k,1));   
         pert = max(pert,p_ratio);
         g.mac.psikd(k,2) = g.mac.psikd(k,1) + pert;
         p_file   % m file of perturbations
         st_name(k,j) = 4;
      end
      
      if k_tra==1 || k_sub==1
         j=j+1;
         pert = p_ratio*abs(g.mac.edprime(k,1));   
         pert = max(pert,p_ratio);
         g.mac.edprime(k,2) = g.mac.edprime(k,1) + pert;
         p_file   % m file of perturbations
         st_name(k,j) = 5;
      end 
      
      if k_sub==1 
         j=j+1;
         pert = p_ratio*abs(g.mac.psikq(k,1));   
         pert = max(pert,p_ratio);
         g.mac.psikq(k,2) = g.mac.psikq(k,1) + pert;
         p_file  % m file of perturbations
         st_name(k,j) = 6;
      end 
      
      % exciters
      if ~isempty(g.exc.exc_con)
         p_exc    
      end
      %pss
      if ~isempty(g.pss.pss_con)
         p_pss
      end
      if ~isempty(dpw_con)
         p_dpw
      end
      % turbine/governor
      if ~isempty(g.tg.tg_con)
         p_tg
      end
      
      % disturb the input variables
      if g.exc.n_exc ~= 0
         c_state = 1;
         exc_number = find(g.mac.mac_int(g.exc.exc_con(:,2)) ==k);
         if ~isempty(exc_number)
            disp('disturb V_ref')
            vr_input = vr_input + 1;  
            pert = p_ratio*abs(g.exc.exc_pot(exc_number,3));
            pert = max(pert,p_ratio);
            nominal = g.exc.exc_pot(exc_number,3);
            g.exc.exc_pot(exc_number,3) = g.exc.exc_pot(exc_number,3) + pert;
            p_file
            g.exc.exc_pot(exc_number,3) = nominal;
         end
      end
      
      if g.tg.n_tg ~=0 || g.tg.n_tgh~=0
         c_state = 2;
         tg_number = find(mac_tg == k);
         if isempty(tg_number)
            tg_number = find(mac_tgh==k);
         end
         if ~isempty(tg_number)
            disp('disturb P_ref')
            pr_input = pr_input + 1;
            pert = p_ratio*abs(g.tg.tg_pot(tg_number,5));
            pert = max(pert,p_ratio);
            nominal = g.tg.tg_pot(tg_number,5);
            g.tg.tg_pot(tg_number,5) = g.tg.tg_pot(tg_number,5) + pert;
            p_file
            g.tg.tg_pot(tg_number,5) = nominal;
         end
      end
      c_state = 0;           
   end
end

% disturb induction motor states
if n_mot~=0
   disp('disturbing induction motors')
   for k = g.mac.n_mac+1:ngm
      j=1;
      k_ind = k - g.mac.n_mac;
      pert = p_ratio*abs(vdp(k_ind,1));   
      pert = max(pert,p_ratio);
      vdp(k_ind,2) = vdp(k_ind,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 26;
      j=j+1;
      pert = p_ratio*abs(vqp(k_ind,1));   
      pert = max(pert,p_ratio);
      vqp(k_ind,2) = vqp(k_ind,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 27;
      j=j+1;
      pert = p_ratio*abs(slip(k_ind,1));   
      pert = max(pert,0.000001);
      slip(k_ind,2) = slip(k_ind,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 28;
   end
end

% disturb induction generator states
if n_ig~=0
   disp('disturbing induction generators')
   for k = ngm+1:ntot
      j=1;
      k_ig = k - ngm;
      pert = p_ratio*abs(vdpig(k_ig,1));   
      pert = max(pert,p_ratio);
      vdpig(k_ig,2) = vdpig(k_ig,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 29;
      j=j+1;
      pert = p_ratio*abs(vqpig(k_ig,1));   
      pert = max(pert,p_ratio);
      vqpig(k_ig,2) = vqpig(k_ig,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 30;
      j=j+1;
      pert = p_ratio*abs(slig(k_ig,1));   
      pert = max(pert,0.000001);
      slig(k_ig,2) = slig(k_ig,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 31;
   end
end

nts = ntot + g.svc.n_svc;
% disturb svc states
if g.svc.n_svc~=0
   disp('disturbing svc')
   for k = ntot+1:nts
      j=1;
      k_svc = k - ntot;
      pert = p_ratio*abs(g.svc.B_cv(k_svc,1));   
      pert = max(pert,p_ratio);
      g.svc.B_cv(k_svc,2) = g.svc.B_cv(k_svc,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 32;
      % disturb B_con
      if ~isempty(g.svc.svcll_idx)
         j = j+1;
         kcon = find(g.svc.svcll_idx==k_svc);
         if ~isempty(kcon)
            pert = p_ratio*abs(g.svc.B_con(k_svc,1));
            pert = max(pert,p_ratio);
            g.svc.B_con(k_svc,2) = g.svc.B_con(k_svc,1) + pert;
            p_file %m-file of perturbations
            st_name(k,j)= 33;
         end
      end
      % disturb the input variable
      disp('disturbing svc_sig') 
      c_state = 3; 
      svc_input = k_svc;
      pert = p_ratio;
      nominal = 0.0;
      g.svc.svc_sig(k_svc,2) = g.svc.svc_sig(k_svc,2) + pert;
      p_file
      g.svc.svc_sig(k_svc,2) = nominal; 
      c_state = 0;
   end
end
nts = ntot + g.svc.n_svc;

ntf = ntot + g.svc.n_svc + g.tcsc.n_tcsc;
% disturb tcsc states
if g.tcsc.n_tcsc~=0
   disp('disturbing tcsc')
   
   pert_scale = 1e2; % scaling factor for minimum perturbation
                     %  (need larger perturbation for B_TCSC)
                     % added by SGG, 03/30/11
   for k = nts+1:ntf
      j=1;
      k_tcsc = k - nts;
      pert = p_ratio*abs(g.tcsc.B_tcsc(k_tcsc,1));   
      pert = max(pert,p_ratio*pert_scale);  % modified by SGG, 03/30/11
      g.tcsc.B_tcsc(k_tcsc,2) = g.tcsc.B_tcsc(k_tcsc,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 34;
   end
   % disturb the input variable
   disp('disturbing tcsc_sig') 
   c_state = 4; 
   tcsc_input = k_tcsc;
   pert = p_ratio;
   nominal = 0.0;
   g.tcsc.tcsc_sig(k_tcsc,2) = g.tcsc.tcsc_sig(k_tcsc,2) + pert;
   p_file
   g.tcsc.tcsc_sig(k_tcsc,2) = nominal; 
   c_state = 0;
end

ntl = ntf + g.lmod.n_lmod;
% disturb lmod states
if g.lmod.n_lmod~=0
   disp('disturbing load modulation')
   for k = ntf+1:ntl
      j=1;
      k_lmod = k - ntf;
      pert = p_ratio*abs(g.lmod.lmod_st(k_lmod,1));   
      pert = max(pert,p_ratio);
      g.lmod.lmod_st(k_lmod,2) = g.lmod.lmod_st(k_lmod,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 35;
      % disturb the input variable
      disp('disturbing lmod_sig') 
      c_state = 5; 
      lmod_input = k_lmod;
      pert = p_ratio;
      nominal = 0.0;
      g.lmod.lmod_sig(k_lmod,2) = g.lmod.lmod_sig(k_lmod,2) + pert;
      p_file
      g.lmod.lmod_sig(k_lmod,2) = nominal;  
      c_state = 0;
   end
end
ntrl = ntl + g.rlmod.n_rlmod;
% disturb rlmod states
if g.rlmod.n_rlmod~=0
   disp('disturbing reactive load modulation')
   for k = ntl+1:ntrl
      j=1;
      k_rlmod = k - ntl;
      pert = p_ratio*abs(g.rlmod.rlmod_st(k_rlmod,1));   
      pert = max(pert,p_ratio);
      g.rlmod.rlmod_st(k_rlmod,2) = g.rlmod.rlmod_st(k_rlmod,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 36;
      % disturb the input variable
      disp('disturbing rlmod_sig') 
      c_state = 6; 
      rlmod_input = k_rlmod;
      pert = p_ratio;
      nominal = 0.0;
      g.rlmod.rlmod_sig(k_rlmod,2) = g.rlmod.rlmod_sig(k_rlmod,2) + pert;
      p_file
      g.rlmod.rlmod_sig(k_rlmod,2) = nominal;  
      c_state = 0;
   end
end
% disturb pwrmod_p states
ntpwr_p = ntrl + g.pwr.n_pwrmod;
if g.pwr.n_pwrmod~=0
   disp('disturbing real power modulation')
   for k = ntrl+1:ntpwr_p
      j=1;
      k_pwrmod_p = k - ntrl;
      pert = p_ratio*abs(g.pwr.pwrmod_p_st(k_pwrmod_p,1));   
      pert = max(pert,p_ratio);
      g.pwr.pwrmod_p_st(k_pwrmod_p,2) = g.pwr.pwrmod_p_st(k_pwrmod_p,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 37;
      % disturb the input variable
      disp('disturbing pwrmod_p_sig') 
      c_state = 7; 
      pwrmod_p_input = k_pwrmod_p;
      pert = p_ratio;
%       nominal = 0.0;
      g.pwr.pwrmod_p_sig(k_pwrmod_p,2) = g.pwr.pwrmod_p_sig(k_pwrmod_p,2) + pert;
      p_file
%       pwrmod_p_sig(k_pwrmod_p,2) = nominal;  
      c_state = 0;
   end
end
% disturb pwrmod_q states
ntpwr_q = ntpwr_p + g.pwr.n_pwrmod;
if g.pwr.n_pwrmod~=0
   disp('disturbing reac power modulation')
   for k = ntpwr_p+1:ntpwr_q
      j=1;
      k_pwrmod_q = k - ntpwr_p;
      pert = p_ratio*abs(g.pwr.pwrmod_q_st(k_pwrmod_q,1));   
      pert = max(pert,p_ratio);
      g.pwr.pwrmod_q_st(k_pwrmod_q,2) = g.pwr.pwrmod_q_st(k_pwrmod_q,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 38;
      % disturb the input variable
      disp('disturbing pwrmod_q_sig') 
      c_state = 8; 
      pwrmod_q_input = k_pwrmod_q;
      pert = p_ratio;
%       nominal = 0.0;
      g.pwr.pwrmod_q_sig(k_pwrmod_q,2) = g.pwr.pwrmod_q_sig(k_pwrmod_q,2) + pert;
      p_file
%       pwrmod_q_sig(k_pwrmod_q,2) = nominal;  
      c_state = 0;
   end
end

ntdc = ntpwr_q + n_dcl;
% disturb the HVDC states
if n_conv~=0
   disp('disturbing HVDC')
   for k = ntpwr_q+1:ntdc
      j = 1;
      k_hvdc = k - ntpwr_q;
      pert = p_ratio*abs(v_conr(k_hvdc,1));
      pert = max(pert,p_ratio);
      v_conr(k_hvdc,2) = v_conr(k_hvdc,1) + pert;
      p_file;
      st_name(k,j) = 39;
      j = j + 1;
      pert = p_ratio*abs(v_coni(k_hvdc,1));
      pert = max(pert,p_ratio);
      v_coni(k_hvdc,2) = v_coni(k_hvdc,1) + pert;
      p_file;
      st_name(k,j) = 40;
      j= j+1;
      pert = p_ratio*abs(i_dcr(k_hvdc,1));
      pert = max(pert,p_ratio);
      i_dcr(k_hvdc,2) = i_dcr(k_hvdc,1) + pert;
      p_file;
      st_name(k,j) = 41;
      if ~isempty(cap_idx)
         k_cap_idx = find(cap_idx == k_hvdc);
         if ~isempty(k_cap_idx)
            pert = p_ratio*abs(i_dci(k_hvdc,1));
            pert = max(pert,p_ratio);
            i_dci(k_hvdc,2) = i_dci(k_hvdc,1) + pert;
            p_file;
            st_name(k,j) = 42;
            j = j + 1;
            pert = p_ratio*abs(v_dcc(k_hvdc,1));
            pert = max(pert,p_ratio);
            v_dcc(k_hvdc,2) = v_dcc(k_hvdc,1) + pert;
            p_file;
            st_name(k,j) = 43;
         end
      end
      disp('disturbing rectifier dc_sig') 
      c_state = 9; 
      dcmod_input = k_hvdc;
      pert = p_ratio;
      nominal = 0.0;
      dc_sig(r_idx(k_hvdc),2) = dc_sig(r_idx(k_hvdc),1) + pert;
      p_file
      dc_sig(r_idx(k_hvdc),2) = nominal;  
      c_state = 0;
      disp('disturbing inverter dc_sig') 
      c_state = 10; 
      dcmod_input = k_hvdc;
      pert = p_ratio;
      nominal = 0.0;
      dc_sig(i_idx(k_hvdc),2) = dc_sig(i_idx(k_hvdc),1) + pert;
      p_file
      dc_sig(i_idx(k_hvdc),2) = nominal;  
      c_state = 0;
   end
end
