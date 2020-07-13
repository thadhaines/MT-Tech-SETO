% This m-file forms a sparse permutation matrix p_mat 
% which converts the vector of rates of change 
% of state perturbations (d_vector) into a column of
% the state matrix a_mat
% 5:16 PM 18/8/97
% called by svm_mgen
% Author: Graham Rogers
% Date: September 1996
% Modified December 1998
% tcsc model added
% Modified: 2015, added pwrmod, D. Trudnowski
% Modified: June 1998
%           hydro turbine added
% Modified: April 1998
%           lead lag added to SVC
% Mofified: August 1997
%           Induction Generator added
% Modified: August 1997
%           Load modulation added
% Modified: April 1997
%           HVDC added
% Modified: November 1996
%           svc added 
% Copyright: Joe Chow/Cherry Tree Scientific Software 1991 to 1997, All rights reserved


k_row = 1;
exc_count = 0;
pss_count = 0;
dpw_count = 0;
tg_count = 0;

for k = 1:g.mac.n_mac
   if state(k)~=0
      if k ~=1
         k_row = 1+sum(state(1:k-1));
      end
      % generators
      k_col = k;
      for kgs = 1:2
         p_mat(k_row+kgs-1,k_col+(kgs-1)*g.mac.n_mac) = 1;%all generators
      end
      % sub or tra models
      if gen_state(k)>2
          p_mat(k_row+2, k_col+2*g.mac.n_mac)=1;
      end
      %tra models
      if gen_state(k)==4
          p_mat(k_row+3,k_col+4*g.mac.n_mac)=1;
      end
      %sub models
      if gen_state(k)==6
         for kgs = 4:6
            p_mat(k_row+kgs-1,k_col+(kgs-1)*g.mac.n_mac) = 1;
         end
      end
      k_row = k_row + gen_state(k)-1;
      k_col = g.mac.n_mac*6;
      
      if mac_exc ~=0
         k_exc = find(mac_exc==k);
         if ~isempty(k_exc)
            exc_count = exc_count +1;
            k_cex = k_col+k_exc;
            %exciters
            if TR_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cex) = 1;
            end
            k_cex = k_cex + g.exc.n_exc;     
            if TB_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cex) = 1;
            end
            k_cex = k_cex + g.exc.n_exc;
            if TA_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cex) = 1;
            end
            k_cex = k_cex + g.exc.n_exc;     
            if Efd_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cex) = 1;
            end
            k_cex = k_cex + g.exc.n_exc;     
            if R_f_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cex) = 1;
            end 
         end
      end
      k_col = 6*g.mac.n_mac+5*g.exc.n_exc;
      
      if mac_pss~=0
         k_pss = find(mac_pss==k);
         if ~isempty(k_pss) 
            pss_count = pss_count+1;   
            %pss
            k_cpss = k_col + k_pss;
            if pss1_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cpss) = 1;  
            end
            k_cpss = k_cpss + g.pss.n_pss;     
            if pss2_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cpss) = 1;
            end
            k_cpss = k_cpss + g.pss.n_pss;     
            if pss3_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_cpss) = 1;
            end 
         end 
      end
      k_col = 6*g.mac.n_mac+5*g.exc.n_exc+3*g.pss.n_pss;

      if mac_dpw~=0
         k_dpw = find(mac_dpw==k);
         if ~isempty(k_dpw) 
            dpw_count = dpw_count+1;   
            %dpw filter
            k_cdpw = k_col + k_dpw;
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1;  
            k_cdpw = k_cdpw + n_dpw;     
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1;
            k_cdpw = k_cdpw + k_dpw;
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1; 
            k_cdpw = k_cdpw + k_dpw;
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1;  
            k_cdpw = k_cdpw + k_dpw;
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1;  
            k_cdpw = k_cdpw + k_dpw;
            k_row = k_row + 1;
            p_mat(k_row,k_cdpw) = 1;  
         end 
      end
      k_col= 6*g.mac.n_mac+5*g.exc.n_exc+3*g.pss.n_pss+6*n_dpw;    
      % governors
      if mac_tg ~=0
         k_tg = find(mac_tg == k);
         if ~isempty(k_tg)
            tg_count = tg_count + 1;
            k_ctg = k_col + g.tg.tg_idx(k_tg);
            if tg_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg+g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg+g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
            end
         end
      end
      
      if mac_tgh ~=0
         k_tgh = find(mac_tgh==k);
         if ~isempty(k_tgh)
            tg_count = tg_count + 1;
            k_ctg = k_col + g.tg.tgh_idx(k_tgh);
            if tgh_state(k)~=0
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg + g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg + g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg + g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
               k_ctg = k_ctg + g.tg.n_tg + g.tg.n_tgh;
               k_row = k_row + 1;
               p_mat(k_row,k_ctg) = 1;
            end
         end
      end
   end
end

% p_mat for induction motors

k_colg = 6*g.mac.n_mac+5*g.exc.n_exc+3*g.pss.n_pss+6*n_dpw...
    +5*(g.tg.n_tg + g.tg.n_tgh);
k_col = k_colg;
if n_mot~=0
   for k=1:n_mot
      k_row = k_row+1;
      k_col = k_col+k;
      p_mat(k_row,k_col)=1;
      k_row=k_row+1;
      k_col = k_col+n_mot;
      p_mat(k_row,k_col)=1;
      k_row = k_row+1;
      k_col = k_col+n_mot;
      p_mat(k_row,k_col) = 1;
      k_col =6*g.mac.n_mac+5*g.exc.n_exc+3*g.pss.n_pss+5*(g.tg.n_tg+g.tg.n_tgh);
   end
end

k_col = k_colg+3*n_mot;

if g.igen.n_ig~=0
   for k=1:g.igen.n_ig
      k_row = k_row+1;
      k_col = k_col+k;
      p_mat(k_row,k_col)=1;
      k_row=k_row+1;
      k_col = k_col+g.igen.n_ig;
      p_mat(k_row,k_col)=1;
      k_row = k_row+1;
      k_col = k_col+g.igen.n_ig;
      p_mat(k_row,k_col) = 1;
      k_col =6*g.mac.n_mac+5*g.exc.n_exc+3*g.pss.n_pss ...
          +5*(g.tg.n_tg+g.tg.n_tgh)+3*n_mot;
   end
end

% p_mat for svcs
k_col = k_colg+3*n_mot+3*g.igen.n_ig;

if g.svc.n_svc ~=0
   for k = 1:g.svc.n_svc
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      if ~isempty(g.svc.svcll_idx)
         kcon = find(k == g.svc.svcll_idx);
         if ~isempty(kcon)
            k_row = k_row + 1;
            k_col = k_col + g.svc.n_svc;
            p_mat(k_row,k_col)=1;
         end
      end
      k_col = k_colg+3*n_mot+3*g.igen.n_ig;
   end
end  
% p_mat for tcscs

k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc;

if g.tcsc.n_tcsc ~=0
   for k = 1:g.tcsc.n_tcsc
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc;
   end
end  

% p_mat for lmod
k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc;
if g.lmod.n_lmod ~=0
   for k = 1:g.lmod.n_lmod
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc;
   end
end

% p_mat for rlmod
k_col = k_colg+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod;
if g.rlmod.n_rlmod ~=0
   for k = 1:g.rlmod.n_rlmod
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod;
   end
end  

% p_mat for pwrmod_p
k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod;
if g.pwr.n_pwrmod ~=0
   for k = 1:g.pwr.n_pwrmod
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod;
   end
end

% p_mat for pwrmod_q
k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod+g.pwr.n_pwrmod;
if g.pwr.n_pwrmod ~=0
   for k = 1:g.pwr.n_pwrmod
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col)=1;
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod+g.pwr.n_pwrmod;
   end
end

% p_mat for hvdc links
k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc+g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod+2*g.pwr.n_pwrmod;
if n_conv~=0
   for k = 1:n_dcl
      %  converter controls 
      k_row = k_row + 1;
      k_col = k_col + k;
      p_mat(k_row,k_col) = 1;
      k_row = k_row + 1;
      k_col = k_col + n_dcl;
      p_mat(k_row,k_col) = 1;
      % hvdc line
      k_row = k_row + 1;
      k_col = k_col + n_dcl;
      p_mat(k_row,k_col) = 1;
      if l_cap~=0
         if ~isempty(cap_idx(k))
            % capicitor in this line model
            k_row = k_row + 1;
            k_col = k_col + n_dcl;
            p_mat(k_row,k_col) = 1;
            k_row = k_row+1;
            k_col = k_col + n_dcl;
            p_mat(k_row, k_col) = 1;
         end
      end
      k_col = k_colg+3*n_mot+3*g.igen.n_ig+2*g.svc.n_svc...
          +g.tcsc.n_tcsc+g.lmod.n_lmod+g.rlmod.n_rlmod...
          +2*g.pwr.n_pwrmod;
   end
end 