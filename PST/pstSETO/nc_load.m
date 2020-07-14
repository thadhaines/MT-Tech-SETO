function V_nc = nc_load(bus,flag,Y22,Y21,psi,V_o,tol,k,kdc)
% Syntax: V_nc = nc_load(bus,flag,Y22,Y21,psi,V_o,tol,k,kdc)
% 4:59 PM 15/08/97
% Purpose: non-conforming load model, for constant power 
%          and constant current loads; the nonlinear
%          equations are solved using a Newton's algorithm 
%
% Input: bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - computation not needed
%        Y22 - reduced self Y matrix of non-conforming loads
%        Y21 - reduced mutual Y matrix of generator internal
%              node to non-conforming loads 
%        psi - machine internal node voltage,
%              not used in initialization
%        V_o - initial non-conforming load bus voltage 
%              values (vector), not used in initialization
%        tol - tolerance for Newton's algorithm convergence, 
%              not used in initialization
%        k   - integer time (only for svc/facts models),
%              not used in initialization
%        kdc - dc integer time
%  
% Output: V_nc - solved non-conforming load bus voltage 
%                values (vector)
%
% Files:
%
% See Also: red_Ybus, svc

% Algorithm: 
%
% Calls: dc_load  dc_cur
%
% Called By:svm_mgen, s_simu

% (c) Copyright 1991-1997 Joe H. Chow/Cherry Tree Scientific Software
%     All Rights Reserved

% History (in reverse chronological order)
%
% Version:      2.3   
% Date:         Feb 2015  
% Author:       D. Trudnowski  
% Purpose:      Add power modulation 
%  
% Version:      2.2
% Date:         August 1997
% Author:       Graham Rogers
% Purpose:      Add load modulation
%
% Version:      2.1
% Date:         March 1997
% Author:       Graham Rogers
% Purpose:      Add dc
% Modification: Added section to determine VHT for dc converters
% Version:      2.0
% Date:         November 1996
% Author:       Graham Rogers
% Purpose:      Improve convergence
% Modification: reduction to constant impedance load for low voltages
%               modification of Jacobian

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1991
%     %% HVDC link variables - 63
%     global  dcsp_con  dcl_con  dcc_con
%     global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
%     global  inv_ac_line  rec_ac_line ac_line dcli_idx
%     global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
%     global  Vdc  i_dc P_dc i_dcinj dc_pot alpha gamma VHT dc_sig  cur_ord dcr_dsig dci_dsig
%     global  ric_idx  rpc_idx Vdc_ref dcc_pot
%     global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
%     global  ndcr_ud ndci_ud dcrud_idx dciud_idx dcrd_sig dcid_sig
% 
%     % States
%     %line
%     global i_dcr i_dci  v_dcc
%     global di_dcr  di_dci  dv_dcc
%     global dc_dsig % added 07/13/20 -thad
%     %rectifier
%     global v_conr dv_conr
%     %inverter
%     global v_coni dv_coni
%     
%     % added to global dc
%     global xdcr_dc dxdcr_dc xdci_dc dxdci_dc angdcr angdci t_dc
%     global dcr_dc dci_dc % damping control
%     global  ldc_idx
% 
% % % dc variables
% % global  i_dci  i_dcr  dcc_pot  alpha  gamma  basmva  r_idx  i_idx
% % global  n_conv n_dcl ldc_idx

global g 



if ~isempty(g.ncl.load_con)
   jay = sqrt(-1);
   if flag == 0; % initialization
       
      g.ncl.nload = size(g.ncl.load_con,1);
      %  set up constant power and current load components in 
      %    load_pot
      %  vectorized computation
      j = g.sys.bus_int(g.ncl.load_con(:,1));
      % no need for special treatment for dc buses on initialization
      V_nc = bus(j,2).*exp(jay*bus(j,3)*pi/180);
      % constant power component
      g.ncl.load_pot(:,1) = bus(j,6).*g.ncl.load_con(:,2) ...
                      + jay*bus(j,7).*g.ncl.load_con(:,3);
      S_cc = bus(j,6).*g.ncl.load_con(:,4) ...
             + jay*bus(j,7).*g.ncl.load_con(:,5);
      % constant current component
      g.ncl.load_pot(:,2) = S_cc./abs(V_nc);
      g.ncl.load_pot(:,3) = g.ncl.load_pot(:,1)./V_nc./conj(V_nc);% const impedance equ of const power load
      g.ncl.load_pot(:,4) = S_cc./V_nc./conj(V_nc);% const impedance equ of constant current load
   end
   if flag == 1 % network interface computation 
      if nargin < 7
         tol = 1e-6;
%         tol = 1e-5;   % JHC August 3, 2018
      end
      
      g.ncl.nload = size(g.ncl.load_con,1);
      
      if nargin < 6
         V_o = ones(g.ncl.nload,1);%set default trial voltages
      end
      V_nc = V_o;
      
      if g.svc.n_svc~=0
         j = g.svc.svc_idx;    
         Y22(j,j) = Y22(j,j)+jay*diag(g.svc.B_cv(:,k));% note that Y22 is a local variable
      end
      if g.tcsc.n_tcsc~=0
         j = g.tcsc.tcscf_idx; 
         jj = g.tcsc.tcsct_idx;
         Y22(j,j)=Y22(j,j) + jay*diag(g.tcsc.B_tcsc(:,k));Y22(j,jj)=Y22(j,jj) - jay*diag(g.tcsc.B_tcsc(:,k));
         Y22(jj,j)=Y22(jj,j) - jay*diag(g.tcsc.B_tcsc(:,k));Y22(jj,jj)=Y22(jj,jj) + jay*diag(g.tcsc.B_tcsc(:,k));
      end
      if g.lmod.n_lmod ~=0 
         j = g.lmod.lmod_idx;
         Y22(j,j) = Y22(j,j) + diag(g.lmod.lmod_st(:,k));
      end
      if g.rlmod.n_rlmod ~=0
         j = g.rlmod.rlmod_idx;
         Y22(j,j) = Y22(j,j) + jay*diag(g.rlmod.rlmod_st(:,k));
      end
      
      %pwrmod % added 06/11/20
      load_pot_mod = g.ncl.load_pot;
      if g.pwr.n_pwrmod~=0
          for index=1:g.pwr.n_pwrmod
            if (g.ncl.load_con(g.pwr.pwrmod_idx(index),2)==1 && g.ncl.load_con(g.pwr.pwrmod_idx(index),3)==1)
                load_pot_mod(g.pwr.pwrmod_idx(index),1) = - (g.pwr.pwrmod_p_st(index,k) + jay*g.pwr.pwrmod_q_st(index,k)); %power modulation
%                 load_pot_mod(pwrmod_idx(index),3) = ...
%                     - (pwrmod_p_st(index,k) + jay*pwrmod_q_st(index,k))/(V_nc(index)*conj(V_nc(index))); %power modulation with v<0.5 - maybe disable?
            elseif (g.ncl.load_con(g.pwr.pwrmod_idx(index),4)==1 && g.ncl.load_con(g.pwr.pwrmod_idx(index),5)==1); 
                load_pot_mod(g.pwr.pwrmod_idx(index),2) = - (g.pwr.pwrmod_p_st(index,k) + jay*g.pwr.pwrmod_q_st(index,k)); %current modulation
%                 load_pot_mod(pwrmod_idx(index),4) = ...
%                     - (pwrmod_p_st(index,k) + jay*pwrmod_q_st(index,k))/abs(V_nc(index)); %current modulation with v<.5 - maybe disable?
            end
          end
      end
      
      lv_idx = find(abs(V_nc)<=0.5);
      
      hv_idx = find(abs(V_nc) > 0.5);
      curr_mis = zeros(g.ncl.nload,1);
      curr_load = Y21*psi + Y22*V_nc;
      if ~isempty(hv_idx)
         curr_nc = -conj((load_pot_mod(hv_idx,1)+load_pot_mod(hv_idx,2)...
            .*abs(V_nc(hv_idx)))./V_nc(hv_idx));
      end     
      
      if g.dc.n_conv~=0
         %modify nc-current to take account of ac current
         i_ac = dc_cur(V_nc(g.dc.ldc_idx),k,kdc);
         curr_nc(g.dc.ldc_idx) = curr_nc(g.dc.ldc_idx)-i_ac;
      end
      if ~isempty(hv_idx)
         curr_mis(hv_idx)=curr_nc - curr_load(hv_idx); % current mismatch
      end
      % converts to constant impedance abs(V_nc)<=0.5
      if ~isempty(lv_idx)
         %l_vl = length(lv_idx);
         curr_mis(lv_idx) = -conj(diag(load_pot_mod(lv_idx,3)+load_pot_mod(lv_idx,4)))*V_nc(lv_idx)...
            - curr_load(lv_idx);
      end
      count = 0;
      Y22_real = real(Y22); Y22_imag = imag(Y22);
      % Newton's algorithm
      while(norm(curr_mis,'inf') > tol)
         if ~isempty(hv_idx)
            % form Jacobian
            v_re = real(V_nc(hv_idx));
            v_im = imag(V_nc(hv_idx));
            v_mag = abs(V_nc(hv_idx));
            v_mag2 = v_mag.*v_mag;
            v1 = (v_re.*v_re - v_im.*v_im);
            v2 = v_re.*v_im;
            vp11 = -v1.*real(load_pot_mod(hv_idx,1))-2*v2.*imag(load_pot_mod(hv_idx,1));
            vp12 = v1.*imag(load_pot_mod(hv_idx,1))-2*v2.*real(load_pot_mod(hv_idx,1));
            vi11 = v_im.*v_im.*real(load_pot_mod(hv_idx,2)) - v2.*imag(load_pot_mod(hv_idx,2));
            vi12 = v_re.*v_re.*imag(load_pot_mod(hv_idx,2)) - v2.*real(load_pot_mod(hv_idx,2));
            vi21 = -v_im.*v_im.*imag(load_pot_mod(hv_idx,2)) -v2.*real(load_pot_mod(hv_idx,2));
            vi22 = v_re.*v_re.*real(load_pot_mod(hv_idx,2)) + v2.*imag(load_pot_mod(hv_idx,2));
            vp11 = vp11./v_mag2./v_mag2;
            vp12 = vp12./v_mag2./v_mag2;
            vi11 = vi11./v_mag./v_mag2;
            vi12 = vi12./v_mag./v_mag2;
            vi21 = vi21./v_mag./v_mag2;
            vi22 = vi22./v_mag./v_mag2;
            v_s11(hv_idx) = (vi11 + vp11);
            v_s12(hv_idx) = (vp12 + vi12);
            v_s21(hv_idx) = (vp12 + vi21);
            v_s22(hv_idx) = -vp11 +vi22;
         end;
         if ~isempty(lv_idx)
            v_s11(lv_idx) = -real(load_pot_mod(lv_idx,3)+load_pot_mod(lv_idx,4));
            v_s12(lv_idx) = imag(load_pot_mod(lv_idx,3)+load_pot_mod(lv_idx,4));
            v_s22(lv_idx) = v_s11(lv_idx);
            v_s21(lv_idx) = -v_s12(lv_idx);
         end
         % modify Jacobian for dc
         if g.dc.n_conv~=0
            V = V_o(g.dc.ldc_idx);
            [y1,y2,y3,y4] = dc_load(V,k,kdc);
            v_s11(g.dc.ldc_idx) = v_s11(g.dc.ldc_idx) + y1';
            v_s12(g.dc.ldc_idx) = v_s12(g.dc.ldc_idx) + y2';
            v_s21(g.dc.ldc_idx) = v_s21(g.dc.ldc_idx) + y3';
            v_s22(g.dc.ldc_idx) = v_s22(g.dc.ldc_idx) + y4';
         end
         Jac_nc = [ Y22_real+diag(v_s11) diag(v_s12)-Y22_imag;
            Y22_imag+diag(v_s21)   Y22_real+diag(v_s22)];
         b = [ real(curr_mis); imag(curr_mis) ];
         x = Jac_nc\b;   % solve for voltage increment
         V_nc = V_nc + x(1:g.ncl.nload,1) + jay*x(g.ncl.nload+1:2*g.ncl.nload,1);
         % update voltage
         count = count + 1;
         lv_idx = find(abs(V_nc)<=0.5);
         hv_idx = find(abs(V_nc) > 0.5);
         curr_load = Y21*psi + Y22*V_nc;
         if ~isempty(hv_idx)
            curr_mis(hv_idx)=-conj((load_pot_mod(hv_idx,1)+load_pot_mod(hv_idx,2)...
               .*abs(V_nc(hv_idx)))./V_nc(hv_idx))...
               -curr_load(hv_idx); % current mismatch
         end
         if ~isempty(lv_idx)
            l_vl = length(lv_idx);
            curr_mis(lv_idx) = -conj(diag(load_pot_mod(lv_idx,3)+load_pot_mod(lv_idx,4)))*V_nc(lv_idx)...
               - curr_load(lv_idx);
         end
         if g.dc.n_conv~=0
            % modify mismatch for dc
            i_ac = dc_cur(V_nc(g.dc.ldc_idx),k,kdc);
            curr_mis(g.dc.ldc_idx) = curr_mis(g.dc.ldc_idx) - i_ac;
         end
         if count > 30
            disp('NC_LOAD: Newton algorithm not converged in 30 iterations')
            fprintf('current mismatch is %g \n', curr_mis) 
            error('executuion terminated')
         end
      end
   end
   if flag == 2 % no dynamics calculation needed
   end
end
