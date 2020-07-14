function  dc_line(i,k,kdc,bus,flag)
%Syntax:  f = dc_line(i,kdc,bus,flag)
% 5:14 PM 15/08/97
%Purpose: Models HVDC line dynamics
% Input: i - 0 for vectorized computation only option
%        k  - integer time
%        kc - integer time for dc
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - hvdc dynamics computation 
%               3 - state matrix building 
%
% Output: f - dummy variable 
%
% Calls:
%
% Called By:

% (c) Copyright 1991-1997 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:  1.0
% Date:     February 1997
% Author:   Graham Rogers

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
% %define global variables
% global dcsp_con  dcl_con  dcc_con  dcc_pot  dc_pot
% global  r_idx  i_idx n_dcl  n_conv  
% global  Vdc  i_dc  
% global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
% 
% % States
% %line
% global i_dcr i_dci  v_dcc
% global di_dcr  di_dci  dv_dcc  

global g
% jay = sqrt(-1);
k=fix(1+(kdc+1)/10);
% check for dcline data
if ~isempty(g.dc.dcsp_con)
  if flag==0
  % initialization
    % check that inductances are specified
    if ~isempty(g.dc.no_ind_idx)
        error('you must specify inductances for dc lines and smoothing reactors')
    end

    if i==0
    % vector computation
      g.dc.i_dcr(:,1) = g.dc.i_dc(g.dc.r_idx,1);
      g.dc.i_dci(:,1) = g.dc.i_dc(g.dc.i_idx,1);
      g.dc.v_dcc(:,1) = g.dc.Vdc(g.dc.r_idx,1) - g.dc.i_dcr(:,1).*g.dc.dcl_con(:,3)/2.0;
      if g.dc.l_no_cap~=0
        g.dc.dc_pot(g.dc.no_cap_idx,3) = zeros(g.dc.l_no_cap,1);
        g.dc.dc_pot(g.dc.no_cap_idx,1) = 1000*ones(g.dc.l_no_cap,1)./(g.dc.dcl_con(g.dc.no_cap_idx,4) ...
                               + g.dc.dcl_con(g.dc.no_cap_idx,6) + g.dc.dcl_con(g.dc.no_cap_idx,7));
        g.dc.dc_pot(g.dc.no_cap_idx,2) = g.dc.dcl_con(g.dc.no_cap_idx,3);
        g.dc.dc_pot(g.dc.no_cap_idx,4) = g.dc.dc_pot(g.dc.no_cap_idx,1);
        g.dc.dc_pot(g.dc.no_cap_idx,5) = g.dc.dc_pot(g.dc.no_cap_idx,2);
      end
      if g.dc.l_cap ~= 0
        g.dc.dc_pot(g.dc.cap_idx,1) = 1000*ones(g.dc.l_cap,1)./(g.dc.dcl_con(g.dc.cap_idx,4)*.5 + g.dc.dcl_con(g.dc.cap_idx,6));
        g.dc.dc_pot(g.dc.cap_idx,2) = g.dc.dcl_con(g.dc.cap_idx,3)/2;
        g.dc.dc_pot(g.dc.cap_idx,3) = 1e6*ones(g.dc.l_cap,1)./g.dc.dcl_con(g.dc.cap_idx,5);
        g.dc.dc_pot(g.dc.cap_idx,4) = 1000*ones(g.dc.l_cap,1)./(g.dc.dcl_con(g.dc.cap_idx,4)*.5 + g.dc.dcl_con(g.dc.cap_idx,7));
        g.dc.dc_pot(g.dc.cap_idx,5) = g.dc.dcl_con(g.dc.cap_idx,3)/2;
      end
    else
      error(' no non-vector computation in HVDC')
    end
  end
  if flag == 1
  % network inter face - no calculation required
  end
  if flag == 2
  % rate of change of states
    if i== 0
    %vector compuation
      if g.dc.l_cap~=0 
        g.dc.di_dcr(g.dc.cap_idx,kdc) = - g.dc.dc_pot(cg.dc.ap_idx,1)...
                           .*(g.dc.dc_pot(g.dc.cap_idx,2).*g.dc.i_dcr(g.dc.cap_idx,kdc) + g.dc.v_dcc(g.dc.cap_idx,kdc) - g.dc.Vdc(g.dc.r_idx(g.dc.cap_idx),kdc));
        g.dc.di_dci(g.dc.cap_idx,kdc) = -g.dc.dc_pot(g.dc.cap_idx,4)...
                           .*(g.dc.dc_pot(g.dc.cap_idx,5).*g.dc.i_dci(g.dc.cap_idx,kdc) - g.dc.v_dcc(g.dc.cap_idx,kdc) + g.dc.Vdc(g.dc.i_idx(g.dc.cap_idx),kdc));
        g.dc.dv_dcc(g.dc.cap_idx,kdc) = g.dc.dc_pot(g.dc.cap_idx,3).*(g.dc.i_dcr(g.dc.cap_idx,kdc) - g.dc.i_dci(g.dc.cap_idx,kdc));
      end
      if g.dc.l_no_cap~=0
        g.dc.di_dcr(g.dc.no_cap_idx,kdc) = (-g.dc.dc_pot(g.dc.no_cap_idx,2).*g.dc.i_dcr(g.dc.no_cap_idx,kdc) + ...
                               g.dc.Vdc(g.dc.r_idx(g.dc.no_cap_idx),kdc)-g.dc.Vdc(g.dc.i_idx(g.dc.no_cap_idx),kdc))... 
                               .*g.dc.dc_pot(g.dc.no_cap_idx,1);
        g.dc.di_dci(g.dc.no_cap_idx,kdc) = g.dc.di_dcr(g.dc.no_cap_idx,kdc);
        g.dc.dv_dcc(g.dc.no_cap_idx) = 0;
      end
    else
      error(' no non-vector calculation in HVDC')
    end
  end
end
