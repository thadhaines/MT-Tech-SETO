function dc_cont(i,k,kdc,bus,flag)
%DC_CONT models hvdc pole controls
% DC_CONT  models hvdc pole controls
%
% Syntax: dc_cont(i,k,kdc,bus,flag)
%
%   NOTES:  
% 
%   Input: 
%   i - 0 vector computaion only for HVDC control
%   k - integer time (data index)
%   kdc - integer time for dc (dc data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - dynamics computation and state state matrix building
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   04/xx/97    XX:XX   Graham Rogers  	Version 1.0
%   (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%   07/15/20    10:34   Thad Haines     Revised format of globals and internal function documentation

global g

% check that dc controls are defined
if ~isempty(g.dc.dcc_con)
   if flag == 0
      % initialize controls
      if i == 0 
         % vector computation
         % rectifier controls
         % calculate current order
         g.dc.cur_ord(g.dc.r_idx,1) = g.dc.i_dc(g.dc.r_idx,1);
         g.dc.cur_ord(g.dc.i_idx,1) = g.dc.i_dc(g.dc.i_idx,1).*(ones(g.dc.n_dcl,1)-g.dc.dcl_con(:,9)/100.0);
         g.dc.Vdc_ref = g.dc.Vdc(g.dc.i_idx,1);
         g.dc.v_conr(:,1) = g.dc.alpha(:,1)./g.dc.dcc_con(g.dc.r_idx,4);
         vrmax=find(g.dc.v_conr(:,1)>g.dc.dcc_con(g.dc.r_idx,5));
         vrmin=find(g.dc.v_conr(:,1)<g.dc.dcc_con(g.dc.r_idx,6));
         % check rectifier integrator limit violation
         if ~isempty(vrmax)
            svrm=int2str(vrmax);
            disp('v_conr greater than maximum limit at the following rectifiers')
            disp(svrm)
            error
         end
         if ~isempty(vrmin)
            svrm=int2str(vrmin);
            disp('v_conr less than minimum limit at the following rectifiers')
            disp(svrm)
            error
         end
         
         g.dc.v_coni(:,1) = g.dc.gamma(:,1)./g.dc.dcc_con(g.dc.i_idx,4);
         vimax=find(g.dc.v_coni(:,1)>g.dc.dcc_con(g.dc.i_idx,5));
         vimin=find(g.dc.v_coni(:,1)<g.dc.dcc_con(g.dc.i_idx,6));
         % check rectifier integrator limit violation
         if ~isempty(vimax)
            svim=int2str(vimax);
            disp('v_coni greater than maximum limit at the following inverters')
            disp(svim)
            error
         end
         if ~isempty(vimin)
            svim=int2str(vimin);
            disp('v_coni less than minimum limit at the following inverters')
            disp(svim)
            error
         end
         
         g.dc.dcc_pot(:,1) = g.dc.gamma(:,1); %gamma reference value
         g.dc.dcc_pot(:,2) = g.dc.dcsp_con(g.dc.r_idx,5)./g.dc.dcsp_con(g.dc.r_idx,6); 
         g.dc.dcc_pot(:,2) = g.dc.dcc_pot(:,2)*g.sys.basmva./bus(g.dc.rec_ac_bus,13)./bus(g.dc.rec_ac_bus,13); % xeqpu rec
         g.dc.dcc_pot(:,3) = g.dc.dcsp_con(g.dc.r_idx,6).*g.dc.dcsp_con(g.dc.r_idx,5)*3/pi;  % Rc rectifiers   
         g.dc.dcc_pot(:,4) = g.dc.dcsp_con(g.dc.i_idx,5)./g.dc.dcsp_con(g.dc.i_idx,6); 
         g.dc.dcc_pot(:,4) = g.dc.dcc_pot(:,4)*g.sys.basmva./bus(g.dc.inv_ac_bus,13)./bus(g.dc.inv_ac_bus,13); %xeqpu inv
         g.dc.dcc_pot(:,5) = g.dc.dcsp_con(g.dc.i_idx,6).*g.dc.dcsp_con(g.dc.i_idx,5)*3/pi;  % Rc inverters  
         g.dc.dcc_pot(:,6) = g.dc.dcc_pot(:,5); 
         % multiplier ideal rec dc voltage 
         g.dc.dcc_pot(:,7) = 3*sqrt(2).*g.dc.dcsp_con(g.dc.r_idx,6).*bus(g.dc.rec_ac_bus,13)/pi;
         % multiplier ideal inv dc voltage
         g.dc.dcc_pot(:,8) = 3*sqrt(2).*g.dc.dcsp_con(g.dc.i_idx,6).*bus(g.dc.inv_ac_bus,13)/pi;
         % multiplier acpu to dc amps
         g.dc.dcc_pot(:,9) = pi*g.sys.basmva/3/sqrt(2)./bus(g.dc.rec_ac_bus,13)./g.dc.dcsp_con(g.dc.r_idx,6);%rectifier
         g.dc.dcc_pot(:,10) = pi*g.sys.basmva/3/sqrt(2)./bus(g.dc.inv_ac_bus,13)./g.dc.dcsp_con(g.dc.i_idx,6);%inverter
         g.dc.dc_dsig(:,1)=zeros(g.dc.n_conv,1); % zero damping control signals
         if g.dc.dc_sig(:,1)~=zeros(g.dc.n_conv,1)
            % reset initial values of alpha and gamma
            g.dc.alpha(:,1) = ((-g.dc.cur_ord(g.dc.r_idx,1) + g.dc.dc_sig(g.dc.r_idx,1) + g.dc.dcr_dsig(:,1)...
               + g.dc.i_dcr(:,1)).*g.dc.dcc_con(g.dc.r_idx,2) + ...    
               g.dc.v_conr(:,1)).*g.dc.dcc_con(g.dc.r_idx,4); 
            g.dc.gamma(:,1) = ((g.dc.Vdc(g.dc.i_idx,1)-g.dc.Vdc_ref)./g.dc.Vdc_ref.*g.dc.dcc_con(g.dc.i_idx,2)...
               + g.dc.v_coni(:,1)).*g.dc.dcc_con(g.dc.i_idx,4);
         end
         
      else
         error('vector computation only in dc controls')
      end
      % end of initialization
   end
   if flag== 1
      % network interface
      if i ~= 0
         error('vector computation only with dc')
      else
         % vector computation 
         % i_dc vcdc and the control states are fixed
         
         % determine firing and extinction angles 
         g.dc.alpha(:,kdc) = ((-g.dc.cur_ord(g.dc.r_idx,k) + g.dc.dc_sig(g.dc.r_idx,k) + g.dc.dcr_dsig(:,k)...
            + g.dc.i_dcr(:,kdc)).*g.dc.dcc_con(g.dc.r_idx,2) + ...    
            g.dc.v_conr(:,kdc)).*g.dc.dcc_con(g.dc.r_idx,4); 
         % check for alpha limits
         g.dc.alpha(:,kdc) = max(g.dc.alpha(:,kdc), g.dc.dcc_con(g.dc.r_idx,8)*pi/180);
         g.dc.alpha(:,kdc) = min(g.dc.alpha(:,kdc), g.dc.dcc_con(g.dc.r_idx,7)*pi/180);
         
         g.dc.gamma(:,kdc) = ((g.dc.Vdc(g.dc.i_idx,kdc)-g.dc.Vdc_ref)./g.dc.Vdc_ref.*g.dc.dcc_con(g.dc.i_idx,2)...
            + g.dc.v_coni(:,kdc)).*g.dc.dcc_con(g.dc.i_idx,4);
        
         cur_error = g.dc.i_dci(:,kdc) - g.dc.cur_ord(g.dc.i_idx,k);
         ce_idx = find(cur_error < 0);
         if ~isempty(ce_idx)
            g.dc.gamma(ce_idx,kdc) = g.dc.gamma(ce_idx,kdc) + cur_error(ce_idx).* ...
               g.dc.dcc_con(g.dc.i_idx(ce_idx),2).*g.dc.dcc_con(g.dc.i_idx(ce_idx),4);
         end
         % check gamma limits
         g.dc.gamma(:,kdc) = max(g.dc.gamma(:,kdc), g.dc.dcc_con(g.dc.i_idx,8)*pi/180);
         g.dc.gamma(:,kdc) = min( g.dc.gamma(:,kdc), g.dc.dcc_con(g.dc.i_idx,7)*pi/180);
      end 
   end
   if flag == 2
      % calculate rates of change of states
      if i==0
         % vector computation 
         % rectifier 
         g.dc.dv_conr(:,kdc) = (-g.dc.cur_ord(g.dc.r_idx,k) + g.dc.i_dcr(:,kdc) ...
             + g.dc.dc_sig(g.dc.r_idx,k) + g.dc.dcr_dsig(:,k))...
            .*g.dc.dcc_con(g.dc.r_idx,3);
         %check for state limits
         recmx  = find(g.dc.v_conr(:,kdc)>g.dc.dcc_con(g.dc.r_idx,5));
         if ~isempty(recmx)
            g.dc.v_conr(recmx,kdc) = g.dc.dcc_con(g.dc.r_idx(recmx),5);
            recdmx = find(g.dc.dv_conr(recmx,k)>0);
            if ~isempty(recdmx)
               g.dc.dv_conr(recmx(recdmx),kdc) = zeros(length(recdmx),1);
            end
         end
         recmn  = find(g.dc.v_conr(:,kdc)<g.dc.dcc_con(g.dc.r_idx,6));
         if ~isempty(recmn)
            g.dc.v_conr(recmn,kdc) = g.dc.dcc_con(g.dc.r_idx(recmn),6);
            recdmn = find(g.dc.dv_conr(recmn,kdc)<0);
            if ~isempty(recdmn)
               g.dc.dv_conr(recmn(recdmn),kdc) = zeros(length(recdmn),1);
            end
         end  
         
         %inverter
         cur_err = g.dc.cur_ord(g.dc.i_idx,k) - g.dc.i_dci(:,kdc);
         n_ce = find(cur_err<0);
         if ~isempty(n_ce)
            cur_err(n_ce) = zeros(length(n_ce),1);
         end
         inv_err = (g.dc.Vdc(g.dc.i_idx,kdc) - g.dc.Vdc_ref)./g.dc.Vdc_ref - cur_err;
         g.dc.dv_coni(:,kdc) = (inv_err + g.dc.dc_sig(g.dc.i_idx,k)+g.dc.dci_dsig(:,k))...
            .*g.dc.dcc_con(g.dc.i_idx,3);
         % check state limits
         invmx  = find(g.dc.v_coni(:,kdc)>g.dc.dcc_con(g.dc.i_idx,5));
         if ~isempty(invmx)
            g.dc.v_coni(invmx,kdc) = g.dc.dcc_con(g.dc.i_idx(invmx),5);
            invdmx = find(g.dc.dv_coni(invmx,k)>0);
            if ~isempty(invdmx)
               g.dc.dv_coni(recmx(invdmx),kdc) = zeros(length(invdmx),1);
            end
         end
         invmn  = find(g.dc.v_coni(:,kdc)<g.dc.dcc_con(g.dc.i_idx,6));
         if ~isempty(invmn)
            g.dc.v_coni(invmn,kdc) = g.dc.dcc_con(g.dc.i_idx(invmn),6);
            invdmn = find(g.dc.dv_coni(invmn,kdc)<0);
            if ~isempty(invdmn)
               g.dc.dv_conr(invmn(invdmn),kdc) = zeros(length(invdmn),1);
            end
         end  
      end
   end
end