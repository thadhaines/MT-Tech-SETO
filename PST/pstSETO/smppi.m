function smppi(i,k,flag)
%SMPPI Simple excitation system with pi avr
% SMPPI is a simple excitation system with pi avr, (exc_con(i,1)=4)
% with vectorized computation option.
%
% Syntax: smppi(i,k,flag) 
%
%   Input: 
%   i - generator number
%           0 for vectorized computation
%   k - integer time (data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation 
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   01/29/99    08:51   -               -
%   09/xx/99    xx:xx   Graham Rogers   Version 1.0
%   (c) Copyright 1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved
%   06/19/20    10:34   Thad Haines     Revised format of globals and internal function documentation
%   07/06/20    14:09   Thad Haines     Completion of global g alterations

global g

if i ~= 0
   if g.exc.exc_con(i,1) ~= 4
      error('SMPPI: inappropriate exciter model')
   end
end

if flag == 0 % initialization
   if i ~= 0  % scalar computation
      n = g.mac.mac_int(g.exc.exc_con(i,2)); % machine number
      g.exc.Efd(i,1) = g.mac.vex(n,1);
      g.exc.V_As(i,1) = g.exc.Efd(i,1); % integrator state
      g.exc.V_TR(i,1)=g.mac.eterm(n,1); % input filter state
      %err = 0; % summing junction error -commented out, not used -thad
      g.exc.exc_pot(i,3) = g.mac.eterm(n,1); % reference voltage
   else  % vectorized computation
      if g.exc.n_smppi ~= 0
         n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2)); % machine number with simple exciters
         g.exc.Efd(g.exc.smppi_idx,1) = g.mac.vex(n,1);
         g.exc.V_As(g.exc.smppi_idx,1) = g.exc.Efd(g.exc.smppi_idx,1); % integrator
         g.exc.V_TR(g.exc.smppi_idx,1) = g.mac.eterm(n,1); % input filter state
         g.exc.exc_pot(g.exc.smppi_idx,3) = g.mac.eterm(n,1); % reference voltage
      end
   end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
      n = g.mac.mac_int(g.exc.exc_con(i,2)); % machine number
      g.mac.vex(n,k) = g.exc.Efd(i,k); % field voltage for machines
   else      % vectorized computation
      if g.exc.n_smppi ~=0 % check for any simple pi exciters
         n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2)); % machine numbers for simple exciters
         g.mac.vex(n,k) = g.exc.Efd(g.exc.smppi_idx,k);
      end
   end
end

if flag == 2 % exciter dynamics calculation
   if i ~= 0 % scalar computation
      n = g.mac.mac_int(g.exc.exc_con(i,2)); % machine number
      if g.exc.exc_con(i,3)== 0 %no input filter
         g.exc.dV_TR(i,k) = 0;
         g.exc.V_TR(i,k)=g.mac.eterm(n,k);
      else
         g.exc.dV_TR(i,k)=(g.mac.eterm(n,k)-g.exc.V_TR(i,k))/g.exc.exc_con(i,3);
      end
      err = g.exc.exc_sig(i,k)+g.exc.exc_pot(i,3)-g.exc.V_TR(i,k)...
         + g.pss.pss_out(i,k);
      g.exc.dV_As(i,k) = err*g.exc.exc_con(i,4);
      g.exc.dEfd(i,k) = (-g.exc.Efd(i,k)+g.exc.V_A(i,k)+g.exc.exc_con(i,6)*err)...
         /g.exc.exc_con(i,5);
      % anti-windup reset
      if g.exc.Efd(i,k) > g.exc.exc_con(i,8)
         g.exc.Efd(i,k) = g.exc.exc_con(i,8);
         if g.exc.dEfd(i,k)>0 % handles windup by setting derivative to zero... -thad 06/17/20
            g.exc.dEfd(i,k) = 0;
         end
      end
      if g.exc.Efd(i,k) < g.exc.exc_con(i,9)
         g.exc.Efd(i,k) = g.exc.exc_con(i,9);
         if g.exc.dEfd(i,k) < 0
            g.exc.dEfd(i,k) =0;
         end
      end
      g.exc.R_f(i,k) = 0; 
      g.exc.dR_f(i,k) = 0;
      g.exc.V_R(i,k) = 0; 
      g.exc.dV_R(i,k) = 0;
  
      
   else % vectorized computation
      
      if g.exc.n_smppi~=0
         % machine numbers
         n = g.mac.mac_int(g.exc.exc_con(g.exc.smppi_idx,2));
         TR = g.exc.smppi_TR_idx;
         no_TR = g.exc.smppi_noTR_idx;
         if ~isempty(no_TR) % some exciters have zero TR
            n_nTR = n(no_TR);
            g.exc.dV_TR(g.exc.smppi_idx(no_TR),k)=zeros(length(no_TR),1);
            g.exc.V_TR(g.exc.smppi_idx(no_TR),k)=g.mac.eterm(n_nTR,k);
         end
         if ~isempty(TR) %some exciters have nonzero TR
            n_TR = n(TR);
            g.exc.dV_TR(g.exc.smppi_idx(TR),k)=(g.mac.eterm(n_TR,k)-g.exc.V_TR(g.exc.smppi_idx(TR),k))...
               ./g.exc.exc_con(g.exc.smppi_idx(TR),3);
         end
         % error defined for all simple exciters
         err = g.exc.exc_sig(g.exc.smppi_idx,k)+g.exc.exc_pot(g.exc.smppi_idx,3)...
             -g.exc.V_TR(g.exc.smppi_idx,k) + g.pss.pss_out(g.exc.smppi_idx,k);
         g.exc.dV_As(g.exc.smppi_idx,k) = err*g.exc.exc_con(g.exc.smppi_idx,4);
         g.exc.dEfd(g.exc.smppi_idx,k) = (-g.exc.Efd(g.exc.smppi_idx,k)...
             +g.exc.V_As(g.exc.smppi_idx,k)+g.exc.exc_con(g.exc.smppi_idx,6)*err)...
            ./g.exc.exc_con(g.exc.smppi_idx,5);
         TA_max = find(g.exc.Efd(g.exc.smppi_idx,k)>g.exc.exc_con(g.exc.smppi_idx,8));
         TA_min = find(g.exc.Efd(g.exc.smppi_idx,k)<g.exc.exc_con(g.exc.smppi_idx,9));
         if ~isempty(TA_max)
            g.exc.Efd(g.exc.smppi_idx(TA_max),k) = g.exc.exc_con(g.exc.smppi_idx(TA_max),8);
            dTA_max = find(g.exc.dEfd(g.exc.smppi_idx(TA_max),k)>0);
            if ~isempty(dTA_max)
               n_dTA = length(dTA_max);
               g.exc.dEfd(g.exc.smppi_idx(TA_max(dTA_max)),k) = zeros(n_dTA,1);
            end
         end
         if ~isempty(TA_min)
            g.exc.Efd(g.exc.smppi_idx(TA_min),k) = g.exc.exc_con(g.exc.smppi_idx(TA_min),9);
            dTA_min = find(g.exc.dEfd(g.exc.smppi_idx(TA_min),k)<0);
            if ~isempty(dTA_min)
               n_dTA = length(dTA_min);
               g.exc.dEfd(g.exc.smppi_idx(TA_min(dTA_min)),k) = zeros(n_dTA,1);
            end
         end
         g.exc.R_f(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
         g.exc.dR_f(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
         g.exc.V_R(g.exc.smppi_idx,k) =zeros(g.exc.n_smppi,1); 
         g.exc.dV_R(g.exc.smppi_idx,k) = zeros(g.exc.n_smppi,1);
      end
   end
end
