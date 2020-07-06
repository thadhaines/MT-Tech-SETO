function smpexc(i,k,flag)
%SMPEXC is a simple excitation system, (exc_con(i,1)=0)
% SMPPI is a simple excitation system, (exc_con(i,1)=0) with vectorized
% computation option.
%
% Syntax: smpexc(i,k,flag) 
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
%   03/xx/91    xx:xx   Joe H. Chow     Version 1.0
%   (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%   01/29/99    08:51   -               -
%   06/xx/99    xx:xx   Graham Rogers   Version 2.0 Modified vectorization to allow different exciters.
%                                       note change in way exc_sig is referenced
%                                       i.e., in terms of the generator number and not the
%                                       exciter number
%   08/XX/97    XX:XX   Graham Rogers   Version 2.1 Reversed change in exc_sig referencing
%                                       Added pss_out to free exc_sig for other control functions
%   01/29/99    08:51   -               Error Correction in non-vector dynamics calculation
%                                       Efd(i,k) = exc_cin(i,9); changed to Efd(i,k) = exc_con(i,9);
%   06/19/20    10:34   Thad Haines     Revised format of globals and internal function documentation
%   07/06/20    14:09   Thad Haines     Completion of global g alterations

global g

if i ~= 0
  if g.exc.exc_con(i,1) ~= 0
    error('SMPEXC: inappropriate exciter model')
  end
end

if flag == 0 % initialization
   if i ~= 0  % scalar computation
     n = g.mac.mac_int(g.exc.exc_con(i,2)); % machine number
     g.exc.Efd(i,1) = g.mac.vex(n,1);
     g.exc.V_A(i,1) = g.exc.Efd(i,1)/g.exc.exc_con(i,4); % laglead
     g.exc.V_As(i,1) = g.exc.V_A(i,1); % leadlag state variable
     g.exc.V_TR(i,1)= g.mac.eterm(n,1); % input filter state
     err = g.exc.V_A(i,1); % summing junction error
     if g.exc.exc_con(i,6)~=0
       g.exc.exc_pot(i,5) = g.exc.exc_con(i,7)/g.exc.exc_con(i,6);
     else
       g.exc.exc_pot(i,5)=1;
     end
     g.exc.exc_pot(i,3) = g.mac.eterm(n,1)+err; % reference voltage
    
   else  % vectorized computation

     if g.exc.n_smp ~= 0
       n = g.mac.mac_int(g.exc.exc_con(g.exc.smp_idx,2)); % machine number with simple exciters
       g.exc.Efd(g.exc.smp_idx,1) = g.mac.vex(n,1);
       g.exc.V_A(g.exc.smp_idx,1) = g.exc.Efd(g.exc.smp_idx,1)./g.exc.exc_con(g.exc.smp_idx,4); % laglead
       g.exc.V_As(g.exc.smp_idx,1) = g.exc.V_A(g.exc.smp_idx,1); % leadlag state variable
       g.exc.V_TR(g.exc.smp_idx,1) = g.mac.eterm(n,1); % input filter state
       err = g.exc.V_A(g.exc.smp_idx,1); % summing junction error
       g.exc.exc_pot(g.exc.smp_idx,5) = ones(g.exc.n_smp,1);
       TB = g.exc.smp_TB_idx;
       if ~isempty(TB)
         g.exc.exc_pot(g.exc.smp_idx(TB),5) = g.exc.exc_con(g.exc.smp_idx(TB),7)./g.exc.exc_con(g.exc.smp_idx(TB),6);
       end
       g.exc.exc_pot(g.exc.smp_idx,3) = g.mac.eterm(n,1)+err; % reference voltage
     end
   end
end

if flag == 1 % network interface computation
   if i ~= 0 % scalar computation
     n = g.mac.mac_int(g.exc.exc_con(i,2)); % machine number
     if g.exc.exc_con(i,5)~=0
       g.mac.vex(n,k) = g.exc.Efd(i,k); % field voltage for machines
     else
       if g.exc.exc_con(i,6)==0
         g.mac.vex(n,k)=(g.exc.exc_sig(i,k)+g.exc.exc_pot(i,3)-g.exc.V_TR(i,k)...
                   + g.pss.pss_out(i,k))*g.exc.exc_con(i,4);
       else
         g.exc.V_A(i,k)=(g.exc.exc_sig(i,k)+g.exc.exc_pot(i,3)-g.exc.V_TR(i,k)+g.pss.pss_out(i,k))...
                  *g.exc.exc_pot(i,5) + (1.-g.exc.exc_pot(i,5))*g.exc.V_As(i,k);
         g.mac.vex(n,k) = g.exc.V_A(i,k)*g.exc.exc_con(i,4);
       end
     end 
   else      % vectorized computation

     if g.exc.n_smp ~=0 % check for any simple exciters
        n = g.mac.mac_int(g.exc.exc_con(g.exc.smp_idx,2)); % machine numbers for simple exciters
        % field voltage for machines having TA and TB  zero
        not_TATB = find((g.exc.smp_TA + g.exc.smp_TB)<0.001);
        if ~isempty(not_TATB)
          n_nTATB = n(not_TATB);% machine numbers for exciters with TA & TB zero
          g.mac.vex(n_nTATB,k) = (g.exc.exc_sig(g.exc.smp_idx(not_TATB),k)...
                           +g.exc.exc_pot(g.exc.smp_idx(not_TATB),3)-...
                           g.exc.V_TR(g.exc.smp_idx(not_TATB),k) + ...
                           g.pss.pss_out(g.exc.smp_idx(not_TATB),k))...
                           .*g.exc.exc_con(g.exc.smp_idx(not_TATB),4);
        end
       
 % field voltage for machines with TB non zero and TA zero      
        TB = find((g.exc.smp_TA<0.001)&(g.exc.smp_TB>0.001));
        if ~isempty(TB)   
           n_TB = n(TB);
           g.exc.V_A(g.exc.smp_idx(TB),k) = (g.exc.exc_sig(g.exc.smp_idx(TB),k)...
                                +g.exc.exc_pot(g.exc.smp_idx(TB),3)...
                                -g.exc.V_TR(g.exc.smp_idx(TB),k)...
                                +g.pss.pss_out(g.exc.smp_idx(TB),k))...
                                .*g.exc.exc_pot(g.exc.smp_idx(TB),5) + ...
                                (ones(length(TB),1)-...
                                g.exc.exc_pot(g.exc.smp_idx(TB),5)).*g.exc.V_As(g.exc.smp_idx(TB),k);
           g.mac.vex(n_TB,k) = g.exc.V_A(g.exc.smp_idx(TB),k).*g.exc.exc_con(g.exc.smp_idx(TB),4);
        end 

        % field voltage for non zero TA
        TA=g.exc.smp_TA_idx;
        if ~isempty(TA)
           n_TA = n(TA); % machine number 
           g.mac.vex(n_TA,k) = g.exc.Efd(g.exc.smp_idx(TA),k); % field voltage for machines with TA
        end 
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
     if g.exc.exc_con(i,6) == 0 % no leadlag
        g.exc.dV_As(i,k) = 0;
        g.exc.V_As(i,k) = err;
        g.exc.V_A(i,k) = err;
     else
        g.exc.dV_As(i,k) = (-g.exc.V_As(i,k)+err)/g.exc.exc_con(i,6);
        g.exc.V_A(i,k) = g.exc.exc_pot(i,5)*err + ...
                   (1-g.exc.exc_pot(i,5))*g.exc.V_As(i,k);
     end
     if g.exc.exc_con(i,5) == 0 % Ta zero
        g.exc.dEfd(i,k) = 0;
        g.exc.Efd(i,k) = g.exc.exc_con(i,4)*g.exc.V_A(i,k);
        g.exc.Efd(i,k) = max(g.exc.exc_con(i,9),...
                   min(g.exc.Efd(i,k),g.exc.exc_con(i,8))); % voltage limit
     else
        g.exc.dEfd(i,k) = (-g.exc.Efd(i,k)+g.exc.exc_con(i,4)*g.exc.V_A(i,k))...
                    /g.exc.exc_con(i,5);
        % anti-windup reset
        if g.exc.Efd(i,k) > g.exc.exc_con(i,8)
          g.exc.Efd(i,k) = g.exc.exc_con(i,8);
          if g.exc.dEfd(i,k)>0
            g.exc.dEfd(i,k) = 0;
          end
        end
        if g.exc.Efd(i,k) < g.exc.exc_con(i,9)
          g.exc.Efd(i,k) = g.exc.exc_con(i,9);
          if g.exc.dEfd(i,k) < 0
            g.exc.dEfd(i,k) =0;
          end
        end
     end
     g.exc.R_f(i,k) = 0; 
     g.exc.dR_f(i,k) = 0;
     g.exc.V_R(i,k) = 0; 
     g.exc.dV_R(i,k) = 0;

  else % vectorized computation
  
    if g.exc.n_smp~=0
         % machine numbers
         n = g.mac.mac_int(g.exc.exc_con(g.exc.smp_idx,2));
         TR = g.exc.smp_TR_idx;
         no_TR = g.exc.smp_noTR_idx;
         if ~isempty(no_TR) % some exciters have zero TR
           n_nTR = n(no_TR);
           g.exc.dV_TR(g.exc.smp_idx(no_TR),k) = zeros(length(no_TR),1);
           g.exc.V_TR(g.exc.smp_idx(no_TR),k) = g.mac.eterm(n_nTR,k);
         end
         if ~isempty(TR) %some exciters have nonzero TR
           n_TR = n(TR);
           g.exc.dV_TR(g.exc.smp_idx(TR),k) = (g.mac.eterm(n_TR,k)-...
                                g.exc.V_TR(g.exc.smp_idx(TR),k))...
                              ./g.exc.exc_con(g.exc.smp_idx(TR),3);
         end
         % error defined for all simple exciters
         err = g.exc.exc_sig(g.exc.smp_idx,k) ...
               + g.exc.exc_pot(g.exc.smp_idx,3)...
               - g.exc.V_TR(g.exc.smp_idx,k)...
               + g.pss.pss_out(g.exc.smp_idx,k);
         no_TB = g.exc.smp_noTB_idx;
         if ~isempty(no_TB)
           g.exc.dV_As(g.exc.smp_idx(no_TB),k) = zeros(length(no_TB),1);
           g.exc.V_As(g.exc.smp_idx(no_TB),k) = err(no_TB);
           g.exc.V_A(g.exc.smp_idx(no_TB),k) = err(no_TB);
         end
         TB = g.exc.smp_TB_idx;
         if ~isempty(TB)
           g.exc.dV_As(g.exc.smp_idx(TB),k) = (-g.exc.V_As(g.exc.smp_idx(TB),k)+err(TB))...
                                 ./g.exc.exc_con(g.exc.smp_idx(TB),6);
           g.exc.V_A(g.exc.smp_idx(TB),k) = g.exc.exc_pot(g.exc.smp_idx(TB),5).*err(TB) ...
                     +(ones(length(TB),1) ...
                     -g.exc.exc_pot(g.exc.smp_idx(TB),5)).*g.exc.V_As(g.exc.smp_idx(TB),k);
         end
         no_TA = g.exc.smp_noTA_idx;
         if ~isempty(no_TA)
           g.exc.dEfd(g.exc.smp_idx(no_TA),k) = zeros(length(no_TA),1);
           g.exc.Efd(g.exc.smp_idx(no_TA),k) = g.exc.exc_con(g.exc.smp_idx(no_TA),4)...
                                  .*g.exc.V_A(g.exc.smp_idx(no_TA),k);
           % apply output limits
           g.exc.Efd(g.exc.smp_idx(no_TA),k) = min(g.exc.exc_con(g.exc.smp_idx(no_TA),8)...
                                  ,max(g.exc.Efd(g.exc.smp_idx(no_TA),k)...
                                  ,g.exc.exc_con(g.exc.smp_idx(no_TA),9)));
         end
         TA = g.exc.smp_TA_idx;
         if ~isempty(TA)
           g.exc.dEfd(g.exc.smp_idx(TA),k) = (-g.exc.Efd(g.exc.smp_idx(TA),k)...
                    + g.exc.exc_con(g.exc.smp_idx(TA),4)...
                               .*g.exc.V_A(g.exc.smp_idx(TA),k))...
                               ./g.exc.exc_con(g.exc.smp_idx(TA),5);
           %apply non-windup limits
           TA_max = find(g.exc.Efd(g.exc.smp_idx(TA),k)> g.exc.exc_con(g.exc.smp_idx(TA),8));
           TA_min = find(g.exc.Efd(g.exc.smp_idx(TA),k)< g.exc.exc_con(g.exc.smp_idx(TA),9));
           if ~isempty(TA_max)
             g.exc.Efd(g.exc.smp_idx(TA(TA_max)),k) = g.exc.exc_con(g.exc.smp_idx(TA(TA_max)),8);
             dTA_max = find(g.exc.dEfd(g.exc.smp_idx(TA(TA_max)),k)>0);
             if ~isempty(dTA_max)
               n_dTA = length(dTA_max);
               g.exc.dEfd(g.exc.smp_idx(TA(TA_max(dTA_max))),k) = zeros(n_dTA,1);
             end
           end
           if ~isempty(TA_min)
             g.exc.Efd(g.exc.smp_idx(TA(TA_min)),k) = g.exc.exc_con(g.exc.smp_idx(TA(TA_min)),9);
             dTA_min = find(g.exc.dEfd(g.exc.smp_idx(TA(TA_min)),k)<0);
             if ~isempty(dTA_min)
               n_dTA = length(dTA_min);
               g.exc.dEfd(g.exc.smp_idx(TA(TA_min(dTA_min))),k) = zeros(n_dTA,1);
             end
           end
         end
         g.exc.R_f(g.exc.smp_idx,k) = zeros(g.exc.n_smp,1); 
         g.exc.dR_f(g.exc.smp_idx,k) = zeros(g.exc.n_smp,1);
         g.exc.V_R(g.exc.smp_idx,k) =zeros(g.exc.n_smp,1); 
         g.exc.dV_R(g.exc.smp_idx,k) = zeros(g.exc.n_smp,1);
    end
  end
end
