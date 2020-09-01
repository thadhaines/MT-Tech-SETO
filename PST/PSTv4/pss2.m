function pss(i,k,flag)
%PSS is a power system stabilization model
% PSS is a power system stabilization model
%
% Syntax: pss(i,k,flag)
%
%   Input:
%   i - generator number
%           0 for vectorized computation
%   k - integer time (data index)
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   04/xx/92    xx:xx   Joe H. Chow     Version 1.0
%   (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%   06/xx/96    xx:xx   Graham Rogers   Version 2.0 modified vector code
%                                       and added negative output limit to
%                                       allow vector caclutaion on only some units.
%   10/xx/97    xx:xx   Graham Rogers   Version 2.1 - changed length to is empty in index checks.
%   07/xx/98    xx:xx   xxx             Version 2.2 - Added interface to deltaP/omega filter
%   04/xx/11    xx:xx   JHC             Version 3.1 - Corrected washout filter calculation
%   07/06/20    13:57   Thad Haines     Revised format of globals and internal function documentation
%   07/09/20    08:55   Thad Haines     'Un-corrected' JCH washout filter calculation so pss is same as model from V2
%   09/01/20    10:08   Thad Haines     Modified non-vector initialize to use k as data index

global g
global dpw_out n_dpw % deals with delta omega filter...

if g.pss.n_pss~=0
    if i ~= 0
        if g.pss.pss_con(i,1) ~= 1 && g.pss.pss_con(i,1) ~= 2
            error('PSS: inappropriate power system stablizer model')
        end
    end
    
    if flag == 0 % initialization
         if i ~= 0  % scalar computation
            n = g.pss.pss_mb_idx(i); % machine number
            if g.pss.pss_con(i,1) == 1 
                g.pss.pss1(i,k) = g.mac.mac_spd(n,k);% speed ref
            else
                g.pss.pss1(i,k) = g.mac.pelect(n,k)*g.sys.basmva/g.mac.mac_con(n,3); % power ref
            end
            if n_dpw ~= 0
                i_dpw = find(dpw_pss_idx==i);
                if ~isempty(i_dpw)
                    g.pss.pss1(i,k)= dpw_out(i_dpw,k);
                end
            end
            g.pss.pss2(i,k) = 0.;
            g.pss.pss3(i,k) = 0.;
            g.pss.pss_out(g.pss.pss_exc_idx(i),k) = 0.;
            % constant gains
            g.pss.pss_pot(i,1) = g.pss.pss_con(i,5)/g.pss.pss_con(i,6);
            g.pss.pss_pot(i,2) = 1.0;
            if g.pss.pss_con(i,8) ~= 0
                g.pss.pss_pot(i,2) = g.pss.pss_con(i,7)/g.pss.pss_con(i,8);
            end
        else
            
            %vector computation
            g.pss.pss_pot=ones(g.pss.n_pss,2);
            n= g.pss.pss_mb_idx;
            if ~isempty(g.pss.pss_sp_idx)
                n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                g.pss.pss1(g.pss.pss_sp_idx,1)=g.mac.mac_spd(n_sp,1);
            end
            if ~isempty(g.pss.pss_p_idx)
                n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                g.pss.pss1(g.pss.pss_p_idx,1)=g.mac.pelect(n_p,1)*g.sys.basmva./g.mac.mac_con(n_p,3);
            end
            if n_dpw ~=0
                g.pss.pss1(dpw_pss_idx,1) = dpw_out(:,1);
            end
            g.pss.pss2(g.pss.pss_idx,1)=zeros(g.pss.n_pss,1);
            g.pss.pss3(g.pss.pss_idx,1)=zeros(g.pss.n_pss,1);
            g.pss.pss_out(g.pss.pss_exc_idx,1)=zeros(g.pss.n_pss,1);
            g.pss.pss_pot(:,1)=g.pss.pss_con(:,5)./g.pss.pss_con(:,6);
            if ~isempty(g.pss.pss_T4_idx)
                g.pss.pss_pot(g.pss.pss_T4_idx,2)=g.pss.pss_con(g.pss.pss_T4_idx,7)./g.pss.pss_T4(g.pss.pss_T4_idx);
            end
        end
    end
    
    if flag == 1 % network interface computation
        if i ~= 0 % scalar computation
            n = g.pss.pss_mb_idx(i); % machine number
            if g.pss.pss_con(i,1) == 1
                var1 = g.mac.mac_spd(i,k)-g.pss.pss1(i,k)/g.pss.pss_con(i,4); % DO divide by pss_con(i,4)
            else
                n = g.mac.mac_int(g.pss.pss_con(i,2)); % machine number
                var1 = g.mac.pelect(i,k)*g.sys.basmva/g.mac.mac_con(n,3)-g.pss.pss1(i,k)/g.pss.pss_con(i,4); % DO divide by pss_con(i,4)
            end
            if n_dpw~=0
                if n_dpw ~= 0
                    i_dpw = find(dpw_pss_idx==i);
                    if ~isempty(i_dpw)
                        var1 = (dpw_out(i_dpw,k)-g.pss.pss1(i,k))/g.pss.pss_con(i,4); % DO divide by pss_con(i,4)
                    end
                end
            end
            var2 = g.pss.pss_pot(i,1)*g.pss.pss_con(i,3)*var1 + g.pss.pss2(i,k);
            
            if g.pss.pss_con(i,8) == 0
                var3 = var2;
            else
                var3 = g.pss.pss_pot(i,2)*var2 + g.pss.pss3(i,k);
            end
            g.pss.pss_out(g.pss.pss_exc_idx(i),k) = min(g.pss.pss_con(i,9),max(var3,-g.pss.pss_con(i,9)));
        else
            % vector computation
            if g.pss.n_pss~=0
                n = g.pss.pss_mb_idx; % machine number vector
                
                var1 = zeros(g.pss.n_pss,1);
                var2 = var1;
                var3 = var1;
                if length(g.pss.pss_sp_idx)~=0
                    n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                    var1(g.pss.pss_sp_idx) = g.mac.mac_spd(n_sp,k)-g.pss.pss1(g.pss.pss_sp_idx,k)./g.pss.pss_con(g.pss.pss_sp_idx,4); % DO divide by pss_con(pss_sp_idx,4)
                end
                if ~isempty(g.pss.pss_p_idx)
                    n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                    var1(g.pss.pss_p_idx) = (g.mac.pelect(n_p,k)*g.sys.basmva./g.mac.mac_con(n_p,3)-g.pss.pss1(g.pss.pss_p_idx,k))./g.pss.pss_con(g.pss.pss_p_idx,4); % DO divide by pss_con(pss_p_idx,4)
                end
                if n_dpw ~= 0
                    var1 = (dpw_out(:,k)-g.pss.pss1(dpw_pss_idx,k))./g.pss.pss_con(i,4); % DO divide by pss_con(i,4)
                end
            end
            
            var2(g.pss.pss_idx) = g.pss.pss_pot(g.pss.pss_idx,1).*(g.pss.pss_con(g.pss.pss_idx,3).*var1) + g.pss.pss2(g.pss.pss_idx,k);
            var3 = var2;
            
            if ~isempty(g.pss.pss_T4_idx)
                var3(g.pss.pss_T4_idx,1) = g.pss.pss_pot(g.pss.pss_T4_idx,2).*var2(g.pss.pss_T4_idx,1)...
                    + g.pss.pss3(g.pss.pss_T4_idx,k);
                
            end
            
            g.pss.pss_out(g.pss.pss_exc_idx,k) = min(g.pss.pss_con(g.pss.pss_idx,9),max(var3,g.pss.pss_con(g.pss.pss_idx,10)));
        end
    end
    
    if flag == 2 % pss dynamics calculation
        if i ~= 0 % scalar computation
            n = g.pss.pss_mb_idx(i); % machine number
            if g.pss.pss_con(i,1) == 1
                var1 = (g.mac.mac_spd(i,k)-g.pss.pss1(i,k))/g.pss.pss_con(i,4);  % DO divide by pss_con(i,4)
            else
                n = g.mac.mac_int(g.pss.pss_con(i,2)); % machine number
                var1 = (g.mac.pelect(i,k)*g.sys.basmva./g.mac.mac_con(n,3)-g.pss.pss1(i,k))/g.pss.pss_con(i,4); % DO divide by pss_con(i,4)
            end
            if n_dpw~=0
                if n_dpw ~= 0
                    i_dpw = find(dpw_pss_idx==i);
                    if ~isempty(i_dpw)
                        var1 = (dpw_out(i_dpw,k)-g.pss.pss1(i,k))/g.pss.pss_con(i,4);
                    end
                end
            end
            
            g.pss.dpss1(i,k) = var1;   % do NOT divide by pss_con(i,4)
            
            var2 = g.pss.pss_pot(i,1)*g.pss.pss_con(i,3)*var1 + g.pss.pss2(i,k);
            g.pss.dpss2(i,k) = ((1-g.pss.pss_pot(i,1))*g.pss.pss_con(i,3)*var1 - g.pss.pss2(i,k))/g.pss.pss_con(i,6);
            
            if g.pss.pss_con(i,8) == 0
                var3 = var2;
                g.pss.dpss3(i,k) = g.pss.dpss2(i,k);
            else
                var3 = g.pss.pss_pot(i,2)*var2 + g.pss.pss3(i,k);
                g.pss.dpss3(i,k) = ((1-g.pss.pss_pot(i,2))*var2 - g.pss.pss3(i,k))/g.pss.pss_con(i,8);
            end
            g.pss.pss_out(g.pss.pss_exc_idx(i),k) = min(g.pss.pss_con(i,9),max(var3,-g.pss.pss_con(i,9)));
            
        else
            
            % vector computation
            if g.pss.n_pss~=0
                n = g.pss.pss_mb_idx; % machine number vector
                var1 = zeros(g.pss.n_pss,1);
                var2 = var1;
                var3 = var1;
                if ~isempty(g.pss.pss_sp_idx)
                    n_sp = g.mac.mac_int(g.pss.pss_con(g.pss.pss_sp_idx,2));
                    var1(g.pss.pss_sp_idx) = (g.mac.mac_spd(n_sp,k)-g.pss.pss1(g.pss.pss_sp_idx,k))./g.pss.pss_con(g.pss.pss_sp_idx,4); % DO divide by pss_con(pss_sp_idx,4)
                end
                if ~isempty(g.pss.pss_p_idx)
                    n_p = g.mac.mac_int(g.pss.pss_con(g.pss.pss_p_idx,2));
                    var1(g.pss.pss_p_idx) = (g.mac.pelect(n_p,k)*g.sys.basmva./g.mac.mac_con(n_p,3)...
                        -g.pss.pss1(g.pss.pss_p_idx,k))./g.pss.pss_con(g.pss.pss_sp_idx,4); % DO divide by pss_con(pss_sp_idx,4)
                end
                if n_dpw ~= 0
                    var1 = (dpw_out(:,k)-g.pss.pss1(dpw_pss_idx,k))./g.pss.pss_con(dpw_pss_idx,4);
                end
            end
            
            g.pss.dpss1(g.pss.pss_idx,k) = var1;%./g.pss.pss_con(g.pss.pss_sp_idx,4); % DO NOT divide by pss_con(pss_sp_idx,4)
            
            var2 = g.pss.pss_pot(g.pss.pss_idx,1).*(g.pss.pss_con(g.pss.pss_idx,3).*var1) + g.pss.pss2(g.pss.pss_idx,k);
            g.pss.dpss2(g.pss.pss_idx,k) = ((ones(g.pss.n_pss,1)-g.pss.pss_pot(g.pss.pss_idx,1))...
                .*(g.pss.pss_con(g.pss.pss_idx,3).*var1 )...
                - g.pss.pss2(g.pss.pss_idx,k))./g.pss.pss_con(g.pss.pss_idx,6);
            
            var3 = var2;
            g.pss.dpss3(:,k) = g.pss.dpss2(:,k);
            if ~isempty(g.pss.pss_T4_idx)
                var3(g.pss.pss_T4_idx) = g.pss.pss_pot(g.pss.pss_T4_idx,2).*var2(g.pss.pss_T4_idx)...
                    + g.pss.pss3(g.pss.pss_T4_idx,k);
                g.pss.dpss3(g.pss.pss_T4_idx,k) = ((ones(length(g.pss.pss_T4_idx),1)...
                    -g.pss.pss_pot(g.pss.pss_T4_idx,2)).*var2(g.pss.pss_T4_idx)...
                    - g.pss.pss3(g.pss.pss_T4_idx,k))./g.pss.pss_T4(g.pss.pss_T4_idx);
            end
            g.pss.pss_out(g.pss.pss_exc_idx,k) = min(g.pss.pss_con(g.pss.pss_idx,9),max(var3,g.pss.pss_con(g.pss.pss_idx,10)));
        end
    end
end
