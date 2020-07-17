function pwrmod_p(i,k,bus,flag)
%PWRMOD_P is a real-power modulation model.
% PWRMOD_P is a is a real-power modulation model with 
% a vectorized computation option. Load modulation bus must be 
% declared as a non-conforming constant-power load bus.
%
% Syntax: pwrmod_p(i,k,bus,flag)
%
%   Input: 
%   i - power modulation number
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
%   02/xx/15    xx:xx   Dan Trudnowski  Version 1.0
%   07/02/20    14:15   Thad Haines     Revised format of globals and internal function documentation

global g

if ~isempty(g.pwr.pwrmod_con)
   if flag == 0; % initialization
      if i~=0
          jj = g.bus.bus_int(g.pwr.pwrmod_con(i,1));
          if (g.ncl.load_con(g.pwr.pwrmod_idx(i),2)==1 && g.ncl.load_con(g.pwr.pwrmod_idx(i),3)==1)
             g.pwr.pwrmod_p_st(i,1) = bus(jj,4);
          else
             g.pwr.pwrmod_p_st(i,1) = bus(jj,4)/abs(bus(jj,2));
          end
          clear jj
      else % vectorized calculation
          for ii=1:g.pwr.n_pwrmod
              jj = g.bus.bus_int(g.pwr.pwrmod_con(ii,1));
              if (g.ncl.load_con(g.pwr.pwrmod_idx(ii),2)==1 && g.ncl.load_con(g.pwr.pwrmod_idx(ii),3)==1)
                 g.pwr.pwrmod_p_st(ii,1) = bus(jj,4);
              else
                 g.pwr.pwrmod_p_st(ii,1) = bus(jj,4)/abs(bus(jj,2));
              end
          end
          clear jj
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end
   
   if flag == 2 %  dynamics calculation
      % for linearization with operating condition at limits,
      % additional code will be needed
      if i ~= 0
         g.pwr.dpwrmod_p_st(i,k) = (-g.pwr.pwrmod_p_st(i,k)+g.pwr.pwrmod_p_sig(i,k))/g.pwr.pwrmod_con(i,2);
         % anti-windup reset
         if g.pwr.pwrmod_p_st(i,k) > g.pwr.pwrmod_con(i,3)
            if g.pwr.dpwrmod_p_st(i,k)>0
               g.pwr.dpwrmod_p_st(i,k) = 0;
            end
         end
         if g.pwr.pwrmod_p_st(i,k) < g.pwr.pwrmod_con(i,4)
            if g.pwr.dpwrmod_p_st(i,k)<0
               g.pwr.dpwrmod_p_st(i,k) = 0;
            end
         end
      else %vectorized computation
         g.pwr.dpwrmod_p_st(:,k) = (-g.pwr.pwrmod_p_st(:,k)+g.pwr.pwrmod_p_sig(:,k))./g.pwr.pwrmod_con(:,2);
         % anti-windup reset
         indmx =find(g.pwr.pwrmod_p_st(:,k) > g.pwr.pwrmod_con(:,3));
         if ~isempty(indmx)
            indrate = find(g.pwr.dpwrmod_p_st(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               g.pwr.dpwrmod_p_st(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(g.pwr.pwrmod_p_st(:,k) < g.pwr.pwrmod_con(:,4));
         if ~isempty(indmn)
            indrate = find(g.pwr.dpwrmod_p_st(indmn)<0);
            if ~isempty(indrate)
               % set rate to zero
               g.pwr.dpwrmod_p_st(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
