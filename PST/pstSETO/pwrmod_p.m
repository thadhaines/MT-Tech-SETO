function f = pwrmod_p(i,k,bus,flag)
% Syntax: f = pwrmod_p(i,k,bus,flag)
% Purpose: real-power modulation model 
%          with vectorized computation option
%          NOTE - load modulation bus must be declared as a
%                 non-conforming constant-power load bus
% Input: i - power modulation number
%            if i= 0, vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%
% Output: f is a dummy variable
%                    

% History (in reverse chronological order)
%
% Version:1
% Date:Feb 2015
% Author: Dan Trudnowski

% power modulation variables
global pwrmod_con n_pwrmod pwrmod_p_st dpwrmod_p_st pwrmod_p_sig pwrmod_idx bus_int load_con

f = 0;% dummy variable

jay = sqrt(-1);
if ~isempty(pwrmod_con)
   if flag == 0; % initialization
      if i~=0
          jj = bus_int(pwrmod_con(i,1));
          if (load_con(pwrmod_idx(i),2)==1 && load_con(pwrmod_idx(i),3)==1)
             pwrmod_p_st(i,1) = bus(jj,4);
          else
             pwrmod_p_st(i,1) = bus(jj,4)/abs(bus(jj,2));
          end
          clear jj
      else % vectorized calculation
          for ii=1:n_pwrmod
              jj = bus_int(pwrmod_con(ii,1));
              if (load_con(pwrmod_idx(ii),2)==1 && load_con(pwrmod_idx(ii),3)==1)
                 pwrmod_p_st(ii,1) = bus(jj,4);
              else
                 pwrmod_p_st(ii,1) = bus(jj,4)/abs(bus(jj,2));
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
         dpwrmod_p_st(i,k) = (-pwrmod_p_st(i,k)+pwrmod_p_sig(i,k))/pwrmod_con(i,2);
         % anti-windup reset
         if pwrmod_p_st(i,k) > pwrmod_con(i,3)
            if dpwrmod_p_st(i,k)>0
               dpwrmod_p_st(i,k) = 0;
            end
         end
         if pwrmod_p_st(i,k) < pwrmod_con(i,4)
            if dpwrmod_p_st(i,k)<0
               dpwrmod_p_st(i,k) = 0;
            end
         end
      else %vectorized computation
         dpwrmod_p_st(:,k) = (-pwrmod_p_st(:,k)+pwrmod_p_sig(:,k))./pwrmod_con(:,2);
         % anti-windup reset
         indmx =find(pwrmod_p_st(:,k) > pwrmod_con(:,3));
         if ~isempty(indmx)
            indrate = find(dpwrmod_p_st(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               dpwrmod_p_st(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(pwrmod_p_st(:,k) < pwrmod_con(:,4));
         if ~isempty(indmn)
            indrate = find(dpwrmod_p_st(indmn)<0);
            if ~isempty(indrate)
               % set rate to zero
               dpwrmod_p_st(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
