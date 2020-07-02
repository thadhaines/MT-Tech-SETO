function rlmod(i,k,flag)
%RLMOD Handles reactive load modulation
% RLMOD Handles reactive load modulation model with vectorized 
% computation option.
%
%   NOTES:
%   The reactive load modulation bus must be declared as a 
%   non-conforming load bus.
%
%   Syntax: 
%   rlmod(i,k,flag)
%
%   Input: 
%   i - reactive load modulation number
%            if i= 0, vectorized computation
%	k - integer time
%	flag -  0 - initialization
%         	1 - network interface computation
%         	2 - generator dynamics computation
%
%   Output:
%   VOID
%                    
%   Define Variables:
%
%   History:
%   Date        Time    Engineer        Description
%   08/27/97    17:20   Graham Rogers   Version 1
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   06/15/20    09:45   Thad Haines     Revised format of globals and internal function documentation

global g

if ~isempty(g.rlmod.rlmod_con)
   if flag == 0 % initialization
      if i~=0
         g.rlmod.rlmod_pot(i,1) = g.rlmod.rlmod_con(i,4)*g.rlmod.rlmod_con(i,3)/g.sys.basmva;
         % max modulation on system base
         g.rlmod.rlmod_pot(i,2) = g.rlmod.rlmod_con(i,5)*g.rlmod.rlmod_con(i,3)/g.sys.basmva;
         % min modulation on system base
         g.rlmod.rlmod_st(i,1) = 0;
      else % vectorized calculation
         g.rlmod.rlmod_pot(:,1) = g.rlmod.rlmod_con(:,4).*g.rlmod.rlmod_con(:,3)/g.sys.basmva;
         % max modulation on system base
         g.rlmod.rlmod_pot(:,2) = g.rlmod.rlmod_con(:,5).*g.rlmod.rlmod_con(:,3)/g.sys.basmva;
         % min modulation on system base
         g.rlmod.rlmod_st(:,1) = zeros(g.rlmod.n_rlmod,1);
      end
   end
   if flag == 1 % network interface computation
      % no interface calculation required - done in nc_load
   end
   
   if flag == 2 %  dynamics calculation
      % for linearization with operating condition at limits,
      % additional code will be needed
      if i ~= 0
         g.rlmod.drlmod_st(i,k) = (-g.rlmod.rlmod_st(i,k)+ ... 
             g.rlmod.rlmod_con(i,6)*g.rlmod.rlmod_sig(i,k))/g.rlmod.rlmod_con(i,7);
         % anti-windup reset
         if g.rlmod.rlmod_st(i,k) > g.rlmod.rlmod_pot(i,1)
            if g.rlmod.drlmod_st(i,k)>0
               g.rlmod.drlmod_st(i,k) = 0;
            end
         end
         if g.rlmod.rlmod_st(i,k) < g.rlmod.rlmod_pot(i,2)
            if g.rlmod.drlmod_st(i,k)<0
               g.rlmod.drlmod_st(i,k) = 0;
            end
         end
      else %vectorized computation
         g.rlmod.drlmod_st(:,k) = (-g.rlmod.rlmod_st(:,k)+ ...
             g.rlmod.rlmod_con(:,6).*g.rlmod.rlmod_sig(:,k))./g.rlmod.rlmod_con(:,7);
         % anti-windup reset
         indmx =find( g.rlmod.rlmod_st(:,k) > g.rlmod.rlmod_pot(:,1));
         if ~isempty(indmx)
            g.rlmod.rlmod_st(indmx,k) = g.rlmod.rlmod_pot(indmx,1);
            indrate = find(g.rlmod.drlmod_st(indmx,k)>0);
            if ~isempty(indrate)
               % set rate to zero
               g.rlmod.drlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
            end
         end
         indmn = find(g.rlmod.rlmod_st(:,k) < g.rlmod.rlmod_pot(:,2));
         if ~isempty(indmn)
            g.rlmod.rlmod_st(indmn,k) = g.rlmod.rlmod_pot(indmn,2); % Fixed rl_mod_st to rlmod_st?
            indrate = find(g.rlmod.drlmod_st(indmn)<0);
            if ~isempty(indrate)
               % set rate to zero
               g.rlmod.drlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
            end
         end
      end
   end
end
