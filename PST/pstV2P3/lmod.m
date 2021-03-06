function lmod(i,k,flag)
%LMOD Handles real load modulation
% LMOD handles the real load modulation model with a vectorized 
% computation option.
%
%   NOTES:
%   Load modulation bus must be declared as a non-conforming load bus.
%   Additional code will be needed for linearization with operating condition at limits.
%
%   Syntax: 
%   lmod(i,k,flag)
%
%   Inputs: 
%   i -     load modulation number. if i= 0, vectorized computation is used.
%   k -     integer time
%   flag -  0 - initialization
%           1 - network interface computation
%         	2 - generator dynamics computation
%
%   Output: 
%   VOID

%   Define Variables:
%
%
%   History:
%   Date        Time    Engineer        Description
%   08/15/97    16:58   Graham Rogers   Version 1
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   05/28/20    09:45   Thad Haines     Revised format internal function documentation

% system variables
global  basmva

% load modulation variables
global  lmod_con n_lmod lmod_idx
global lmod_pot lmod_st dlmod_st
global  lmod_sig

% new combined global
% global g



%f = 0;% dummy variable
%jay = sqrt(-1);

% modified g code
% if ~isempty(g.lmod.lmod_con)
%     if flag == 0 % initialization
%         if i~=0
%             % max modulation on system base
%             g.lmod.lmod_pot(i,1) = g.lmod.lmod_con(i,4)*g.lmod.lmod_con(i,3)/basmva;
%             % min modulation on system base
%             g.lmod.lmod_pot(i,2) = g.lmod.lmod_con(i,5)*g.lmod.lmod_con(i,3)/basmva;
%             g.lmod.lmod_st(i,1) = 0;
%         else % vectorized calculation
%             % max modulation on system base
%             g.lmod.lmod_pot(:,1) = g.lmod.lmod_con(:,4).*g.lmod.lmod_con(:,3)/basmva;
%             % min modulation on system base
%             g.lmod.lmod_pot(:,2) = g.lmod.lmod_con(:,5).*g.lmod.lmod_con(:,3)/basmva;
%             g.lmod.lmod_st(:,1) = zeros(g.lmod.n_lmod,1);
%         end
%     end
%     if flag == 1 % network interface computation
%         % no interface calculation required - done in nc_load
%     end
%     
%     if flag == 2 %  dynamics calculation
%         % for linearization with operating condition at limits,
%         % additional code will be needed
%         if i ~= 0
%             g.lmod.dlmod_st(i,k) = (-g.lmod.lmod_st(i,k)+g.lmod.lmod_con(i,6)*g.lmod.lmod_sig(i,k))/g.lmod.lmod_con(i,7);
%             % anti-windup reset
%             if g.lmod.lmod_st(i,k) > g.lmod.lmod_pot(i,1)
%                 if g.lmod.dlmod_st(i,k)>0
%                     g.lmod.dlmod_st(i,k) = 0;
%                 end
%             end
%             if g.lmod.lmod_st(i,k) < g.lmod.lmod_pot(i,2)
%                 if g.lmod.dlmod_st(i,k)<0
%                     g.lmod.dlmod_st(i,k) = 0; %corrected dmod_st 5/28/20 thad
%                 end
%             end
%         else %vectorized computation
%             g.lmod.dlmod_st(:,k) = (-g.lmod.lmod_st(:,k)+g.lmod.lmod_con(:,6).*g.lmod.lmod_sig(:,k))./g.lmod.lmod_con(:,7);
%             % anti-windup reset
%             indmx = find( g.lmod.lmod_st(:,k) > g.lmod.lmod_pot(:,1));
%             if ~isempty(indmx)
%                 indrate = find(g.lmod.dlmod_st(indmx,k)>0);
%                 if ~isempty(indrate)
%                     % set rate to zero
%                     g.lmod.dlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
%                 end
%             end
%             indmn = find(g.lmod.lmod_st(:,k) < g.lmod.lmod_pot(:,2));
%             if ~isempty(indmn)
%                 indrate = find(g.lmod.dlmod_st(indmn)<0);
%                 if ~isempty(indrate)
%                     % set rate to zero
%                     g.lmod.dlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
%                 end
%             end
%         end
%     end
% end

% Original code
f = 0;% dummy variable

jay = sqrt(-1);
if ~isempty(lmod_con)
    if flag == 0; % initialization
        if i~=0
            lmod_pot(i,1) = lmod_con(i,4)*lmod_con(i,3)/basmva;
            % max modulation on system base
            lmod_pot(i,2) = lmod_con(i,5)*lmod_con(i,3)/basmva;
            % min modulation on system base
            lmod_st(i,1) = 0;
        else % vectorized calculation
            lmod_pot(:,1) = lmod_con(:,4).*lmod_con(:,3)/basmva;
            % max modulation on system base
            lmod_pot(:,2) = lmod_con(:,5).*lmod_con(:,3)/basmva;
            % min modulation on system base
            lmod_st(:,1) = zeros(n_lmod,1);
        end
    end
    if flag == 1 % network interface computation
        % no interface calculation required - done in nc_load
    end
    
    if flag == 2 %  dynamics calculation
        % for linearization with operating condition at limits,
        % additional code will be needed
        if i ~= 0
            dlmod_st(i,k) = (-lmod_st(i,k)+lmod_con(i,6)*lmod_sig(i,k))/lmod_con(i,7);
            % anti-windup reset
            if lmod_st(i,k) > lmod_pot(i,1)
                if dlmod_st(i,k)>0
                    dlmod_st(i,k) = 0;
                end
            end
            if lmod_st(i,k) < lmod_pot(i,2)
                if dlmod_st(i,k)<0
                    dmod_st(i,k) = 0;
                end
            end
        else %vectorized computation
            dlmod_st(:,k) = (-lmod_st(:,k)+lmod_con(:,6).*lmod_sig(:,k))./lmod_con(:,7);
            % anti-windup reset
            indmx =find( lmod_st(:,k) > lmod_pot(:,1));
            if ~isempty(indmx)
                indrate = find(dlmod_st(indmx,k)>0);
                if ~isempty(indrate)
                    % set rate to zero
                    dlmod_st(indmx(indrate),k) = zeros(length(indrate),1);
                end
            end
            indmn = find(lmod_st(:,k) < lmod_pot(:,2));
            if ~isempty(indmn)
                indrate = find(dlmod_st(indmn)<0);
                if ~isempty(indrate)
                    % set rate to zero
                    dlmod_st(indmn(indrate),k) = zeros(length(indrate),1);
                end
            end
        end
    end
end
