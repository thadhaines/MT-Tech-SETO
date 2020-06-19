function mac_indx()
%MAC_INDX Forms indexes for the machine models
% MAC_INDX Forms indexes for the machine models to allow vector computation
% with different generator types
%
% Syntax: mac_indx()
%
%   Input:
%   VOID
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   12/16/97    09:38   -               Version 1
%   06/19/20    09:56   Thad Haines     Revised format of globals and internal function documentation


%global mac_pot mac_con mac_int ibus_con n_mac n_em n_tra n_sub n_ib
%global mac_em_idx mac_tra_idx mac_sub_idx mac_ib_idx not_ib_idx
%global mac_ib_em mac_ib_tra mac_ib_sub n_ib_em n_ib_tra n_ib_sub

global n_ivm mac_ivm_idx
global ind_int ind_con n_mot
global igen_int igen_con n_ig

global g

% insert the default values of the %percentage p and q
[g.mac.n_mac, n_par] = size(g.mac.mac_con);
g.mac.mac_pot = zeros(g.mac.n_mac,15);
if n_par < 22
    g.mac.mac_con(:,22:23) = ones(g.mac.n_mac,2);
end
pqpc_idx = find(g.mac.mac_con(:,22)==0&g.mac.mac_con(:,23)==0);
if ~isempty(pqpc_idx)
    g.mac.mac_con(pqpc_idx,22:23)=ones(length(pqpc_idx),2);
end
% set up internal machine list
macmax = max(g.mac.mac_con(:,1));
g.mac.mac_int = zeros(macmax,1);
g.mac.mac_int(round(g.mac.mac_con(:,1))) = 1:g.mac.n_mac;
n_tot = g.mac.n_mac;
ngm = g.mac.n_mac;
n_mot = 0;
n_ig = 0;
if ~isempty(ind_con)
    n_mot = length(ind_con(:,1));
    n_tot = g.mac.n_mac + n_mot;
    ngm = n_tot;
    motmax= max(ind_con(:,1));
    ind_int = zeros(motmax,1);
    ind_int(round(ind_con(:,1)))= g.mac.n_mac+1:n_tot;
end

if ~isempty(igen_con)
    n_ig = length(igen_con(:,1));
    n_tot = n_tot + n_ig;
    igmax= max(igen_con(:,1));
    igen_int = zeros(igmax,1);
    igen_int(round(igen_con(:,1)))=ngm+1:n_tot;
end


% check for types of generators
% infinite buses
g.mac.n_ib = 0;
g.mac.n_ib_em = 0;
g.mac.n_ib_tra = 0;
g.mac.n_ib_sub = 0;
g.mac.not_ib_idx = (1:g.mac.n_mac)';% sets default to all generators not infinite buses
if ~isempty(g.mac.ibus_con)
    g.mac.mac_ib_idx = find(g.mac.ibus_con==1);
    g.mac.not_ib_idx = find(g.mac.ibus_con==0);
    g.mac.n_ib = length(g.mac.mac_ib_idx);
else
    g.mac.mac_ib_idx = [];
end
%em has no xd or xq
g.mac.mac_em_idx = find((g.mac.mac_con(:,6)==0) & (g.mac.mac_con(:,16)~=0));
if ~isempty(g.mac.mac_em_idx)
    g.mac.n_em = length(g.mac.mac_em_idx);
else
    g.mac.n_em = 0;
end
%tra has no xdpp
g.mac.mac_tra_idx = find((g.mac.mac_con(:,6)~=0) & ...
    (g.mac.mac_con(:,8)==0) & (g.mac.mac_con(:,16)~=0));
if ~isempty(g.mac.mac_tra_idx)
    g.mac.n_tra = length(g.mac.mac_tra_idx);
else
    g.mac.n_tra = 0;
end
%sub has xdpp
g.mac.mac_sub_idx = find((g.mac.mac_con(:,8)~=0) & (g.mac.mac_con(:,16)~=0));
if ~isempty(g.mac.mac_sub_idx)
    g.mac.n_sub = length(g.mac.mac_sub_idx);
else
    g.mac.n_sub = 0;
end
%IVM generators
mac_ivm_idx = find(g.mac.mac_con(:,16)==0);
if ~isempty(mac_ivm_idx)
    n_ivm = length(mac_ivm_idx);
else
    n_ivm = 0;
end
if g.mac.n_ib~=0
    ib_em = find(g.mac.mac_con(g.mac.mac_ib_idx,6)==0);
    if ~isempty(ib_em)
        g.mac.mac_ib_em = g.mac.mac_ib_idx(ib_em);
        g.mac.n_ib_em = length(ib_em);
    end
    ib_tra = find((g.mac.mac_con(g.mac.mac_ib_idx,6)~=0)&(g.mac.mac_con(g.mac.mac_ib_idx,8)==0));
    if ~isempty(ib_tra)
        g.mac.mac_ib_tra = g.mac.mac_ib_idx(ib_tra);
        g.mac.n_ib_tra = length(ib_tra);
    end
    ib_sub = find(g.mac.mac_con(g.mac.mac_ib_idx,8)~=0);
    if ~isempty(ib_sub)
        g.mac.mac_ib_sub = g.mac.mac_ib_idx(ib_sub);
        ib_sub = length(g.mac.mac_ib_sub);
    end
end