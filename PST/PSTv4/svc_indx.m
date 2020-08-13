function svc_indx()
% syntax: f = svc_indx
% 12:00 PM 8/8/97
% determines the relationship between svc and nc loads
% checks for svc
% determines number of SVCs
% checks for user defined damping controls

global g

g.svc.n_svc = 0;
g.svc.svc_idx = [];
g.svc.dcud_idx = [];
g.svc.n_dcud = 0;

if ~isempty(g.svc.svc_con)
    [g.svc.n_svc, npar] = size(g.svc.svc_con);
    g.svc.svc_idx = zeros(g.svc.n_svc,1);
    % set defaults for lead lag
    if npar<9
        g.svc.svc_con(:,8:9) = zeros(g.svc.n_svc,2);
    end
    g.svc.svcll_idx = find(g.svc.svc_con(:,9)~=0);
    for j = 1:g.svc.n_svc
        index = find(g.svc.svc_con(j,2)== g.ncl.load_con(:,1));
        if ~isempty(index)
            g.svc.svc_idx(j) = index;
        else
            error('you must have the svc bus declared as a non-conforming load')
        end
    end
    % check for user defined controls
    if ~isempty(g.svc.svc_dc)
        g.svc.n_dcud = size(g.svc.svc_dc,1);
        for j = 1:g.svc.n_dcud
            g.svc.dcud_idx(j) = g.svc.svc_dc{j,2};
        end
    end
end

