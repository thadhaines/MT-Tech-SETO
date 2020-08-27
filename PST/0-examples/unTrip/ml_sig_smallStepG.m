function ml_sig(k)
% Syntax: f = ml_sig(k)

% defines modulation signal for lmod control

global g
if g.lmod.n_lmod~=0
    if g.sys.t(k) > 2
        % load step
        g.lmod.lmod_sig(1,k) = 0.1; % modify first load only
    end
end
return