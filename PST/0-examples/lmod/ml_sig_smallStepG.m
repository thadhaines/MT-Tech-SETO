function ml_sig(k)
% Syntax: ml_sig(k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control

global g

if g.lmod.n_lmod~=0
    if g.sys.t(k) > 1
        % load step
        g.lmod.lmod_sig(1,k) = 0.01; % modify first load only
    end
    if g.sys.t(k) > 10
        g.lmod.lmod_sig(1,k) = 0.0; % modify first load only
    end
end
return