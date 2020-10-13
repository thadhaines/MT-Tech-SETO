function ml_sig(k)
% Syntax: ml_sig(k)
global g

if g.lmod.n_lmod~=0
    if g.sys.t(k) >= 1
        % load step
        g.lmod.lmod_sig(1,k) = 0.25; % modify first load only
    end
    
end
return