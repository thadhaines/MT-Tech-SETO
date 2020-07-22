function ml_sig(k)
% system blows up if timestep becomes larger than 0.008...
%% global G
global g

if g.lmod.n_lmod~=0
    %if n_lmod~=0
    if g.sys.t(k) > 1 && g.sys.t(k) < 5
        % load ramp
        g.lmod.lmod_sig(1,k) = g.lmod.lmod_sig(1,k-1) + ...
            0.125 * abs(g.sys.t(k)-g.sys.t(k-1));
    end
    if g.sys.t(k) >= 5
        g.lmod.lmod_sig(1,k) = 0.5;
    end
end
return