function ml_sig(t, k)
% Syntax: ml_sig(k)
% for pre structured global runs

global n_lmod lmod_sig

if n_lmod~=0
    if t(k) >= 1
        % load step
        lmod_sig(1,k) = 0.25; % modify first load only
    end
    
end
return