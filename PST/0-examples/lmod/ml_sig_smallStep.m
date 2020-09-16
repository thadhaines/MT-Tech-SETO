function ml_sig(t, k)
% Syntax: ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global lmod_sig n_lmod

if n_lmod~=0
    
    if t(k) > 1
        % load step
        lmod_sig(1,k) = 0.01; % modify first load only
    end
    if t(k) > 10
        lmod_sig(1,k) = 0.0; % modify first load only
    end
end
return