function ml_sig(t, k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
%% legacy
global lmod_sig n_lmod

if n_lmod~=0
    if t(k) > 1
        % load step
        lmod_sig(1,k) = 0.1; % modify first load only
    end
    if t(k) > 2
        % load step
        lmod_sig(1,k) = 0.0; % modify first load only
    end

end
return