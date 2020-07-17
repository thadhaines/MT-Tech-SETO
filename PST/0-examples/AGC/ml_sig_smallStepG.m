function ml_sig(k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
%% legacy
%global lmod_sig n_lmod

%% global G
global g

if g.lmod.n_lmod~=0
%if n_lmod~=0
    if g.sys.t(k) > 1
        % load step
        g.lmod.lmod_sig(1,k) = 5.0; % modify first load only
        %lmod_sig(1,k) = 0.01; % modify first load only
    end
end
return