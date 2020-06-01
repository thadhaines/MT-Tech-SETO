function ml_sig(k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global lmod_sig n_lmod
global g
%f=0; %dummy variable
%% Function to modulate load - meant to replace ml_sig in main PST dir
% lmod_con must be specified in the data file
% and the load bus must be in the nonconforming load list.
%lmod_sig(:,k) = zeros(size(lmod_sig(:,k))); %Initialized in s_simu

%fprintf('%4.4f \t %d\n', t(k), k); % DEBUG
%if g.lmod.n_lmod~=0
if n_lmod~=0
    if g.sys.t(k) > 1
        % load step
        %g.lmod.lmod_sig(1,k) = 0.25; % modify first load only
        
        lmod_sig(1,k) = 0.25; % modify first load only
    end
    
%     if (t(k) > 12) && (t(k) < 18)
%         % random noise
%         lmod_sig(1,k) = 0.25+ 0.1*rand();
%     end
        
end
return