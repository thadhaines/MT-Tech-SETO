function f = rml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control

global g
f=0; %dummy variable
%% Function to modulate reactive load - meant to replace ml_sig in main PST dir
% rlmod_con must be specified in the data file
% and the load bus must be in the nonconforming load list.

%% V3 code
global rlmod_sig
global n_rlmod
if n_rlmod~=0
    if t> 1
        rlmod_sig(1,k) = 0.01;
    end

end