function rml_sig(k)
% Syntax: rml_sig(k)
%4:40 PM 15/08/97
% defines modulation signal for rlmod control

global g
%% Function to modulate reactive load - meant to replace ml_sig in main PST dir
% rlmod_con must be specified in the data file
% and the load bus must be in the nonconforming load list.

if g.rlmod.n_rlmod~=0
    if g.sys.t(k) > 1
        g.rlmod.rlmod_sig(1,k) = 0.01;
    end
    if g.sys.t(k) > 8
        g.rlmod.rlmod_sig(1,k) = 0.0;
    end
    
end