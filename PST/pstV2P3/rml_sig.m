function f = rml_sig(t,k)
% Syntax: f = rml_sig(t,k)
%5:43 PM 27/8/97
% defines modulation signal for rlmod control
global rlmod_sig n_rlmod bus_v mac_spd
f=0; %dummy variable
if n_rlmod~=0
    rlmod_sig(:,k)=zeros(n_rlmod,1);
end
% if t(k)>=6 & t(k)<6.5
%     rlmod_sig(46:64,k) = 0.025;
% end
