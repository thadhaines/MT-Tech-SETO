function f = mexc_sig(t,k)
% Syntax: f = mexc_sig(t,k)
% 1:20 PM 15/08/97
% defines modulation signal for exciter control
global exc_sig n_exc
f=0; %dummy variable
if n_exc~=0
    exc_sig(:,k) = zeros(n_exc,1);
end