function f = mtg_sig(t,k)
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig n_tg n_tgh
f=0; %dummy variable
% kTg = 20; %turnbine gove number
          %Gen 34 -- kTg=21
          %Gen 33 -- kTg=20
if n_tg~=0|n_tgh~=0
  tg_sig(:,k) = zeros(n_tg+n_tgh,1);
%   if t>=10;
%       tg_sig(kTg,k) = 0.01*sin(0.3719*2*pi*t); %
%   end
%   if t>=1.0 & t<1.2
%       tg_sig(4,k) = 0.;
%   end
end
% if k==8; size(tg_sig), end