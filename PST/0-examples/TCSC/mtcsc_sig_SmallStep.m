function mtcsc_sig(t,k)
% Syntax: f = mtcsc_sig(t,k)
% 4:39 PM 15/08/97
% defines modulation signal for tcsc control
global tcsc_sig n_tcsc
% f=0; %dummy variable
if n_tcsc ~=0
  if t(k) >= 0.5 && t(k)<1.5
    tcsc_sig(1,k) = 0.01;
  else
     tcsc_sig(:,k) = 0;
  end
end
return