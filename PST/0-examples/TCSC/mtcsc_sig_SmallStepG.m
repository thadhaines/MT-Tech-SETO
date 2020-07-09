function mtcsc_sig(k)
% Syntax: f = mtcsc_sig(k)
% 4:39 PM 15/08/97
% defines modulation signal for tcsc control
global g
% f=0; %dummy variable
if g.tcsc.n_tcsc ~=0
  if g.sys.t(k) >= 0.5 && g.sys.t(k) <1.5
    g.tcsc.tcsc_sig(1,k) = 0.01;
  else
     g.tcsc.tcsc_sig(:,k) = 0;
  end
end
return