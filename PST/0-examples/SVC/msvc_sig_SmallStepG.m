function msvc_sig(k)
% Syntax: f = msvc_sig(t,k)
% 4:39 PM 15/08/97
% defines modulation signal for svc control
global g
if g.svc.n_svc ~=0
  if g.sys.t(k)>=0.5 && g.sys.t(k) < 1.5
     g.svc.svc_sig(1,k) = 0.01;
  else
     g.svc.svc_sig(1,k) = 0;
  end
end
return