function msvc_sig(t,k)
% Syntax: f = msvc_sig(t,k)
% 4:39 PM 15/08/97
% defines modulation signal for svc control
global svc_sig n_svc
% f=0; %dummy variable
if n_svc ~=0
  if t>=0.5 && t < 1.5
     svc_sig(1,k) = 0.01;
  else
     svc_sig(1,k) = 0;
  end
end
return