function [Gov] = funGov_TGOV1(Gov,kTime,dt,Flag)

if Flag==0 %Initialize
    %Parameters
    if ~isfield(Gov.param,'k')
        Gov.param.k(1) = 1-Gov.param.T2/Gov.param.T3;
        Gov.param.k(2) = Gov.param.k(1)/Gov.param.T3;
        Gov.param.k(3) = Gov.param.T2/Gov.param.T3;
    end
    Gov.x1(kTime) = Gov.Pref(kTime);
    Gov.x2(kTime) = Gov.param.k(1)*Gov.Pref(kTime);
elseif Flag==1; %Pmech
    Gov.Pmech(kTime) = Gov.x2(kTime) + Gov.param.k(3)*Gov.x1(kTime) - Gov.param.Dt*(Gov.w(kTime)-1);
elseif Flag==2 %Update states
    dw = Gov.w(kTime)-1;
    d = Gov.Pref(kTime)-dw/Gov.param.R;
    if d>=Gov.param.Vmax
        d = Gov.param.Vmax;
    elseif d<=Gov.param.Vmin
        d = Gov.param.Vmin;
    end
    Gov.x1dot(kTime) = (-Gov.x1(kTime) + d)/Gov.param.T1;
    Gov.x2dot(kTime) = -Gov.x2(kTime)/Gov.param.T3 + Gov.param.k(2)*Gov.x1(kTime);  
    Gov.x1(kTime+1) = Gov.x1(kTime) + dt*(1.5*Gov.x1dot(kTime) - 0.5*Gov.x1dot(kTime-1));
    Gov.x2(kTime+1) = Gov.x2(kTime) + dt*(1.5*Gov.x2dot(kTime) - 0.5*Gov.x2dot(kTime-1));
else
    error('Invalid Flag')
end
end

