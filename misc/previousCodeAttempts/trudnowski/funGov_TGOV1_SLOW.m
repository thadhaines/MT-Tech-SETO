function [Gov] = funGov_TGOV1_SLOW(Gov,kTime,dt,dtStartFlag,Flag)

if Flag==0 %Initialize with non-zero states
    %Parameters
%     if ~isfield(Gov.param,'k')
%         Gov.param.k(1) = 1-Gov.param.T2/Gov.param.T3;
%         Gov.param.k(2) = Gov.param.k(1)/Gov.param.T3;
%         Gov.param.k(3) = Gov.param.T2/Gov.param.T3;
%     end
    dw = Gov.w(kTime)-1;
    Gov.x2(kTime) = Gov.Pmech(kTime) - Gov.param.k(3)*(Gov.Pref(kTime) - dw/Gov.param.R) - Gov.param.Dt*dw;
elseif Flag==1; %Pmech
    dw = Gov.w(kTime)-1;
    Gov.Pmech(kTime) = Gov.x2(kTime) + Gov.param.k(3)*(Gov.Pref(kTime) - dw/Gov.param.R) - Gov.param.Dt*dw;
elseif Flag==2 %Update states
    dw = Gov.w(kTime)-1;
    d = Gov.Pref(kTime)-dw/Gov.param.R;
    if d>=Gov.param.Vmax
        d = Gov.param.Vmax;
    elseif d<=Gov.param.Vmin
        d = Gov.param.Vmin;
    end
    Gov.x2dot(kTime) = -Gov.x2(kTime)/Gov.param.T3 + Gov.param.k(2)*d;  
    if dtStartFlag
        Gov.x2(kTime+1) = Gov.x2(kTime) + dt*Gov.x2dot(kTime); %Euler to get started
    else
        Gov.x2(kTime+1) = Gov.x2(kTime) + dt*(1.5*Gov.x2dot(kTime) - 0.5*Gov.x2dot(kTime-1));
    end
else
    error('Invalid Flag')
end
end

