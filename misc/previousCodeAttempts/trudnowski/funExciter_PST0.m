function [Exc] = funExciter_PST0(Exc,kTime,dt,Flag)

if Flag==0 %Initialize
    %Parameters
    if ~isfield(Exc.param,'k')
        if Exc.param.TB==0
            if Exc.param.TC~=0; error(' '); end
            Exc.param.k(1) = 1;
            Exc.param.k(2) = 0;
        else
            Exc.param.k(1) = Exc.param.TC/Exc.param.TB;
            Exc.param.k(2) = 1 - Exc.param.k(1);
        end
    end
    Exc.Vref(kTime) = abs(Exc.VT(kTime)) + Exc.Efd(kTime)/Exc.param.KA;
    e = Exc.Vref(kTime) - abs(Exc.VT(kTime));
    Exc.x1(kTime) = Exc.param.k(2)*e;
%     Exc.Efd(kTime) = Exc.param.KA*(Exc.param.k(1)*e + Exc.x1(kTime)) + Exc.param.EfdOffSet;
    clear e
elseif Flag==1 %Network interface
    Exc.Efd(kTime) = Exc.Efd(kTime);
elseif Flag==2 %State eqns
    e = Exc.Vref(kTime) - abs(Exc.VT(kTime));
    Exc.x1dot(kTime) = (-Exc.x1(kTime) + Exc.param.k(2)*e)/Exc.param.TB;
    if Exc.param.TB==0 
        Exc.x1dot(kTime) = 0;
%         Exc.Efddot(kTime) = (-Exc.Efd(kTime) + Exc.param.KA*e + Exc.param.EfdOffSet)/Exc.param.TA;
        Exc.Efddot(kTime) = (-Exc.Efd(kTime) + Exc.param.KA*e)/Exc.param.TA;
    else
        Exc.x1dot(kTime) = (-Exc.x1(kTime) + Exc.param.k(2)*e)/Exc.param.TB;
%         Exc.Efddot(kTime) = (-Exc.Efd(kTime) + Exc.param.KA*(Exc.param.k(1)*e + Exc.x1(kTime)) + Exc.param.EfdOffSet)/Exc.param.TA;
        Exc.Efddot(kTime) = (-Exc.Efd(kTime) + Exc.param.KA*(Exc.param.k(1)*e + Exc.x1(kTime)))/Exc.param.TA;
    end
    
    %Integrate
    Exc.x1(kTime+1) = Exc.x1(kTime) + dt*(1.5*Exc.x1dot(kTime) - 0.5*Exc.x1dot(kTime-1));
    Exc.Efd(kTime+1) = Exc.Efd(kTime) + dt*(1.5*Exc.Efddot(kTime) - 0.5*Exc.Efddot(kTime-1));
    if Exc.Efd(kTime+1)>Exc.param.Efdmax; 
        Exc.Efddot(kTime) = 0;
        Exc.Efd(kTime+1) = Exc.param.Efdmax;
    elseif Exc.Efd(kTime+1)<Exc.param.Efdmin; 
        Exc.Efddot(kTime) = 0;
        Exc.Efd(kTime+1) = Exc.param.Efdmin;
    end
    clear e
else
    error('Invalid Flag')
end
    
end
    