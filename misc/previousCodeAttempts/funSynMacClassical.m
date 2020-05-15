function Gen = funSynMacClassical(Gen,kTime,dt,wbase,Flag)

if Flag==0
    e = Gen.VT(kTime) + (1i*Gen.param.xdp)*Gen.IT(kTime);
    Gen.delta(kTime) = angle(e);
    Gen.Efd(kTime) = abs(e);
elseif Flag==1 %Internal voltage
    Gen.E(kTime) = Gen.Efd(kTime).*exp(1i*Gen.delta(kTime));
elseif Flag==2 %states
    Pe = real(Gen.E(kTime)*conj(Gen.IT(kTime)));
    Gen.deltadot(kTime) = wbase*(Gen.w(kTime) - 1);
%     Gen.wdot(kTime) = (Gen.Pmech(kTime)/Gen.w(kTime) - Pe/Gen.w(kTime) - Gen.param.D*(Gen.w(kTime)-1))/(2*Gen.param.H);
    Gen.wdot(kTime) = (Gen.Pmech(kTime) - Pe - Gen.param.D*(Gen.w(kTime)-1))/(2*Gen.param.H);
    
    %Integrate
    Gen.w(kTime+1) = Gen.w(kTime) + dt*(1.5*Gen.wdot(kTime) - 0.5*Gen.wdot(kTime-1));
    Gen.delta(kTime+1) = Gen.delta(kTime) + dt*(1.5*Gen.deltadot(kTime) - 0.5*Gen.deltadot(kTime-1));
end

