function [Gen] = funSynMacSub_PSTN2_Slow(Gen,kTime,dt,dtStartFlag,wbase,Flag)

if Flag==0 %Initialize
    %Parameters
%     if ~isfield(Gen.param,'k')
%         Gen.param.k(1) = Gen.param.xdp - Gen.param.xl;
%         Gen.param.k(2) = Gen.param.xd - Gen.param.xdp;
%         Gen.param.k(3) = (Gen.param.xdpp - Gen.param.xl)/Gen.param.k(1);
%         Gen.param.k(4) = (Gen.param.xdp - Gen.param.xdpp)/Gen.param.k(1);
%         Gen.param.k(5) = Gen.param.k(4)/Gen.param.k(1);
% 
%         Gen.param.k(6) = Gen.param.xqp - Gen.param.xl;
%         Gen.param.k(7) = Gen.param.xq - Gen.param.xqp;
%         Gen.param.k(8) = (Gen.param.xqpp - Gen.param.xl)/Gen.param.k(6);
%         Gen.param.k(9) = (Gen.param.xqp - Gen.param.xqpp)/Gen.param.k(6);
%         Gen.param.k(10) = Gen.param.k(9)/Gen.param.k(6);
%     end
    
    %Rotor angle
    e = Gen.VT(kTime) + (Gen.param.ra + 1i*Gen.param.xq)*Gen.IT(kTime);
    Gen.delta(kTime) = atan2(imag(e),real(e));
    
    %Field voltage
    rot = 1i*exp(-1i*Gen.delta(kTime)); %rotates v and i from system frame to dq frame
    vdq = rot.*Gen.VT(kTime);
    vd = real(vdq);
    vq = imag(vdq);
    idq = rot.*Gen.IT(kTime);
    id = real(idq);
    iq = imag(idq);
    Gen.Efd(kTime) = vq + Gen.param.xd*id + Gen.param.ra*iq;
    
    %States
    Gen.psikd(kTime) = vq + Gen.param.ra*iq + Gen.param.xl*id;
    Gen.psikq(kTime) = -vd - Gen.param.ra*id + Gen.param.xl*iq;
    Gen.edp(kTime) = -Gen.param.k(7)*iq;
    Gen.eqp(kTime) = Gen.psikd(kTime) + Gen.param.k(1)*id;
    clear e rot vdq idq vd vq id iq
elseif Flag==1 %Network interface
    rot = -1i*exp(1i*Gen.delta(kTime)); %rotates v and i from dq frame to system frame
    psidpp = Gen.param.k(3)*Gen.eqp(kTime) + Gen.param.k(4)*Gen.psikd(kTime);
    psiqpp = Gen.param.k(8)*Gen.edp(kTime) + Gen.param.k(9)*Gen.psikq(kTime);
    edq = (-psiqpp + 1i*psidpp)*Gen.w(kTime);
    Gen.E(kTime) = edq*rot; %gen internal voltage on system frame
    clear rot edq psidpp psiqpp
elseif Flag==2 %State eqns
    rot = 1i*exp(-1i*Gen.delta(kTime)); %rotates v and i from system frame to dq frame
    idq = Gen.IT(kTime)*rot;
    id = real(idq);
    iq = imag(idq);
    psidpp = Gen.param.k(3)*Gen.eqp(kTime) + Gen.param.k(4)*Gen.psikd(kTime);
    psiqpp = Gen.param.k(8)*Gen.edp(kTime) + Gen.param.k(9)*Gen.psikq(kTime);
    
    %Circuit state derivatives
    e = Gen.eqp(kTime) - Gen.psikd(kTime) - Gen.param.k(1)*id;
    Gen.eqpdot(kTime) = (Gen.Efd(kTime) - Gen.eqp(kTime) - Gen.param.k(2)*(id + Gen.param.k(5)*e))/Gen.param.Td0p;
    Gen.psikddot(kTime) = 0;%e/Gen.param.Td0pp;
    e = Gen.edp(kTime) - Gen.psikq(kTime) - Gen.param.k(6)*iq;
    Gen.edpdot(kTime) = (-Gen.edp(kTime) + Gen.param.k(7)*(-iq - Gen.param.k(10)*e))/Gen.param.Tq0p;
    Gen.psikqdot(kTime) = 0;%e/Gen.param.Tq0pp;
    
    %Swing state derivatives
%     Te = psidpp*iq - psiqpp*id;
%     Gen.deltadot(kTime) = wbase*(Gen.w(kTime) - 1);
%     Gen.wdot(kTime) = (Gen.Pmech(kTime)/Gen.w(kTime) - Te - Gen.param.D*(Gen.w(kTime)-1))/(2*Gen.param.H);
    
    %Integrate
    if dtStartFlag
        Gen.eqp(kTime+1) = Gen.eqp(kTime) + dt*Gen.eqpdot(kTime);
        Gen.edp(kTime+1) = Gen.edp(kTime) + dt*Gen.edpdot(kTime);
    else
        Gen.eqp(kTime+1) = Gen.eqp(kTime) + dt*(1.5*Gen.eqpdot(kTime) - 0.5*Gen.eqpdot(kTime-1));
        Gen.edp(kTime+1) = Gen.edp(kTime) + dt*(1.5*Gen.edpdot(kTime) - 0.5*Gen.edpdot(kTime-1));
    end
    Gen.psikd(kTime+1) = Gen.psikd(kTime); %Gen.psikd(kTime) + dt*(1.5*Gen.psikddot(kTime) - 0.5*Gen.psikddot(kTime-1));
    Gen.psikq(kTime+1) = Gen.psikq(kTime); %Gen.psikq(kTime) + dt*(1.5*Gen.psikqdot(kTime) - 0.5*Gen.psikqdot(kTime-1));
%     Gen.w(kTime+1) = Gen.w(kTime) + dt*(1.5*Gen.wdot(kTime) - 0.5*Gen.wdot(kTime-1));
%     Gen.delta(kTime+1) = Gen.delta(kTime) + dt*(1.5*Gen.deltadot(kTime) - 0.5*Gen.deltadot(kTime-1));
    clear e Te psidpp psiqpp id iq idq rot
else
    error('Invalid Flag')
end
    
end
    