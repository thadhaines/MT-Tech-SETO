function [Gen] = funSynMacSub_DansPSTN2(Gen,ii,dt,Flag)

for jj  = 1:length(Gen)

if Flag==0 %Initialize
    %Parameters
    if 1==1%~isfield(Gen(jj),'k')
        Gen(jj).k(1) = Gen(jj).Xd_p - Gen(jj).Xl;
        Gen(jj).k(2) = Gen(jj).Xd - Gen(jj).Xd_p;
        Gen(jj).k(3) = (Gen(jj).Xd_pp - Gen(jj).Xl)/Gen(jj).k(1);
        Gen(jj).k(4) = (Gen(jj).Xd_p - Gen(jj).Xd_pp)/Gen(jj).k(1);
        Gen(jj).k(5) = Gen(jj).k(4)/Gen(jj).k(1);

        Gen(jj).k(6) = Gen(jj).Xq_p - Gen(jj).Xl;
        Gen(jj).k(7) = Gen(jj).Xq - Gen(jj).Xq_p;
        Gen(jj).k(8) = (Gen(jj).Xq_pp - Gen(jj).Xl)/Gen(jj).k(6);
        Gen(jj).k(9) = (Gen(jj).Xq_p - Gen(jj).Xq_pp)/Gen(jj).k(6);
        Gen(jj).k(10) = Gen(jj).k(9)/Gen(jj).k(6);
        
    end
    
    %Rotor angle
    e = Gen(jj).Vt(ii,1) + (Gen(jj).R + 1i*Gen(jj).Xq)*Gen(jj).I(ii,1);
    Gen(jj).RotAng(ii,1) = atan2(imag(e),real(e));
    
    %Field voltage
    rot = 1i*exp(-1i*Gen(jj).RotAng(ii,1)); %rotates v and i from system frame to dq frame
    vdq = rot.*Gen(jj).Vt(ii,1);
    vd = real(vdq);
    vq = imag(vdq);
    idq = rot.*Gen(jj).I(ii,1);
    id = real(idq);
    iq = imag(idq);
    Gen(jj).Efd(ii,1) = vq + Gen(jj).Xd*id + Gen(jj).R*iq;
    
    %States
    Gen(jj).Psikd(ii,1) = vq + Gen(jj).R*iq + Gen(jj).Xl*id;
    Gen(jj).Psikq(ii,1) = -vd - Gen(jj).R*id + Gen(jj).Xl*iq;
    Gen(jj).Edp(ii,1) = -Gen(jj).k(7)*iq;
    Gen(jj).Eqp(ii,1) = Gen(jj).Psikd(ii,1) + Gen(jj).k(1)*id;
    clear e rot vdq idq vd vq id iq
elseif Flag==1 %Network interface
    rot = -1i*exp(1i*Gen(jj).RotAng(ii,1)); %rotates v and i from dq frame to system frame
    psidpp = Gen(jj).k(3)*Gen(jj).Eqp(ii,1) + Gen(jj).k(4)*Gen(jj).Psikd(ii,1);
    psiqpp = Gen(jj).k(8)*Gen(jj).Edp(ii,1) + Gen(jj).k(9)*Gen(jj).Psikq(ii,1);
    edq = (-psiqpp + 1i*psidpp)*Gen(jj).Wr(ii,1);
    Gen(jj).E(ii,1) = edq*rot; %gen internal voltage on system frame
    clear rot edq psidpp psiqpp
elseif Flag==2 %State eqns
    rot = 1i*exp(-1i*Gen(jj).RotAng(ii,1)); %rotates v and i from system frame to dq frame
    idq = Gen(jj).I(ii,1)*rot;
    id = real(idq);
    iq = imag(idq);
    psidpp = Gen(jj).k(3)*Gen(jj).Eqp(ii,1) + Gen(jj).k(4)*Gen(jj).Psikd(ii,1);
    psiqpp = Gen(jj).k(8)*Gen(jj).Edp(ii,1) + Gen(jj).k(9)*Gen(jj).Psikq(ii,1);
    
    %Circuit state derivatives
    e = Gen(jj).Eqp(ii,1) - Gen(jj).Psikd(ii,1) - Gen(jj).k(1)*id;
    Gen(jj).dotEqp(ii,1) = (Gen(jj).Efd(ii,1) - Gen(jj).Eqp(ii,1) - Gen(jj).k(2)*(id + Gen(jj).k(5)*e))/Gen(jj).Td0_p;
    Gen(jj).dotPsikd(ii,1) = e/Gen(jj).Td0_pp;
    e = Gen(jj).Edp(ii,1) - Gen(jj).Psikq(ii,1) - Gen(jj).k(6)*iq;
    Gen(jj).dotEdp(ii,1) = (-Gen(jj).Edp(ii,1) + Gen(jj).k(7)*(-iq - Gen(jj).k(10)*e))/Gen(jj).Tq0_p;
    Gen(jj).dotPsikq(ii,1) = e/Gen(jj).Tq0_pp;
    
    %Swing state derivatives

    Te = psidpp*iq - psiqpp*id;
    Gen(jj).dotRotAng(ii,1) = Gen(jj).Ws*(Gen(jj).Wr(ii,1) - 1);
    Gen(jj).dotWr(ii,1) = (Gen(jj).Pm(ii,1)/Gen(jj).Wr(ii,1) - Te - Gen(jj).D*(Gen(jj).Wr(ii,1)-1))/(2*Gen(jj).H);
    
    %Integrate
    Gen(jj).Eqp(ii+1,1) =   Gen(jj).Eqp(ii,1) +     dt*(1.5*Gen(jj).dotEqp(ii,1) - 0.5*Gen(jj).dotEqp(ii-1,1));
    Gen(jj).Edp(ii+1,1) =   Gen(jj).Edp(ii,1) +     dt*(1.5*Gen(jj).dotEdp(ii,1) - 0.5*Gen(jj).dotEdp(ii-1,1));
    Gen(jj).Psikd(ii+1,1) = Gen(jj).Psikd(ii,1) +   dt*(1.5*Gen(jj).dotPsikd(ii,1) - 0.5*Gen(jj).dotPsikd(ii-1,1));
    Gen(jj).Psikq(ii+1,1) = Gen(jj).Psikq(ii,1) +   dt*(1.5*Gen(jj).dotPsikq(ii,1) - 0.5*Gen(jj).dotPsikq(ii-1,1));
    Gen(jj).Wr(ii+1,1) =    Gen(jj).Wr(ii,1) +      dt*(1.5*Gen(jj).dotWr(ii,1) - 0.5*Gen(jj).dotWr(ii-1,1));
    Gen(jj).RotAng(ii+1,1) =Gen(jj).RotAng(ii,1) +  dt*(1.5*Gen(jj).dotRotAng(ii,1) - 0.5*Gen(jj).dotRotAng(ii-1,1));
    clear e Te psidpp psiqpp id iq idq rot
else
    error('Invalid Flag')
end

end
    
end
    