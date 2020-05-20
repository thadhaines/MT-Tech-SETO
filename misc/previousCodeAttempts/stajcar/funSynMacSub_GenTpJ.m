function [Gen] = funSynMacSub_GenTpJ(Gen,ii,dt,Flag)

for jj  = 1:length(Gen)

if Flag==0 %Initialize
    
    %Parameters
    if 1==1%~isfield(Gen(jj),'k')
        Gen(jj).k(1) = Gen(jj).Xd - Gen(jj).Xd_p;
        Gen(jj).k(2) = Gen(jj).Xd_p - Gen(jj).Xd_pp;
        Gen(jj).k(3) = Gen(jj).Xd - Gen(jj).Xd_pp;
        Gen(jj).k(4) = Gen(jj).k(1)/Gen(jj).k(2);
        Gen(jj).k(5) = Gen(jj).k(3)/Gen(jj).k(2);

        Gen(jj).k(6) = Gen(jj).Xq - Gen(jj).Xq_p;
        Gen(jj).k(7) = Gen(jj).Xq_p - Gen(jj).Xq_pp;
        Gen(jj).k(8) = Gen(jj).Xq - Gen(jj).Xq_pp;
        Gen(jj).k(9) = Gen(jj).k(6)/Gen(jj).k(7);
        Gen(jj).k(10) = Gen(jj).k(8)/Gen(jj).k(7);
        
        % Assume no saturation.
        Gen(jj).Se = 1;
        
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
    Epp = vdq + idq * (Gen(jj).R + 1j * Gen(jj).Xd_pp);
    Gen(jj).Psid_pp(ii,1) = imag(Epp);
    Gen(jj).Psiq_pp(ii,1) = -real(Epp);
    Gen(jj).Ed_p(ii,1) = 1 / Gen(jj).Se * (...
        Gen(jj).Psiq_pp(ii,1) * Gen(jj).Se + iq * Gen(jj).k(7) );
    Gen(jj).Eq_p(ii,1) = 1 / Gen(jj).Se * (...
        Gen(jj).Psid_pp(ii,1) + id * Gen(jj).k(2) );
    clear e rot vdq idq vd vq id iq Epp
    
elseif Flag==1 %Network interface
    
    rot = -1i*exp(1i*Gen(jj).RotAng(ii,1)); %rotates v and i from dq frame to system frame
    edq = (-Gen(jj).Psiq_pp(ii,1) + 1i * Gen(jj).Psid_pp(ii,1)) * Gen(jj).Wr(ii,1);
    Gen(jj).E(ii,1) = edq*rot; %gen internal voltage on system frame
    clear rot edq psidpp psiqpp
    
elseif Flag==2 %State eqns
    
    rot = 1i*exp(-1i*Gen(jj).RotAng(ii,1)); %rotates v and i from system frame to dq frame
    idq = Gen(jj).I(ii,1)*rot;
    id = real(idq);
    iq = imag(idq);
    
    %Circuit state derivatives
    Gen(jj).dotPsid_pp(ii,1) = 1 / Gen(jj).Td0_pp * (...
        Gen(jj).Eq_p(ii,1) * Gen(jj).Se - ...
        Gen(jj).Psid_pp(ii,1) * Gen(jj).Se - ...
        id * Gen(jj).k(2) );        
    Gen(jj).dotPsiq_pp(ii,1) = 1 / Gen(jj).Tq0_pp * (...
        Gen(jj).Ed_p(ii,1) * Gen(jj).Se - ...
        Gen(jj).Psiq_pp(ii,1) * Gen(jj).Se - ...
        iq * Gen(jj).k(7) );
    Gen(jj).dotEq_p(ii,1) = 1 / Gen(jj).Td0_p * (...
        Gen(jj).Efd(ii,1) + ...
        Gen(jj).Psid_pp(ii,1) * Gen(jj).Se * Gen(jj).k(4) - ...
        Gen(jj).Eq_p(ii,1) * Gen(jj).Se * Gen(jj).k(5) );
    Gen(jj).dotEd_p(ii,1) = 1 / Gen(jj).Tq0_p * (...
        Gen(jj).Psiq_pp(ii,1) * Gen(jj).Se * Gen(jj).k(9) - ...
        Gen(jj).Ed_p(ii,1) * Gen(jj).Se * Gen(jj).k(10) );
    
    %Swing state derivatives

    Te = Gen(jj).Psid_pp(ii,1) * iq - Gen(jj).Psiq_pp(ii,1) * id;
    Gen(jj).dotRotAng(ii,1) = Gen(jj).Ws*(Gen(jj).Wr(ii,1) - 1);
    Gen(jj).dotWr(ii,1) = (Gen(jj).Pm(ii,1)/Gen(jj).Wr(ii,1) - Te - Gen(jj).D*(Gen(jj).Wr(ii,1)-1))/(2*Gen(jj).H);
    
    %Integrate
    Gen(jj).Psid_pp(ii+1,1) =   Gen(jj).Psid_pp(ii,1) +     dt*(1.5*Gen(jj).dotPsid_pp(ii,1) -  0.5*Gen(jj).dotPsid_pp(ii-1,1));
    Gen(jj).Psiq_pp(ii+1,1) =   Gen(jj).Psiq_pp(ii,1) +     dt*(1.5*Gen(jj).dotPsiq_pp(ii,1) -  0.5*Gen(jj).dotPsiq_pp(ii-1,1));
    Gen(jj).Eq_p(ii+1,1) =      Gen(jj).Eq_p(ii,1) +        dt*(1.5*Gen(jj).dotEq_p(ii,1) -     0.5*Gen(jj).dotEq_p(ii-1,1));
    Gen(jj).Ed_p(ii+1,1) =      Gen(jj).Ed_p(ii,1) +        dt*(1.5*Gen(jj).dotEd_p(ii,1) -     0.5*Gen(jj).dotEd_p(ii-1,1));
    Gen(jj).Wr(ii+1,1) =        Gen(jj).Wr(ii,1) +          dt*(1.5*Gen(jj).dotWr(ii,1) -       0.5*Gen(jj).dotWr(ii-1,1));
    Gen(jj).RotAng(ii+1,1) =    Gen(jj).RotAng(ii,1) +      dt*(1.5*Gen(jj).dotRotAng(ii,1) -   0.5*Gen(jj).dotRotAng(ii-1,1));
    clear e Te id iq idq rot
else
    error('Invalid Flag')
end

end
    
end
    