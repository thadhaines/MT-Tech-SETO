function [Gen] = funSynMacSub_Genrou(Gen,ii,dt,Flag)

for jj  = 1:length(Gen)

if Flag==0 %Initialize
    
    %Parameters
    if 1==1%~isfield(Gen(jj),'k')
        Gen(jj).k(1) = Gen(jj).Xd_p - Gen(jj).Xl;
        Gen(jj).k(2) = Gen(jj).Xd - Gen(jj).Xd_p;
        Gen(jj).k(3) = (Gen(jj).Xd - Gen(jj).Xd_p)*(Gen(jj).Xd_p - Gen(jj).Xd_pp)/(Gen(jj).k(1) ^ 2);
        Gen(jj).k(4) = Gen(jj).Xd_pp - Gen(jj).Xl;
        Gen(jj).k(5) = Gen(jj).Xd_p - Gen(jj).Xd_pp;
        Gen(jj).k(6) = Gen(jj).k(4)/Gen(jj).k(1);
        Gen(jj).k(7) = Gen(jj).k(5)/Gen(jj).k(1);
        
        Gen(jj).k(8) = Gen(jj).Xq_p - Gen(jj).Xl;
        Gen(jj).k(9) = Gen(jj).Xq - Gen(jj).Xq_p;
        Gen(jj).k(10) = (Gen(jj).Xq - Gen(jj).Xq_p)*(Gen(jj).Xq_p - Gen(jj).Xq_pp)/(Gen(jj).k(8) ^ 2);
        Gen(jj).k(11) = Gen(jj).Xq_pp - Gen(jj).Xl;
        Gen(jj).k(12) = Gen(jj).Xq_p - Gen(jj).Xq_pp;
        Gen(jj).k(13) = Gen(jj).k(11)/Gen(jj).k(8);
        Gen(jj).k(14) = Gen(jj).k(12)/Gen(jj).k(8);
        Gen(jj).k(15) = (Gen(jj).Xq - Gen(jj).Xl)/(Gen(jj).Xd - Gen(jj).Xl);
        
        % Assume no saturation.
        Gen(jj).Se = 0;
        
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
    Gen(jj).Psi_pp(ii,1) = sqrt( Gen(jj).Psid_pp(ii,1)^2 + Gen(jj).Psiq_pp(ii,1)^2 );
    num = Gen(jj).Psid_pp(ii,1) * inv(Gen(jj).k(7)) + id * Gen(jj).k(1);
    den = 1 + Gen(jj).k(4) / Gen(jj).k(5);
    Gen(jj).Psifd(ii,1) = num/den;
    clear num den
    Gen(jj).Psikd(ii,1) = Gen(jj).Psifd(ii,1) - id * Gen(jj).k(1);
    num = Gen(jj).Psiq_pp(ii,1) * inv(Gen(jj).k(14)) + iq * Gen(jj).k(8);
    den = 1 + Gen(jj).k(11) / Gen(jj).k(12);
    Gen(jj).Psi1q(ii,1) = num/den;
    Gen(jj).Psi2q(ii,1) = Gen(jj).Psi1q(ii,1) - iq * Gen(jj).k(8);
    clear e rot vdq idq vd vq id iq Epp num den
    
elseif Flag==1 %Network interface
    
    rot = -1i*exp(1i*Gen(jj).RotAng(ii,1)); %rotates v and i from dq frame to system frame
    Gen(jj).Psid_pp(ii,1) = Gen(jj).Psifd(ii,1) * Gen(jj).k(6) + Gen(jj).Psikd(ii,1) * Gen(jj).k(7);
    Gen(jj).Psiq_pp(ii,1) = Gen(jj).Psi1q(ii,1) * Gen(jj).k(13) + Gen(jj).Psi2q(ii,1) * Gen(jj).k(14);
    edq = (-Gen(jj).Psiq_pp(ii,1) + 1i * Gen(jj).Psid_pp(ii,1)) * Gen(jj).Wr(ii,1);
    Gen(jj).E(ii,1) = edq*rot; %gen internal voltage on system frame
    Gen(jj).Psi_pp(ii,1) = sqrt( Gen(jj).Psid_pp(ii,1)^2 + Gen(jj).Psiq_pp(ii,1)^2 );
    clear rot edq psidpp psiqpp
    
elseif Flag==2 %State eqns
    
    rot = 1i*exp(-1i*Gen(jj).RotAng(ii,1)); %rotates v and i from system frame to dq frame
    idq = Gen(jj).I(ii,1)*rot;
    id = real(idq);
    iq = imag(idq);
    
    %Circuit state derivatives
    
    temp = Gen(jj).Psifd(ii,1) - Gen(jj).Psikd(ii,1) - id * Gen(jj).k(1);
    Gen(jj).dotPsikd(ii,1) = temp / Gen(jj).Td0_pp;
    
    Gen(jj).dotPsifd(ii,1) = 1 / Gen(jj).Td0_p * (...
        Gen(jj).Efd(ii,1) - ...
        Gen(jj).Psifd(ii,1) - ...
        id * Gen(jj).k(2) - ...
        temp * Gen(jj).k(3) - ...
        Gen(jj).Se );
    clear temp
    
    temp = Gen(jj).Psi1q(ii,1) - Gen(jj).Psi2q(ii,1) - iq * Gen(jj).k(8);
    Gen(jj).dotPsi2q(ii,1) = temp / Gen(jj).Tq0_pp;
    
    Gen(jj).dotPsi1q(ii,1) = 1 / Gen(jj).Tq0_p * (...
        -Gen(jj).Psi1q(ii,1) - ...
        iq * Gen(jj).k(9) - ...
        temp * Gen(jj).k(10)  -...
        Gen(jj).Se );
    clear temp
    
    %Swing state derivatives

    Te = Gen(jj).Psid_pp(ii,1) * iq - Gen(jj).Psiq_pp(ii,1) * id;
    Gen(jj).dotRotAng(ii,1) = Gen(jj).Ws * (Gen(jj).Wr(ii,1) - 1);
    Gen(jj).dotWr(ii,1) = (Gen(jj).Pm(ii,1)/Gen(jj).Wr(ii,1) - Te - Gen(jj).D*(Gen(jj).Wr(ii,1)-1))/(2*Gen(jj).H);
    
    %Integrate
    Gen(jj).Psifd(ii+1,1) = Gen(jj).Psifd(ii,1) +	dt*(1.5*Gen(jj).dotPsifd(ii,1) - 	0.5*Gen(jj).dotPsifd(ii-1,1));
    Gen(jj).Psikd(ii+1,1) = Gen(jj).Psikd(ii,1) + 	dt*(1.5*Gen(jj).dotPsikd(ii,1) -	0.5*Gen(jj).dotPsikd(ii-1,1));
    Gen(jj).Psi1q(ii+1,1) = Gen(jj).Psi1q(ii,1) +  	dt*(1.5*Gen(jj).dotPsi1q(ii,1) -  	0.5*Gen(jj).dotPsi1q(ii-1,1));
    Gen(jj).Psi2q(ii+1,1) = Gen(jj).Psi2q(ii,1) + 	dt*(1.5*Gen(jj).dotPsi2q(ii,1) -  	0.5*Gen(jj).dotPsi2q(ii-1,1));
    Gen(jj).Wr(ii+1,1) = 	Gen(jj).Wr(ii,1) +     	dt*(1.5*Gen(jj).dotWr(ii,1) -       0.5*Gen(jj).dotWr(ii-1,1));
    Gen(jj).RotAng(ii+1,1) = Gen(jj).RotAng(ii,1) +	dt*(1.5*Gen(jj).dotRotAng(ii,1) -   0.5*Gen(jj).dotRotAng(ii-1,1));
    clear e Te id iq idq rot
else
    error('Invalid Flag')
end

end
    
end
    