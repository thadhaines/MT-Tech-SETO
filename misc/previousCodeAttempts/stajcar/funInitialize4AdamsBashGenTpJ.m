function [Gen,Bus] = funInitialize4AdamsBashGenTpJ(Gen,Bus,ii)
    
    % The Adams-Bashforth integration method uses the derivative caculation
    % from the previous time step. Therefore, the initial values calculated
    % from the power flow are set equal to the next time step,
    % Time.SimTime(ii+1), and one is added to the itteration step, ii+1.

    for jj = 1:length(Gen)
        Gen(jj).Vt(ii+1,1) = Gen(jj).Vt(ii);
        Gen(jj).I(ii+1,1) = Gen(jj).I(ii);
        Gen(jj).E(ii+1,1) = Gen(jj).E(ii,1);
        Gen(jj).RotAng(ii+1,1) = Gen(jj).RotAng(ii,1);
        Gen(jj).Pe(ii+1,1) = Gen(jj).Pe(ii,1);
        Gen(jj).Pm(ii+1,1) = Gen(jj).Pm(ii,1);
        Gen(jj).RotAng(ii+1,1) = Gen(jj).RotAng(ii,1);
        Gen(jj).Wr(ii+1,1) = Gen(jj).Wr(ii,1);
        
        Gen(jj).dotPsid_pp(ii+1,1) = 0;
        Gen(jj).dotPsiq_pp(ii+1,1) = 0;
        Gen(jj).dotEd_p(ii+1,1) = 0;
        Gen(jj).dotEq_p(ii+1,1) = 0;
        Gen(jj).dotWr(ii+1,1) = 0;
        Gen(jj).dotRotAng(ii+1,1) = 0;
        
        Gen(jj).Psid_pp(ii+1,1) = Gen(jj).Psid_pp(ii,1);
        Gen(jj).Psiq_pp(ii+1,1) = Gen(jj).Psiq_pp(ii,1);
        Gen(jj).Eq_p(ii+1,1) = Gen(jj).Eq_p(ii,1);
        Gen(jj).Ed_p(ii+1,1) = Gen(jj).Ed_p(ii,1);
        Gen(jj).Wr(ii+1,1) = Gen(jj).Wr(ii,1);
        Gen(jj).RotAng(ii+1,1) = Gen(jj).RotAng(ii,1);
    end
    clear jj
    
    for jj = 1:length(Bus)
        Bus(jj).V(ii+1,1) = Bus(jj).V(ii,1); 
    end
    clear jj
    
end