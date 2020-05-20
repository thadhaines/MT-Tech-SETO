function [Gen,Bus] = funInitialize4AdamsBashDansPSTN2(Gen,Bus,ii)
    
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
        
        Gen(jj).dotPsikd(ii+1,1) = 0;
        Gen(jj).dotPsikq(ii+1,1) = 0;
        Gen(jj).dotEdp(ii+1,1) = 0;
        Gen(jj).dotEqp(ii+1,1) = 0;
        Gen(jj).dotWr(ii+1,1) = 0;
        Gen(jj).dotRotAng(ii+1,1) = 0;
        
        Gen(jj).Eqp(ii+1,1) = Gen(jj).Eqp(ii,1);
        Gen(jj).Edp(ii+1,1) = Gen(jj).Edp(ii,1);
        Gen(jj).Psikd(ii+1,1) = Gen(jj).Psikd(ii,1);
        Gen(jj).Psikq(ii+1,1) = Gen(jj).Psikq(ii,1);
        Gen(jj).Wr(ii+1,1) = Gen(jj).Wr(ii,1);
        Gen(jj).RotAng(ii+1,1) = Gen(jj).RotAng(ii,1);
    end
    clear jj
    
    for jj = 1:length(Bus)
        Bus(jj).V(ii+1,1) = Bus(jj).V(ii,1); 
    end
    clear jj
    
end