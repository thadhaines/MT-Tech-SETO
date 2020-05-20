function [Gen,Bus] = funCalculate_V_I_Pe(Gen,Bus,Ybus,ii)

    % Calculate voltages.
    E = nan(length(Gen),1);
    for jj = 1:length(Gen)
        E(jj) = Gen(jj).E(ii,1);
    end
    clear jj

    V = Ybus.Yrecov * E;
    
    % Calculate bus voltages.
    for jj = 1:length(Bus)
        Bus(jj).V(ii,1) = V(jj); 
    end
    clear jj V    

    % Set gen voltages.
    for jj = 1:length(Gen)
        for kk = 1:length(Bus)
            if Bus(kk).Num == Gen(jj).BusNum
                Gen(jj).Vt(ii,1) = Bus(kk).V(ii,1);
            end
        end
    end
    clear jj kk

    % Calculate gen currents.
    I = Ybus.Yred * E;
    
    % Set gen currents.
    for jj = 1:length(Gen)
        Gen(jj).I(ii,1) = I(jj);
    end
    clear jj E I

    % Calculate electrical power, Pe.
    for jj = 1:length(Gen)
        Gen(jj).Pe(ii,1) = real( Gen(jj).Vt(ii,1) * conj(Gen(jj).I(ii,1)) ); 
    end
    
end