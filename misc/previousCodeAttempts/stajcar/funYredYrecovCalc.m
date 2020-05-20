function [Ybus] = funYredYrecovCalc(Ybus)

    Ybus.Yrecov = -Ybus.Y11\Ybus.Y12;
    Ybus.Yred = Ybus.Y22 - Ybus.Y12.'*(Ybus.Y11\Ybus.Y12);  % When 
    % calculating the transpose of a complex matrix A, the expression A.'
    % must be used instead of A' in order to take the transpose WITHOUT
    % calculating the complex conjugate.
    
end