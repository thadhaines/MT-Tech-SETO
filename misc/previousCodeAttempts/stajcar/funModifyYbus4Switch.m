function [Ybus] = funModifyYbus4Switch(Switch,Ybus)

    switch Switch(1).Type
        
        case 1 % Bus fault
            
            Ybus.N = Switch.BusNumNear(1);
            Ybus.PostSwitchNN = Ybus.Y11(Ybus.N,Ybus.N);
            Ybus.Y11(Ybus.N,Ybus.N) = 1e6; 
        
        case 2 % Open branch.
            
            Ybus.N = Switch.BusNumNear(1);
            Ybus.PostSwitchNN = Ybus.Y11(Ybus.N,Ybus.N);
            Ybus.Y11(Ybus.N,Ybus.N) = 1e6;  
            
    end
    
    [Ybus] = funYredYrecovCalc(Ybus);

end