function [Switch] = funSwitching()   
    
    Switch.TstartSwitch(1) = 1; 
    Switch.Type(1) = 1;     % Bus fault.
    Switch.TclearNear(1) = 1/60 + Switch.TstartSwitch(1);
    Switch.TclearFar(1) = 1/60 + Switch.TstartSwitch(1);
    Switch.BusNumNear(1) = 3;
    Switch.BusNumFar(1) = 3;

end