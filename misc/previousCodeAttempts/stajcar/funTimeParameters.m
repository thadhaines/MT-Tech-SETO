function [Time] = funTimeParameters()

    Time.Tstart = 0;
    Time.SimEnd = 100-1/600;
    
    Time.dtSubTrans = 1/600;
    Time.dtSubTrans = 1/240;
%     Time.dtClassical = 1/600;
    Time.dtSlow = 1;
    
    Time.SimTime = 0;   
    
end