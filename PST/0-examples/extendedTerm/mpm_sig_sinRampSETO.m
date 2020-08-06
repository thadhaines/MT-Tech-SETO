function mpm_sig(k)
% MPM_SIG defines modulation signal for generator mechanical power
% Syntax: mpm_sig(k)
%
%   NOTE: This is a user defined file that can be empty....
%   A history section seems a bit much...
%
%   History:
%   Date        Time    Engineer        Description
%   07/xx/98    12:37   Graham Rogers   Version 1
%   06/05/20    16:21   Thad Haines     V2 - using global g, no t passed
%   in, no dummy varibale f passed out. Default behavior commented out.
%
%   Purpose:    Ramp generator 3's mechanical power to simulate to a solar ramp
%               Meant to give an estimate of run times between versions/methods
% 4 minutes (240 seconds)
% 0-30      - no action
% 30-90     - ramp up
% 90-150    - hold peak
% 150-210   - ramp down
% 210-240   - no action

% cloud cover events
% 45-55 - 20% max gen
% 120-140 - 30% cover
% 180-190 - 15% cover

global g

if g.mac.n_pm~=0

    % Handle basic generation curve    
    maxGen = 0.5;
    
    if g.sys.t(k) >= 30 && g.sys.t(k) < 90
        % up ramp
        % use first 1/4 cycle of a sin to ramp up
        g.mac.pm_sig(3,k) = maxGen*sin( (g.sys.t(k)-30) * 2*pi/240);
    elseif g.sys.t(k) >= 150 && g.sys.t(k) < 210
        % down ramp
        % use second 1/4 cycle of a sin to ramp down
        g.mac.pm_sig(3,k) = maxGen*sin( (g.sys.t(k)-150) * 2*pi/240 + pi/2);
    elseif g.sys.t(k) >= 90 && g.sys.t(k) < 150
        % held peak
        g.mac.pm_sig(3,k) = maxGen;
    else
        % no generation
        g.mac.pm_sig(3,k) = 0;
    end
    
    % Handle Cloud cover generation reductions
    if g.sys.t(k) >= 45 && g.sys.t(k) < 55
        g.mac.pm_sig(3,k) = maxGen*.2;
    elseif g.sys.t(k) >= 120 && g.sys.t(k) < 140
        g.mac.pm_sig(3,k) = g.mac.pm_sig(3,k)*0.7;
    elseif g.sys.t(k) >= 180 && g.sys.t(k) < 190
        g.mac.pm_sig(3,k) = g.mac.pm_sig(3,k)*0.85;
    end
    
    
end

return
