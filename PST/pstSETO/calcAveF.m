function calcAveF(k,flag)
%CALCAVEF calculates weighted average frequency for system and areas.
% CALCAVEF calculates weighted average frequency for system and areas
% based on generator speed and inertia.
% Total inertia for the system and each area is also calculated.
%
% Syntax: calcAveF(k,flag)
%
%   NOTES:  Flag option included incase more initialization operations
%           are desired in the future...
%           Does not account for tripped gens (for now).
%
%   Input:
%   k - data index
%   flag  - dictate what operation to perform
%       0 - Handle future intialization operations?
%       1 - calculate frequencies
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/18/20    11:11   Thad Haines     Version 1

global g

if flag == 0
    % all machines have bus number as col 2, MVA base as col 3
    % mac_em, mac_sub, and mac_tra have H as col 16
    
    % Induction models have an H value, but no speed calculations...?
    % mac_igen has H as col 9... (col 16 blank...)
    % mac_ind has H as col 9... (col 16 blank...)
    
    % Generator trips will affect system inertia.
    % this flag could be used to account for that and then the totH(k)
    % calculations in flag == 1 could be removed.
    
end

if flag == 1
    % Calculate and log system and area average weighted frequency according to index k
    
    if g.area.n_area == 0
        % calculate average system frequency only
        % calculate total system H NOTE: does not account for tripped gens
        g.sys.totH(k) = sum(g.mac.mac_con(:,16).*g.mac.mac_con(:,3));
        % calculate ave weighted f = sum(mac_speed.*MVA.*H) / totH
        g.sys.aveF(k) = sum(g.mac.mac_spd(:,k).*g.mac.mac_con(:,3).*g.mac.mac_con(:,16))/g.sys.totH(k);
    else
        % calculate weighted average frequency for each area, sum for system
        runningSysF = 0;
        for areaN=1:g.area.n_area
            % calculate area total inertia
            mNdx = g.area.area(areaN).macBusNdx;
            g.area.area(areaN).totH(k) = sum(g.mac.mac_con(mNdx,16).*g.mac.mac_con(mNdx,3));
            % calculate ave weighted f = sum(mac_speed.*MVA.*H) / totH
            g.area.area(areaN).aveF(k) = sum(g.mac.mac_spd(mNdx,k).*g.mac.mac_con(mNdx,3).*g.mac.mac_con(mNdx,16)) ...
                /g.area.area(areaN).totH(k);
            
            runningSysF = runningSysF + g.area.area(areaN).aveF(k);
        end
        
        g.sys.totH(k) = sum(g.mac.mac_con(:,16).*g.mac.mac_con(:,3));
        g.sys.aveF(k) = runningSysF/g.area.n_area;
    end
end
end% end function