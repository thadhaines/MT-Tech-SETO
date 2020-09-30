function [tripOut,mac_trip_states] = mac_trip_logic(tripStatus,mac_trip_states,t,kT)
% Purpose: trip generators.
%
% Inputs:
%   tripStatus = n_mac x 1 bool vector of current trip status.  If
%       tripStatus(n) is true, then the generator corresponding to the nth
%       row of mac_con is already tripped.  Else, it is false.
%   mac_trip_states = storage matrix defined by user.
%   t = vector of simulation time (sec.).
%   kT = current integer time (sample).  Corresponds to t(kT)
%
% Output:
%   tripOut = n_mac x 1 bool vector of desired trips.  If
%       tripOut(n)==1, then the generator corresponding to the nth
%       row of mac_con is will be tripped.  Note that each element of
%       tripOut must be either 0 or 1.

% Version 1.0
% Author:   Dan Trudnowski
% Date:   Jan 2017

% Modified 09/30/20 to use global g - Thad Haines

%% define global variables
global g
persistent t1 t2 % used to only display trip notifications once.

if kT<2
    tripOut = false(g.mac.n_mac,1);
    mac_trip_states = [0 0;0 0]; % to store two generators trip data...
    t1 = 0;
    t2 = 0;
else
    tripOut = tripStatus;
    if abs(t(kT)-5)<1e-6
        tripOut(1) = true; %trip gen 1 at t=5 sec.
        mac_trip_states(1,:) = [1; t(kT)]; %keep track of when things trip
        if ~t1
            disp(['Tripping gen 1 at t = ' num2str(t(kT))])
            t1 = 1;
        end
    end
    if abs(t(kT)-8)<1e-6
        tripOut(13) = true; %trip gen 13 at t=8 sec.
        mac_trip_states(2,:) = [13; t(kT)]; %keep track of when things trip
        if ~t2
            disp(['Tripping gen 13 at t = ' num2str(t(kT))])
            t2 = 1;
        end
    end
end

end
