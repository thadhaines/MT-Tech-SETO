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

%% define global variables
global g

if kT<2
    tripOut = false(g.mac.n_mac,1);
    mac_trip_states = [0 0;0 0]; % to store two generators trip data...
else
    tripOut = tripStatus;
    if abs(t(kT)-2)<1e-5
        tripOut(3) = true; %trip gen 1 at t=5 sec.
        mac_trip_states(3,:) = [3; t(kT)]; %keep track of when things trip
        disp(['Tripping gen 3 at t = ' num2str(t(kT))])
    end
    
    %untrip gen
    if abs(t(kT)-25)<1e-5
        tripOut(3) = false; %trip gen 1 at t=5 sec.
        mac_trip_states(3,:) = [3; t(kT)]; %keep track of when things trip
        disp(['Un-Tripping gen 3 at t = ' num2str(t(kT))])
        g.mac.mac_trip_flags = zeros(size(g.mac.mac_con,1),1); % set global flags to zero.
        mac_sub(3,kT,g.bus.bus,0) % re-init single gen
%         % testing of handle states
%         for n=0:1
%         g.mac.mac_spd(3,kT+n) = g.mac.mac_spd(1,kT-1)*1.03; % 3 percent faster...
%         
%         %g.mac.pelect(3,kT+n) = 0;
%         %g.mac.qelect(3,kT+n) = 0;
%         %g.mac.pmech(3,kT+n) = 0;
%         end
    end
end

end
