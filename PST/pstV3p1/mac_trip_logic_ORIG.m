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
global bus_v %pu bus voltages; bus_v(n,kT) = bus n voltage at time t(kT)
global mac_spd %Gen pu speeds; mac_spd(n,kT) = gen n speed at time t(kT)
global pelect qelect %gen pu powers; pelect(n,kT) = gen n real power at time t(kT)
global cur_re cur_im %gen pu currents
global mac_con n_mac

%% Initialize
tripOut = false(n_mac,1);
mac_trip_states = 0;


end
