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

% 08/28/20  12:35   Thad Haines     Trip a generator, then bring it back online
% exciter and governor ramp at same time

%% define global variables
global g

persistent excVrefNEW excVrefOLD % variables for exciter ramping
persistent wRef0 wRef1 wDelta r0 % variables for governor ramping

if kT<2
    tripOut = false(g.mac.n_mac,1);
    mac_trip_states = [0 0;0 0]; % to store two generators trip data...
else
    tripOut = tripStatus;
    
    %% Trip generator
    if abs(t(kT)-5)<1e-5
        tripOut(3) = true; %trip gen 1 at t=5 sec.
        mac_trip_states(3,:) = [3; t(kT)]; %keep track of when things trip
        disp(['MAC_TRIP_LOGIC:  Tripping gen 3 at t = ' num2str(t(kT))])
        for n=0:1
            g.mac.pmech(3,kT+n) = 0; % set pmech to zero
        end
        
        % bypass governor
        g.tg.tg_pot(3,5) = 0.0;     % set Pref to zero
        r0 = g.tg.tg_con(3,4);      % store orginal 1/R
        g.tg.tg_con(3,4) = 0.0;     % set 1/R = 0
        reInitGov(3,kT)             % reset governor states
    end
    
    %% untrip gen
    if abs(t(kT)-15.0)<1e-5 %
        disp(['MAC_TRIP_LOGIC:  "Un-Tripping" gen 3 at t = ' num2str(t(kT))])
        tripOut(3) = false;
        mac_trip_states(3,:) = [3; t(kT)];  % keep track of when things trip
        g.mac.mac_trip_flags(3) = 0;        % set global flag to zero.
        
        % bypass exciter (and pss)
        g.exc.exc_bypass(3) = 1;            % set bypass flag
        excVrefOLD = g.exc.exc_pot(3,3);   	% save initial voltage reference
        reInitSub(3,kT)                     % init machine states and voltage to connected bus at index kT
    end
    
    %% ramp wref to wref0
    if abs(g.sys.t(kT)-20) < 1e-5
        disp(['MAC_TRIP_LOGIC:  reinit gov, start ramping wref at t = ', num2str(t(kT))])
        g.tg.tg_con(3,4) = r0;       % restore original 1/R value
        reInitGov(3,kT)              % re-init gov states
        wRef0 = g.tg.tg_con(3,3);    % wref0
        wRef1 = g.mac.mac_spd(3,kT); % current machine speed
        g.tg.tg_con(3,3) = wRef1;    % set reference to current speed
        wDelta = wRef0 - wRef1;      % amount to ramp in
    end
    
    if g.sys.t(kT)>= 20 && g.sys.t(kT)< 35      % ramp w ref to original value
        g.tg.tg_con(3,3) = wRef1+ wDelta*(1 - exp( 20-g.sys.t(kT) ) ); % concave down
    end
    
    if abs(t(kT)-35.0)<1e-5 % Reset governor w ref
        g.tg.tg_con(3,3) = wRef0;
        disp(['MAC_TRIP_LOGIC:  wref ramp in complete at t = ', num2str(t(kT)) ])
    end
    
    %% Re-connect exciter
    if abs(t(kT)-35.0)<1e-5     % remove bypass on exciter
        disp(['MAC_TRIP_LOGIC:  connecting exciter at t = ', num2str(t(kT))])
        reInitSmpExc(3,kT)      % re-init single exciter
        pss(3,kT,0)             % re-init pss
        g.exc.exc_bypass(3) = 0;% remove exciter bypass
    end
    
    %% Ramp exciter
    if abs(t(kT)-45.0)<1e-5 % ramp exciter reference voltage
        disp(['MAC_TRIP_LOGIC:  ramping exciter to original Vref at t = ', num2str(t(kT))])
        excVrefNEW = excVrefOLD - g.exc.exc_pot(3,3); % calculate difference to make up
        excVrefOLD = g.exc.exc_pot(3,3);
    end
    if t(kT)>=45 && t(kT) <65
        g.exc.exc_pot(3,3) =  excVrefOLD + (t(kT)-45)*excVrefNEW/20;
    end
    
end
end

