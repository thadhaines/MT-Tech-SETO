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
% exciter and governor ramp over 40-65

%% define global variables
global g

persistent excVrefNEW excVrefOLD

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
        g.tg.tg_pot(3,5) = 0.0; % set Pref to zero
        g.tg.tg_con(3,4) = 0.0; % set 1/R = 0
        reInitGov(3,kT) % reset governor states
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
    
    
    %% ramp R in
    if abs(g.sys.t(kT)-20) < 1e-5
        disp(['MAC_TRIP_LOGIC:  reinit gov, start ramping R in at t = ', num2str(t(kT))])
       reInitGov(3,kT) 
    end
    if g.sys.t(kT)>= 20 && g.sys.t(kT)< 25 
        g.tg.tg_con(3,4) = 20*(1 - exp( 20-g.sys.t(kT) ) ); % concave down
        %g.tg.tg_con(3,4) = (g.sys.t(kT)-20)*20/5; % 5 second ramp up linear ramp
    end
    
    if abs(t(kT)-25.0)<1e-5 % Reset governor delta w gain (keep Pref = 0)
       % Remove bypass of governor R
       g.tg.tg_con(3,4) = 20.0; % restore 1/R value
        disp(['MAC_TRIP_LOGIC:  R ramp in complete, allow governor to account for frequency deviation at t = ', num2str(t(kT)) ])
    end
    
    %% Re-connect exciter
    if abs(t(kT)-35.0)<1e-5 % remove bypass on exciter
        disp(['MAC_TRIP_LOGIC:  connecting exciter at t = ', num2str(t(kT))])
        reInitSmpExc(3,kT) % re-init single exciter
        g.exc.exc_bypass(3) = 0; % remove exciter bypass
        
    end
    
    if abs(t(kT)-50.0)<1e-5 % ramp exciter reference voltage
        disp(['MAC_TRIP_LOGIC:  ramping exciter to original ref voltage at t = ', num2str(t(kT))])
        excVrefNEW = excVrefOLD - g.exc.exc_pot(3,3); % calculate difference to make up
        excVrefOLD = g.exc.exc_pot(3,3);
    end
    
    if t(kT)>=50 && t(kT) <75
         g.exc.exc_pot(3,3) =  excVrefOLD + (t(kT)-50)*excVrefNEW/25;
    end
    
end
end

