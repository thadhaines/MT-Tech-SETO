function [delta,E,ddelta_states,dE_states,delta_statesIni,E_statesIni] = ivmmod_dyn(delta_states,E_states,bus,Time,kSim,Flag)
% Implement state or output variables to model power injection
% 
% Inputs:
%   Time = vector of simulation time
%   kSim = current simulation index.  Current time = Time(kSim).
%   Flag:
%       If Flag==0, Initialize d_statesIni, E_statesIni at t = 0.
%       If Flag==1, Calculate d, E at Time(kSim)
%       If Flag==2, Calculate dd_states, dE_states at Time(kSim)
%   d_state = cell array of mac_ang injection states.  Used to set d.
%       d_state{k} = column vector of states corresponding to the kth IVM in
%       mac_con.
%   E_state = cell array of E injection states.  Same format as d_states.
%   bus = initial bus matrix from the solved power flow.
%   
% Outputs:
%   d_statesIni = = cell array of initial of d_states
%   E_statesIni = = cell array of initial of E_states
%   dd_state = cell array of d/dt of d_states.
%   dE_state = cell array of d/dt of E_states.
%   d = n_ivm by 1 column vector of ivmmod_d_sig commands at t = Time(kSim).
%       ivmmod_d_sig set the mac_ang for an IVM generator in mac_ivm.m
%       ivmmod_d_sig(k,kSim) = d(k).  k corresponds to the kth IVM in mac_con.
%   E = n_ivm by 1 column vector of ivmmod_e_sig commands at t = Time(kSim).
%       ivmmod_e_sig set the edprime for an IVM generator in mac_ivm.m
%       ivmmod_e_sig(k,kSim) = E(k).  k corresponds to the kth IVM in mac_con.
%
% Global:
%   ivmmod_data = general variable for storing data when necessary.
%
% D. Trudnowski, 2019

global ivmmod_data
global bus_v bus_int edprime mac_ang pelect qelect
global n_ivm mac_ivm_idx mac_con

%% Parameters
busnum = bus_int(mac_con(mac_ivm_idx,2)); % bus numbers where ivm's are connected
nOrderd = 1; %order of state equations for d modulation
nOrderE = 1; %order of state equations for E modulation
Te = 0.05*ones(n_ivm,1); %Time constant voltage controller
Ke = 10*ones(n_ivm,1); %gain of voltage controller
R = 0.05*ones(n_ivm,1); %droop gain

%% Initialize output variables
delta = zeros(n_ivm,1);
E = zeros(n_ivm,1);
ddelta_states = cell(n_ivm,1);
dE_states = cell(n_ivm,1);
delta_statesIni = cell(n_ivm,1);
E_statesIni = cell(n_ivm,1);

%% Define and initialize state derivatives at t = 0.
if Flag==0
    ivmmod_data = zeros(n_ivm,2); %Store Vref and Pref in ivmmod_data
    Vt = bus(busnum,2).*exp(1j*bus(busnum,3).*pi./180); %Terminal voltage
    %It = conj((pelect(mac_ivm_idx,1) + 1j*qelect(mac_ivm_idx,1))./Vt); %Current out of terminal
    for k=1:n_ivm
        E_statesIni{k}(1) = edprime(mac_ivm_idx(k),1);
        ivmmod_data(k,1) = abs(Vt(k)) + E_statesIni{k}./Ke(k); %initial Vref
        delta_statesIni{k}(1) = mac_ang(mac_ivm_idx(k),1);
        %ivmmod_data(k,2) =  pelect(mac_ivm_idx(k),1)+ sign(pelect(mac_ivm_idx(k),1))*mac_con(mac_ivm_idx(k),5)*abs(It(k))^2; %initial Pref
        ivmmod_data(k,2) =  pelect(mac_ivm_idx(k),1); %initial Pref
    end
    clear k Vt It

%% Calculate delta and E
elseif Flag==1
    for k=1:n_ivm
        E(k) = E_states{k}(1);
        delta(k) = delta_states{k}(1);
    end
    %Not used

%% Calculate derivatives
elseif Flag==2
    %Voltage controller
    Vt = abs(bus_v(busnum,kSim));
    for k=1:n_ivm
        dE_states{k}(1) = (-E_states{k}(1) + Ke(k)*(ivmmod_data(k,1) - Vt(k))) / Te(k);
    end
    
    %droop controller
    wr = 1 - R.*(pelect(mac_ivm_idx,kSim) - ivmmod_data(:,2));
    if kSim>1; delT = Time(kSim)-Time(kSim-1);
    else delT = Time(kSim);
    end
    if abs(delT)>1e-8 && kSim>2
        w = angle(bus_v(busnum,kSim)./bus_v(busnum,kSim-1))./delT;
        w = 1 + w./(2*pi*60); %freq in pu
    else
        w = wr;
    end
    ddelta_states{k}(1) = wr - w;
end


end

