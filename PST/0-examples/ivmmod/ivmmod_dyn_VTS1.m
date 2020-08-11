function [dc,Ec,ddelta_states,dE_states,delta_statesIni,E_statesIni] = ivmmod_dyn(delta_states,E_states,bus,Time,kSim,Flag)
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
%   dc = n_ivm by 1 column vector of ivmmod_d_sig commands at t = Time(kSim).
%       ivmmod_d_sig set the mac_ang for an IVM generator in mac_ivm.m
%       ivmmod_d_sig(k,kSim) = d(k).  k corresponds to the kth IVM in mac_con.
%   Ec = n_ivm by 1 column vector of ivmmod_e_sig commands at t = Time(kSim).
%       ivmmod_e_sig set the edprime for an IVM generator in mac_ivm.m
%       ivmmod_e_sig(k,kSim) = E(k).  k corresponds to the kth IVM in mac_con.
%
% Global:
%   ivmmod_data = general variable for storing data when necessary.
%
% D. Trudnowski, 2020

% global ivmmod_data
% global n_ivm mac_ivm_idx

global g
%% Parameters
busnum = g.bus.bus_int(g.mac.mac_con(g.ivm.mac_ivm_idx,2)); % bus numbers where ivm's are connected

%% Initialize output variables
dc = zeros(g.ivm.n_ivm,1);
Ec = zeros(g.ivm.n_ivm,1);
ddelta_states = cell(g.ivm.n_ivm,1);
dE_states = cell(g.ivm.n_ivm,1);
delta_statesIni = cell(g.ivm.n_ivm,1);
E_statesIni = cell(g.ivm.n_ivm,1);

%% Define and initialize state derivatives at t = 0.
if Flag==0
    %Note, there are no differential equations.  Just make the DFE zero.
    g.ivm.ivmmod_data = nan; %not used
    for k=1:g.ivm.n_ivm
        delta_statesIni{k}(1,1) = 0;
        E_statesIni{k}(1,1) = 0;
    end
    
    
    %% Calculate dc and Ec
elseif Flag==1
    for k=1:g.ivm.n_ivm
        if Time(kSim) >=1 && Time(kSim) <2
            dE = 0.1 *sin( (Time(kSim)-1) * 2*pi/4); % amount to step voltage (pu)
        elseif Time(kSim) >= 2
            dE = 0.1;
        else
            dE = 0;
        end
        Ec(k) = g.mac.edprime(g.ivm.mac_ivm_idx(k),1) + dE; %edprime(*,1) is initial internal voltage E
        if Time(kSim) >= 10
            % sine ramp signal up over 5 seconds.
            if Time(kSim) >=10 && Time(kSim) <15
                da = 0.5 *sin( (Time(kSim)-10) * 2*pi/20); % signal to step angle by (radians)
            else
                da = 0.5;
            end
        else
            da = 0;
        end
        dc(k) = g.mac.mac_ang(g.ivm.mac_ivm_idx(k),1) + da; %mac_ang(*,1) is initial internal voltage d
    end
    
    %% Calculate derivatives
elseif Flag==2
    %No differential equations for this example.  set derivatives to zero.
    for k=1:g.ivm.n_ivm
        ddelta_states{k}(1,1) = 0;
        dE_states{k}(1,1) = 0;
    end
    
end

end

