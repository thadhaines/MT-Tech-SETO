function [f] = mac_ivm(i,k,bus,flag)
% Syntax: [f] = mac_ivm(i,k,bus,flag)
%
% Purpose: Internal Voltage Model type generator
% 
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation and state
%                   state matrix building
%
% Output: f - dummy variable 
%
% Files:
%
% See Also:

% History (in reverse chronological order)
%
% Version: 1.0
% Date:    June 2019
% Author:  D. Trudnowski

% system variables
global basmva bus_v psi_re psi_im bus_int
global n_ivm mac_ivm_idx mac_pot mac_con eterm theta pelect qelect pmech

% IVM states and inputs
global mac_ang mach_ref edprime dmac_ang dedprime
global ivmmod_d_sig ivmmod_e_sig

f = 0;
if n_ivm ~=0
    if flag == 0 % initialization
        if i~=0
          % non-vector calculation
          error(' ')
        else
          % vectorized computation
          busnum = bus_int(mac_con(mac_ivm_idx,2)); % bus numbers where ivm's are connected 
          mac_pot(mac_ivm_idx,1) = basmva*ones(n_ivm,1)./mac_con(mac_ivm_idx,3); % factor to convert between MVA basis
          mac_pot(mac_ivm_idx,2) = 1.0*ones(n_ivm,1); % base kv?
          eterm(mac_ivm_idx,1) = bus(busnum,2);  % terminal bus voltage magnitude
          theta(busnum,1) = bus(busnum,3)*pi/180;% terminal bus angle in radians
          pelect(mac_ivm_idx,1) = bus(busnum,4).*mac_con(mac_ivm_idx,22); %electrical output power, active
          qelect(mac_ivm_idx,1) = bus(busnum,5).*mac_con(mac_ivm_idx,23); %electrical output power, reactive
          Vt = eterm(mac_ivm_idx,1).*exp(1j*theta(busnum,1)); %Terminal voltage
          It = conj((pelect(mac_ivm_idx,1) + 1j*qelect(mac_ivm_idx,1))./Vt); %Current out of terminal
          E = Vt + (mac_con(mac_ivm_idx,5) + 1j*mac_con(mac_ivm_idx,7))*It; %Internal voltage
          edprime(mac_ivm_idx,1) = abs(E); %Internal voltage magnitude
          mac_ang(mac_ivm_idx,1) = atan2(imag(E),real(E)); %Internal voltage angle (rad.)
          psi_re(mac_ivm_idx,1) = real(E); %real-part of internal voltage used for network interface
          psi_im(mac_ivm_idx,1) = imag(E); %imag-part of internal voltage used for network interface
          %pmech(mac_ivm_idx,1) = pelect(mac_ivm_idx,1) + sign(pelect(mac_ivm_idx,1)).*(abs(It).^2).*mac_con(mac_ivm_idx,5);
        end
    end

    if flag == 1 % network interface computation 
        if i ~= 0
          % non-vector calculation
          error(' ')
        else
          % vectorized computation
          mac_ang(mac_ivm_idx,k) = mac_ang(mac_ivm_idx,k)-mach_ref(k)*ones(n_ivm,1); % wrt machine reference
          E = edprime(mac_ivm_idx,k).*exp(1j*mac_ang(mac_ivm_idx,k));
          psi_re(mac_ivm_idx,k) = real(E); %real part of internal voltage
          psi_im(mac_ivm_idx,k) = imag(E); %imag part of internal voltage
        end
    end
 
    if flag == 2 % generator dynamics calculation
        if i ~= 0
          % non-vector calculation
          error(' ')
        else
          % vectorized computation
          dmac_ang(mac_ivm_idx,k) = (-mac_ang(mac_ivm_idx,k) + ivmmod_d_sig(:,k))./mac_con(mac_ivm_idx,9);
          dedprime(mac_ivm_idx,k) = (-edprime(mac_ivm_idx,k) + ivmmod_e_sig(:,k))./mac_con(mac_ivm_idx,10);
          busnum = bus_int(mac_con(mac_ivm_idx,2)); % bus numbers where ivm's are connected
          Vt = bus_v(busnum,k);
          E = edprime(mac_ivm_idx,k).*exp(1j*mac_ang(mac_ivm_idx,k));
          It = (E - Vt)./(mac_con(mac_ivm_idx,5) + 1j*mac_con(mac_ivm_idx,7));
          S = Vt.*conj(It);
          pelect(mac_ivm_idx,k) = real(S);
          qelect(mac_ivm_idx,k) = imag(S);
          eterm(mac_ivm_idx,k) = abs(Vt);
        end
    end
 
end
end
