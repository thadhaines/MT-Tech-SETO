function mac_ivm(i,k,bus,flag)
%MAC_IVM Internal Voltage Model type generator
% MAC_IVM Internal Voltage Model type generator
%
% Syntax: mac_ivm(i,k,bus,flag)
%
%   NOTES:
%
%   Input:
%   i - generator number
%     	- 0 for vectorized computation
%   k - integer time
%   bus - solved loadflow bus data
% 	flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation and state matrix building
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   06/xx/19    XX:XX   D. Trudnowski  	Version 1.0
%   08/11/20    10:33   Thad Haines     Revised format of globals and internal function documentation

global g

if g.ivm.n_ivm ~=0
    if flag == 0 % initialization
        if i~=0
            % non-vector calculation
            error(' ')
        else
            % vectorized computation
            busnum = g.bus.bus_int(g.mac.mac_con(g.ivm.mac_ivm_idx,2)); % bus numbers where ivm's are connected
            g.mac.mac_pot(g.ivm.mac_ivm_idx,1) = g.sys.basmva*ones(g.ivm.n_ivm,1)./g.mac.mac_con(g.ivm.mac_ivm_idx,3); % factor to convert between MVA basis
            g.mac.mac_pot(g.ivm.mac_ivm_idx,2) = 1.0*ones(g.ivm.n_ivm,1); % base kv?
            g.mac.eterm(g.ivm.mac_ivm_idx,1) = bus(busnum,2);  % terminal bus voltage magnitude
            g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;% terminal bus angle in radians
            g.mac.pelect(g.ivm.mac_ivm_idx,1) = bus(busnum,4).*g.mac.mac_con(g.ivm.mac_ivm_idx,22); %electrical output power, active
            g.mac.qelect(g.ivm.mac_ivm_idx,1) = bus(busnum,5).*g.mac.mac_con(g.ivm.mac_ivm_idx,23); %electrical output power, reactive
            Vt = g.mac.eterm(g.ivm.mac_ivm_idx,1).*exp(1j*g.bus.theta(busnum,1)); %Terminal voltage
            It = conj((g.mac.pelect(g.ivm.mac_ivm_idx,1) + 1j*g.mac.qelect(g.ivm.mac_ivm_idx,1))./Vt); %Current out of terminal
            E = Vt + (g.mac.mac_con(g.ivm.mac_ivm_idx,5) + 1j*g.mac.mac_con(g.ivm.mac_ivm_idx,7))*It; %Internal voltage
            g.mac.edprime(g.ivm.mac_ivm_idx,1) = abs(E); %Internal voltage magnitude
            g.mac.mac_ang(g.ivm.mac_ivm_idx,1) = atan2(imag(E),real(E)); %Internal voltage angle (rad.)
            g.mac.psi_re(g.ivm.mac_ivm_idx,1) = real(E); %real-part of internal voltage used for network interface
            g.mac.psi_im(g.ivm.mac_ivm_idx,1) = imag(E); %imag-part of internal voltage used for network interface
            %pmech(mac_ivm_idx,1) = pelect(mac_ivm_idx,1) + sign(pelect(mac_ivm_idx,1)).*(abs(It).^2).*mac_con(mac_ivm_idx,5);
        end
    end
    
    if flag == 1 % network interface computation
        if i ~= 0
            % non-vector calculation
            error(' ')
        else
            % vectorized computation
            g.mac.mac_ang(g.ivm.mac_ivm_idx,k) = g.mac.mac_ang(g.ivm.mac_ivm_idx,k)-g.sys.mach_ref(k)*ones(g.ivm.n_ivm,1); % wrt machine reference
            E = g.mac.edprime(g.ivm.mac_ivm_idx,k).*exp(1j*g.mac.mac_ang(g.ivm.mac_ivm_idx,k));
            g.mac.psi_re(g.ivm.mac_ivm_idx,k) = real(E); %real part of internal voltage
            g.mac.psi_im(g.ivm.mac_ivm_idx,k) = imag(E); %imag part of internal voltage
        end
    end
    
    if flag == 2 % generator dynamics calculation
        if i ~= 0
            % non-vector calculation
            error(' ')
        else
            % vectorized computation
            g.mac.dmac_ang(g.ivm.mac_ivm_idx,k) = (-g.mac.mac_ang(g.ivm.mac_ivm_idx,k) + g.ivm.ivmmod_d_sig(:,k))./g.mac.mac_con(g.ivm.mac_ivm_idx,9);
            g.mac.dedprime(g.ivm.mac_ivm_idx,k) = (-g.mac.edprime(g.ivm.mac_ivm_idx,k) + g.ivm.ivmmod_e_sig(:,k))./g.mac.mac_con(g.ivm.mac_ivm_idx,10);
            busnum = g.bus.bus_int(g.mac.mac_con(g.ivm.mac_ivm_idx,2)); % bus numbers where ivm's are connected
            Vt = g.bus.bus_v(busnum,k);
            E = g.mac.edprime(g.ivm.mac_ivm_idx,k).*exp(1j*g.mac.mac_ang(g.ivm.mac_ivm_idx,k));
            It = (E - Vt)./(g.mac.mac_con(g.ivm.mac_ivm_idx,5) + 1j*g.mac.mac_con(g.ivm.mac_ivm_idx,7));
            S = Vt.*conj(It);
            g.mac.pelect(g.ivm.mac_ivm_idx,k) = real(S);
            g.mac.qelect(g.ivm.mac_ivm_idx,k) = imag(S);
            g.mac.eterm(g.ivm.mac_ivm_idx,k) = abs(Vt);
        end
    end
    
end
end
