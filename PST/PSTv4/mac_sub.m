function mac_sub(i,k,bus,flag)
%MAC_SUB Voltage-behind-subtransient-reactance generator model.
% MAC_SUB is a voltage-behind-subtransient-reactance generator model,
% with vectorized computation option.
%
% Syntax: mac_sub(i,k,bus,flag)
%
%   NOTES: State variables are: mac_ang, mac_spd, eqprime, psikd, edprime,
%   psikq in the g.mac. field
%   Algorithm: R. P. Schulz, "Synchronous machine modeling".
%
%   Input:
%   i - generator number
%          - 0 for vectorized computation
%   k - integer time (data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation and state state matrix building
%           3 - state matrix building
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   03/xx/91    XX:XX   Joe H. Chow     Version 1
%   05/06/95    xx:xx   GJR             Correction to saturation
%   06/xx/96    XX:XX   Graham Rogers   Version 2 - change to allow multple generator model types
%   09/xx/97    xx:xx   -               Version 2.1 - modified for low Eqprime
%   (c) Copyright 1991-1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved
%   06/19/20    10:16   Thad Haines     Revised format of globals and internal function documentation
%   07/29/20    15:20   Thad Haines     jay -> 1j

global g

if g.mac.n_sub~=0
    if flag == 0 % initialization
        if i ~= 0
            % check that subtransient machine
            sub = find(g.mac.mac_sub_idx==i);
            if ~isempy(sub)
                % check data
                % ensure xd" = xq"
                if g.mac.mac_con(i,8)~=g.mac.mac_con(i,13)
                    
                    disp(['making xdpp = xqpp for generator ',num2str(i)])
                    g.mac.mac_con(i,13) = g.mac.mac_con(i,8);
                end
                if g.mac.mac_con(i,14) == 0
                    g.mac.mac_con(i,14) = 999.0;
                end %default transient time constant
                if g.mac.mac_con(i,15) == 0
                    g.mac.mac_con(i,15) = 999.0;
                end %default subtransient time constant
                
                busnum = g.bus.bus_int(g.mac.mac_con(i,2)); % bus number
                g.mac.mac_pot(i,1)=g.sys.basmva/g.mac.mac_con(i,3); % scaled MVA base
                g.mac.mac_pot(i,2)=1.0; % base kv
                g.mac.mac_pot(i,8)=g.mac.mac_con(i,7)-g.mac.mac_con(i,4);
                g.mac.mac_pot(i,9)=(g.mac.mac_con(i,8)-g.mac.mac_con(i,4))/g.mac.mac_pot(i,8);
                g.mac.mac_pot(i,7)=(g.mac.mac_con(i,6)-g.mac.mac_con(i,7))*g.mac.mac_pot(i,9);
                g.mac.mac_pot(i,10)=(g.mac.mac_con(i,7)-g.mac.mac_con(i,8))...
                    /g.mac.mac_pot(i,8);
                g.mac.mac_pot(i,6)=(g.mac.mac_con(i,6)-g.mac.mac_con(i,7))...
                    /g.mac.mac_pot(i,8)*g.mac.mac_pot(i,10);
                g.mac.mac_pot(i,13)=g.mac.mac_con(i,12)-g.mac.mac_con(i,4);
                g.mac.mac_pot(i,14)=(g.mac.mac_con(i,13)-g.mac.mac_con(i,4))...
                    /g.mac.mac_pot(i,13);
                g.mac.mac_pot(i,12)=(g.mac.mac_con(i,11)-g.mac.mac_con(i,12))...
                    *g.mac.mac_pot(i,14);
                g.mac.mac_pot(i,15)=(g.mac.mac_con(i,12)-g.mac.mac_con(i,13))...
                    /g.mac.mac_pot(i,13);
                g.mac.mac_pot(i,11)=(g.mac.mac_con(i,11)-g.mac.mac_con(i,12))...
                    /g.mac.mac_pot(i,13)*g.mac.mac_pot(i,15);
                % extract bus information
                g.mac.eterm(i,1) = bus(busnum,2);  % terminal bus voltage
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;
                % terminal bus angle in radians
                g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22); % system base
                % electrical output power, active
                g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);
                % electrical output power, reactive
                curr = sqrt(g.mac.pelect(i,1)^2+g.mac.qelect(i,1)^2) ...
                    /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);  % current magnitude
                %on genenearor base
                phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));
                % power factor angle
                v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));
                % complex voltage
                % in system reference frame
                curr = curr*exp(1j*(g.bus.theta(busnum,1)-phi)); % complex current
                % in system reference frame
                ei = v + (g.mac.mac_con(i,5)+1j*g.mac.mac_con(i,11))*curr;
                g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));
                % machine angle (delta)
                g.mac.mac_spd(i,1) = 1; % machine speed at steady state
                rot = 1j*exp(-1j*g.mac.mac_ang(i,1));  % system reference frame rotation
                curr = curr*rot;
                g.mac.curdg(i,1) = real(curr);
                g.mac.curqg(i,1) = imag(curr);% current in Park's frame
                g.mac.curd(i,1) = real(curr)/g.mac.mac_pot(i,1);% current on system base
                g.mac.curq(i,1) = imag(curr)/g.mac.mac_pot(i,1);
                mcurmag = abs(curr); % current magnitude on machine base
                g.mac.pmech(i,1) = g.mac.pelect(i,1)*g.mac.mac_pot(i,1) + g.mac.mac_con(i,5)*(mcurmag*mcurmag);
                %pmech = g.mac.pelect + losses on machine base
                v = v*rot;
                g.mac.ed(i,1) = real(v);
                g.mac.eq(i,1) = imag(v);
                eqra = g.mac.eq(i,1)+g.mac.mac_con(i,5)*g.mac.curqg(i,1);
                g.mac.psidpp = eqra + g.mac.mac_con(i,8)*g.mac.curdg(i,1);
                g.mac.psikd(i,1) = eqra + g.mac.mac_con(i,4)*g.mac.curdg(i,1);
                g.mac.eqprime(i,1) = eqra + g.mac.mac_con(i,7)*g.mac.curdg(i,1);
                edra = -g.mac.ed(i,1)-g.mac.mac_con(i,5)*g.mac.curdg(i,1);
                g.mac.psiqpp = edra + g.mac.mac_con(i,13)*g.mac.curqg(i,1);
                g.mac.psikq(i,1) = edra + g.mac.mac_con(i,4)*g.mac.curqg(i,1);
                g.mac.edprime(i,1) = edra + g.mac.mac_con(i,12)*g.mac.curqg(i,1);
                % compute saturation
                inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
                b = [0.8 1+g.mac.mac_con(i,20) 1.2*(1+g.mac.mac_con(i,21))];
                g.mac.mac_pot(i,3) = b*inv_sat(1,:)';
                g.mac.mac_pot(i,4) = b*inv_sat(2,:)';
                g.mac.mac_pot(i,5) = b*inv_sat(3,:)';
                E_Isat = g.mac.mac_pot(i,3)*g.mac.eqprime(i,1)^2 ...
                    + g.mac.mac_pot(i,4)*g.mac.eqprime(i,1) + g.mac.mac_pot(i,5);
                if g.mac.eqprime(i,1)<0.8
                    E_Isat=g.mac.eqprime(i,1);
                end
                g.mac.vex(i,1) = E_Isat + g.mac.mac_pot(i,6)*(g.mac.eqprime(i,1)-...
                    g.mac.psikd(i,1))+g.mac.mac_pot(i,7)*g.mac.curdg(i,1);
                g.mac.fldcur(i,1) = g.mac.vex(i,1);
                g.mac.psi_re(i,1) = sin(g.mac.mac_ang(i,1)).*(-g.mac.psiqpp) + ...
                    cos(g.mac.mac_ang(i,1)).*g.mac.psidpp; % real part of psi
                g.mac.psi_im(i,1) = -cos(g.mac.mac_ang(i,1)).*(-g.mac.psiqpp) + ...
                    sin(g.mac.mac_ang(i,1)).*g.mac.psidpp; % imag part of psi
            end
        else
            % vectorized computation
            % check parameters
            uets_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,8)~=g.mac.mac_con(g.mac.mac_sub_idx,13));
            if ~isempty(uets_idx)
                g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),13)=g.mac.mac_con(g.mac.mac_sub_idx(uets_idx),8);
                disp('xqpp made equal to xdpp at generators  '); disp((g.mac.mac_sub_idx(uets_idx))')
            end
            notp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,14)==0);
            if ~isempty(notp_idx)
                g.mac.mac_con(g.mac.mac_sub_idx(notp_idx),14) = 999.0*ones(length(notp_idx),1);
            end
            notpp_idx = find(g.mac.mac_con(g.mac.mac_sub_idx,15)==0);
            if ~isempty(notpp_idx)
                g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),15) = 999.0*ones(length(notpp_idx),1);
                % set x'q = x"q
                g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),12) =...
                    g.mac.mac_con(g.mac.mac_sub_idx(notpp_idx),13);
            end
            busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_sub_idx,2)); % bus number
            g.mac.mac_pot(g.mac.mac_sub_idx,1) = g.sys.basmva*ones(g.mac.n_sub,1)./g.mac.mac_con(g.mac.mac_sub_idx,3);
            % scaled MVA base
            g.mac.mac_pot(g.mac.mac_sub_idx,2) = ones(g.mac.n_sub,1); % base kv
            g.mac.mac_pot(g.mac.mac_sub_idx,8)=g.mac.mac_con(g.mac.mac_sub_idx,7)-g.mac.mac_con(g.mac.mac_sub_idx,4);
            g.mac.mac_pot(g.mac.mac_sub_idx,9)=(g.mac.mac_con(g.mac.mac_sub_idx,8)-g.mac.mac_con(g.mac.mac_sub_idx,4))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,8);
            g.mac.mac_pot(g.mac.mac_sub_idx,7)=(g.mac.mac_con(g.mac.mac_sub_idx,6)-g.mac.mac_con(g.mac.mac_sub_idx,7))...
                .*g.mac.mac_pot(g.mac.mac_sub_idx,9);
            g.mac.mac_pot(g.mac.mac_sub_idx,10)=(g.mac.mac_con(g.mac.mac_sub_idx,7)-g.mac.mac_con(g.mac.mac_sub_idx,8))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,8);
            g.mac.mac_pot(g.mac.mac_sub_idx,6)=(g.mac.mac_con(g.mac.mac_sub_idx,6)-g.mac.mac_con(g.mac.mac_sub_idx,7))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,8).*g.mac.mac_pot(g.mac.mac_sub_idx,10);
            g.mac.mac_pot(g.mac.mac_sub_idx,13)=g.mac.mac_con(g.mac.mac_sub_idx,12)-g.mac.mac_con(g.mac.mac_sub_idx,4);
            g.mac.mac_pot(g.mac.mac_sub_idx,14)=(g.mac.mac_con(g.mac.mac_sub_idx,13)-g.mac.mac_con(g.mac.mac_sub_idx,4))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,13);
            g.mac.mac_pot(g.mac.mac_sub_idx,12)=(g.mac.mac_con(g.mac.mac_sub_idx,11)-g.mac.mac_con(g.mac.mac_sub_idx,12))...
                .*g.mac.mac_pot(g.mac.mac_sub_idx,14);
            g.mac.mac_pot(g.mac.mac_sub_idx,15)=(g.mac.mac_con(g.mac.mac_sub_idx,12)-g.mac.mac_con(g.mac.mac_sub_idx,13))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,13);
            g.mac.mac_pot(g.mac.mac_sub_idx,11)=(g.mac.mac_con(g.mac.mac_sub_idx,11)-g.mac.mac_con(g.mac.mac_sub_idx,12))...
                ./g.mac.mac_pot(g.mac.mac_sub_idx,13).*g.mac.mac_pot(g.mac.mac_sub_idx,15);
            % extract bus information
            g.mac.eterm(g.mac.mac_sub_idx,1) = bus(busnum,2);  % terminal bus voltage
            g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;
            % terminal bus angle in radians
            g.mac.pelect(g.mac.mac_sub_idx,1) = bus(busnum,4).*g.mac.mac_con(g.mac.mac_sub_idx,22);
            % electrical output power, active
            g.mac.qelect(g.mac.mac_sub_idx,1) = bus(busnum,5).*g.mac.mac_con(g.mac.mac_sub_idx,23);
            % electrical output power, reactive
            curr = sqrt(g.mac.pelect(g.mac.mac_sub_idx,1).^2+g.mac.qelect(g.mac.mac_sub_idx,1).^2) ...
                ./g.mac.eterm(g.mac.mac_sub_idx,1).*g.mac.mac_pot(g.mac.mac_sub_idx,1);
            % current magnitude on generator base
            phi = atan2(g.mac.qelect(g.mac.mac_sub_idx,1),g.mac.pelect(g.mac.mac_sub_idx,1));
            % power factor angle
            v = g.mac.eterm(g.mac.mac_sub_idx,1).*exp(1j*g.bus.theta(busnum,1));
            % voltage in real and imaginary parts
            % in system reference frame
            curr = curr.*exp(1j*(g.bus.theta(busnum,1)-phi));
            % complex current in system reference frame
            ei = v + (g.mac.mac_con(g.mac.mac_sub_idx,5)+1j*g.mac.mac_con(g.mac.mac_sub_idx,11)).*curr;
            % voltage behind sub-transient reactance in system frame
            g.mac.mac_ang(g.mac.mac_sub_idx,1) = atan2(imag(ei),real(ei));
            % machine angle (delta)
            g.mac.mac_spd(g.mac.mac_sub_idx,1) = ones(g.mac.n_sub,1);
            % machine speed at steady state
            rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_sub_idx,1));
            % system reference frame rotation to Park's frame
            curr = curr.*rot;
            % current on generator base in Park's frame
            mcurmag = abs(curr);
            g.mac.pmech(g.mac.mac_sub_idx,1) = g.mac.pelect(g.mac.mac_sub_idx,1).*g.mac.mac_pot(g.mac.mac_sub_idx,1)...
                + g.mac.mac_con(g.mac.mac_sub_idx,5).*(mcurmag.*mcurmag);
            % mechanical power = electrical power + losses on generator base
            g.mac.curdg(g.mac.mac_sub_idx,1) = real(curr);
            g.mac.curqg(g.mac.mac_sub_idx,1) = imag(curr);
            % d and q axis current on generator base
            g.mac.curd(g.mac.mac_sub_idx,1) = real(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1);
            g.mac.curq(g.mac.mac_sub_idx,1) = imag(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1);
            % d and q axis currents on system base
            v = v.*rot;% voltage in Park's frame
            g.mac.ed(g.mac.mac_sub_idx,1) = real(v);
            g.mac.eq(g.mac.mac_sub_idx,1) = imag(v);
            % d and q axis voltages in Park's frame
            eqra = g.mac.eq(g.mac.mac_sub_idx,1)+g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,1);
            % q axis voltage behind resistance
            g.mac.psidpp = eqra + g.mac.mac_con(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,1);
            g.mac.psikd(g.mac.mac_sub_idx,1) = eqra + g.mac.mac_con(g.mac.mac_sub_idx,4).*g.mac.curdg(g.mac.mac_sub_idx,1);
            g.mac.eqprime(g.mac.mac_sub_idx,1) = eqra + g.mac.mac_con(g.mac.mac_sub_idx,7).*g.mac.curdg(g.mac.mac_sub_idx,1);
            edra = -g.mac.ed(g.mac.mac_sub_idx,1)-g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curdg(g.mac.mac_sub_idx,1);
            g.mac.psiqpp = edra + g.mac.mac_con(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,1);
            g.mac.psikq(g.mac.mac_sub_idx,1) = edra + g.mac.mac_con(g.mac.mac_sub_idx,4).*g.mac.curqg(g.mac.mac_sub_idx,1);
            g.mac.edprime(g.mac.mac_sub_idx,1) = edra + g.mac.mac_con(g.mac.mac_sub_idx,12).*g.mac.curqg(g.mac.mac_sub_idx,1);
            % this is the negative of Edprime in block diagram
            % compute saturation
            inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
            b = [0.8*ones(g.mac.n_sub,1) ones(g.mac.n_sub,1)+g.mac.mac_con(g.mac.mac_sub_idx,20)...
                1.2*(ones(g.mac.n_sub,1)+g.mac.mac_con(g.mac.mac_sub_idx,21))];
            g.mac.mac_pot(g.mac.mac_sub_idx,3) = b*inv_sat(1,:)';
            g.mac.mac_pot(g.mac.mac_sub_idx,4) = b*inv_sat(2,:)';
            g.mac.mac_pot(g.mac.mac_sub_idx,5) = b*inv_sat(3,:)';
            E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3).*g.mac.eqprime(g.mac.mac_sub_idx,1).^2 ...
                + g.mac.mac_pot(g.mac.mac_sub_idx,4).*g.mac.eqprime(g.mac.mac_sub_idx,1) + g.mac.mac_pot(g.mac.mac_sub_idx,5);
            nosat_idx=find(g.mac.eqprime(g.mac.mac_sub_idx,1)<.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx)=g.mac.eqprime(g.mac.mac_sub_idx(nosat_idx),1);
            end
            g.mac.vex(g.mac.mac_sub_idx,1) = E_Isat + g.mac.mac_pot(g.mac.mac_sub_idx,6).*(g.mac.eqprime(g.mac.mac_sub_idx,1)-...
                g.mac.psikd(g.mac.mac_sub_idx,1))+g.mac.mac_pot(g.mac.mac_sub_idx,7).*g.mac.curdg(g.mac.mac_sub_idx,1);
            g.mac.fldcur(g.mac.mac_sub_idx,1) = g.mac.vex(g.mac.mac_sub_idx,1);
            g.mac.psi_re(g.mac.mac_sub_idx,1) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) + ...
                cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp; % real part of psi
            g.mac.psi_im(g.mac.mac_sub_idx,1) = -cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) + ...
                sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp; % imag part of psi
            % psi is in system base and is the voltage behind xpp
        end
        %end initialization
    end
    if flag == 1 % network interface computation
        if i ~= 0
            % check for subsynchronous machine
            sub = find(g.mac.mac_sub_idx == i,1);
            if ~isempty(sub)
                g.mac.mac_ang(i,k) = g.mac.mac_ang(i,k) - g.sys.mach_ref(k);
                % wrt machine referencek
                g.mac.psidpp = g.mac.mac_pot(i,9)*g.mac.eqprime(i,k) + ...
                    g.mac.mac_pot(i,10)*g.mac.psikd(i,k);
                g.mac.psiqpp = g.mac.mac_pot(i,14)*g.mac.edprime(i,k) + ...
                    g.mac.mac_pot(i,15)*g.mac.psikq(i,k);
                g.mac.psi_re(i,k) = sin(g.mac.mac_ang(i,k))*(-g.mac.psiqpp) + ...
                    cos(g.mac.mac_ang(i,k))*g.mac.psidpp; % real part of psi
                g.mac.psi_im(i,k) = -cos(g.mac.mac_ang(i,k))*(-g.mac.psiqpp) + ...
                    sin(g.mac.mac_ang(i,k))*g.mac.psidpp; % imag part of psi
            end
        else
            % vectorized computation
            
            g.mac.mac_ang(g.mac.mac_sub_idx,k) = g.mac.mac_ang(g.mac.mac_sub_idx,k)-g.sys.mach_ref(k)*ones(g.mac.n_sub,1);
            % wrt machine reference
            g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9).*g.mac.eqprime(g.mac.mac_sub_idx,k) + ...
                g.mac.mac_pot(g.mac.mac_sub_idx,10).*g.mac.psikd(g.mac.mac_sub_idx,k);
            g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14).*g.mac.edprime(g.mac.mac_sub_idx,k) + ...
                g.mac.mac_pot(g.mac.mac_sub_idx,15).*g.mac.psikq(g.mac.mac_sub_idx,k);
            g.mac.psi_re(g.mac.mac_sub_idx,k) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) + ...
                cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp; % real part of psi
            g.mac.psi_im(g.mac.mac_sub_idx,k) = -cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) + ...
                sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp; % imag part of psi
        end
        % end of interface
    end
    
    if flag == 2 || flag == 3 % generator dynamics calculation
        if i ~= 0
            %check for subsynchronous machine
            sub = find(g.mac.mac_sub_idx==i,1);
            if ~isempty(sub)
                g.mac.psiqpp = g.mac.mac_pot(i,14)*g.mac.edprime(i,k) + ...
                    g.mac.mac_pot(i,15)*g.mac.psikq(i,k);
                g.mac.psidpp = g.mac.mac_pot(i,9)*g.mac.eqprime(i,k) + ...
                    g.mac.mac_pot(i,10)*g.mac.psikd(i,k);
                g.mac.curd(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) - ...
                    cos(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k); % d-axis current
                g.mac.curq(i,k) = cos(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) + ...
                    sin(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k); % q-axis current
                g.mac.curdg(i,k) = g.mac.curd(i,k)*g.mac.mac_pot(i,1);
                g.mac.curqg(i,k) = g.mac.curq(i,k)*g.mac.mac_pot(i,1);
                mcurmag = abs(g.mac.curdg(i,k) + 1j*g.mac.curqg(i,k));
                E_Isat = g.mac.mac_pot(i,3)*g.mac.eqprime(i,k)^2 ...
                    + g.mac.mac_pot(i,4)*g.mac.eqprime(i,k) + g.mac.mac_pot(i,5);
                if g.mac.eqprime(k,1)<0.8
                    E_Isat=g.mac.eqprime(k,1);
                end
                g.mac.fldcur(i,k) = E_Isat + g.mac.mac_pot(i,6)...
                    *(g.mac.eqprime(i,k)-g.mac.psikd(i,k)) + g.mac.mac_pot(i,7)...
                    *g.mac.curdg(i,k);
                g.mac.deqprime(i,k) = (g.mac.vex(i,k)-g.mac.fldcur(i,k))/g.mac.mac_con(i,9);
                g.mac.dpsikd(i,k) = (-g.mac.psikd(i,k)+g.mac.eqprime(i,k)-g.mac.mac_pot(i,8)...
                    *g.mac.curdg(i,k))/g.mac.mac_con(i,10);
                g.mac.dedprime(i,k) = (-g.mac.edprime(i,k) - g.mac.mac_pot(i,11)...
                    *(g.mac.edprime(i,k)-g.mac.psikq(i,k)) - g.mac.mac_pot(i,12)...
                    *g.mac.curqg(i,k))/g.mac.mac_con(i,14);
                g.mac.dpsikq(i,k) = (g.mac.edprime(i,k)-g.mac.psikq(i,k)-g.mac.mac_pot(i,13)...
                    *g.mac.curqg(i,k))/g.mac.mac_con(i,15);
                g.mac.ed(i,k) = -g.mac.mac_con(i,5)*g.mac.curdg(i,k) - (g.mac.psiqpp...
                    -g.mac.mac_con(i,13)*g.mac.curqg(i,k));
                g.mac.eq(i,k) = -g.mac.mac_con(i,5)*g.mac.curqg(i,k) + (g.mac.psidpp...
                    -g.mac.mac_con(i,8)*g.mac.curdg(i,k));
                g.mac.eterm(i,k) = sqrt(g.mac.ed(i,k)^2+g.mac.eq(i,k)^2);
                g.mac.pelect(i,k) = g.mac.eq(i,k)*g.mac.curq(i,k) + g.mac.ed(i,k)*g.mac.curd(i,k);
                g.mac.qelect(i,k) = g.mac.eq(i,k)*g.mac.curd(i,k) - g.mac.ed(i,k)*g.mac.curq(i,k);
                g.mac.dmac_ang(i,k) = g.sys.basrad*(g.mac.mac_spd(i,k)-1.);
                Te = g.mac.pelect(i,k)*g.mac.mac_pot(i,1) + g.mac.mac_con(i,5)*mcurmag*mcurmag;
                g.mac.dmac_spd(i,k) = (g.mac.pmech(i,k)+g.mac.pm_sig(i,k)-Te...
                    -g.mac.mac_con(i,17)*(g.mac.mac_spd(i,k)-1)...
                    -g.mac.mac_con(i,18)*(g.mac.mac_spd(i,k)-g.sys.sys_freq(k)))... % sys freq always one...
                    /(2*g.mac.mac_con(i,16));
            end
        else
            % vectorized computation
            
            g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14).*g.mac.edprime(g.mac.mac_sub_idx,k) + ...
                g.mac.mac_pot(g.mac.mac_sub_idx,15).*g.mac.psikq(g.mac.mac_sub_idx,k);
            g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9).*g.mac.eqprime(g.mac.mac_sub_idx,k) + ...
                g.mac.mac_pot(g.mac.mac_sub_idx,10).*g.mac.psikd(g.mac.mac_sub_idx,k);
            g.mac.curd(g.mac.mac_sub_idx,k) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.cur_re(g.mac.mac_sub_idx,k) - ...
                cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.cur_im(g.mac.mac_sub_idx,k); % d-axis current
            g.mac.curq(g.mac.mac_sub_idx,k) = cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.cur_re(g.mac.mac_sub_idx,k) + ...
                sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.cur_im(g.mac.mac_sub_idx,k); % q-axis current
            g.mac.curdg(g.mac.mac_sub_idx,k) = g.mac.curd(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1);
            g.mac.curqg(g.mac.mac_sub_idx,k) = g.mac.curq(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1);
            mcurmag = abs(g.mac.curdg(g.mac.mac_sub_idx,k)+1j*g.mac.curqg(g.mac.mac_sub_idx,k));
            E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3).*g.mac.eqprime(g.mac.mac_sub_idx,k).^2 ...
                + g.mac.mac_pot(g.mac.mac_sub_idx,4).*g.mac.eqprime(g.mac.mac_sub_idx,k) + g.mac.mac_pot(g.mac.mac_sub_idx,5);
            nosat_idx=find(g.mac.eqprime(g.mac.mac_sub_idx,1)<.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx)=g.mac.eqprime(g.mac.mac_sub_idx(nosat_idx),k);
            end
            g.mac.fldcur(g.mac.mac_sub_idx,k) = E_Isat + g.mac.mac_pot(g.mac.mac_sub_idx,6)...
                .*(g.mac.eqprime(g.mac.mac_sub_idx,k)-g.mac.psikd(g.mac.mac_sub_idx,k)) + g.mac.mac_pot(g.mac.mac_sub_idx,7)...
                .*g.mac.curdg(g.mac.mac_sub_idx,k);
            g.mac.deqprime(g.mac.mac_sub_idx,k) = (g.mac.vex(g.mac.mac_sub_idx,k)-g.mac.fldcur(g.mac.mac_sub_idx,k))./g.mac.mac_con(g.mac.mac_sub_idx,9);
            g.mac.dpsikd(g.mac.mac_sub_idx,k) = (-g.mac.psikd(g.mac.mac_sub_idx,k)+g.mac.eqprime(g.mac.mac_sub_idx,k)-g.mac.mac_pot(g.mac.mac_sub_idx,8)...
                .*g.mac.curdg(g.mac.mac_sub_idx,k))./g.mac.mac_con(g.mac.mac_sub_idx,10);
            g.mac.dedprime(g.mac.mac_sub_idx,k) = (-g.mac.edprime(g.mac.mac_sub_idx,k) - g.mac.mac_pot(g.mac.mac_sub_idx,11)...
                .*(g.mac.edprime(g.mac.mac_sub_idx,k)-g.mac.psikq(g.mac.mac_sub_idx,k)) - g.mac.mac_pot(g.mac.mac_sub_idx,12)...
                .*g.mac.curqg(g.mac.mac_sub_idx,k))./g.mac.mac_con(g.mac.mac_sub_idx,14);
            g.mac.dpsikq(g.mac.mac_sub_idx,k) = (g.mac.edprime(g.mac.mac_sub_idx,k)-g.mac.psikq(g.mac.mac_sub_idx,k)-g.mac.mac_pot(g.mac.mac_sub_idx,13)...
                .*g.mac.curqg(g.mac.mac_sub_idx,k))./g.mac.mac_con(g.mac.mac_sub_idx,15);
            g.mac.ed(g.mac.mac_sub_idx,k) = -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curdg(g.mac.mac_sub_idx,k) - (g.mac.psiqpp...
                -g.mac.mac_con(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,k));
            g.mac.eq(g.mac.mac_sub_idx,k) = -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,k) + (g.mac.psidpp...
                -g.mac.mac_con(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,k));
            g.mac.eterm(g.mac.mac_sub_idx,k) = sqrt(g.mac.ed(g.mac.mac_sub_idx,k).^2+g.mac.eq(g.mac.mac_sub_idx,k).^2);
            g.mac.pelect(g.mac.mac_sub_idx,k) = g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k) + g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k);
            g.mac.qelect(g.mac.mac_sub_idx,k) = g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k) - g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k);
            g.mac.dmac_ang(g.mac.mac_sub_idx,k) = g.sys.basrad*(g.mac.mac_spd(g.mac.mac_sub_idx,k)-ones(g.mac.n_sub,1));
            Te = g.mac.pelect(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1) + g.mac.mac_con(g.mac.mac_sub_idx,5).*mcurmag.*mcurmag;
            g.mac.dmac_spd(g.mac.mac_sub_idx,k) = (g.mac.pmech(g.mac.mac_sub_idx,k)+ g.mac.pm_sig(g.mac.mac_sub_idx,k)-Te...
                -g.mac.mac_con(g.mac.mac_sub_idx,17).*(g.mac.mac_spd(g.mac.mac_sub_idx,k)-ones(g.mac.n_sub,1))...
                -g.mac.mac_con(g.mac.mac_sub_idx,18).*(g.mac.mac_spd(g.mac.mac_sub_idx,k)-g.sys.sys_freq(k)...
                *ones(g.mac.n_sub,1)))./(2*g.mac.mac_con(g.mac.mac_sub_idx,16));
        end
        %end rate calculation
    end
end
