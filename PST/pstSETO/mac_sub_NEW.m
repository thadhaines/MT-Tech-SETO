function mac_sub(i,k,bus,flag)
%MAC_SUB Voltage-behind-subtransient-reactance generator model.
% MAC_SUB is a voltage-behind-subtransient-reactance generator model, 
% with vectorized computation option.
%
% Syntax: mac_sub(i,k,bus,flag)
%
%   NOTES: State variables are: mac_ang, mac_spd, eqprime, psikd, edprime,
%   psikq in the g.mac. field
%   Algorithm: PSLF model from John Undrill without saturation.
% 
%   Input: 
%   i - generator number
%          - 0 for vectorized computation
%   k - integer time (data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation and state state matrix building
%           3 - Same as 2
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
%   07/xx/08    XX:XX   Dan Trudnowski  Version 2.x? Implemented genrou model of PSLF
%   06/19/20    10:13   Thad Haines     Revised format of globals and internal function documentation

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang  bus_int
global  psi_re psi_im cur_re cur_im

global g

jay = sqrt(-1);
if g.mac.n_sub~=0
  if flag == 0 % initialization
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
      busnum = bus_int(g.mac.mac_con(g.mac.mac_sub_idx,2)); % bus number 
      g.mac.mac_pot(g.mac.mac_sub_idx,1) = basmva*ones(g.mac.n_sub,1)./g.mac.mac_con(g.mac.mac_sub_idx,3); 
                          % scaled MVA base
      g.mac.mac_pot(g.mac.mac_sub_idx,2) = ones(g.mac.n_sub,1); % base kv
      g.mac.mac_pot(g.mac.mac_sub_idx,8)=g.mac.mac_con(g.mac.mac_sub_idx,7)-g.mac.mac_con(g.mac.mac_sub_idx,4);
      g.mac.mac_pot(g.mac.mac_sub_idx,9)=(g.mac.mac_con(g.mac.mac_sub_idx,8)-g.mac.mac_con(g.mac.mac_sub_idx,4))...
                   ./g.mac.mac_pot(g.mac.mac_sub_idx,8);
      g.mac.mac_pot(g.mac.mac_sub_idx,7)=g.mac.mac_con(g.mac.mac_sub_idx,6)-g.mac.mac_con(g.mac.mac_sub_idx,7);
      g.mac.mac_pot(g.mac.mac_sub_idx,10)=(g.mac.mac_con(g.mac.mac_sub_idx,7)-g.mac.mac_con(g.mac.mac_sub_idx,8))...
                    ./g.mac.mac_pot(g.mac.mac_sub_idx,8);
      g.mac.mac_pot(g.mac.mac_sub_idx,6)=g.mac.mac_pot(g.mac.mac_sub_idx,10)./g.mac.mac_pot(g.mac.mac_sub_idx,8);
      g.mac.mac_pot(g.mac.mac_sub_idx,13)=g.mac.mac_con(g.mac.mac_sub_idx,12)-g.mac.mac_con(g.mac.mac_sub_idx,4);
      g.mac.mac_pot(g.mac.mac_sub_idx,14)=(g.mac.mac_con(g.mac.mac_sub_idx,13)-g.mac.mac_con(g.mac.mac_sub_idx,4))...
                     ./g.mac.mac_pot(g.mac.mac_sub_idx,13);
      g.mac.mac_pot(g.mac.mac_sub_idx,12)=g.mac.mac_con(g.mac.mac_sub_idx,11)-g.mac.mac_con(g.mac.mac_sub_idx,12);
      g.mac.mac_pot(g.mac.mac_sub_idx,15)=(g.mac.mac_con(g.mac.mac_sub_idx,12)-g.mac.mac_con(g.mac.mac_sub_idx,13))...
                     ./g.mac.mac_pot(g.mac.mac_sub_idx,13);
      g.mac.mac_pot(g.mac.mac_sub_idx,11)=g.mac.mac_pot(g.mac.mac_sub_idx,15)./g.mac.mac_pot(g.mac.mac_sub_idx,13);
      
      % extract bus information
      g.mac.eterm(g.mac.mac_sub_idx,1) = bus(busnum,2);  % terminal bus voltage
      g.mac.theta(busnum,1) = bus(busnum,3)*pi/180;  % terminal bus angle in radians
      g.mac.pelect(g.mac.mac_sub_idx,1) = bus(busnum,4).*g.mac.mac_con(g.mac.mac_sub_idx,22);  % electrical output power, active
      g.mac.qelect(g.mac.mac_sub_idx,1) = bus(busnum,5).*g.mac.mac_con(g.mac.mac_sub_idx,23);  % electrical output power, reactive
      curr = sqrt(g.mac.pelect(g.mac.mac_sub_idx,1).^2+g.mac.qelect(g.mac.mac_sub_idx,1).^2) ...
            ./g.mac.eterm(g.mac.mac_sub_idx,1).*g.mac.mac_pot(g.mac.mac_sub_idx,1);  % current magnitude on generator base
      phi = atan2(g.mac.qelect(g.mac.mac_sub_idx,1),g.mac.pelect(g.mac.mac_sub_idx,1)); % power factor angle
      v = g.mac.eterm(g.mac.mac_sub_idx,1).*exp(jay*g.mac.theta(busnum,1)); % voltage in real and imaginary parts in system reference frame 
      curr = curr.*exp(jay*(g.mac.theta(busnum,1)-phi));  % complex current in system reference frame 
      ei = v + (g.mac.mac_con(g.mac.mac_sub_idx,5)+jay*g.mac.mac_con(g.mac.mac_sub_idx,11)).*curr; % voltage behind sub-transient reactance in system frame
      g.mac.mac_ang(g.mac.mac_sub_idx,1) = atan2(imag(ei),real(ei)); % machine angle (delta)
      g.mac.mac_spd(g.mac.mac_sub_idx,1) = ones(g.mac.n_sub,1); % machine speed at steady state
      rot = jay*exp(-jay*g.mac.mac_ang(g.mac.mac_sub_idx,1)); % system reference frame rotation to Park's frame
      curr = curr.*rot;% current on generator base in Park's frame
      mcurmag = abs(curr);
      g.mac.pmech(g.mac.mac_sub_idx,1) = g.mac.pelect(g.mac.mac_sub_idx,1).*g.mac.mac_pot(g.mac.mac_sub_idx,1)...
                             + g.mac.mac_con(g.mac.mac_sub_idx,5).*(mcurmag.*mcurmag);% mechanical power = electrical power + losses on generator base
      g.mac.curdg(g.mac.mac_sub_idx,1) = real(curr); 
      g.mac.curqg(g.mac.mac_sub_idx,1) = imag(curr); % d and q axis current on generator base
      g.mac.curd(g.mac.mac_sub_idx,1) = real(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1); 
      g.mac.curq(g.mac.mac_sub_idx,1) = imag(curr)./g.mac.mac_pot(g.mac.mac_sub_idx,1);% d and q axis currents on system base
      v = v.*rot;% voltage in Park's frame
      g.mac.ed(g.mac.mac_sub_idx,1) = real(v); 
      g.mac.eq(g.mac.mac_sub_idx,1) = imag(v);% d and q axis voltages in Park's frame
      eqra = g.mac.eq(g.mac.mac_sub_idx,1)+g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,1);% q axis voltage behind resistance 
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
      
      %No saturation for now
      E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,1);
%       E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3).*eqprime(g.mac.mac_sub_idx,1).^2 ...
%          + g.mac.mac_pot(g.mac.mac_sub_idx,4).*eqprime(g.mac.mac_sub_idx,1) + g.mac.mac_pot(g.mac.mac_sub_idx,5);
%       nosat_idx=find(eqprime(g.mac.mac_sub_idx,1)<.8);
%       if ~isempty(nosat_idx)
%          E_Isat(nosat_idx)=eqprime(g.mac.mac_sub_idx(nosat_idx),1);
%       end
      Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,1) - g.mac.psikd(g.mac.mac_sub_idx,1) - g.mac.mac_pot(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,1);
      g.mac.vex(g.mac.mac_sub_idx,1) = g.mac.eqprime(g.mac.mac_sub_idx,1) + ...
                              g.mac.mac_pot(g.mac.mac_sub_idx,7).*(g.mac.curdg(g.mac.mac_sub_idx,1)+g.mac.mac_pot(g.mac.mac_sub_idx,6).*Eqpe);
      g.mac.fldcur(g.mac.mac_sub_idx,1) = g.mac.vex(g.mac.mac_sub_idx,1);   
      psi_re(g.mac.mac_sub_idx,1) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) + ...
                              cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp; % real part of psi
      psi_im(g.mac.mac_sub_idx,1) = -cos(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*(-g.mac.psiqpp) + ...
                              sin(g.mac.mac_ang(g.mac.mac_sub_idx,1)).*g.mac.psidpp; % imag part of psi
      % psi is in system base and is the voltage behind xpp
  %end initialization
  end
  if flag == 1 % network interface computation 
     % vectorized computation
      
      g.mac.mac_ang(g.mac.mac_sub_idx,k) = g.mac.mac_ang(g.mac.mac_sub_idx,k)-mach_ref(k)*ones(g.mac.n_sub,1); % wrt machine reference
      g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9).*g.mac.eqprime(g.mac.mac_sub_idx,k) + ...
               g.mac.mac_pot(g.mac.mac_sub_idx,10).*g.mac.psikd(g.mac.mac_sub_idx,k);
      g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14).*g.mac.edprime(g.mac.mac_sub_idx,k) + ...
               g.mac.mac_pot(g.mac.mac_sub_idx,15).*g.mac.psikq(g.mac.mac_sub_idx,k);
      psi_re(g.mac.mac_sub_idx,k) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) + ...
                              cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp; % real part of psi
      psi_im(g.mac.mac_sub_idx,k) = -cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*(-g.mac.psiqpp) + ...
                              sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*g.mac.psidpp; % imag part of psi   
  % end of interface
  end

  if flag == 2 || flag == 3 % generator dynamics calculation
    % vectorized computation
      
      g.mac.psiqpp = g.mac.mac_pot(g.mac.mac_sub_idx,14).*g.mac.edprime(g.mac.mac_sub_idx,k) + ...
               g.mac.mac_pot(g.mac.mac_sub_idx,15).*g.mac.psikq(g.mac.mac_sub_idx,k); 
      g.mac.psidpp = g.mac.mac_pot(g.mac.mac_sub_idx,9).*g.mac.eqprime(g.mac.mac_sub_idx,k) + ...
               g.mac.mac_pot(g.mac.mac_sub_idx,10).*g.mac.psikd(g.mac.mac_sub_idx,k);
      g.mac.curd(g.mac.mac_sub_idx,k) = sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*cur_re(g.mac.mac_sub_idx,k) - ...
                            cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*cur_im(g.mac.mac_sub_idx,k); % d-axis current
      g.mac.curq(g.mac.mac_sub_idx,k) = cos(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*cur_re(g.mac.mac_sub_idx,k) + ...
                            sin(g.mac.mac_ang(g.mac.mac_sub_idx,k)).*cur_im(g.mac.mac_sub_idx,k); % q-axis current
      g.mac.curdg(g.mac.mac_sub_idx,k) = g.mac.curd(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1);
      g.mac.curqg(g.mac.mac_sub_idx,k) = g.mac.curq(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1);
      mcurmag = abs(g.mac.curdg(g.mac.mac_sub_idx,k)+jay*g.mac.curqg(g.mac.mac_sub_idx,k));
      
      %No saturation for now
      E_Isat = g.mac.eqprime(g.mac.mac_sub_idx,k);      
%       E_Isat = g.mac.mac_pot(g.mac.mac_sub_idx,3).*eqprime(g.mac.mac_sub_idx,k).^2 ...
%               + g.mac.mac_pot(g.mac.mac_sub_idx,4).*eqprime(g.mac.mac_sub_idx,k) + g.mac.mac_pot(g.mac.mac_sub_idx,5); 
%       nosat_idx=find(eqprime(g.mac.mac_sub_idx,1)<.8);
%       if ~isempty(nosat_idx)
%          E_Isat(nosat_idx)=eqprime(g.mac.mac_sub_idx(nosat_idx),k);
%       end
      Eqpe = g.mac.eqprime(g.mac.mac_sub_idx,k) - g.mac.psikd(g.mac.mac_sub_idx,k) - g.mac.mac_pot(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,k);
      g.mac.fldcur(g.mac.mac_sub_idx,k) = g.mac.eqprime(g.mac.mac_sub_idx,k) + ...
                              g.mac.mac_pot(g.mac.mac_sub_idx,7).*(g.mac.curdg(g.mac.mac_sub_idx,k)+g.mac.mac_pot(g.mac.mac_sub_idx,6).*Eqpe);      
      g.mac.deqprime(g.mac.mac_sub_idx,k) = (g.mac.vex(g.mac.mac_sub_idx,k)-g.mac.fldcur(g.mac.mac_sub_idx,k))./g.mac.mac_con(g.mac.mac_sub_idx,9);
      g.mac.dpsikd(g.mac.mac_sub_idx,k) = Eqpe./g.mac.mac_con(g.mac.mac_sub_idx,10);
      
      Edpe = g.mac.edprime(g.mac.mac_sub_idx,k) - g.mac.psikq(g.mac.mac_sub_idx,k) - g.mac.mac_pot(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,k);
      Hold = -g.mac.edprime(g.mac.mac_sub_idx,k) + ...
                              g.mac.mac_pot(g.mac.mac_sub_idx,12).*(-g.mac.curqg(g.mac.mac_sub_idx,k)-g.mac.mac_pot(g.mac.mac_sub_idx,11).*Edpe); 
      g.mac.dedprime(g.mac.mac_sub_idx,k) = (Hold)./g.mac.mac_con(g.mac.mac_sub_idx,14);
      g.mac.dpsikq(g.mac.mac_sub_idx,k) = Edpe./g.mac.mac_con(g.mac.mac_sub_idx,15);

      g.mac.ed(g.mac.mac_sub_idx,k) = -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curdg(g.mac.mac_sub_idx,k) - ...
                (g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psiqpp... %Multiply g.mac.psiqpp by gen spd
                -g.mac.mac_con(g.mac.mac_sub_idx,13).*g.mac.curqg(g.mac.mac_sub_idx,k));
      g.mac.eq(g.mac.mac_sub_idx,k) = -g.mac.mac_con(g.mac.mac_sub_idx,5).*g.mac.curqg(g.mac.mac_sub_idx,k) + ...
                (g.mac.mac_spd(g.mac.mac_sub_idx,k).*g.mac.psidpp... %Multiply g.mac.psidpp by gen spd
                -g.mac.mac_con(g.mac.mac_sub_idx,8).*g.mac.curdg(g.mac.mac_sub_idx,k));
      g.mac.eterm(g.mac.mac_sub_idx,k) = sqrt(g.mac.ed(g.mac.mac_sub_idx,k).^2+g.mac.eq(g.mac.mac_sub_idx,k).^2);
      g.mac.pelect(g.mac.mac_sub_idx,k) = g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k) + g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k);
      g.mac.qelect(g.mac.mac_sub_idx,k) = g.mac.eq(g.mac.mac_sub_idx,k).*g.mac.curd(g.mac.mac_sub_idx,k) - g.mac.ed(g.mac.mac_sub_idx,k).*g.mac.curq(g.mac.mac_sub_idx,k);
      g.mac.dmac_ang(g.mac.mac_sub_idx,k) = basrad*(g.mac.mac_spd(g.mac.mac_sub_idx,k)-ones(g.mac.n_sub,1));
      Te = g.mac.pelect(g.mac.mac_sub_idx,k).*g.mac.mac_pot(g.mac.mac_sub_idx,1) + g.mac.mac_con(g.mac.mac_sub_idx,5).*mcurmag.*mcurmag;
      g.mac.dmac_spd(g.mac.mac_sub_idx,k) = (g.mac.pmech(g.mac.mac_sub_idx,k)./g.mac.mac_spd(g.mac.mac_sub_idx,k) ...
            + g.mac.pm_sig(g.mac.mac_sub_idx,k) - Te...
            - g.mac.mac_con(g.mac.mac_sub_idx,17).*(g.mac.mac_spd(g.mac.mac_sub_idx,k)-ones(g.mac.n_sub,1)))...
            ./(2*g.mac.mac_con(g.mac.mac_sub_idx,16)); 

  %end rate calculation
  end
end
