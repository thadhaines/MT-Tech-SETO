function [f] = mac_sub(i,k,bus,flag)
% Syntax: [f] = mac_sub(i,k,bus,flag)
%
% Purpose: voltage-behind-subtransient-reactance generator
%            model, with vectorized computation option
%          state variables are: mac_ang, mac_spd, eqprime,
%                               psikd, edprime, psikq
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation 
%               3 - state matrix building 
%
% Output: f - dummy variable 
%
% Files:
%
% See Also: pst_var, mac_em, mac_tra

% Algorithm: PSLF model from John Undrill without saturation.
%
% Calls:
%
% Called By:

% (c) Copyright 1991-1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved
% modified by Dan Trudnowski, 2008

% History (in reverse chronological order)
%
% Version 2.x?
% IMplemented genrou model of PSLF
% Date: July 2008
% Author: Dan Trudnowski
%
% Version 2.1
% Date:   September 1997
% Saturation modified for low Eqprime
% Version:  2.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  change to allow multple generator model types
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
% Modification: Correction to saturation GJR May 6 1995

% define global variables
% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang  bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime psikd psikq
global  psi_re psi_im cur_re cur_im
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq
global  n_mac n_sub mac_sub_idx
global  pm_sig  n_pm
f = 0;

jay = sqrt(-1);
if n_sub~=0
  if flag == 0; % initialization
      % vectorized computation
      % check parameters
      uets_idx = find(mac_con(mac_sub_idx,8)~=mac_con(mac_sub_idx,13));
      if ~isempty(uets_idx)
         mac_con(mac_sub_idx(uets_idx),13)=mac_con(mac_sub_idx(uets_idx),8);
         disp('xqpp made equal to xdpp at generators  '); disp((mac_sub_idx(uets_idx))')
      end
      notp_idx = find(mac_con(mac_sub_idx,14)==0);
      if ~isempty(notp_idx)
         mac_con(mac_sub_idx(notp_idx),14) = 999.0*ones(length(notp_idx),1);
      end
      notpp_idx = find(mac_con(mac_sub_idx,15)==0);
      if ~isempty(notpp_idx)
         mac_con(mac_sub_idx(notpp_idx),15) = 999.0*ones(length(notpp_idx),1);
         % set x'q = x"q
         mac_con(mac_sub_idx(notpp_idx),12) =...
                                 mac_con(mac_sub_idx(notpp_idx),13);
      end
      busnum = bus_int(mac_con(mac_sub_idx,2)); % bus number 
      mac_pot(mac_sub_idx,1) = basmva*ones(n_sub,1)./mac_con(mac_sub_idx,3); 
                          % scaled MVA base
      mac_pot(mac_sub_idx,2) = ones(n_sub,1); % base kv
      mac_pot(mac_sub_idx,8)=mac_con(mac_sub_idx,7)-mac_con(mac_sub_idx,4);
      mac_pot(mac_sub_idx,9)=(mac_con(mac_sub_idx,8)-mac_con(mac_sub_idx,4))...
                   ./mac_pot(mac_sub_idx,8);
      mac_pot(mac_sub_idx,7)=mac_con(mac_sub_idx,6)-mac_con(mac_sub_idx,7);
      mac_pot(mac_sub_idx,10)=(mac_con(mac_sub_idx,7)-mac_con(mac_sub_idx,8))...
                    ./mac_pot(mac_sub_idx,8);
      mac_pot(mac_sub_idx,6)=mac_pot(mac_sub_idx,10)./mac_pot(mac_sub_idx,8);
      mac_pot(mac_sub_idx,13)=mac_con(mac_sub_idx,12)-mac_con(mac_sub_idx,4);
      mac_pot(mac_sub_idx,14)=(mac_con(mac_sub_idx,13)-mac_con(mac_sub_idx,4))...
                     ./mac_pot(mac_sub_idx,13);
      mac_pot(mac_sub_idx,12)=mac_con(mac_sub_idx,11)-mac_con(mac_sub_idx,12);
      mac_pot(mac_sub_idx,15)=(mac_con(mac_sub_idx,12)-mac_con(mac_sub_idx,13))...
                     ./mac_pot(mac_sub_idx,13);
      mac_pot(mac_sub_idx,11)=mac_pot(mac_sub_idx,15)./mac_pot(mac_sub_idx,13);
      
      % extract bus information
      eterm(mac_sub_idx,1) = bus(busnum,2);  % terminal bus voltage
      theta(busnum,1) = bus(busnum,3)*pi/180;  % terminal bus angle in radians
      pelect(mac_sub_idx,1) = bus(busnum,4).*mac_con(mac_sub_idx,22);  % electrical output power, active
      qelect(mac_sub_idx,1) = bus(busnum,5).*mac_con(mac_sub_idx,23);  % electrical output power, reactive
      curr = sqrt(pelect(mac_sub_idx,1).^2+qelect(mac_sub_idx,1).^2) ...
            ./eterm(mac_sub_idx,1).*mac_pot(mac_sub_idx,1);  % current magnitude on generator base
      phi = atan2(qelect(mac_sub_idx,1),pelect(mac_sub_idx,1)); % power factor angle
      v = eterm(mac_sub_idx,1).*exp(jay*theta(busnum,1)); % voltage in real and imaginary parts in system reference frame 
      curr = curr.*exp(jay*(theta(busnum,1)-phi));  % complex current in system reference frame 
      ei = v + (mac_con(mac_sub_idx,5)+jay*mac_con(mac_sub_idx,11)).*curr; % voltage behind sub-transient reactance in system frame
      mac_ang(mac_sub_idx,1) = atan2(imag(ei),real(ei)); % machine angle (delta)
      mac_spd(mac_sub_idx,1) = ones(n_sub,1); % machine speed at steady state
      rot = jay*exp(-jay*mac_ang(mac_sub_idx,1)); % system reference frame rotation to Park's frame
      curr = curr.*rot;% current on generator base in Park's frame
      mcurmag = abs(curr);
      pmech(mac_sub_idx,1) = pelect(mac_sub_idx,1).*mac_pot(mac_sub_idx,1)...
                             + mac_con(mac_sub_idx,5).*(mcurmag.*mcurmag);% mechanical power = electrical power + losses on generator base
      curdg(mac_sub_idx,1) = real(curr); 
      curqg(mac_sub_idx,1) = imag(curr); % d and q axis current on generator base
      curd(mac_sub_idx,1) = real(curr)./mac_pot(mac_sub_idx,1); 
      curq(mac_sub_idx,1) = imag(curr)./mac_pot(mac_sub_idx,1);% d and q axis currents on system base
      v = v.*rot;% voltage in Park's frame
      ed(mac_sub_idx,1) = real(v); 
      eq(mac_sub_idx,1) = imag(v);% d and q axis voltages in Park's frame
      eqra = eq(mac_sub_idx,1)+mac_con(mac_sub_idx,5).*curqg(mac_sub_idx,1);% q axis voltage behind resistance 
      psidpp = eqra + mac_con(mac_sub_idx,8).*curdg(mac_sub_idx,1);
      psikd(mac_sub_idx,1) = eqra + mac_con(mac_sub_idx,4).*curdg(mac_sub_idx,1);
      eqprime(mac_sub_idx,1) = eqra + mac_con(mac_sub_idx,7).*curdg(mac_sub_idx,1);
      edra = -ed(mac_sub_idx,1)-mac_con(mac_sub_idx,5).*curdg(mac_sub_idx,1);
      psiqpp = edra + mac_con(mac_sub_idx,13).*curqg(mac_sub_idx,1);
      psikq(mac_sub_idx,1) = edra + mac_con(mac_sub_idx,4).*curqg(mac_sub_idx,1);
      
      edprime(mac_sub_idx,1) = edra + mac_con(mac_sub_idx,12).*curqg(mac_sub_idx,1);
      % this is the negative of Edprime in block diagram
      % compute saturation
      inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
      b = [0.8*ones(n_sub,1) ones(n_sub,1)+mac_con(mac_sub_idx,20)...
           1.2*(ones(n_sub,1)+mac_con(mac_sub_idx,21))];
      mac_pot(mac_sub_idx,3) = b*inv_sat(1,:)';
      mac_pot(mac_sub_idx,4) = b*inv_sat(2,:)';
      mac_pot(mac_sub_idx,5) = b*inv_sat(3,:)';
      
      %No saturation for now
      E_Isat = eqprime(mac_sub_idx,1);
%       E_Isat = mac_pot(mac_sub_idx,3).*eqprime(mac_sub_idx,1).^2 ...
%          + mac_pot(mac_sub_idx,4).*eqprime(mac_sub_idx,1) + mac_pot(mac_sub_idx,5);
%       nosat_idx=find(eqprime(mac_sub_idx,1)<.8);
%       if ~isempty(nosat_idx)
%          E_Isat(nosat_idx)=eqprime(mac_sub_idx(nosat_idx),1);
%       end
      Eqpe = eqprime(mac_sub_idx,1) - psikd(mac_sub_idx,1) - mac_pot(mac_sub_idx,8).*curdg(mac_sub_idx,1);
      vex(mac_sub_idx,1) = eqprime(mac_sub_idx,1) + ...
                              mac_pot(mac_sub_idx,7).*(curdg(mac_sub_idx,1)+mac_pot(mac_sub_idx,6).*Eqpe);
      fldcur(mac_sub_idx,1) = vex(mac_sub_idx,1);   
      psi_re(mac_sub_idx,1) = sin(mac_ang(mac_sub_idx,1)).*(-psiqpp) + ...
                              cos(mac_ang(mac_sub_idx,1)).*psidpp; % real part of psi
      psi_im(mac_sub_idx,1) = -cos(mac_ang(mac_sub_idx,1)).*(-psiqpp) + ...
                              sin(mac_ang(mac_sub_idx,1)).*psidpp; % imag part of psi
      % psi is in system base and is the voltage behind xpp
  %end initialization
  end
  if flag == 1 % network interface computation 
     % vectorized computation
      
      mac_ang(mac_sub_idx,k) = mac_ang(mac_sub_idx,k)-mach_ref(k)*ones(n_sub,1); % wrt machine reference
      psidpp = mac_pot(mac_sub_idx,9).*eqprime(mac_sub_idx,k) + ...
               mac_pot(mac_sub_idx,10).*psikd(mac_sub_idx,k);
      psiqpp = mac_pot(mac_sub_idx,14).*edprime(mac_sub_idx,k) + ...
               mac_pot(mac_sub_idx,15).*psikq(mac_sub_idx,k);
      psi_re(mac_sub_idx,k) = mac_spd(mac_sub_idx,k).*(...
                              sin(mac_ang(mac_sub_idx,k)).*(-psiqpp) + ...
                              cos(mac_ang(mac_sub_idx,k)).*psidpp); % real part of psi
      psi_im(mac_sub_idx,k) = mac_spd(mac_sub_idx,k).*(...
                              -cos(mac_ang(mac_sub_idx,k)).*(-psiqpp) + ...
                              sin(mac_ang(mac_sub_idx,k)).*psidpp); % imag part of psi   
  % end of interface
  end

  if flag == 2 | flag == 3 % generator dynamics calculation
    % vectorized computation
      
      psiqpp = mac_pot(mac_sub_idx,14).*edprime(mac_sub_idx,k) + ...
               mac_pot(mac_sub_idx,15).*psikq(mac_sub_idx,k); 
      psidpp = mac_pot(mac_sub_idx,9).*eqprime(mac_sub_idx,k) + ...
               mac_pot(mac_sub_idx,10).*psikd(mac_sub_idx,k);
      curd(mac_sub_idx,k) = sin(mac_ang(mac_sub_idx,k)).*cur_re(mac_sub_idx,k) - ...
                            cos(mac_ang(mac_sub_idx,k)).*cur_im(mac_sub_idx,k); % d-axis current
      curq(mac_sub_idx,k) = cos(mac_ang(mac_sub_idx,k)).*cur_re(mac_sub_idx,k) + ...
                            sin(mac_ang(mac_sub_idx,k)).*cur_im(mac_sub_idx,k); % q-axis current
      curdg(mac_sub_idx,k) = curd(mac_sub_idx,k).*mac_pot(mac_sub_idx,1);
      curqg(mac_sub_idx,k) = curq(mac_sub_idx,k).*mac_pot(mac_sub_idx,1);
      mcurmag = abs(curdg(mac_sub_idx,k)+jay*curqg(mac_sub_idx,k));
      
      %No saturation for now
      E_Isat = eqprime(mac_sub_idx,k);      
%       E_Isat = mac_pot(mac_sub_idx,3).*eqprime(mac_sub_idx,k).^2 ...
%               + mac_pot(mac_sub_idx,4).*eqprime(mac_sub_idx,k) + mac_pot(mac_sub_idx,5); 
%       nosat_idx=find(eqprime(mac_sub_idx,1)<.8);
%       if ~isempty(nosat_idx)
%          E_Isat(nosat_idx)=eqprime(mac_sub_idx(nosat_idx),k);
%       end
      Eqpe = eqprime(mac_sub_idx,k) - psikd(mac_sub_idx,k) - mac_pot(mac_sub_idx,8).*curdg(mac_sub_idx,k);
      fldcur(mac_sub_idx,k) = eqprime(mac_sub_idx,k) + ...
                              mac_pot(mac_sub_idx,7).*(curdg(mac_sub_idx,k)+mac_pot(mac_sub_idx,6).*Eqpe);      
      deqprime(mac_sub_idx,k) = (vex(mac_sub_idx,k)-fldcur(mac_sub_idx,k))./mac_con(mac_sub_idx,9);
      dpsikd(mac_sub_idx,k) = Eqpe./mac_con(mac_sub_idx,10);
      
      Edpe = edprime(mac_sub_idx,k) - psikq(mac_sub_idx,k) - mac_pot(mac_sub_idx,13).*curqg(mac_sub_idx,k);
      Hold = -edprime(mac_sub_idx,k) + ...
                              mac_pot(mac_sub_idx,12).*(-curqg(mac_sub_idx,k)-mac_pot(mac_sub_idx,11).*Edpe); 
      dedprime(mac_sub_idx,k) = (Hold)./mac_con(mac_sub_idx,14);
      dpsikq(mac_sub_idx,k) = Edpe./mac_con(mac_sub_idx,15);

      ed(mac_sub_idx,k) = -mac_con(mac_sub_idx,5).*curdg(mac_sub_idx,k) - ...
                (mac_spd(mac_sub_idx,k).*psiqpp... %Multiply psiqpp by gen spd
                -mac_con(mac_sub_idx,13).*curqg(mac_sub_idx,k));
      eq(mac_sub_idx,k) = -mac_con(mac_sub_idx,5).*curqg(mac_sub_idx,k) + ...
                (mac_spd(mac_sub_idx,k).*psidpp... %Multiply psidpp by gen spd
                -mac_con(mac_sub_idx,8).*curdg(mac_sub_idx,k));
      eterm(mac_sub_idx,k) = sqrt(ed(mac_sub_idx,k).^2+eq(mac_sub_idx,k).^2);
      pelect(mac_sub_idx,k) = eq(mac_sub_idx,k).*curq(mac_sub_idx,k) + ed(mac_sub_idx,k).*curd(mac_sub_idx,k); %System base
      qelect(mac_sub_idx,k) = eq(mac_sub_idx,k).*curd(mac_sub_idx,k) - ed(mac_sub_idx,k).*curq(mac_sub_idx,k);
      dmac_ang(mac_sub_idx,k) = basrad*(mac_spd(mac_sub_idx,k)-ones(n_sub,1));
%        Te = pelect(mac_sub_idx,k).*mac_pot(mac_sub_idx,1) + mac_con(mac_sub_idx,5).*mcurmag.*mcurmag;%gen base
      Te = psidpp.*curqg(mac_sub_idx,k) - psiqpp.*curdg(mac_sub_idx,k);%gen base   
      dmac_spd(mac_sub_idx,k) = (pmech(mac_sub_idx,k)./mac_spd(mac_sub_idx,k) ...
            + pm_sig(mac_sub_idx,k) - Te...
            - mac_con(mac_sub_idx,17).*(mac_spd(mac_sub_idx,k)-ones(n_sub,1)))...
            ./(2*mac_con(mac_sub_idx,16)); 
 
  %end rate calculation
  end
end
