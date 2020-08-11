function mac_em(i,k,bus,flag)
%MAC_EM Classical generator electromechanical model
% MAC_EM is a classical generator electromechanical model, with
% vectorized computation option.
%
% Syntax: mac_em(i,k,bus,flag)
%
%   NOTES: State variables are: mac_ang, mac_spd in the g.mac. field
% 
%   Input: 
%   i - generator number
%          - 0 for vectorized computation
%   k - integer time (data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation and state state matrix building
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   01/xx/91    XX:XX   Joe H. Chow     Version 1
%   (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%   06/xx/96    xx:xx   Graham Rogers   Version 2 added facility to allow different machine models in vector run
%   (c) Copyright 1991-1997 Joe H. Chow/ Cherry Tree Scientific Software - All Rights Reserved
%   06/19/20    09:52   Thad Haines     Revised format of globals and internal function documentation
%   07/07/20    14:32   Thad Haines     Completion of global g alteration
%   07/29/20    15:20   Thad Haines     jay -> 1j
%   08/11/20    11:48   Thad Haines     added ivm to global

global g

if g.mac.n_em ~=0
 if flag == 0 % initialization
   if i~=0
      % non-vector calculation
      % check for em model
      em = find(g.mac.mac_em_idx==i);
      if ~isempty(em)
        busnum = g.bus.bus_int(g.mac.mac_con(i,2)); % bus number 
        g.mac.mac_pot(i,1) = g.sys.basmva/g.mac.mac_con(i,3); % scaled MVA base
        g.mac.mac_pot(i,2) = 1.0; % base kv
        % extract bus information
        g.mac.eterm(i,1) = bus(busnum,2);  % terminal bus voltage
        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
        g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22);  
                        % electrical output power, active
        g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);  
                        % electrical output power, reactive
        curr = sqrt(g.mac.pelect(i,1)^2+g.mac.qelect(i,1)^2) ...
              /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);  % current magnitude
                                         % on generator base
        phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1)); 
                                        % power factor angle
        v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
        curr = curr*exp(1j*(g.bus.theta(busnum,1)-phi)); % complex current  
                                                    % in system reference frame 
        eprime = v + 1j*g.mac.mac_con(i,7)*curr; 
        % voltage behind transient reactance
        ei = eprime;
        g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
        g.mac.mac_spd(i,1) = 1; % machine speed at steady state
        rot = 1j*exp(-1j*g.mac.mac_ang(i,1)); 
                          % system reference frame rotation
        g.mac.psi_re(i,1) = real(eprime);
        g.mac.psi_im(i,1) = imag(eprime);
        eprime = eprime*rot;
        g.mac.edprime(i,1) = real(eprime); 
        g.mac.eqprime(i,1) = imag(eprime); 
        curr = curr*rot;%current on Park's frame
        g.mac.curdg(i,1) = real(curr); 
        g.mac.curqg(i,1) = imag(curr);
        g.mac.curd(i,1) = real(curr)/g.mac.mac_pot(i,1); 
        g.mac.curq(i,1) = imag(curr)/g.mac.mac_pot(i,1);
        % convert to system base
        v = v*rot;
        g.mac.ed(i,1) = real(v); 
        g.mac.eq(i,1) = imag(v);% in Park's frame
        g.mac.vex(i,1) = g.mac.eqprime(i,1);
        g.mac.pmech(i,1) = g.mac.pelect(i,1)*g.mac.mac_pot(i,1); % set input
         % mechanical power equal to electrical output power
         % since losses are zero for em model. On generator base
      end
   else
      % vectorized computation
      busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_em_idx,2)); % bus number 
      g.mac.mac_pot(g.mac.mac_em_idx,1) = g.sys.basmva*ones(g.mac.n_em,1)./g.mac.mac_con(g.mac.mac_em_idx,3); 
                          % scaled MVA base
      g.mac.mac_pot(g.mac.mac_em_idx,2) = 1.0*ones(g.mac.n_em,1); % base kv
      % extract bus information
      g.mac.eterm(g.mac.mac_em_idx,1) = bus(busnum,2);  % terminal bus voltage
      g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
      g.mac.pelect(g.mac.mac_em_idx,1) = bus(busnum,4).*g.mac.mac_con(g.mac.mac_em_idx,22);  
                        % electrical output power, active
      g.mac.qelect(g.mac.mac_em_idx,1) = bus(busnum,5).*g.mac.mac_con(g.mac.mac_em_idx,23);  
                        % electrical output power, reactive
      curr = sqrt(g.mac.pelect(g.mac.mac_em_idx,1).^2+g.mac.qelect(g.mac.mac_em_idx,1).^2)...
            ./g.mac.eterm(g.mac.mac_em_idx,1).*g.mac.mac_pot(g.mac.mac_em_idx,1);  % current magnitude
                                                           % on generator base
      phi = atan2(g.mac.qelect(g.mac.mac_em_idx,1),g.mac.pelect(g.mac.mac_em_idx,1)); 
                                        % power factor angle
      v = g.mac.eterm(g.mac.mac_em_idx,1).*exp(1j*g.bus.theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
      curr = curr.*exp(1j*(g.bus.theta(busnum,1)-phi)); % current in real and 
                 % imaginary parts on system reference frame 
      eprime = v + 1j*g.mac.mac_con(g.mac.mac_em_idx,7).*curr; 
      ei = eprime;
      g.mac.mac_ang(g.mac.mac_em_idx,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
      g.mac.mac_spd(g.mac.mac_em_idx,1) = ones(g.mac.n_em,1); 
                            % machine speed at steady state
      rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_em_idx,1)); 
                          % system reference frame rotation
      g.mac.psi_re(g.mac.mac_em_idx,1) = real(eprime);
      g.mac.psi_im(g.mac.mac_em_idx,1) = imag(eprime);
      eprime = eprime.*rot;% in Park's frame
      g.mac.edprime(g.mac.mac_em_idx,1) = real(eprime); 
      g.mac.eqprime(g.mac.mac_em_idx,1) = imag(eprime);
      curr = curr.*rot;%in Park's frame
      g.mac.curdg(g.mac.mac_em_idx,1) = real(curr); 
      g.mac.curqg(g.mac.mac_em_idx,1) = imag(curr);
      g.mac.curd(g.mac.mac_em_idx,1) = real(curr)./g.mac.mac_pot(g.mac.mac_em_idx,1); 
      g.mac.curq(g.mac.mac_em_idx,1) = imag(curr)./g.mac.mac_pot(g.mac.mac_em_idx,1);
      v = v.*rot;%in Park's frame
      g.mac.ed(g.mac.mac_em_idx,1) = real(v); 
      g.mac.eq(g.mac.mac_em_idx,1) = imag(v);
      g.mac.vex(g.mac.mac_em_idx,1) = g.mac.eqprime(g.mac.mac_em_idx,1);
      g.mac.pmech(g.mac.mac_em_idx,1) = g.mac.pelect(g.mac.mac_em_idx,1)...
          .*g.mac.mac_pot(g.mac.mac_em_idx,1); % set input
         % mechanical power equal to electrical output power
         % since losses are zero in em model. On generator base
   end
   %end initialization
 end

 if flag == 1 % network interface computation 
  if i ~= 0
      % check for em machine
      em = find(g.mac.mac_em_idx==i);
      if ~isempty(em)
        g.mac.mac_ang(i,k) = g.mac.mac_ang(i,k) - g.sys.mach_ref(k);   
                     % wrt machine reference
        g.mac.psi_re(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) + ...
                      cos(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k); % real part of psi
        g.mac.psi_im(i,k) = -cos(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) + ...
                      sin(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k); % imag part of psi
      end
  else
      % vectorized computation
      g.mac.mac_ang(g.mac.mac_em_idx,k) = g.mac.mac_ang(g.mac.mac_em_idx,k)...
          -g.sys.mach_ref(k)*ones(g.mac.n_em,1);  % wrt machine reference
      g.mac.psi_re(g.mac.mac_em_idx,k) = sin(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.edprime(g.mac.mac_em_idx,k) + ...
         cos(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.eqprime(g.mac.mac_em_idx,k); % real part of psi
      g.mac.psi_im(g.mac.mac_em_idx,k) = -cos(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.edprime(g.mac.mac_em_idx,k) + ...
         sin(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.eqprime(g.mac.mac_em_idx,k); % imag part of psi
  end
  % end interface
 end
 if flag == 2 % generator dynamics calculation
  if i ~= 0
    % check for em machine
    em = find(g.mac.mac_em_idx==i);
    if ~isempty(em)
      g.mac.curd(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) - ...
            cos(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k); % d-axis current
      g.mac.curq(i,k) = cos(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) + ...
            sin(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k); % q-axis current
      g.mac.curdg(i,k) = g.mac.curd(i,k)*g.mac.mac_pot(i,1);
      g.mac.curqg(i,k) = g.mac.curq(i,k)*g.mac.mac_pot(i,1);
      g.mac.dedprime(i,k) = 0;
      g.mac.deqprime(i,k) = 0;
      g.mac.ed(i,k) = g.mac.edprime(i,k) + g.mac.mac_con(i,7)*g.mac.curqg(i,k);
      g.mac.eq(i,k) = g.mac.eqprime(i,k) - g.mac.mac_con(i,7)*g.mac.curdg(i,k);
      g.mac.eterm(i,k) = sqrt(g.mac.ed(i,k)^2+g.mac.eq(i,k)^2);
      g.mac.pelect(i,k) = g.mac.eq(i,k)*g.mac.curq(i,k) + g.mac.ed(i,k)*g.mac.curd(i,k);
      g.mac.qelect(i,k) = g.mac.eq(i,k)*g.mac.curd(i,k) - g.mac.ed(i,k)*g.mac.curq(i,k);
      g.mac.dmac_ang(i,k) = g.sys.basrad*(g.mac.mac_spd(i,k)-1.);
      g.mac.dmac_spd(i,k) = (g.mac.pmech(i,k)-g.mac.pelect(i,k)*g.mac.mac_pot(i,1)... 
        -g.mac.mac_con(i,17)*(g.mac.mac_spd(i,k)-1))/(2.*g.mac.mac_con(i,16));
    end
  else
      % vectorized computation
      g.mac.curd(g.mac.mac_em_idx,k) = sin(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.cur_re(g.mac.mac_em_idx,k) - ...
            cos(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.cur_im(g.mac.mac_em_idx,k); % d-axis current
      g.mac.curq(g.mac.mac_em_idx,k) = cos(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.cur_re(g.mac.mac_em_idx,k) + ...
            sin(g.mac.mac_ang(g.mac.mac_em_idx,k)).*g.mac.cur_im(g.mac.mac_em_idx,k); % q-axis current
      g.mac.curdg(g.mac.mac_em_idx,k) = g.mac.curd(g.mac.mac_em_idx,k).*g.mac.mac_pot(g.mac.mac_em_idx,1);
      g.mac.curqg(g.mac.mac_em_idx,k) = g.mac.curq(g.mac.mac_em_idx,k).*g.mac.mac_pot(g.mac.mac_em_idx,1);
      g.mac.dedprime(g.mac.mac_em_idx,k) = zeros(g.mac.n_em,1);
      g.mac.deqprime(g.mac.mac_em_idx,k) = zeros(g.mac.n_em,1);
      g.mac.ed(g.mac.mac_em_idx,k) = g.mac.edprime(g.mac.mac_em_idx,k)...
                        + g.mac.mac_con(g.mac.mac_em_idx,7).*g.mac.curqg(g.mac.mac_em_idx,k);
      g.mac.eq(g.mac.mac_em_idx,k) = g.mac.eqprime(g.mac.mac_em_idx,k)...
                        - g.mac.mac_con(g.mac.mac_em_idx,7).*g.mac.curdg(g.mac.mac_em_idx,k);
      g.mac.eterm(g.mac.mac_em_idx,k) = sqrt(g.mac.ed(g.mac.mac_em_idx,k).^2+g.mac.eq(g.mac.mac_em_idx,k).^2);
      g.mac.pelect(g.mac.mac_em_idx,k) = g.mac.eq(g.mac.mac_em_idx,k).*g.mac.curq(g.mac.mac_em_idx,k)...
                           + g.mac.ed(g.mac.mac_em_idx,k).*g.mac.curd(g.mac.mac_em_idx,k);
      g.mac.qelect(g.mac.mac_em_idx,k) = g.mac.eq(g.mac.mac_em_idx,k).*g.mac.curd(g.mac.mac_em_idx,k)...
                           - g.mac.ed(g.mac.mac_em_idx,k).*g.mac.curq(g.mac.mac_em_idx,k);
      g.mac.dmac_ang(g.mac.mac_em_idx,k) = g.sys.basrad*(g.mac.mac_spd(g.mac.mac_em_idx,k)-ones(g.mac.n_em,1));
      g.mac.dmac_spd(g.mac.mac_em_idx,k) =(g.mac.pmech(g.mac.mac_em_idx,k)+ g.mac.pm_sig(g.mac.mac_em_idx,k) ...
                      -g.mac.pelect(g.mac.mac_em_idx,k).*g.mac.mac_pot(g.mac.mac_em_idx,1)...
                      -g.mac.mac_con(g.mac.mac_em_idx,17).*(g.mac.mac_spd(g.mac.mac_em_idx,k)...
                      -ones(g.mac.n_em,1)))./(2*g.mac.mac_con(g.mac.mac_em_idx,16));
  end
  % end rate calculation
 end
end
