function [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc)
%LFDCS Solves load flow with one or more dc lines
% LFDCS     Solves load flow with one or more dc lines by iterating between
%           AC and DC solutions until converged solution is reached.
%
% Syntax: [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc)  
%
%   NOTES:  dci_dc and dcr_dc are the same as g.dc.dci_dc g.dcr_dc, but appear unused?
%           Called by s_simu
%           Calls loadflow, dc_lf
% 
%   Input: 
%   bus - ac bus specification matrix
% 	line - ac line specification matrix
%   dcr_dc - user defined damping control at rectifier cell
%   dci_dc - user defined damping control at inverter cell
%
%   Output: 
%   bus_sol - solved ac bus specification file
%  	line_sol - solved ac line specification file
% 	line_flow - calculated flow on ac lines
%  	rec_par - rectifier parameters
%           	- rec_par(:,1) firing angle degrees
%           	- rec_par(:,2) dc voltage kV
%               - rec_par(:,3) dc power MW
%               - rec_par(:,4) equi HT bus voltage kV
%   inv_par - inverterer parameters
%             	- inv_par(:,1) extinction angle degrees
%             	- inv_par(:,2) dc voltage kV
%             	- inv_par(:,3) dc power MW
%             	- inv_par(:,4) equi HT bus voltage kV
% 	line_par - dc line current kA
%
%   History:
%   Date        Time    Engineer        Description
%   02/xx/97    XX:XX   Graham Rogers  	Version 1.0
%   (c) copyright Joe Chow 1991-1997  - All rights reserved
%   07/15/20    13:43   Thad Haines     Revised format of globals and internal function documentation
%   07/29/20    15:20   Thad Haines     jay -> 1j

global g

disp('load flow with HVDC')


% perform load flow iterations
errv = 0;
itermax = 40;
iter = 0;
bus_old = bus;

while (errv == 0&&iter<=itermax)
  iter = iter + 1;
  [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30, ...
                                 1.0,'n',2);
  if iter==1
     dc_indx(bus,line,g.dc.dci_dc,g.dc.dcr_dc);
  end

    % perform dc load flow
    [rec_par,inv_par,line_par,g.dc.tap,Sr,Si] = dc_lf(bus_sol,line_sol,g.dc.dci_dc,g.dc.dcr_dc);

    % set taps in load flow data
    line_sol(g.dc.rec_ac_line,6) = g.dc.tap(g.dc.r_idx);
    line_sol(g.dc.inv_ac_line,6) = g.dc.tap(g.dc.i_idx);
    % set loads at rectifier and inverter
    bus_sol(g.dc.rec_ac_bus,6) = real(Sr);
    bus_sol(g.dc.inv_ac_bus,6) = real(Si);
    bus_sol(g.dc.rec_ac_bus,7) = imag(Sr);
    bus_sol(g.dc.inv_ac_bus,7) = imag(Si);
    Sr_old = bus(g.dc.rec_ac_bus,6)+1j*bus(g.dc.rec_ac_bus,7);
    Si_old = bus(g.dc.inv_ac_bus,6)+1j*bus(g.dc.inv_ac_bus,7);
    % check convergence
    errdc = max(abs([(Sr_old-Sr), (Si_old-Si)]));
    if errdc>1e-5
      bus = bus_sol;
      line = line_sol;
    else
      errv = 1;
     % perform ac load flow
     [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,30, ...
        1.0,'n',2);
     % reperform dc load flow
     [rec_par,inv_par,line_par,g.dc.tap,Sr,Si] = dc_lf(bus_sol,line_sol,g.dc.dci_dc,g.dc.dcr_dc);
     
     % set dc quantities based on final ac load flow
     P = bus_sol(g.dc.ac_bus,6);
     Q = bus_sol(g.dc.ac_bus,7);
     Vac = bus_sol(g.dc.ac_bus,2);
     Vang = bus_sol(g.dc.ac_bus,3)*pi/180;
     i_ac = (P - 1j*Q)./Vac./exp(-1j*Vang);
     Vac = Vac.*bus_sol(g.dc.ac_bus,13); %convert to kV
     % convert ac currents to dc
     idceq = abs(i_ac)*pi*g.sys.basmva/3/sqrt(2)...
             ./bus_sol(g.dc.ac_bus,13)./g.dc.dcsp_con(:,6);
     % convert ac currents to kA
     i_ac = i_ac*g.sys.basmva/sqrt(3)./bus_sol(g.dc.ac_bus,13);
     %calculate equivalent HT bus voltage
     xequ = sqrt(3)*g.dc.dcsp_con(:,5)./g.dc.dcsp_con(:,6);% eq transformer reactance
     g.dc.VHT = abs(Vac.*exp(1j*Vang) + 1j*xequ.*i_ac);
     Vdo = 3*sqrt(2)*g.dc.VHT.*g.dc.dcsp_con(:,6)/pi;% ideal dc voltages
     %calculate dc voltages
     dc_ang(g.dc.r_idx,1)  = rec_par(:,1)*pi/180;
     dc_ang(g.dc.i_idx,1) = inv_par(:,1)*pi/180;
     Rc = g.dc.dcsp_con(:,6).*g.dc.dcsp_con(:,5)*3/pi;% in series as far as dc is concerned
     g.dc.Vdc = Vdo.*cos(dc_ang) - Rc.*idceq;  
   end
end
if iter >= itermax
   imstr = int2str(itermax);
   disp(['dc load flow not converged in',imstr,' iterations'])
   error('stop')
else
   % dc and ac converged
   itstr= int2str(iter);
   disp([itstr,' dc load flow iterations'])
   disp('ac/dc solution converged')
end
return