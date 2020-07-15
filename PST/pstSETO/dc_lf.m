function [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus,line,dci_dc,dcr_dc)
%DC_LF Dc initialization model
% DC_LF Dc initialization model
%
% Syntax: [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus,line,dci_dc,dcr_dc)
%
%   NOTES:  Output tap is not global during function computation.
%           dci_dc and dcr_dc are the same as g.dc.dci_dc g.dcr_dc, though not used?
%           Output documentation is copied from original - accuracy unknown.
% 
%   Input: 
%   bus - ac bus matrix (from data or as modified in loadflow)
%   line - ac line matrix
%   dcr_dc - user defined damping control at rectifier cell
%   dci_dc - user defined damping control at inverter cell
%
%   Output: 
%   rectifier firing angle
%   inverter extinction angle
%   converter transformer tap settings
%   complex ac load at LT terminals on system base
%
%   History:
%   Date        Time    Engineer        Description
%   10/xx/96    XX:XX   Graham Rogers  	Version 1
%   (c) Copyright Joe Chow 1996 - All rights reserved
%   07/15/20    10:55   Thad Haines     Revised format of globals and internal function documentation

% NOTE: output tap NOT global. ? -thad 07/14/20

global g
jay = sqrt(-1);
% determine dc indexes
dc_indx(bus,line,g.dc.dci_dc,g.dc.dcr_dc);

if g.dc.n_conv~=0 && g.dc.n_dcl~=0
  g.dc.Vdc=zeros(g.dc.n_conv,1);% initialize to zero vector
  dc_bus = g.dc.dcsp_con(:,1);
  Vac = bus(g.dc.ac_bus,2);
  Vang = bus(g.dc.ac_bus,3)*pi/180;
  al_min = g.dc.dcsp_con(g.dc.r_idx,7)*pi/180;
  al_max = g.dc.dcsp_con(g.dc.r_idx,8)*pi/180;
  rdc_bus = dc_bus(g.dc.r_idx);
  ga_min = g.dc.dcsp_con(g.dc.i_idx,7)*pi/180;
  ga_max = g.dc.dcsp_con(g.dc.i_idx,8)*pi/180;
  idc_bus = dc_bus(g.dc.i_idx);
  % get transformer taps from load flow
  tap = line(g.dc.ac_line,6);
  g.dc.tapr = line(g.dc.rec_ac_line,6);
  g.dc.tapi = line(g.dc.inv_ac_line,6);
  g.dc.tmax = line(g.dc.ac_line,8);
  g.dc.tmaxr = line(g.dc.rec_ac_line,8);
  g.dc.tmaxi = line(g.dc.inv_ac_line,8);
  g.dc.tmin = line(g.dc.ac_line,9);
  g.dc.tmini = line(g.dc.inv_ac_line,9);
  g.dc.tminr = line(g.dc.rec_ac_line,9);
  g.dc.tstep = line(g.dc.ac_line,10);
  g.dc.tstepi = line(g.dc.inv_ac_line,10);
  g.dc.tstepr = line(g.dc.rec_ac_line,10);
  % set Vdc limits
  Vdc_rated = g.dc.dcsp_con(g.dc.i_idx,4);
  % set max and min values of dc voltage at 1% from nominal
  Vdc_max = 1.01*Vdc_rated;
  Vdc_min = 0.99*Vdc_rated;
  % in load flow the rectifier and inverter are represented by 
  % P and Q loads
  % Get P and Q in MW
  P = bus(g.dc.ac_bus,6)*g.sys.basmva;
  Q = bus(g.dc.ac_bus,7)*g.sys.basmva;
  % convert ac voltages to kV line-to-line
  Vac = Vac.*bus(g.dc.ac_bus,13);
  % calculate commutating voltage
  xequ = sqrt(3)*g.dc.dcsp_con(:,5)./g.dc.dcsp_con(:,6);% eq transformer reactance
  Rc = g.dc.dcsp_con(:,6).*g.dc.dcsp_con(:,5)*3/pi;% in series as far as dc is concerned
  iac = (P-jay*Q)./(Vac.*exp(-jay*Vang))/sqrt(3);
  ang_iac = angle(iac);
  % get specified idc in kA
  idc = g.dc.dcl_con(g.dc.dcli_idx,8)./g.dc.dcsp_con(g.dc.i_idx,4);%dc power/dc voltage
  % adjust ac currents to have same angle as load flow but magnitude
  % corresponding to desired dc current
  iaci = idc.*g.dc.dcsp_con(g.dc.i_idx,6).*exp(jay*ang_iac(g.dc.i_idx))*sqrt(6)/pi;
  iacr = idc.*g.dc.dcsp_con(g.dc.r_idx,6).*exp(jay*ang_iac(g.dc.r_idx))*sqrt(6)/pi;
  % calculate the equivalent HT bus voltage (commutating voltage)
  g.dc.VHT(g.dc.i_idx,1) = abs(Vac(g.dc.i_idx).*exp(jay*Vang(g.dc.i_idx)) + jay*xequ(g.dc.i_idx).*iaci);
  g.dc.VHT(g.dc.r_idx,1) = abs(Vac(g.dc.r_idx).*exp(jay*Vang(g.dc.r_idx)) + jay*xequ(g.dc.r_idx).*iacr);
  Vdo = 3*sqrt(2)*g.dc.VHT.*g.dc.dcsp_con(:,6)/pi; % ideal dc voltages
  rdc = g.dc.dcl_con(g.dc.dcli_idx,3);% dc line resistance
  cm = ones(g.dc.n_dcl,1)-g.dc.dcl_con(g.dc.dcli_idx,9)/100;

  % Assume that inverters are operating at gamm_min
  mode = ones(g.dc.n_dcl,1);
  [g.dc.gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo(g.dc.i_idx),Rc(g.dc.i_idx));
  
  % calculate rectifier dc voltage
  g.dc.Vdc(g.dc.r_idx) = g.dc.Vdc(g.dc.i_idx) + idc.*rdc;

  % calculate alpha

  [g.dc.alpha,mode_new,idc] = rec_lf(idc,al_min,al_max,Vdo(g.dc.r_idx),Rc(g.dc.r_idx),cm);
  
  % check for mode change
  mode_ch = sum(mode-mode_new);
  if mode_ch ~=0

    % mode has changed - need to recalculate the inverter conditions
    [g.dc.gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo(g.dc.i_idx),Rc(g.dc.i_idx));
    
  end
%   alpha = alpha; % oh? -thad 07/13/20
%   gamma = gamma;
  % recalculate Vdo based on the modified firing and extinction angles
  Vdo(g.dc.r_idx) = (g.dc.Vdc(g.dc.r_idx)+Rc(g.dc.r_idx).*idc)./cos(g.dc.alpha*pi/180);
  Vdo(g.dc.i_idx) = (g.dc.Vdc(g.dc.i_idx)+Rc(g.dc.i_idx).*idc)./cos(g.dc.gamma*pi/180);
  % calculate ac power factor at LT bus 
  cphi = g.dc.Vdc./Vdo;
  oor_idx = find(abs(cphi)>=1);
  if ~isempty(oor_idx)
    bad_num = num2str(oor_idx');
    disp(['dc_lf failed to converge at converter ',bad_num])
    error('stop')
  end
  % converter power in per unit on ac system base 
  Pr = g.dc.Vdc(g.dc.r_idx).*idc/g.sys.basmva;
  Pi = -g.dc.Vdc(g.dc.i_idx).*idc/g.sys.basmva;
  tphi = sqrt(ones(g.dc.n_conv,1)-cphi.*cphi)./cphi;
  Qr = Pr.*tphi(g.dc.r_idx);
  Qi = -Pi.*tphi(g.dc.i_idx);
  % convert to load at the converter LT bus
  iacr = idc.*g.dc.dcsp_con(g.dc.r_idx,6)*sqrt(6)/pi;
  %convert to per unit
  iacr = iacr*sqrt(3).*bus(g.dc.rec_ac_bus,13)/g.sys.basmva;
  % convert xequ to per unit on ac system base
  xequ = xequ*g.sys.basmva./bus(g.dc.ac_bus,13)./bus(g.dc.ac_bus,13)/sqrt(3);
  iaci = idc.*g.dc.dcsp_con(g.dc.i_idx,6)*sqrt(6)/pi;
  % convert to per unit
  iaci = iaci*sqrt(3).*bus(g.dc.inv_ac_bus,13)/g.sys.basmva;
  % calculate LT reactive power
  Qr = Qr - xequ(g.dc.r_idx).*iacr.*iacr;
  Qi = Qi - xequ(g.dc.i_idx).*iaci.*iaci;
  Sr = (Pr + jay*Qr);
  Si = (Pi + jay*Qi);
  rec_par = zeros(g.dc.n_dcl,3);
  inv_par = rec_par;
  rec_par(:,1) = g.dc.alpha;
  inv_par(:,1) = g.dc.gamma;
  rec_par(:,2) = g.dc.Vdc(g.dc.r_idx);
  inv_par(:,2) = g.dc.Vdc(g.dc.i_idx);
  rec_par(:,3) = g.dc.Vdc(g.dc.r_idx).*idc;
  inv_par(:,3) = g.dc.Vdc(g.dc.i_idx).*idc;
  line_par = idc;
  tap(g.dc.i_idx) = g.dc.tapi;
  tap(g.dc.r_idx) = g.dc.tapr;
end
return
