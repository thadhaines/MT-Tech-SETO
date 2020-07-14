function [xr,dxr,xi,dxi] = ...
   dc_sim(k,kk,dcr,dci,xr,xi,bus_sim,hdc_sol)
%     %% HVDC link variables - 63
%     global  dcsp_con  dcl_con  dcc_con
%     global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
%     global  inv_ac_line  rec_ac_line ac_line dcli_idx
%     global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
%     global  Vdc  i_dc P_dc i_dcinj dc_pot alpha gamma VHT dc_sig  cur_ord dcr_dsig dci_dsig
%     global  ric_idx  rpc_idx Vdc_ref dcc_pot
%     global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
%     global  ndcr_ud ndci_ud dcrud_idx dciud_idx dcrd_sig dcid_sig
% 
%     % States
%     %line
%     global i_dcr i_dci  v_dcc
%     global di_dcr  di_dci  dv_dcc
%     global dc_dsig % added 07/13/20 -thad
%     %rectifier
%     global v_conr dv_conr
%     %inverter
%     global v_coni dv_coni
%     
%     % added to global dc
%     global xdcr_dc dxdcr_dc xdci_dc dxdci_dc angdcr angdci t_dc
%     global dcr_dc dci_dc % damping control
%     global  ldc_idx


% global  dcsp_con  dcl_con dcc_con 
% global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
% global  Vdc  i_dc  dcc_pot  alpha  gamma Vdc_ref
% global  cur_ord dc_sig dcc_pot dcr_dsig dci_dsig ndcr_ud ndci_ud dcrd_sig dcid_sig
% global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
% 
% % States
% %line
% global i_dcr i_dci  v_dcc
% global di_dcr  di_dci  dv_dcc  
% %rectifier
% global v_conr  
% global dv_conr  
% %inverter
% global v_coni
% global dv_coni 


global g


% predictor
kdc = 10*(k-1)+ kk; 
jdc=kdc+1;

if g.dc.n_conv~=0
   if g.dc.ndcr_ud~=0
      tot_states=0;
      for jj = 1:g.dc.ndcr_ud
         ydcrmx = dcr{jj,5};ydcrmn = dcr{jj,6};
         rec_num = dcr{jj,2};
         st_state = tot_states+1; dcr_states = dcr{jj,7}; tot_states = tot_states+dcr_states; 
         [g.dc.dcr_dsig(rec_num,k),xr(st_state:tot_states,1),dxr(st_state:tot_states,1)] =...
            dcr_sud(jj,kdc,2,dcr{jj,1},g.dc.dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,1));
      end
   else
      dxr(1,1)=0;
   end
   if g.dc.ndci_ud~=0
      tot_states=0;
      for jj = 1:g.dc.ndci_ud
         ydcimx = dci{jj,5};ydcimn = dci{jj,6};
         inv_num = dci{jj,2};
         st_state = tot_states+1; dci_states = dci{jj,7}; tot_states = tot_states+dci_states; 
         [g.dc.dci_dsig(inv_num,k),xi(st_state:tot_states,1),dxi(st_state:tot_states,1)] =...
            dci_sud(jj,kdc,2,dci{jj,1},g.dc.dcid_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,1));
      end
   else
      dxi(1,1)=0;
   end
   dc_cont(0,k,kdc,bus_sim,2);
   dc_line(0,k,kdc,bus_sim,2);
end
g.dc.v_conr(:,jdc) = g.dc.v_conr(:,kdc) + hdc_sol*g.dc.dv_conr(:,kdc);
g.dc.v_coni(:,jdc) = g.dc.v_coni(:,kdc) + hdc_sol*g.dc.dv_coni(:,kdc);
g.dc.i_dcr(:,jdc) = g.dc.i_dcr(:,kdc) + hdc_sol*g.dc.di_dcr(:,kdc);
g.dc.i_dci(:,jdc) = g.dc.i_dci(:,kdc) + hdc_sol*g.dc.di_dci(:,kdc);
g.dc.v_dcc(:,jdc) = g.dc.v_dcc(:,kdc) + hdc_sol*g.dc.dv_dcc(:,kdc);

% function output...
xr(:,2) = xr(:,1) + hdc_sol* dxr(:,1);
xi(:,2) = xi(:,1) + hdc_sol* dxi(:,1); 

% begining of corrector? -thad 07/14/20
% dc flag = 1?
dc_cont(0,k,jdc,bus_sim,1);% recalculate alpha and gamma 
dc_vidc(k,kdc); % update Vdc and i_dc
if g.dc.ndcr_ud~=0
   tot_states=0;
   for jj = 1:g.dc.ndcr_ud
      ydcrmx = dcr{jj,5};ydcrmn = dcr{jj,6};
      rec_num = dcr{jj,2};
      st_state = tot_states+1; dcr_states = dcr{jj,7}; tot_states = tot_states+dcr_states; 
      [g.dc.dcr_dsig(rec_num,k),xr(st_state:tot_states,2),dxr(st_state:tot_states,2)] =...
         dcr_sud(jj,jdc,2,dcr{jj,1},g.dc.dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,2));
   end
else
   dxr(1,2)=0;
end
if g.dc.ndci_ud~=0
   tot_states=0;
   for jj = 1:g.dc.ndci_ud
      ydcimx = dci{jj,5};ydcimn = dci{jj,6};
      inv_num = dci{jj,2};
      st_state = tot_states+1; dci_states = dci{jj,7}; tot_states = tot_states+dci_states; 
      [g.dc.dci_dsig(inv_num,j),xi(st_state:tot_states,2),dxi(st_state:tot_states,2)] =...
         dci_sud(jj,jdc,flag,dci{jj,1},dci_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,2));
   end
else
   dxi(1,2)=0;
end

% dc flag = 2

dc_cont(0,k,jdc,bus_sim,2);
dc_line(0,k,jdc,bus_sim,2);

g.dc.v_conr(:,jdc) = g.dc.v_conr(:,kdc) + 0.5*hdc_sol*(g.dc.dv_conr(:,kdc)+g.dc.dv_conr(:,jdc));
g.dc.v_coni(:,jdc) = g.dc.v_coni(:,kdc) + 0.5*hdc_sol*(g.dc.dv_coni(:,kdc)+g.dc.dv_coni(:,jdc));
g.dc.i_dcr(:,jdc) = g.dc.i_dcr(:,kdc) + 0.5*hdc_sol*(g.dc.di_dcr(:,kdc)+g.dc.di_dcr(:,jdc));
g.dc.i_dci(:,jdc) = g.dc.i_dci(:,kdc) + 0.5*hdc_sol*(g.dc.di_dci(:,kdc)+g.dc.di_dci(:,jdc));
g.dc.v_dcc(:,jdc) = g.dc.v_dcc(:,kdc) + 0.5*hdc_sol*(g.dc.dv_dcc(:,kdc)+g.dc.dv_dcc(:,jdc));

xr(:,2) = xr(:,1) + 0.5*hdc_sol* (dxr(:,1)+dxr(:,2));
xi(:,2) = xi(:,1) + 0.5*hdc_sol* (dxi(:,1)+dxi(:,2)); 

dc_cont(0,k,jdc,bus_sim,1);%recalculate alpha and gamma
dc_vidc(k,jdc);% update Vdc and i_dc
