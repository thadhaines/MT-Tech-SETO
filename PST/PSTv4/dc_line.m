function  dc_line(i,k,kdc,bus,flag)
%DC_LINE Models HVDC line dynamics
% DC_LINE Models HVDC line dynamics
%
% Syntax: dc_line(i,k,kdc,bus,flag)
%
%   NOTES:
%
%   Input:
%   i - 0 vector computaion only for HVDC control
%   k - integer time (data index)
%   kdc - integer time for dc (dc data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - dynamics computation and state state matrix building
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   02/xx/97    XX:XX   Graham Rogers  	Version 1.0
%   (c) Copyright 1991-1997 Joe H. Chow - All Rights Reserved
%   07/15/20    11:04   Thad Haines     Revised format of globals and internal function documentation

global g

k = fix(1+(kdc+1)/10);
% check for dcline data
if ~isempty(g.dc.dcsp_con)
    if flag==0
        % initialization
        % check that inductances are specified
        if ~isempty(g.dc.no_ind_idx)
            error('you must specify inductances for dc lines and smoothing reactors')
        end
        
        if i==0
            % vector computation
            g.dc.i_dcr(:,1) = g.dc.i_dc(g.dc.r_idx,1);
            g.dc.i_dci(:,1) = g.dc.i_dc(g.dc.i_idx,1);
            g.dc.v_dcc(:,1) = g.dc.Vdc(g.dc.r_idx,1) - g.dc.i_dcr(:,1).*g.dc.dcl_con(:,3)/2.0;
            if g.dc.l_no_cap~=0
                g.dc.dc_pot(g.dc.no_cap_idx,3) = zeros(g.dc.l_no_cap,1);
                g.dc.dc_pot(g.dc.no_cap_idx,1) = 1000*ones(g.dc.l_no_cap,1)./(g.dc.dcl_con(g.dc.no_cap_idx,4) ...
                    + g.dc.dcl_con(g.dc.no_cap_idx,6) + g.dc.dcl_con(g.dc.no_cap_idx,7));
                g.dc.dc_pot(g.dc.no_cap_idx,2) = g.dc.dcl_con(g.dc.no_cap_idx,3);
                g.dc.dc_pot(g.dc.no_cap_idx,4) = g.dc.dc_pot(g.dc.no_cap_idx,1);
                g.dc.dc_pot(g.dc.no_cap_idx,5) = g.dc.dc_pot(g.dc.no_cap_idx,2);
            end
            if g.dc.l_cap ~= 0
                g.dc.dc_pot(g.dc.cap_idx,1) = 1000*ones(g.dc.l_cap,1)./(g.dc.dcl_con(g.dc.cap_idx,4)*.5 + g.dc.dcl_con(g.dc.cap_idx,6));
                g.dc.dc_pot(g.dc.cap_idx,2) = g.dc.dcl_con(g.dc.cap_idx,3)/2;
                g.dc.dc_pot(g.dc.cap_idx,3) = 1e6*ones(g.dc.l_cap,1)./g.dc.dcl_con(g.dc.cap_idx,5);
                g.dc.dc_pot(g.dc.cap_idx,4) = 1000*ones(g.dc.l_cap,1)./(g.dc.dcl_con(g.dc.cap_idx,4)*.5 + g.dc.dcl_con(g.dc.cap_idx,7));
                g.dc.dc_pot(g.dc.cap_idx,5) = g.dc.dcl_con(g.dc.cap_idx,3)/2;
            end
        else
            error(' no non-vector computation in HVDC')
        end
    end
    if flag == 1
        % network inter face - no calculation required
    end
    if flag == 2
        % rate of change of states
        if i== 0
            %vector compuation
            if g.dc.l_cap~=0
                g.dc.di_dcr(g.dc.cap_idx,kdc) = - g.dc.dc_pot(cg.dc.ap_idx,1)...
                    .*(g.dc.dc_pot(g.dc.cap_idx,2).*g.dc.i_dcr(g.dc.cap_idx,kdc) + g.dc.v_dcc(g.dc.cap_idx,kdc) - g.dc.Vdc(g.dc.r_idx(g.dc.cap_idx),kdc));
                g.dc.di_dci(g.dc.cap_idx,kdc) = -g.dc.dc_pot(g.dc.cap_idx,4)...
                    .*(g.dc.dc_pot(g.dc.cap_idx,5).*g.dc.i_dci(g.dc.cap_idx,kdc) - g.dc.v_dcc(g.dc.cap_idx,kdc) + g.dc.Vdc(g.dc.i_idx(g.dc.cap_idx),kdc));
                g.dc.dv_dcc(g.dc.cap_idx,kdc) = g.dc.dc_pot(g.dc.cap_idx,3).*(g.dc.i_dcr(g.dc.cap_idx,kdc) - g.dc.i_dci(g.dc.cap_idx,kdc));
            end
            if g.dc.l_no_cap~=0
                g.dc.di_dcr(g.dc.no_cap_idx,kdc) = (-g.dc.dc_pot(g.dc.no_cap_idx,2).*g.dc.i_dcr(g.dc.no_cap_idx,kdc) + ...
                    g.dc.Vdc(g.dc.r_idx(g.dc.no_cap_idx),kdc)-g.dc.Vdc(g.dc.i_idx(g.dc.no_cap_idx),kdc))...
                    .*g.dc.dc_pot(g.dc.no_cap_idx,1);
                g.dc.di_dci(g.dc.no_cap_idx,kdc) = g.dc.di_dcr(g.dc.no_cap_idx,kdc);
                g.dc.dv_dcc(g.dc.no_cap_idx) = 0;
            end
        else
            error(' no non-vector calculation in HVDC')
        end
    end
end
