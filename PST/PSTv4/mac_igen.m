function [bus_new] = mac_igen(i,k,bus,flag)
%MAC_IGEN is a induction generator model.
% MAC_IGEN is a single cage induction generator model with
% no leakage inductance saturation.
%
% Syntax: [bus_new] = mac_igen(i,k,bus,flag)
%
%   NOTES:  Induction generator pick up power from negative load
%           data format igen_con
%           1 - induction generator number
%           2 - busnumber
%           3 - base MVA
%           4 - rs
%           5 - xs -stator leakage reactance
%           6 - Xm - magnetizing reactance
%           7 - rr
%           8 - xr - rotor leakage reactance
%           9 - H  - inertia constant generator + turbine in sec
%           15 - fraction of bus load power taken by generator
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
%   bus_new     - is a modified bus matrix with the induction generator
%               active and reactive powers subtracted from the original
%               loads at the generator bus.
%
%   History:
%   Date        Time    Engineer        Description
%   08/18/97    18:28   Graham Rogers   Version 1.0
%   (c) Copyright Cherry Tree Scientific Software 1997 _ All rights reserved
%   07/13/20    10:26   Thad Haines     Revised format of globals and internal function documentation
%   07/29/20    15:20   Thad Haines     jay -> 1j

global g

bus_new = bus;
if ~isempty(g.igen.igen_con)
    if flag == 0
        % initialisation
        disp(' initializing')
        if i == 0
            %vector computation
            g.igen.n_ig = length(g.igen.igen_con(:,1));
            g.igen.igbus=g.bus.bus_int(g.igen.igen_con(:,2));
            g.igen.igen_pot = zeros(g.igen.n_ig,7);
            g.igen.igen_pot(:,1)=g.sys.basmva./g.igen.igen_con(:,3); %scaled mva base
            g.igen.igen_pot(:,2)=ones(g.igen.n_ig,1); %base kv
            ig_vm(:,1)=bus(g.igen.igbus,2); %ind gen terminal voltage mag
            ig_ang(:,1)=bus(g.igen.igbus,3)*pi/180; %ind gen term voltage angle
            v=ig_vm(:,1).*exp(1j*ig_ang(:,1));
            g.igen.vdig(:,1)=real(v);
            g.igen.vqig(:,1)=imag(v);
            g.igen.pig(:,1)= bus(g.igen.igbus,6).*g.igen.igen_con(:,15);%ind generator power
            %modify bus load power
            bus_new(g.igen.igbus,6)=bus(g.igen.igbus,6)-g.igen.pig(:,1);
            g.igen.igen_pot(:,3)=g.igen.igen_con(:,5)+g.igen.igen_con(:,6);%Xs
            g.igen.igen_pot(:,4)=g.igen.igen_con(:,8)+g.igen.igen_con(:,6);%Xr
            g.igen.igen_pot(:,5)=g.igen.igen_con(:,5)+g.igen.igen_con(:,6).*...
                g.igen.igen_con(:,8)./g.igen.igen_pot(:,4);%Xsp
            g.igen.igen_pot(:,6)=g.igen.igen_pot(:,3)-g.igen.igen_pot(:,5);%(Xs-Xsp)
            g.igen.igen_pot(:,7)=g.sys.basrad*g.igen.igen_con(:,7)./g.igen.igen_pot(:,4); %1/Tr
            
            rs=g.igen.igen_con(:,4);
            xs=g.igen.igen_con(:,5);
            Xm=g.igen.igen_con(:,6);
            rr=g.igen.igen_con(:,7);
            xr=g.igen.igen_con(:,8);
            % find initial slip
            
            slip_old=zeros(g.igen.n_ig,1);
            slip_new=ones(g.igen.n_ig,1);
            %Newton-Raphson iteration to determine initial slip for
            %induction generators
            iter = 0;
            err=max(abs(slip_new-slip_old));
            while err>=1e-8 && iter<30
                iter=iter+1;
                y=g.sys.basrad.*slip_old./g.igen.igen_pot(:,7);
                denom = ones(g.igen.n_ig,1)+y.*y;
                zr=rs + y.*g.igen.igen_pot(:,6)./denom;
                zi=g.igen.igen_pot(:,5)+g.igen.igen_pot(:,6)./denom;
                dzr=g.igen.igen_pot(:,6).*(ones(g.igen.n_ig,1)-...
                    y.*y)./denom./denom;
                dzi=-2*g.igen.igen_pot(:,6).*y./denom./denom;
                zmod2=zr.*zr+zi.*zi;
                dp=v.*conj(v).*(dzr.*zmod2-...
                    2*zr.*(dzr.*zr+dzi.*zi));
                dp=dp./zmod2./zmod2;
                peig =v.*conj(v).*zr./zmod2;
                ynew=y-(peig - g.igen.pig(:,1).*g.igen.igen_pot(:,1))./dp;
                slip_new = ynew.*g.igen.igen_pot(:,7)/g.sys.basrad;
                err = max(abs(slip_new-slip_old));
                slip_old=slip_new;
            end
            if iter >=30
                error('induction generator slip calculation failed to converge')
            end
            g.igen.slig(:,1)=slip_new;
            y=g.sys.basrad*g.igen.slig(:,1)./g.igen.igen_pot(:,7);
            denom= ones(g.igen.n_ig,1)+y.*y;
            zr=rs+y.*g.igen.igen_pot(:,6)./denom;
            zi=g.igen.igen_pot(:,5)+g.igen.igen_pot(:,6)./denom;
            iig =v./(zr+1j*zi);
            sig = v.*conj(iig);
            peig = real(sig);
            qeig = imag(sig);
            %complex initial rotor states
            vp = v - (rs+ 1j* g.igen.igen_pot(:,5)).*iig;
            g.igen.vdpig = real(vp);
            g.igen.vqpig =imag(vp);
            % initial prime mover torque
            g.igen.tmig(:,1) = real(vp.*conj(iig));
            g.igen.idig(:,1)=real(iig)./g.igen.igen_pot(:,1);
            g.igen.iqig(:,1)=imag(iig)./g.igen.igen_pot(:,1);
            % modify qload
            g.igen.qig(:,1) = qeig./g.igen.igen_pot(:,1);
            bus_new(g.igen.igbus,7)=bus(g.igen.igbus,7)-g.igen.qig(:,1);
        else
            % generator by generator initialization
            error('Only vector computation allowed in induction generators')
        end
    end
    if flag == 1
        %network interface
        %no interface required for induction generators
    end
    if flag == 2
        %induction generator dynamics calculation
        if i == 0
            %vector calculation
            idigm=g.igen.idig(:,k).*g.igen.igen_pot(:,1);%convert to machine base
            iqigm=g.igen.iqig(:,k).*g.igen.igen_pot(:,1);
            %Brereton, Lewis and Young motor model
            g.igen.dvdpig(:,k)=-(iqigm.*g.igen.igen_pot(:,6)+g.igen.vdpig(:,k)).*...
                g.igen.igen_pot(:,7)+g.igen.vqpig(:,k).*g.igen.slig(:,k)*g.sys.basrad;
            g.igen.dvqpig(:,k)=(idigm.*g.igen.igen_pot(:,6)-g.igen.vqpig(:,k)).*...
                g.igen.igen_pot(:,7)-g.igen.vdpig(:,k).*g.igen.slig(:,k)*g.sys.basrad;
            g.igen.dslig(:,k)=(g.igen.tmig(:,k)-g.igen.vdpig(:,k).*...
                idigm-g.igen.vqpig(:,k).*iqigm)/2./g.igen.igen_con(:,9);
        else
            error(' vector computation only for induction generator')
        end
    end
    if flag == 3
        %linearize
        %add code later
    end
end