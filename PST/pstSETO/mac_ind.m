function bus_new = mac_ind(i,k,bus,flag)
%MAC_IND is an induction motor model.
% MAC_IND is an induction motor model with deep bar,double cage 
% and leakage inductance saturation.
%
% Syntax: [bus_new] = mac_ind(i,k,bus,flag)
%
%   NOTES:  Data format ind_con
%           1 - motor number
%           2 - busnumber
%           3 - base MVA
%           4 - rs
%           5 - xs -stator leakage reactance
%           6 - Xm - magnetizing reactance
%           7 - rr
%           8 - xr - rotor leakage reactance
%           9 - H  - inertia constant motor + load in sec
%           10 - r2 - double cage resistance
%           11 - x2 - intercage reactance
%           12 - dbf - deep bar factor
%           13 - isat - current at which leakage inductance starts to saturate
%           15 -    fraction of bus load power taken by motor 
%                   if entry 15 is zero, it is assumed that the motor is 
%                   to be started on the specified bus
% 
%   Input: 
%   i - motor number
%          - 0 for vectorized computation
%   k - integer time (data index)
%   bus - solved loadflow bus data
%   flag -  0 - initialization
%          	1 - network interface computation
%          	2 - generator dynamics computation and state state matrix building
%
%   Output: 
%   bus_new -   is bus with the power and reactive power loads
%               modified to subtract the motor loads
%               modification is made only when the motors are initialized
%               i.e. flag = 0
%
%   History:
%   Date        Time    Engineer        Description
%   11/xx/95    xx:xx   Graham Rogers   Version 1.0
%   06/xx/97    xx:xx   Graham Rogers   Version 2.0 added deep bar, 
%                                       double cage and leakage inductance saturation
%   07/13/20    12:22   Thad Haines     Revised format of globals and internal function documentation

global g

jay=sqrt(-1);
bus_new=bus;

if ~isempty(g.ind.ind_con)
    if flag == 0;
        % initialisation
        if i == 0;
            %vector computation
            motnum=length(g.ind.ind_con(:,1));
            g.ind.motbus=g.bus.bus_int(g.ind.ind_con(:,2));
            g.ind.ind_pot(:,1)=g.sys.basmva./g.ind.ind_con(:,3); %scaled mva base
            g.ind.ind_pot(:,2)=ones(motnum,1); %base kv
            mot_vm(:,1)=bus(g.ind.motbus,2); %motor terminal voltage mag
            mot_ang(:,1)=bus(g.ind.motbus,3)*pi/180; %motor term voltage angle
            v=mot_vm(:,1).*exp(jay*mot_ang(:,1));
            g.ind.vdmot(:,1)=real(v);
            g.ind.vqmot(:,1)=imag(v);
            g.ind.p_mot(:,1)=bus(g.ind.motbus,6).*g.ind.ind_con(:,15);%motor power demand
            %modify bus load power
            bus_new(g.ind.motbus,6)=bus(g.ind.motbus,6)-g.ind.p_mot(:,1);
            
            rs=g.ind.ind_con(:,4);
            xs=g.ind.ind_con(:,5);
            Xm=g.ind.ind_con(:,6);
            rr=g.ind.ind_con(:,7);
            xr=g.ind.ind_con(:,8);
            rr2 = g.ind.ind_con(:,10);
            xr2 = g.ind.ind_con(:,11);
            dbf = g.ind.ind_con(:,12);
            isat = g.ind.ind_con(:,13);
            
            g.ind.dbc_idx = find(rr2~=0);% motors with double cage
            g.ind.db_idx = find(dbf~=0);% motors with deep bars
            g.ind.sat_idx = find(isat~=0);% motors with leakage inductance saturation
            g.ind.ind_pot(:,3)=xs+Xm;%Xs
            g.ind.ind_pot(:,4)=xr+Xm;%Xr
            g.ind.ind_pot(:,5)=xs+Xm.*xr./g.ind.ind_pot(:,4);%Xsp
            g.ind.ind_pot(:,6)=g.ind.ind_pot(:,3)-g.ind.ind_pot(:,5);%(Xs-Xsp)
            g.ind.ind_pot(:,7)=g.sys.basrad*rr./g.ind.ind_pot(:,4); %1/Tr
            
            % index of motors to be initialized for running 
            run_ind = find(g.ind.ind_con(:,15)~=0);
            motrun=length(run_ind);%number of running motors
            start_ind = find(g.ind.ind_con(:,15)==0);
            motstart = length(start_ind);% number of starting motors
            %assumes motor starting if power fraction zero
            % find initial slip
            slip_old=zeros(motnum,1);
            slip_new=ones(motnum,1);
            % reset ind_pot for double cage and deepbar rotor machines
            s = 0.01*ones(motnum,1);s(start_ind)=ones(motstart,1);
            if ~isempty(g.ind.dbc_idx)
                [rdc,xdc]=dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx),rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx),s(g.ind.dbc_idx)); 
                g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx)+xdc;
                g.ind.ind_pot(g.ind.dbc_idx,5)=xs(g.ind.dbc_idx)+Xm(g.ind.dbc_idx).*xdc./g.ind.ind_pot(g.ind.dbc_idx,4);%Xsp
                g.ind.ind_pot(g.ind.dbc_idx,6) = g.ind.ind_pot(g.ind.dbc_idx,3)-g.ind.ind_pot(g.ind.dbc_idx,5);
                g.ind.ind_pot(g.ind.dbc_idx,7)=g.sys.basrad*rdc./g.ind.ind_pot(g.ind.dbc_idx,4); %1/Tr
            end
            if ~isempty(g.ind.db_idx)
                [rdb,xdb]=deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx),s(g.ind.db_idx));
                g.ind.ind_pot(g.ind.db_idx,4) = Xm(g.ind.db_idx)+xr(g.ind.db_idx)+xdb;
                g.ind.ind_pot(g.ind.db_idx,5) = xs(g.ind.db_idx)+Xm(g.ind.db_idx).*(xr(g.ind.db_idx)+xdb)./g.ind.ind_pot(g.ind.db_idx,4);%Xsp
                g.ind.ind_pot(g.ind.db_idx,6) = g.ind.ind_pot(g.ind.db_idx,3)-g.ind.ind_pot(g.ind.db_idx,5);
                g.ind.ind_pot(g.ind.db_idx,7) = g.sys.basrad*rdb./g.ind.ind_pot(g.ind.db_idx,4); %1/Tr
            end
            
            
            % Set defaults for motor starting
            imot=zeros(motnum,1);
            pem = zeros(motnum,1);
            qem = zeros(motnum,1); % should be globals? -thad 07/13/20
            g.ind.vdp(:,1)=zeros(motnum,1);
            g.ind.vqp(:,1)=zeros(motnum,1);
            vp = imot;
            g.ind.t_init = ones(motnum,1); % default for motor starting
            %Newton-Raphson iteration to determine initial slip for
            %running motors
            
            if motrun~=0 %check that some motors are running
                iter = 0;
                err=max(abs(slip_new-slip_old));
                while (err>=1e-8) && (iter<30)
                    iter=iter+1;
                    y=g.sys.basrad.*slip_old(run_ind)./g.ind.ind_pot(run_ind,7);
                    denom = ones(motrun,1)+y.*y;
                    zr=rs(run_ind)+y.*g.ind.ind_pot(run_ind,6)./denom;
                    zi=g.ind.ind_pot(run_ind,5)+g.ind.ind_pot(run_ind,6)./denom;
                    dzr=g.ind.ind_pot(run_ind,6).*(ones(motrun,1)-...
                        y.*y)./denom./denom;
                    dzi=-2*g.ind.ind_pot(run_ind,6).*y./denom./denom;
                    zmod2=zr.*zr+zi.*zi;
                    dp=v(run_ind).*conj(v(run_ind)).*(dzr.*zmod2-...
                        2*zr.*(dzr.*zr+dzi.*zi));
                    dp=dp./zmod2./zmod2;
                    pem(run_ind)=v(run_ind).*conj(v(run_ind)).*zr./zmod2;
                    ynew=y-(pem(run_ind)- ...
                        g.ind.p_mot(run_ind,1).*g.ind.ind_pot(run_ind,1))./dp;
                    slip_new(run_ind)=ynew.*g.ind.ind_pot(run_ind,7)/g.sys.basrad;
                    err = max(abs(slip_new-slip_old));
                    slip_old=slip_new;
                    
                    if ~isempty(g.ind.dbc_idx)
                        [rdc,xdc]=dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx),rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx),slip_new(g.ind.dbc_idx)); 
                        g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx)+xdc;
                        g.ind.ind_pot(g.ind.dbc_idx,5)=xs(g.ind.dbc_idx)+Xm(g.ind.dbc_idx).*xdc./g.ind.ind_pot(g.ind.dbc_idx,4);%Xsp
                        g.ind.ind_pot(g.ind.dbc_idx,6) = g.ind.ind_pot(g.ind.dbc_idx,3)-g.ind.ind_pot(g.ind.dbc_idx,5);
                        g.ind.ind_pot(g.ind.dbc_idx,7)=g.sys.basrad*rdc./g.ind.ind_pot(g.ind.dbc_idx,4); %1/Tr
                    end
                    if ~isempty(g.ind.db_idx)
                        [rdb,xdb]=deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx),slip_new(g.ind.db_idx));
                        g.ind.ind_pot(g.ind.db_idx,4) = Xm(g.ind.db_idx)+xr(g.ind.db_idx)+xdb;
                        g.ind.ind_pot(g.ind.db_idx,5) = xs(g.ind.db_idx)+Xm(g.ind.db_idx).*(xr(g.ind.db_idx)+xdb)./g.ind.ind_pot(g.ind.db_idx,4);%Xsp
                        g.ind.ind_pot(g.ind.db_idx,6) = g.ind.ind_pot(g.ind.db_idx,3)-g.ind.ind_pot(g.ind.db_idx,5);
                        g.ind.ind_pot(g.ind.db_idx,7) = g.sys.basrad*rdb./g.ind.ind_pot(g.ind.db_idx,4); %1/Tr
                    end
                end
                if iter >=30
                    uiwait(msgbox('induction motor slip calculation failed to converge','mac_ind error','modal'))
                    return
                end
            end
            g.ind.slip(:,1)=slip_new;
            ind_ldto(0,1);
            y=g.sys.basrad*g.ind.slip(:,1)./g.ind.ind_pot(:,7);
            denom= ones(motnum,1)+y.*y;
            zr=rs+y.*g.ind.ind_pot(:,6)./denom;
            zi=g.ind.ind_pot(:,5)+g.ind.ind_pot(:,6)./denom;
            if ~isempty(run_ind)
                imot(run_ind)=v(run_ind)./(zr(run_ind)+jay*zi(run_ind));
                sm(run_ind)=v(run_ind).*conj(imot(run_ind));
                pem(run_ind)=real(sm(run_ind));
                qem(run_ind)=imag(sm(run_ind));
                %complex initial rotor states
                vp(run_ind) = v(run_ind) - (rs(run_ind)+ jay* g.ind.ind_pot(run_ind,5)).*imot(run_ind); 
                g.ind.vdp(run_ind,1)=real(vp(run_ind));
                g.ind.vqp(run_ind,1)=imag(vp(run_ind));
            end
            g.ind.idmot(:,1)=real(imot)./g.ind.ind_pot(:,1);
            g.ind.iqmot(:,1)=imag(imot)./g.ind.ind_pot(:,1);
            % modify qload 
            bus_new(g.ind.motbus,7)=bus(g.ind.motbus,7)-qem./g.ind.ind_pot(:,1);
            tlm = g.ind.vdp(:,k).*real(imot)+g.ind.vqp(:,k).*imag(imot);
            trat = tlm./g.ind.tload(:,1);
            % modify load specification to get initial load correct
            
            g.ind.mld_con(run_ind,[3 5])=diag(trat(run_ind))*g.ind.mld_con(run_ind,[3 5]);
        else
            error('motor by motor initialization not supported')  
        end
    end
    if flag == 1
        v = g.bus.bus_v(g.ind.motbus,k);
        g.ind.vdmot(:,k)=real(v);
        g.ind.vqmot(:,k)=imag(v);
    end
    if flag == 2
        %motor dynamics calculation
        if i == 0
            %vector calculation
            
            ind_ldto(0,k);
            idm=g.ind.idmot(:,k).*g.ind.ind_pot(:,1);%convert to machine base
            iqm=g.ind.iqmot(:,k).*g.ind.ind_pot(:,1);
            
            rs=g.ind.ind_con(:,4);
            xs=g.ind.ind_con(:,5);
            Xm=g.ind.ind_con(:,6);
            rr=g.ind.ind_con(:,7);
            xr=g.ind.ind_con(:,8);
            rr2 = g.ind.ind_con(:,10);
            xr2 = g.ind.ind_con(:,11);
            dbf = g.ind.ind_con(:,12);
           
            imot = abs(idm+jay*iqm);
            if ~isempty(g.ind.sat_idx)
                % saturation of leakage inductance
                ism = imot(g.ind.sat_idx);isat = g.ind.ind_con(g.ind.sat_idx,13);
                ir = jay*Xm(g.ind.sat_idx).*(idm(g.ind.sat_idx)+jay*iqm(g.ind.sat_idx))./(rr(g.ind.sat_idx)+jay*g.ind.ind_pot(g.ind.sat_idx,4)); 
                gs = dessat(ism,isat);
                gr = dessat(abs(ir),isat);
                xs(g.ind.sat_idx) = xs(g.ind.sat_idx).*(1+gs)/2;
                xr(g.ind.sat_idx) = xr(g.ind.sat_idx).*(1+gr)/2;
                g.ind.ind_pot(g.ind.sat_idx,3) = Xm(g.ind.sat_idx)+xs(g.ind.sat_idx);
                g.ind.ind_pot(g.ind.sat_idx,4) = Xm(g.ind.sat_idx)+xr(g.ind.sat_idx);
                g.ind.ind_pot(g.ind.sat_idx,5)=xs(g.ind.sat_idx)+Xm(g.ind.sat_idx).*xr(g.ind.sat_idx)./g.ind.ind_pot(g.ind.sat_idx,4);%Xsp
                g.ind.ind_pot(g.ind.sat_idx,6) = g.ind.ind_pot(g.ind.sat_idx,3)-g.ind.ind_pot(g.ind.sat_idx,5);
                g.ind.ind_pot(g.ind.sat_idx,7)=g.sys.basrad*rr(g.ind.sat_idx)./g.ind.ind_pot(g.ind.sat_idx,4); %1/Tr
            end
            if ~isempty(g.ind.dbc_idx)
                % reset double cage
                [rdc,xdc]=dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx),rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx),g.ind.slip(g.ind.dbc_idx,k)); 
                g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx)+xdc;
                g.ind.ind_pot(g.ind.dbc_idx,5)=xs(g.ind.dbc_idx)+Xm(g.ind.dbc_idx).*xdc./g.ind.ind_pot(g.ind.dbc_idx,4);%Xsp
                g.ind.ind_pot(g.ind.dbc_idx,6) = g.ind.ind_pot(g.ind.dbc_idx,3)-g.ind.ind_pot(g.ind.dbc_idx,5);
                g.ind.ind_pot(g.ind.dbc_idx,7)=g.sys.basrad*rdc./g.ind.ind_pot(g.ind.dbc_idx,4); %1/Tr
            end
            if ~isempty(g.ind.db_idx)
                % reset deepbar
                [rdb,xdb]=deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx),g.ind.slip(g.ind.db_idx,k));
                g.ind.ind_pot(g.ind.db_idx,4) = Xm(g.ind.db_idx)+xr(g.ind.db_idx)+xdb;
                g.ind.ind_pot(g.ind.db_idx,5) = xs(g.ind.db_idx)+Xm(g.ind.db_idx).*(xr(g.ind.db_idx)+xdb)./g.ind.ind_pot(g.ind.db_idx,4);%Xsp
                g.ind.ind_pot(g.ind.db_idx,6) = g.ind.ind_pot(g.ind.db_idx,3)-g.ind.ind_pot(g.ind.db_idx,5);
                g.ind.ind_pot(g.ind.db_idx,7) = g.sys.basrad*rdb./g.ind.ind_pot(g.ind.db_idx,4); %1/Tr
            end
            %Brereton, Lewis and Young motor model
            g.ind.dvdp(:,k)=-(iqm.*g.ind.ind_pot(:,6)+g.ind.vdp(:,k)).*...
                g.ind.ind_pot(:,7)+g.ind.vqp(:,k).*g.ind.slip(:,k)*g.sys.basrad;
            g.ind.dvqp(:,k)=(idm.*g.ind.ind_pot(:,6)-g.ind.vqp(:,k)).*...
                g.ind.ind_pot(:,7)-g.ind.vdp(:,k).*g.ind.slip(:,k)*g.sys.basrad;
            g.ind.t_mot(:,k) = g.ind.vdp(:,k).*idm+g.ind.vqp(:,k).*iqm;
            g.ind.dslip(:,k)=(g.ind.tload(:,k)-g.ind.t_mot(:,k))/2./g.ind.ind_con(:,9);
        else
            error('motor by motor mode is not supported')  
        end
    end
    if flag == 3
        %linearize
        %add code later
    end
end