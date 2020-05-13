function [bus_new]=mac_ind(i,k,bus,flag)
%syntax [bus_new]=mac_ind(i,k,bus,flag)
%bus_new is bus with the power and reactive power loads 
%modified to subtract the motor loads
%modification is made only when the motors are initialized
%i.e. flag = 0
%i is the motor number, 0 for vectorized computation
%k is the time step
%bus is the solved load flow bus data
%flag is 0 for initialization
%        1 for network interface
%        2 for dertermination of rates of change of states
%        3 for formation of linearized state matrix
%Purpose
%Simple Induction Motor Model
%Single Cage
%no leakage inductance saturation
%data format ind_con
%1 - motor number
%2 - busnumber
%3 - base MVA
%4 - rs
%5 - xs -stator leakage reactance
%6 - Xm - magnetizing reactance
%7 - rr
%8 - xr - rotor leakage reactance
%9 - H  - inertia constant motor + load in sec
%15 - fraction of bus load power taken by motor 
% if entry 15 is zero, it is assumed that the motor is to be started on the specified bus
% 
%Author Graham Rogers
%Date November 1995
%
global basmva basrad 
global  tload t_init pmot qmot vdmot vqmot  idmot iqmot ind_con ind_pot
global bus_int ind_int motbus
global vdp vqp slip dvdp dvqp dslip
jay=sqrt(-1);
bus_new=bus;

if ~isempty(ind_con)
  if flag == 0;
  % initialisation
	if i == 0;
		%vector computation
            	motnum=length(ind_con(:,1));
		motbus=bus_int(ind_con(:,2));
		ind_pot(:,1)=basmva./ind_con(:,3); %scaled mva base
		ind_pot(:,2)=ones(motnum,1); %base kv
		mot_vm(:,1)=bus(motbus,2); %motor terminal voltage mag
		mot_ang(:,1)=bus(motbus,3)*pi/180; %motor term voltage angle
		v=mot_vm(:,1).*exp(jay*mot_ang(:,1));
		vdmot(:,1)=real(v);
		vqmot(:,1)=imag(v);
		pmot(:,1)=bus(motbus,6).*ind_con(:,15);%motor power demand
                %modify bus load power
                bus_new(motbus,6)=bus(motbus,6)-pmot(:,1);
	        % index of motors to be initialized for running 
               	run_ind = find(ind_con(:,15)~=0);
               	%assumes motor starting if power fraction zero
		ind_pot(:,3)=ind_con(:,5)+ind_con(:,6);%Xs
		ind_pot(:,4)=ind_con(:,8)+ind_con(:,6);%Xr
	        ind_pot(:,5)=ind_con(:,5)+ind_con(:,6).*...
                        ind_con(:,8)./ind_pot(:,4);%Xsp
		ind_pot(:,6)=ind_pot(:,3)-ind_pot(:,5);%(Xs-Xsp)
		ind_pot(:,7)=basrad*ind_con(:,7)./ind_pot(:,4); %1/Tr
		rs=ind_con(:,4);
		xs=ind_con(:,5);
		Xm=ind_con(:,6);
		rr=ind_con(:,7);
		xr=ind_con(:,8);
		% find initial slip
		slip_old=zeros(motnum,1);
		slip_new=ones(motnum,1);
		% Set defaults for motor starting
                imot=zeros(motnum,1);
              	pem = zeros(motnum,1);
                qem = zeros(motnum,1);
                vdp(:,1)=zeros(motnum,1);
                vqp(:,1)=zeros(motnum,1);
                t_init = ones(motnum,1); % default for motor starting
		%Newton-Raphson iteration to determine initial slip for
		%running motors
                motrun=length(run_ind);%number of running motors
               	if motrun~=0 %check that some motors are running
                  rsrun=ind_con(run_ind,4);
		  xsrun=ind_con(run_ind,5);
		  Xmrun=ind_con(run_ind,6);
		  rrrun=ind_con(run_ind,7);
		  xrrun=ind_con(run_ind,8);
                  iter = 0;
                  err=max(abs(slip_new-slip_old));
                  while err>=1e-8 & iter<30
                    iter=iter+1;
		    y=basrad.*slip_old(run_ind)./ind_pot(run_ind,7);
                    denom = ones(motrun,1)+y.*y;
		    zr=rsrun+y.*ind_pot(run_ind,6)./denom;
		    zi=ind_pot(run_ind,5)+ind_pot(run_ind,6)./denom;
		    dzr=ind_pot(run_ind,6).*(ones(motrun,1)-...
				y.*y)./denom./denom;
		    dzi=-2*ind_pot(run_ind,6).*y./denom./denom;
                    zmod2=zr.*zr+zi.*zi;
		    dp=v(run_ind).*conj(v(run_ind)).*(dzr.*zmod2-...
                                 2*zr.*(dzr.*zr+dzi.*zi));
		    dp=dp./zmod2./zmod2;
		    pem(run_ind)=v(run_ind).*conj(v(run_ind)).*zr./zmod2;
		    ynew=y-(pem(run_ind)- ...
                             pmot(run_ind).*ind_pot(run_ind,1))./dp;
                    slip_new(run_ind)=ynew.*ind_pot(run_ind,7)/basrad;
                    err = max(abs(slip_new-slip_old));
                    slip_old=slip_new;
		 end
                 if iter >=30
                       disp('slip calculation failed to converge')
                       return
                 end
	      end
	      slip(:,1)=slip_new;
              y=basrad*slip(:,1)./ind_pot(:,7);
              denom= ones(motnum,1)+y.*y;
	      zr=rs+y.*ind_pot(:,6)./denom;
	      zi=ind_pot(:,5)+ind_pot(:,6)./denom;
              imot(run_ind)=v(run_ind)./(zr(run_ind)+jay*zi(run_ind));
	      sm(run_ind)=v(run_ind).*conj(imot(run_ind));
              pem(run_ind)=real(sm(run_ind));
	      qem(run_ind)=imag(sm(run_ind));
              %complex initial rotor states
	      vp = v - (rs+ jay* ind_pot(:,5)).*imot; 
   	      vdp(run_ind,1)=real(vp(run_ind));
              vqp(run_ind,1)=imag(vp(run_ind));
              % initial motor torque
              t_init(run_ind) = real(vp(run_ind).*conj(imot(run_ind)));
              f=ind_ldto(0,1);
       	      lzero_ind=find(tload(:,1)==0);
              if length(lzero_ind)~=0
                 tload(lzero_ind,1) = ones(length(lzero_ind),1);
              end
              t_init(run_ind) = t_init(run_ind)./tload(run_ind,1);
	      idmot(:,1)=real(imot)./ind_pot(:,1);
	      iqmot(:,1)=imag(imot)./ind_pot(:,1);
	      % modify qload 
	      bus_new(motbus,7)=bus(motbus,7)-qem./ind_pot(:,1);
	else
 	      % motor by motor initialization
	      nmot=size(ind_con,1);
	      motbus=bus_int(ind_con(i,2));
	      ind_pot(i,1)=basmva/ind_con(i,3); %scaled mva base
	      ind_pot(i,2)=1.; %base kv
	      mot_vm(i,1)=bus(motbus,2); %motor terminal voltage mag
	      mot_ang(i,1)=bus(motbus,3)*pi/180.; %motor term voltage angle
	      v=mot_vm(i,1)*exp(jay*mot_ang(i,1));
	      vdmot(i,1)=real(v);
	      vqmot(i,1)=imag(v);
	      pmot(i,1)=bus(motbus,6)*ind_con(i,15);%motor power demand
              %modify motor bus load
              bus_new(motbus,6)=bus(motbus,6)-pmot(i,1); 
	      ind_pot(i,3)=ind_con(i,5)+ind_con(i,6);%Xs
	      ind_pot(i,4)=ind_con(i,8)+ind_con(i,6);%Xr
	      ind_pot(i,5)=ind_con(i,5)+ind_con(i,6)*...
                           ind_con(i,8)/ind_pot(i,4);%Xsp
	      ind_pot(i,6)=ind_pot(i,3)-ind_pot(i,5);%(Xs-Xsp)
	      ind_pot(i,7)=basrad*ind_con(i,7)/ind_pot(i,4); %1/Tr
	      rs=ind_con(i,4);
	      xs=ind_con(i,5);
	      Xm=ind_con(i,6);
	      rr=ind_con(i,7);
	      xr=ind_con(i,8);
	      % find initial slip
	      slip_old=0.005;
	      slip_new=1;
              slip(i,1)=1;
              vdp(i,1)=0;
              vqp(i,1)=0;
              t_init(i)=1;
              idmot(i,1)=0.;
              iqmot(i,1)=0.;
              err = abs(slip_old-slip_new);
              if  ind_con(i,15)>1e-6 % motor running
		%Newton-Raphson iteration to find initial slip
		iter=0;
		while err>=1e-8 & iter<30
                   iter=iter+1;
		   y=basrad*slip_old/ind_pot(i,7); %s omega Tr
                   denom = 1+y*y;
		   zr=rs+y.*ind_pot(i,6)/denom;
		   zi=ind_pot(i,5)+ind_pot(i,6)/denom;
		   dzr=ind_pot(i,6)*(1-y*y)/denom/denom;
		   dzi=-2*ind_pot(i,6)*y/denom/denom;
                   modz2=zr*zr+zi*zi;
		   dp=v*conj(v)*(dzr*modz2 -...
			 2*zr*(dzr*zr+dzi*zi));
		   dp=dp/modz2/modz2;
		   pem=v*conj(v)*zr/modz2;
		   ynew=y-(pem-pmot*ind_pot(i,1))/dp;
                   slip_new=ynew/basrad*ind_pot(i,7);
                   err=abs(slip_old-slip_new);
                   slip_old=slip_new;
		end;
		if iter>=30
		   disp('slip calculation failed to converge')
                   return
                end
	        slip(i,1)=slip_new;
                y=basrad*slip(i,1)/ind_pot(i,7);
                denom= 1+y*y;
	        zr=rs+y*ind_pot(i,6)/denom;
	        zi=ind_pot(i,5) + ind_pot(i,6)/denom;
		imot=v/(zr+jay*zi); %compex motor current
	        smot=v*conj(imot);
		pem=real(smot);
	        qem=imag(smot);
	        vp = v-(rs+jay*ind_pot(i,5))*imot; %complex initial rotor states
                vdp(i,1)=real(vp);
                vqp(i,1)=imag(vp);
                t_init(i)=real(vp*conj(imot));
                f=ind_ldto(i,1);% find the motor load torque
                if tload(i,1)==0
                    error('you must define a load speed characteristic')
                end
                t_init(i)=t_init(i)/tload(i,1);%set the load multiplyer
		%motor currents on system base
	        idmot(i,1)=real(imot)/ind_pot(i,1);
	        iqmot(i,1)=imag(imot)/ind_pot(i,1);
	        % modify qload 
		bus_new(motbus,7)=bus(motbus,7)-qem/ind_pot(i,1);
	     end
	end
  end
  if flag == 1
  %network interface
  %no interface required for induction motor
  end
  if flag == 2
  %motor dynamics calculation
	if i == 0
	%vector calculation
  
             f=ind_ldto(0,k);
		idm=idmot(:,k).*ind_pot(:,1);%convert to machine base
		iqm=iqmot(:,k).*ind_pot(:,1);
		%Brereton, Lewis and Young motor model
		dvdp(:,k)=-(iqm.*ind_pot(:,6)+vdp(:,k)).*...
			ind_pot(:,7)+vqp(:,k).*slip(:,k)*basrad;
		dvqp(:,k)=(idm.*ind_pot(:,6)-vqp(:,k)).*...
			ind_pot(:,7)-vdp(:,k).*slip(:,k)*basrad;	
		dslip(:,k)=(tload(:,k).*t_init(:)-vdp(:,k).*...
			idm-vqp(:,k).*iqm)/2./ind_con(:,9);
	else
		 f=ind_ldto(i,k);
		 idm=idmot(i,k)*ind_pot(i,1);
		 iqm=iqmot(i,k)*ind_pot(i,1);
		 dvdp(i,k)=-(iqm*ind_pot(i,6)+vdp(i,k))*...
			ind_pot(i,7)+vqp(i,k)*slip(i,k)*basrad;
		 dvqp(i,k)=(idm*ind_pot(i,6)-vqp(i,k))*...
			ind_pot(i,7)-vdp(i,k)*slip(i,k)*basrad;	
		 dslip(i,k)=(tload(i,k)*t_init(i)-vdp(i,k)*...
			idm-vqp(i,k)*iqm)/2/ind_con(i,9);
	end
  end
  if flag == 3
  %linearize
  %add code later
  end
end