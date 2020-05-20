%% Simulation Runner

% Matt Stajcar

clear  
close all; 
clc;

tic

%% Load Network Parameters.
[Bus,Branch,Gen] = funBusBranchGenInputs();

for ii = 1:length(Branch)
    Branch(ii).Y = 1/Branch(ii).Zpu;
end
clear ii
for ii = 1:length(Bus)
    Bus(ii).AngR = Bus(ii).AngD * pi/180;
end
clear ii

%% Pre-disturbance Power Flow.

% The variable type of bus, branch, and gen are structures. Struct2Array is
% a function that converts the variables bus, branch, and gen from
% type struct to type array. As type array variables they can be inputs to
% PST's load flow solver.
n = length(Bus);
m = length(Branch);
[PSTbus,PSTbranch] = funStruct2Array(Bus,Branch,Gen,n,m);

% Input parameters for PST's loadflow solver.
tol = 1e-9;     % tolerance for convergence.
iter_max = 10;  % maximum number of iterations.
acc = 1;        % acceleration factor.
display = 'n';  % 'y', generate load-flow study report else, no load-flow 
                %  study report.
flag = 2;       % 1 - form new Jacobian every iteration; 2 - form new 
                % Jacobian every other iteration.
addpath('D:\Users\mstajcar\Documents\2014-2015\Summer\Simulation Project\pstV2')
[bus_sol,~,~] = funLoadFlow(PSTbus,PSTbranch,tol,iter_max,acc,display,flag);
rmpath('D:\Users\mstajcar\Documents\2014-2015\Summer\Simulation Project\pstV2')

% The outputs of PSTs load flow are of type array. The function
% Array2Struct converts the results from type array to type structure.
[Bus] = funArray2Struct(bus_sol,Bus);
clear flag acc bus_sol display iter_max m n PSTbranch PSTbus tol

ii = 1; % Index for simulation.

% Make additional calculations from solved load flow.
for jj = 1:length(Bus)
    Bus(jj).AngR = Bus(jj).AngD*pi/180;
    Bus(jj).V(ii,1) = Bus(jj).Vmag * (cos(Bus(jj).AngR) + 1j*sin(Bus(jj).AngR));
    Bus(jj).Sgen = Bus(jj).Pgen + 1j*Bus(jj).Qgen;
    Bus(jj).Iinj = conj(Bus(jj).Sgen/Bus(jj).V);
    Bus(jj).Sload = Bus(jj).Pload + 1j*Bus(jj).Qload;
    Bus(jj).Yload = conj(Bus(jj).Sload)/Bus(jj).Vmag^2; 
end
clear jj
for jj = 1:length(Gen)
    for kk = 1:length(Bus)
        if Bus(kk).Num == Gen(jj).BusNum
            Gen(jj).Vt(ii,1) = Bus(kk).V(ii,1);
            Gen(jj).I(ii,1) = Bus(kk).Iinj;
        end
    end
    if strcmpi(Gen(jj).Type,'Classical')
    	Gen(jj).Y = 1/(1j*Gen(jj).Xd_p);
    elseif strcmpi(Gen(jj).Type,'SubTrans_Genrou')
        Gen(jj).Y = 1/(Gen(jj).R + 1j*Gen(jj).Xd_pp);
    elseif strcmpi(Gen(jj).Type,'SubTrans_GenTpJ')
        Gen(jj).Y = 1/(Gen(jj).R + 1j*Gen(jj).Xd_pp);
    end
    Gen(jj).Wr(ii,1) = 1;
    Gen(jj).Pe(ii,1) = real(Gen(jj).Vt(ii,1) * conj(Gen(jj).I(ii,1)));
    Gen(jj).Pm(ii,1) = Gen(jj).Pe(ii,1) + abs(Gen(jj).I(ii,1))^2 * Gen(jj).R;    
end
clear jj kk

%% Build Y bus(es).

[Ybus] = funYbusCalc(Bus,Branch,Gen);

[Ybus] = funYredYrecovCalc(Ybus);

%% Initialize.

[Switch] = funSwitching();

[Time] = funTimeParameters();

InitializeSubTrans = 1;

%% Start Simulation.

Flag = 'subtransient';

while Time.SimTime(ii) < Time.SimEnd
    
    % The function, roundn(), has to be used in the following 2 if
    % statements due to the fact that 1/600 will eventually lead to a
    % rounding error in which Time.SimTime will not be an exact number of
    % integer values of 1/600.
    
    % Check to see if a switching case has occured.
    if roundn(Time.SimTime(ii,1),-6) == roundn(Switch.TstartSwitch(1),-6)
        
        % Based on event type, make appropriate changes to system
        % topology.
        [Ybus] = funModifyYbus4Switch(Switch,Ybus);    
    end
    
    % Check to see if fault has cleared. If so, make appropriate 
    % adjustments.
    if roundn(Time.SimTime(ii,1),-6) == roundn(Switch.TclearNear,-6)
        [Ybus,Branch] = funPostFaultYbusCalcs(Bus,Branch,Gen,Ybus,Switch);
    end
    
    % Select appropriate model based on variable, Flag.
    switch Flag
        
        case 'subtransient'
            if InitializeSubTrans == 1
                
                % Calc RotAng, Efd, Psikd, Psikq, Edp, and Eqp @ ii.
                [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,0);
                % Calc E @ ii.
                [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,1);
                
                %
                [Gen,Bus] = funInitialize4AdamsBashGenTpJ(Gen,Bus,ii);
                % Increment time and step ii in order to be able to use Adams Bashforth
                % integration method.
                Time.SimTime(ii+1,1) = Time.SimTime(ii,1) + Time.dtSubTrans;
                ii = ii + 1;

                % Calc RotAng, Efd, Psikd, Psikq, Edp, and Eqp @ ii.
                [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,0);
                % Calc E @ ii.
                [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,1);
                
                InitializeSubTrans = 0;
            end
            
            % Calc E @ ii.
            [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,1);
            
            [Gen,Bus] = funCalculate_V_I_Pe(Gen,Bus,Ybus,ii);
            
            % Check is simulation is finished. Otherwise, update.
            if Time.SimEnd > (Time.SimTime(ii) + Time.dtSubTrans)

                % Integrate Edp, Eqp, Psikd, Psikq, Wr, and RotAng.
                [Gen] = funSynMacSub_GenTpJ(Gen,ii,Time.dtSubTrans,2);

                % Update Pm  and Efd --> delete after turbine model  and exciter is added.
                for jj = 1:length(Gen)
                    Gen(jj).Pm(ii+1,1) = Gen(jj).Pm(ii,1);
                    Gen(jj).Efd(ii+1,1) = Gen(jj).Efd(ii,1);
                end

                Time.SimTime(ii+1,1) = Time.SimTime(ii,1) + Time.dtSubTrans;
            else
                break
            end
            
        case 'classical'
            
            for dumby = 1:1
            
            [Gen,Bus] = funCalculate_V_E_I_Pe(Gen,Bus,Ybus,ii);
            
            % Check is simulation is finished. Otherwise, update.
            if Time.SimEnd > (Time.SimTime(ii) + Time.dtClassical)
                
                [Gen] = funIntegrateAdamsBashClassical(Gen,Time,ii);

                % Update internal machine voltage, E.
                for jj = 1:length(Gen)
                    Gen(jj).E(ii+1,1) = abs(Gen(jj).E(ii,1)) * exp(1j * Gen(jj).RotAng(ii+1,1));
                end
                clear jj

                % Update Pm --> delete after turbine model is added.
                Gen(1).Pm(ii+1,1) = Gen(1).Pm(ii,1);
                Gen(2).Pm(ii+1,1) = Gen(2).Pm(ii,1);

                Time.SimTime(ii+1,1) = Time.SimTime(ii,1) + Time.dtClassical;
            else
                break
            end
            
            end
            
        case 'slow'
            disp('slow dynamics')
            
            Time.SimTime(ii+1,1) = Time.SimTime(ii,1) + Time.dtSlow;
    end

    ii = ii + 1;

end
clear ii

%% Results

% Plot bus voltages.
figure
subplot(221)
plot(Time.SimTime,abs(Bus(1).V),'b')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Bus 1')

subplot(222)
plot(Time.SimTime,abs(Bus(2).V),'b')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Bus 2')

subplot(223)
plot(Time.SimTime,abs(Bus(3).V),'b')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Bus 3')

subplot(224)
plot(Time.SimTime,abs(Bus(4).V),'b')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Bus 4')

% Plot machine internal voltage and terminal voltage.
figure
subplot(211)
plot(Time.SimTime,abs(Gen(1).Vt),'b',Time.SimTime,abs(Gen(1).E),'r')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Gen @ B4')
legend('Vt','E','location','southeast')

subplot(212)
plot(Time.SimTime,abs(Gen(2).Vt),'b',Time.SimTime,abs(Gen(1).E),'r')
xlabel('time (s)')
ylabel('voltage (pu)')
title('Gen @ B1')
legend('Vt','E','location','southeast')

% Plot machine speeds.
figure
plot(Time.SimTime,60.*Gen(1).Wr,'b',Time.SimTime,60.*Gen(2).Wr,'r')
xlabel('Time (s)')
ylabel('Speed (pu)')
title('Machine Speeds')
legend('Gen @ B4','Gen @ B1')
grid on

% Plot machine angles.
figure
plot(Time.SimTime,180/pi*Gen(1).RotAng,'b',Time.SimTime,180/pi*Gen(2).RotAng,'r')
xlabel('Time (s)')
ylabel('Angle (degrees)')
title('Rotor Angle')
legend('Gen @ B4','Gen @ B1','location','east')

% Plot electrical power.
figure
plot(Time.SimTime,Gen(1).Pe,'b',Time.SimTime,Gen(2).Pe,'r')
xlabel('Time (s)')
ylabel('Power (pu)')
title('Electrical Power')
legend('Gen @ B4','Gen @ B1','location','east')

% Plot generator current.
figure
plot(Time.SimTime,abs(Gen(1).I),'b',Time.SimTime,-abs(Gen(2).I),'r')
xlabel('Time (s)')
ylabel('current (pu)')
title('Generator Current')
legend('Gen @ B4','Gen @ B1')

% Play sounds when finished.
beep

toc

Gen1E_PSTN2 = Gen(1).E;
Gen2E_PSTN2 = Gen(2).E;
Gen1Wr_PSTN2 = Gen(1).Wr;
Gen2Wr_PSTN2 = Gen(2).Wr;
Gen1RotAng_PSTN2 = Gen(1).RotAng;
Gen2RotAng_PSTN2 = Gen(2).RotAng;

t = Time.SimTime;

save PSTN2 Gen1E_PSTN2 Gen2E_PSTN2 Gen1Wr_PSTN2 Gen2Wr_PSTN2 Gen1RotAng_PSTN2 Gen2RotAng_PSTN2 t

