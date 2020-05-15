% Two-machine simulation example.
% Uses both classical and sub transient models
% Includes govs
% Compare with PST
% This one does the swap by having one swing equation, getting rid of the
% gen model and making calculating electrical power from
% Pe=Pmech-2*H*w*wdot for each indiviual generator.
% Oct 2015.
%
%  Bus1        Bus3          Bus2
%   |---xline1--|----xline2---|
%  Gen1        Load          Gen2
%
%
% SEE TimeSequenceBlockDiagram3.jpg
clear all; close all; clc

%% Settings
wbase = 2*pi*60;
Sbase = 100;
xLine = [0.65; 0.35]; %Line reactances
Pload = 1.6;

%Gen parameters
Gen(1).param.type = 'Sub_PSTN2';%'Classical';%
Gen(1).param.bus = 1;
Gen(1).param.xl = 0.15;
Gen(1).param.ra = 0*0.001;
Gen(1).param.xd = 1.7;
Gen(1).param.xdp = 0.245;
Gen(1).param.xdpp = 0.185;
Gen(1).param.Td0p = 5.9;
Gen(1).param.Td0pp = 0.03;
Gen(1).param.xq = 1.64;
Gen(1).param.xqp = 1.64;
Gen(1).param.xqpp = 0.185;
Gen(1).param.Tq0p = 0.082;
Gen(1).param.Tq0pp = 0.082;
Gen(1).param.H = 2.37;
if strcmpi(Gen(1).param.type,'Classical')
    Gen(1).param.D = 20; %Classsicals don't have damper coils modelled
else
    Gen(1).param.D = 0;
end

Gen(2).param.type = 'Sub_PSTN2';%'Classical';%
Gen(2).param.bus = 2;
Gen(2).param.xl = 0.15;
Gen(2).param.ra = 0*0.001;
Gen(2).param.xd = 1.7;
Gen(2).param.xdp = 0.245;
Gen(2).param.xdpp = 0.185;
Gen(2).param.Td0p = 5.9;
Gen(2).param.Td0pp = 0.03;
Gen(2).param.xq = 1.64;
Gen(2).param.xqp = 1.64;
Gen(2).param.xqpp = 0.185;
Gen(2).param.Tq0p = 0.082;
Gen(2).param.Tq0pp = 0.082;
Gen(2).param.H = 5;
if strcmpi(Gen(2).param.type,'Classical')
    Gen(2).param.D = 20; %Classsicals don't have damper coils modelled
else
    Gen(2).param.D = 0;
end

%System generator
GenSys.param.H = 0;
GenSys.param.D = 0;
for k=1:length(Gen)
    GenSys.param.H = GenSys.param.H + Gen(k).param.H;
    GenSys.param.D = GenSys.param.D + Gen(k).param.D;
end

%gov parameters
Gov(1).param.type = 'tgov1';
Gov(1).param.nGen = 1;
Gov(1).param.R = 1/20; %Droop
Gov(1).param.T1 = 0.5; %servo valve/Steam bowl time const
Gov(1).param.T2 = 3;
Gov(1).param.T3 = 15; %Tw for hydro, reheater for steam
Gov(1).param.Dt = 0;
Gov(1).param.Vmax = 1.2;
Gov(1).param.Vmin = 0;
Gov(2).param.type = 'tgov1';
Gov(2).param.nGen = 2;
Gov(2).param.R = 1/20;
Gov(2).param.T1 = 0.5;
Gov(2).param.T2 = 3;
Gov(2).param.T3 = 10;
Gov(2).param.Dt = 0;
Gov(2).param.Vmax = 1.0;
Gov(2).param.Vmin = 0;

%Exciter parameters
Exc(1).param.type = 'PST0';
Exc(1).param.nGen = 1;
Exc(1).param.TA = 0.05;
Exc(1).param.KA = 200;
Exc(1).param.TB = 0;
Exc(1).param.TC = 0;
Exc(1).param.Efdmax = 7;
Exc(1).param.Efdmin = -6;
Exc(2).param.type = 'PST0';
Exc(2).param.nGen = 2;
Exc(2).param.TA = 0.05;
Exc(2).param.KA = 200;
Exc(2).param.TB = 0;
Exc(2).param.TC = 0;
Exc(2).param.Efdmax = 7;
Exc(2).param.Efdmin = -6;

%Timing
Time.FaultStart = 0.5; %fault time
Time.FaultEnd = Time.FaultStart + 1/60; %fault time
Time.dtF = 1/600; %Fast-dyamic simulation time step
Time.dtS = 1; %Slow-dynamic simulation time step
Time.End = 30;

%Model Swapping
ModelSwap.Flag = 1; %If true, conduct model swapping
ModelSwap.lim.delw = 1e-4; %Switch to slow dynamics if relative w is below this for on sec
ModelSwap.lim.delV = 1e-3; %PU limit for voltages
ModelSwap.lim.PslackTol = 0.001;
ModelSwap.nGen = 1; %Reference generator (must switch if that gen is tripped)

%% Run PST for comparison
% addpath('C:\Users\dtrudnowski\Documents\2013_14\EELE5550\LabHandouts\pstV2')
% bus = [ ...
% % num volt angle p_gen     q_gen p_load  q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
%   1   1.1  0.0   0.6*Pload 0     0       0      0       0       2    100   -100  20      1.5   0.5;
%   2   1.0  0.0   0.4*Pload 0     0       0      0       0       1    100   -100  20      1.5   0.5;
%   3   1.0  0.0   0         0     Pload/2 0      0       0       3    0     0     20      1.5   0.5;
%   4   1.0  0.0   0         0     Pload/2 0      0       0       3    0     0     20      1.5   0.5];
% 
% line = [ ...
% % bus bus r    x         y    tapratio tapphase tapmax tapmin tapsize
%   1   3   0    xLine(1)        0.0  1.0      0        0      0      0;
%   2   3   0    xLine(2)        0.0  1.0      0        0      0      0;
%   3   4   0    0.01*min(xLine) 0.0  1.0      0        0      0      0;];
% mac_con = zeros(length(Gen),19);
% for k=1:length(Gen)
%     if strcmpi(Gen(k).param.type,'Classical')
%         g = Gen(k).param;
%         mac_con(k,:) = [...
%             % num bus   base  x_l  r_a  x_d x'_d  x"_d  T'_do T"_do x_q  x'_q  x"_q  T'_qo T"_qo H   d_0  d_1  bus
%               k   g.bus 100   0    0    0   g.xdp 0     0     0     0    0     0     0     0     g.H g.D  0    g.bus];
%     elseif strcmpi(Gen(k).param.type,'Sub_PSTN2')
%         g = Gen(k).param;
%         mac_con(k,:) = [...
%             % num bus   base  x_l  r_a  x_d  x'_d  x"_d   T'_do  T"_do   x_q  x'_q  x"_q   T'_qo  T"_qo   H   d_0 d_1  bus
%               k   g.bus 100   g.xl g.ra g.xd g.xdp g.xdpp g.Td0p g.Td0pp g.xq g.xqp g.xqpp g.Tq0p g.Tq0pp g.H g.D 0    g.bus];
%     else
%         error(' ')
%     end
% end
% 
% if exist('Gov','var')
%     tg_con = zeros(length(Gov),10);
%     for k=1:length(Gov)
%         g = Gov(k).param;
%         tg_con(k,:) = [...
%             % type mach   wf 1/R   Tmax   Ts   Tc   T3   T4     T5
%               1    g.nGen 1  1/g.R g.Vmax g.T1 0.01 0    g.T2   g.T3]; %Must have non-zero time constants for T's
%     end
% end
% 
% if exist('Exc','var')
%     exc_con = zeros(length(Exc),9);
%     for k=1:length(Gov)
%         g = Exc(k).param;
%         exc_con(k,:) = [...
%             %  type mach   T_R  K_A T_A  T_B  T_C  Vmax     Vmin
%                 0   g.nGen 0   g.KA g.TA g.TB g.TC g.Efdmax g.Efdmin];
%     end
% end
% 
% sw_con = [...
% 0               0    0    0    0    0    1/600;%sets intitial time step
% Time.FaultStart 3    0    0    0    5    1/600; %loss of load
% Time.FaultEnd   0    0    0    0    0    1/600; %clear near end of fault
% Time.End              0    0    0    0    0    1/600]; % end simulation
% if exist('Gov','var') && exist('Exc','var'); save delme Time Gov Gen xLine wbase Pload Exc
% elseif exist('Gov','var'); save delme Time Gov Gen xLine wbase Pload
% elseif exist('Exc','var'); save delme Time Exc Gen xLine wbase Pload
% else save delme Time Gen xLine wbase Pload
% end
% if exist('tg_con','var') && exist('exc_con','var'); save DataFile bus line mac_con sw_con tg_con exc_con
% elseif exist('tg_con','var'); save DataFile bus line mac_con sw_con tg_con
% elseif exist('exc_con','var'); save DataFile bus line mac_con sw_con exc_con
% else save DataFile bus line mac_con sw_con
% end
% s_simu_Batch2
% clear Time Gov Gen xLine wbase Pload Exc
% save delmePST
% clear all; close all; clc
% load delme
% delete delme.mat
% delete DataFile.mat
% rmpath('C:\Users\dtrudnowski\Documents\2013_14\EELE5550\LabHandouts\pstV2')

%% Power flow
addpath('C:\Users\dtrudnowski\Documents\2013_14\EELE5550\LabHandouts\pstV2')
bus = [ ...
% num volt angle p_gen     q_gen p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.1  0.0   0.6*Pload 0     0      0      0       0       2    100   -100  20      1.5   0.5;
  2   1.0  0.0   0.4*Pload 0     0      0      0       0       1    100   -100  20      1.5   0.5;
  3   1.0  0.0   0         0     Pload  0      0       0       3    0     0     20      1.5   0.5];
line = [ ...
% bus bus r    x         y    tapratio tapphase tapmax tapmin tapsize
  1   3   0    xLine(1)        0.0  1.0      0        0      0      0;
  2   3   0    xLine(2)        0.0  1.0      0        0      0      0];
Bus = loadflow(bus,line,1e-10,10,1,'n',1);
BusV0 = Bus(:,2).*exp(1i*Bus(:,3).*(pi/180));
Yload = (Bus(3,6) - 1i*Bus(3,7))/(Bus(3,2)^2); %Load admittance
Ig0 = (Bus(1:2,4)-1i*Bus(1:2,5))./conj(BusV0(1:2)); %Initial generator currents
rmpath('C:\Users\dtrudnowski\Documents\2013_14\EELE5550\LabHandouts\pstV2')
clear bus

%% Set up Admittance matrices
Ybus = -1i*[1/xLine(1), 0, -1/xLine(1);
            0, 1/xLine(2), -1/xLine(2);
            -1/xLine(1), -1/xLine(2), 1/xLine(1)+1/xLine(2)];
yload = [0;
         0;
         Yload];
ygen = zeros(3,1);
YbusGen = zeros(3,2);
Ygen = zeros(2,2);
for k=1:length(Gen)
    if strcmpi(Gen(k).param.type,'Classical')
        y = -1i/Gen(k).param.xdp;
    elseif strcmpi(Gen(k).param.type,'Sub_PSTN2')
        y = 1/(Gen(k).param.ra+1i*Gen(k).param.xdpp);
    end
    ygen(Gen(k).param.bus) = y;
    Ygen(k,k) = y;
    YbusGen(k,Gen(k).param.bus) = -y;
end

%Calculated admittances 
YbusHat = Ybus + diag(ygen) + diag(yload);
Yred = Ygen - transpose(YbusGen)*(YbusHat\YbusGen);
Yrecov = -(YbusHat\YbusGen);

%% Initialize variables
%Time variables
iTime = 1; %counter
nGens = length(Gen);
nBuses = size(Bus,1);
Time.End = round(Time.End/Time.dtF)*Time.dtF;
N = round(Time.End/Time.dtF);
tOut = [0;Time.dtF;NaN(N-2,1)];
ModelSwap.Type = 1; %If 1, use with fast dynamics, 2 = slow. (always start fast)

Bus_v = NaN(3,N); %Output bus voltages

%Gen variables
Eg0 = Yred\Ig0; %initial internal gen voltages on system frame
Hconst = zeros(nGens,1);
for k=1:nGens
    Hconst(k) = Gen(k).param.H;
    n = find(Gen(k).param.bus==Bus(:,1));
    Gen(k).VT = NaN(N,1); %Terminal voltage
    Gen(k).IT = NaN(N,1); %Output current
    Gen(k).VT(1:2) = (Bus(n,2)*exp((1i*pi/180)*Bus(n,3)))*[1;1];
    Gen(k).IT(1:2) = ((Bus(n,4) - 1i*Bus(n,5))/conj(Gen(k).VT(1)))*[1;1];
    Gen(k).w = [1;1;NaN(N-2,1)]; %PU speed
    Gen(k).wdot = [0;NaN(N-1,1)];
    Gen(k).delta = NaN(N,1);
    Gen(k).deltadot = [0;NaN(N-1,1)];
    p = (Bus(k,4) + Gen(k).param.ra*abs(Gen(k).IT(1))^2);
    Gen(k).Pmech = [p;p;NaN(N-2,1)];
    Gen(k).Pe = [Bus(k,4);NaN(N-1,1)];
    Gen(k).Te = [p;NaN(N-1,1)];
    Gen(k).Efd = NaN(N,1);
    Gen(k).E = NaN(N,1);
    if strcmpi(Gen(k).param.type,'Sub_PSTN2')
        Gen(k).eqp = NaN(N,1);
        Gen(k).eqpdot = [0;NaN(N,1)];
        Gen(k).edp = NaN(N,1);
        Gen(k).edpdot = [0;NaN(N-1,1)];
        Gen(k).psikq = NaN(N,1);
        Gen(k).psikqdot = [0;NaN(N-1,1)];
        Gen(k).psikd = NaN(N,1);
        Gen(k).psikddot = [0;NaN(N-1,1)];
        Gen(k) = funSynMacSub_PSTN2(Gen(k),[1:2],Time.dtF,wbase,0);
    elseif strcmpi(Gen(k).param.type,'Classical')
        Gen(k) = funSynMacClassical(Gen(k),[1:2],Time.dtF,wbase,0);
    else
        error('No gen model')
    end
end
GenSys.wdot = NaN(N,1);
GenSys.w = NaN(N,1);
GenSys.Pe = NaN(N,1);
GenSys.Pmech = NaN(N,1);
clear k p

%Gov variables
nGenGov = NaN(nGens,1);
if exist('Gov','var')
    for n=1:length(Gov)
        if strcmpi(Gov(n).param.type,'tgov1')
            Gov(n).x1dot = [0;NaN(N-1,1)];
            Gov(n).x1 = NaN(N,1);
            Gov(n).x2dot = [0;NaN(N-1,1)];
            Gov(n).x2 = NaN(N,1);
            Gov(n).Pref = NaN(N,1);
            Gov(n).Pref(1:2) = Gen(Gov(n).param.nGen).Pmech(1:2);
            Gov(n).Pmech = NaN(N,1);
            Gov(n).Pmech(1:2) = Gov(n).Pref(1:2);
            Gov(n).w = [1;1;NaN(N-2,1)];
            Gov(n) = funGov_TGOV1(Gov(n),[1;2],Time.dtF,0);
        else
            error(' ')
        end
        nGenGov(Gov(n).param.nGen) = n; %Cross indices for generators
    end
end

%Exc variables
nGenExc = NaN(nGens,1);
if exist('Exc','var')
    for n=1:length(Exc)
        if strcmpi(Exc(n).param.type,'PST0')
            Exc(n).x1dot = [0; NaN(N-1,1)];
            Exc(n).x1 = NaN(N,1);
            Exc(n).Efddot = [0; NaN(N-1,1)];
            Exc(n).Efd = [Gen(Exc(n).param.nGen).Efd(1)*[1;1];NaN(N-2,1)];
            Exc(n).Vref = NaN(N,1);
            Exc(n).VT = NaN(N,1);
            Exc(n).VT(1:2) = Gen(Exc(n).param.nGen).VT(1:2);
            Exc(n) = funExciter_PST0(Exc(n),[1;2],Time.dtF,0);
        else
            error(' ')
        end
        nGenExc(Exc(n).param.nGen) = n; %Cross indices for generators
    end
end
clear N

%Buses
nGenBus = NaN(nGens,1); %Generator bus numbers
for n=1:nGens
    nGenBus(n) = Gen(n).param.bus;
end
if sum(isnan(nGenBus))>1e-4; error(' '); end
if nGens<nBuses
    nNoGenBus = NaN(nBuses-nGens,1); %Non-gen bus numbers
    n = 0;
    for k=1:nBuses
        if min(abs(Bus(k,1)-nGenBus))>1e-4
            n = n + 1;
            nNoGenBus(n) = Bus(k,1);
        end
    end
    if sum(isnan(nNoGenBus))>1e-4; error(' '); end
    clear k n
end

%% Initialize some model swapping variables
if ModelSwap.Flag==1
    ModelSwap.SwitchFlag = 0; %Set true on the integration step used to switch
end

%% Simulate
FaultFlag = 1;
iTimeRef = 1;
iTime = 1;
while tOut(iTime)<=Time.End
    iTime = iTime + 1;
    if abs(tOut(iTime)/5-round(tOut(iTime)/5))<1e-8; disp(['t = ' num2str(tOut(iTime))]); end
    
    %% Check for model swapping
    if ModelSwap.Flag==1
        %Test if we switch to slow (NOTE: DON'T INCLUDE ANY GENERATOR THAT IS TRIPPED)
        if ModelSwap.Type==1
            if tOut(iTime)>Time.FaultEnd && (abs(tOut(iTime)-round(tOut(iTime)))/Time.dtF)<1e-8 %Only check on the sec.
                f = 1;
                n = [iTime:-1:iTime-round(1/Time.dtF)+1];
                for k=1:nGens
                    e = NaN(2,1);
                    e(1) = max(abs(Gen(ModelSwap.nGen).w(n)-Gen(k).w(n)));
                    e(2) = max(abs(diff(abs(Gen(k).VT(n)))));
                    if sum(isnan(e))>0.5; error(' '); end
                    if e(1)>ModelSwap.lim.delw || e(2)>ModelSwap.lim.delV
                        f = 0;
                        break
                    end
                end
                if f
                    disp(['Switching to slow dynamics at t = ' num2str(tOut(iTime))])
                    ModelSwap.SwitchFlag = 1;
                end
                clear f n e
            end
        end
    else
        %Test to switch back to fast
        %TO BE ADDED LATER
    end
    
    %% Initialize switching variables
    if ModelSwap.Flag==1
        if ModelSwap.SwitchFlag==1
            ModelSwap.tSwitch = tOut(iTime);
            ModelSwap.Type = 2; %Go to slow dynamics
            iTimeRef = iTime; %Index where switching occured.

            %Initialize system generator w and Pe
            GenSys.w(iTime) = 0;
            GenSys.Pe(iTime-1) = 0;
            for k=1:nGens 
                GenSys.w(iTime) = GenSys.w(iTime) + Gen(k).w(iTime);
                GenSys.Pe(iTime-1) = GenSys.Pe(iTime-1) + Gen(k).Pe(iTime-1);
            end
            GenSys.w(iTime) = GenSys.w(iTime)/nGens;
            clear k
            
            %Initialize reduced-order governor
            if exist('Gov','var')
                for n=1:length(Gov)
                    %Get gov inputs
                    Gov(n).Pref(iTime) = Gen(Gov(n).param.nGen).Pmech(1); %Input Pref never changes
                    Gov(n).w(iTime) = Gen(Gov(n).param.nGen).w(iTime);
                    Gov(n) = funGov_TGOV1(Gov(n),iTime,Time.dtF,1); %Get Pmech
                    Gov(n).w(iTime) = GenSys.w(iTime); %Switch to new speed
                    Gov(n) = funGov_TGOV1_SLOW(Gov(n),iTime,Time.dtF,1,0); %Initialize gov states
                end
                clear n
            end 
        end
    end
    
    %% Admittance matrices
    if tOut(iTime)>=Time.FaultStart && FaultFlag %Drop load
        FaultFlag = 0;
        yload(3) = 0.5*yload(3); %Have to leave a small load in order of Yred not to be illconditioned.  Don't know why?        
        YbusHat = Ybus + diag(ygen) + diag(yload);
        Yred = Ygen - transpose(YbusGen)*(YbusHat\YbusGen);
        Yrecov = -(YbusHat\YbusGen);
    end
    
    %Set up Yred2 = reduced Ybus for time sequencing (THIS IS WHERE YOU CAN PUT LOAD RAMPING)
    if ModelSwap.Type==2 && ModelSwap.SwitchFlag==1 %Build for the 1st time
       Y = Ybus + diag(yload);
       if nGens==nBuses
           Yred2 = Y(nGenBus,nGenBus);
       else
           n = [nGenBus;nNoGenBus];
           Y = Y(n,n);
           Y11 = Y(1:nGens,1:nGens);
           Y12 = Y(1:nGens,nGens+1:nBuses);
           Y21 = Y(nGens+1:nBuses,1:nGens);
           Y22 = Y(nGens+1:nBuses,nGens+1:nBuses);
           Yred2 = Y11 - Y12*(Y22\Y21);
           clear n Y11 Y12 Y21 Y22
       end
       clear Y    
    end
    
    %% Step time
    if ModelSwap.Type==1
        tOut(iTime+1) = tOut(iTime) + Time.dtF;
    else
        tOut(iTime+1) = tOut(iTime) + Time.dtS;
    end
    
    %% Calculate state eqn inputs and outputs
    %Gov inputs Pref and w; Gov output Pmech
    if exist('Gov','var')
        for n=1:length(Gov)
            %Get inputs
            Gov(n).Pref(iTime) = Gen(Gov(n).param.nGen).Pmech(1); %Input Pref never changes
            if ModelSwap.Type==1
                Gov(n).w(iTime) = Gen(Gov(n).param.nGen).w(iTime);
                Gov(n) = funGov_TGOV1(Gov(n),iTime,Time.dtF,1); %Updates Pmech
            elseif ModelSwap.Type==2
                Gov(n).w(iTime) = GenSys.w(iTime);
                Gov(n) = funGov_TGOV1_SLOW(Gov(n),iTime,Time.dtS,ModelSwap.SwitchFlag,1); %Updates Pmech
            end 
        end
        clear n
    end 
    
    %Detailed model:  inputs Efd, Pmech, and IT; outputs E (also update bus voltages, and gen Pe);
    %For reduced: Pmech, Pe.
    if ModelSwap.Type==1
        e = NaN(nGens,1);
        for k=1:nGens
            if ~isnan(nGenExc(k))
                Gen(k).Efd(iTime) = Exc(nGenExc(k)).Efd(iTime);
            else
                Gen(k).Efd(iTime) = Gen(k).Efd(1); %No exciter
            end
            if strcmpi(Gen(k).param.type,'Sub_PSTN2')
                Gen(k) = funSynMacSub_PSTN2(Gen(k),iTime,Time.dtF,wbase,1); %Get gen internal voltage E
            elseif strcmpi(Gen(k).param.type,'Classical')
                Gen(k) = funSynMacClassical(Gen(k),iTime,Time.dtF,wbase,1); %Get gen internal voltage E
            end
            e(k) = Gen(k).E(iTime);
            if ~isnan(nGenGov(k))
                Gen(k).Pmech(iTime) = Gov(nGenGov(k)).Pmech(iTime);
            else
                Gen(k).Pmech(iTime) = Gen(k).Pmech(1);
            end
        end
        Ig = Yred*e; %Network soln
        V = Yrecov*e;
        Bus_v(:,iTime) = V;
        for k=1:nGens %HOW DO I GET A VECTOR INTO A STRUCT?
            Gen(k).IT(iTime) = Ig(k);
            Gen(k).VT(iTime) = V(Gen(k).param.bus);
            Gen(k).Pe(iTime) = real(V(Gen(k).param.bus)*conj(Ig(k)));
        end
        clear e k Ig V
    else
        %Pmech, Pe
        Pe = zeros(nGens,1);
        GenSys.Pmech(iTime) = 0;
        GenSys.Pe(iTime) = 0;
        for k=1:nGens
            if ~isnan(nGenGov(k))
                Gen(k).Pmech(iTime) = Gov(nGenGov(k)).Pmech(iTime);
            else
                Gen(k).Pmech(iTime) = Gen(k).Pmech(1);
            end
            GenSys.Pmech(iTime) = GenSys.Pmech(iTime) + Gen(k).Pmech(iTime);
        end
        Pacc = GenSys.Pmech(iTime) - GenSys.Pe(iTime-1);
        for k=1:nGens
            Pe(k) = Gen(k).Pmech(iTime) - Hconst(k)*Pacc/GenSys.param.H;
        end
        v = [Gen.VT];
        a = angle(v(iTime-1,:))';
        v = abs(v(iTimeRef-1,:))'; %VOLTAGE MAGNITUDE IS CONSTANT
        [a,Pe] = funGenPowerFlow(Yred2,Pe,v,a,Hconst,ModelSwap.lim.PslackTol);
        for k=1:nGens
            Gen(k).Pe(iTime) = Pe(k);
            Gen(k).VT(iTime) = v(k)*exp(1i*a(k));
            GenSys.Pe(iTime) = GenSys.Pe(iTime) + Pe(k);
        end
        clear v a Pe k
    end
    
    %Exc inputs Vref and VT; Exc output Efd
    if ModelSwap.Type==1
        if exist('Exc','var')
            for n=1:length(Exc)
                %Get inputs
                Exc(n).Vref(iTime) = Exc(n).Vref(1); %Input Vref never changes
                Exc(n).VT(iTime) = Gen(Exc(n).param.nGen).VT(iTime);
            end
            clear n
        end 
    end
    
    %% Update Generator states
    if ModelSwap.Type==1 %Fast dynamics
        for k=1:nGens
            if strcmpi(Gen(k).param.type,'Sub_PSTN2')
                Gen(k) = funSynMacSub_PSTN2(Gen(k),iTime,Time.dtF,wbase,2);
            elseif strcmpi(Gen(k).param.type,'Classical')
                Gen(k) = funSynMacClassical(Gen(k),iTime,Time.dtF,wbase,2);
            end     
        end
    elseif ModelSwap.Type==2 %Slow dynamics
        GenSys.wdot(iTime) = ((GenSys.Pmech(iTime) - GenSys.Pe(iTime))/GenSys.w(iTime) ...
            - GenSys.param.D*(GenSys.w(iTime)-1))/(2*GenSys.param.H);
        if ModelSwap.SwitchFlag==1
            GenSys.w(iTime+1) = GenSys.w(iTime) + Time.dtS*GenSys.wdot(iTime); %Euler 1st time
        else
            GenSys.w(iTime+1) = GenSys.w(iTime) + Time.dtS*(1.5*GenSys.wdot(iTime)-0.5*GenSys.wdot(iTime-1));
        end
        for k=1:nGens
            Gen(k).wdot(iTime) = GenSys.wdot(iTime);
            Gen(k).w(iTime+1) = GenSys.w(iTime+1);
        end
    end
    clear k
    
    %% Governor states
    if exist('Gov','var')
        if ModelSwap.Type==1 %Fast dynamics
            for n=1:length(Gov)
                if strcmpi(Gov(n).param.type,'tgov1')
                    Gov(n) = funGov_TGOV1(Gov(n),iTime,Time.dtF,2);
                else
                    error('No gov model')
                end
            end
            clear n
        elseif ModelSwap.Type==2 %Slow dynamics
            for n=1:length(Gov)
                if strcmpi(Gov(n).param.type,'tgov1')
                    Gov(n) = funGov_TGOV1_SLOW(Gov(n),iTime,Time.dtS,ModelSwap.SwitchFlag,2);
                else
                    error('No gov model')
                end
            end
            clear n
        end

    end 
    
    %% Exciters
    if exist('Exc','var')
        if ModelSwap.Type==1
            for n=1:length(Exc)
                if strcmpi(Exc(n).param.type,'PST0')
                    Exc(n) = funExciter_PST0(Exc(n),iTime,Time.dtF,2);
                else
                    error('No exciter model')
                end
            end
            clear n
        end
    end 
    
    %% Reset ModelSwap.SwitchFlag
    if ModelSwap.Flag
        if ModelSwap.SwitchFlag==1
            ModelSwap.SwitchFlag = 0;
        end
    end
end

%% Truncate
tOut = tOut(1:iTime);
N = length(tOut);

%% Plot results
load delmePST t mac_spd pmech bus_v pelect bus_v

figure
plot(t,mac_spd(1,:)-mac_spd(2,:),'k','LineWidth',2)
hold on
plot(tOut,Gen(1).w(1:N)-Gen(2).w(1:N),'r.','MarkerSize',10);
hold off
legend('PST','Dan')
y=ylim;
hold on
plot(tOut,ModelSwap.lim.delw*ones(N,1),'b',...
     tOut,-ModelSwap.lim.delw*ones(N,1),'b',...
     tOut(iTimeRef-1)*[1;1],y,'b','LineWidth',2,'MarkerSize',10);
hold off
ylabel('Gen 1 - Gen 2 spd (pu)')
xlabel('Time (sec.)')

figure
for k=1:2
    subplot(2,1,k)
    plot(t,mac_spd(k,:),'k','LineWidth',2)
    hold on
    plot(tOut,Gen(k).w(1:N),'r.','MarkerSize',10);
    hold off
    if k==1; legend('PST','Dan'); end
    ylabel(['Gen ' num2str(k) ' spd (pu)'])
    xlabel('Time (sec.)')
end

figure
for k=1:2
    subplot(2,1,k)
    plot(t,pmech(k,:),'k','LineWidth',2)
    hold on
    plot(tOut,Gen(k).Pmech(1:N),'r.','MarkerSize',10);
    hold off
    if k==1; legend('PST','Dan'); end
    ylabel(['Gen ' num2str(k) ' Pmech (pu)'])
    xlabel('Time (sec.)')
end

figure
for k=1:2
    subplot(2,1,k)
    plot(t,pelect(k,:),'k','LineWidth',2)
    hold on
    plot(tOut,Gen(k).Pe(1:N),'r.','MarkerSize',10);
    hold off
    if k==1; legend('PST','Dan'); end
    ylabel(['Gen ' num2str(k) ' Pe (pu)'])
    xlabel('Time (sec.)')
end

figure
for k=1:2
    subplot(2,1,k)
    plot(t,abs(bus_v(k,:)),'k','LineWidth',2)
    hold on
    plot(tOut,abs(Gen(k).VT(1:N)),'r.','MarkerSize',10);
    hold off
    if k==1; legend('PST','Dan'); end
    ylabel(['Gen ' num2str(k) ' VT (pu)'])
    xlabel('Time (sec.)')
end

figure
for k=1:2
    subplot(2,1,k)
    plot(tOut,Gen(k).wdot(1:N),'r.','LineWidth',2,'MarkerSize',10);
    ylabel(['Gen ' num2str(k) ' wdot (pu)'])
    xlabel('Time (sec.)')
end
    






