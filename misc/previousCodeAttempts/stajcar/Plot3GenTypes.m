% Plot E, Wr, and RotAng of Dan's PST, Genrou, and GenTpJ.

clear
close all
clc

load Genrou
load GenTpJ
load PSTN2
load PST

figure
plot(t,abs(Gen1E_PSTN2),'b',t,abs(Gen1E_Genrou),'r',t,abs(Gen1E_GenTpJ),'k')
title('Gen 1 E')
ylabel('Voltage (volts) [V]')
legend('PST','Genrou','GenTpJ')

figure
plot(t,abs(Gen2E_PSTN2),'b',t,abs(Gen2E_Genrou),'r',t,abs(Gen2E_GenTpJ),'k')
title('Gen 2 E')
ylabel('Voltage (volts) [V]')
legend('PST','Genrou','GenTpJ')

figure
plot(t,Gen1Wr_PSTN2,'b',t,Gen1Wr_Genrou,'r',t,Gen1Wr_GenTpJ,'k',tPST,w2PST,'g')
title('Gen 1 Wr')
ylabel('Speed (per unit) [pu]')
legend('PST Model','Genrou','GenTpJ','PST code')

figure
plot(t,Gen2Wr_PSTN2,'b',t,Gen2Wr_Genrou,'r',t,Gen2Wr_GenTpJ,'k',tPST,w1PST,'g')
title('Gen 2 Wr')
ylabel('Speed (per unit) [pu]')
legend('PST Model','Genrou','GenTpJ','PST code')

figure
plot(t,Gen1RotAng_PSTN2*180/pi,'b',t,Gen1RotAng_Genrou*180/pi,'r',t,Gen1RotAng_GenTpJ*180/pi,'k',tPST,ang2PST*180/pi,'g')
title('Gen 1 RotAng')
ylabel('Angle (degrees) [deg]')
legend('PST Model','Genrou','GenTpJ','PST code')

figure
plot(t,Gen2RotAng_PSTN2*180/pi,'b',t,Gen2RotAng_Genrou*180/pi,'r',t,Gen2RotAng_GenTpJ*180/pi,'k',tPST,ang1PST*180/pi,'g')
title('Gen 2 RotAng')
ylabel('Angle (degrees) [deg]')
legend('PST Model','Genrou','GenTpJ','PST code')