% close all

jj = 1;

figure
E = compass(Gen(jj).E,'b');
hold on
Vt = compass(Gen(jj).Vt,'r');
I = compass(Gen(jj).I,'g');
D = compass(Gen(jj).E*exp(-1j*pi/2),'m');
set(E,'LineWidth',2)
set(Vt,'LineWidth',2)
set(I,'LineWidth',2)
set(D,'LineWidth',2)

Vd = abs(Gen(1).Vt)*sin(Gen(jj).RotAng - angle(Gen(jj).Vt));
Vq = abs(Gen(1).Vt)*cos(Gen(jj).RotAng - angle(Gen(jj).Vt));
VD = compass(Vd*exp(1j*(Gen(jj).RotAng - pi/2)),'k');
VQ = compass(Vq*exp(1j*Gen(jj).RotAng),'k');
set(VD,'LineWidth',2)
set(VQ,'LineWidth',2)

Id = abs(Gen(jj).I)*sin(Gen(jj).RotAng - angle(Gen(jj).I));
Iq = abs(Gen(jj).I)*cos(Gen(jj).RotAng - angle(Gen(jj).I));
ID = compass(Id*exp(1j*(Gen(jj).RotAng - pi/2)),'c');
IQ = compass(Iq*exp(1j*Gen(jj).RotAng),'c');
set(ID,'LineWidth',2)
set(IQ,'LineWidth',2)

legend('E','Vt','I','D axis','Vd','Vq','Id','Iq')







