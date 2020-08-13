function [r,x]=dbcage(r1,x1,r2,x2,s)
% calculates the equivalent rotor resistance and reactance of a double cage rotor
% with slip (s)
% called by mac_ind
z = 1j*x1 + (r1./s).*(r2./s + 1j*x2)./((r1+r2)./s + 1j*x2);
r = s.*real(z);
x = imag(z);
