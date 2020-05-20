clear all
close all
clc

% Integrator Example

y(1) = 0;
h = pi/2;
t(1) = 0;

%Euler Trapezoidal
ii = 1;
while t(ii) < 50
    y(ii+1) = y(ii) + h/2*(cos(t(ii)+h)+.02 + cos(t(ii))+.02);    
    t(ii+1) = t(ii) + h;
    ii = ii + 1;
end

%Euler Cauchy
% for ii = 1:20000
%     pred = y(ii) + h/2 * (cos(t(ii))+2);
%     y(ii+1) = y(ii) + h*(cos(t(ii)+h/2)+2);    
%     t(ii+1) = t(ii) + h;
% end

Exact = sin(t) + .02*t;

figure
plot(t,Exact,'b',t,y,'r')

for ii = 1:length(t)
    Delta(ii) = Exact(ii) - y(ii); 
end

figure
plot(t,Delta)
