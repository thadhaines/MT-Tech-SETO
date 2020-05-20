clear all
close all
clc

% Integrator Example

y(1) = 0;
h = pi/10000;
t(1) = 0;

for ii = 1:20000
    pred = y(ii) + h * cos(t(ii));
    y(ii+1) = y(ii) + .5*h*(cos(t(ii)+h) + cos(t(ii)));    
    t(ii+1) = t(ii) + h;
end

Exact = sin(t);

figure
plot(t,Exact,'b',t,y,'r')

for ii = 1:length(t)
    Delta(ii) = Exact(ii) - y(ii); 
end

figure
plot(t,Delta)
