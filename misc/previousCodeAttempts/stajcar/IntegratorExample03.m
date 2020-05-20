% Integrator Example 03


clear all
close all
clc

y(1) = 1;
yy(1) = 1;
t(1) = 0;
h = .1;
f = @(y) -y^2;

for ii = 1:300
    % Euler Trapezoidal.
    pred = yy(ii) + h*(-yy(ii)^2);
    yy(ii+1) = yy(ii) + h/2*(-pred^2 + -yy(ii)^2);
    
    t(ii+1) = t(ii) + h;
    
end

for ii = 1:300

    y(ii+1) = y(ii) + h*(-y(ii)^2);
    
    t(ii+1) = t(ii) + h;
    
end

Exact = 1./(1+t);

figure
plot(t,Exact,'r',t,y,'b',t,yy,'g')

figure
subplot(211)
plot(t,(Exact-y))
subplot(212)
plot(t,(Exact-yy))

