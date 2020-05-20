% Integrator Example 04


clear all
close all
clc

euler(1) = 1;
predcor(1) = 1;
t(1) = 0;
ii = 1;
h = .05;

while t(ii) <= 10
    % Euler Cauchy.
    pred = predcor(ii) + h/2*(-predcor(ii)^2);
    predcor(ii+1) = predcor(ii) + h*(-pred^2);
    
    t(ii+1) = t(ii) + h;
    ii = ii + 1;
    
end
ii = 1;
tt(1) = 0;

while t(ii) <= 10

    euler(ii+1) = euler(ii) + h*(-euler(ii)^2);
    
    tt(ii+1) = tt(ii) + h;
    ii = ii + 1;
    
end

Exact = 1./(1+t);

figure
plot(t,Exact,'r',tt,euler,'b*',t,predcor,'g')

figure
subplot(211)
plot(t,(Exact-euler))
subplot(212)
plot(t,(Exact-predcor))