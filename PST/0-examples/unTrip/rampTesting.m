% script to test manual value ramping over a set interval

clear; close all
ts = 1e-3;
addRmp = 2; % amount to add during a ramp

rmpStart = 2;
rmpMid = 5;
rmpEnd = 8;
rmpDur = rmpEnd-rmpStart;
rmpHalf = rmpMid-rmpStart;

t = 0:ts:10-ts;

u = ones(size(t));
mod = zeros(size(u));

y = zeros(size(u));


a = 30;

for n = 1:length(t)
    
%     
%         % linear ramping
%         if t(n)>= rmpStart && t(n) < rmpEnd
%             mod(n) = addRmp * (t(n)-rmpStart)/(rmpEnd-rmpStart);
%         elseif t(n) >= rmpEnd
%                 mod(n) = addRmp;
%         end
    
    % exponential ramping
    if t(n)>= rmpStart && t(n) < rmpEnd
        mod(n) = addRmp*(1 - exp( rmpStart-t(n) ) ); % concave down, growth
        %mod(n) = addRmp*( exp( t(n)-rmpEnd ) ); % concave up, growth
        
        % non working
        % logistic curve - not really
%         mod(n) = addRmp*( 1/(1 + a*exp(rmpStart-t(n))) - 1/(a+1)  ); 
        % gauss - not really
       % mod(n) = addRmp * exp( (t(n)-rmpEnd) / rmpDur*4)^2 ;
    elseif t(n) >= rmpEnd
        mod(n) = addRmp;
    end

    y(n) = u(n) + mod(n);
    
end

figure
plot(t,u,'k')
hold on
plot(t,mod,'--')
plot(t,y,'m')

