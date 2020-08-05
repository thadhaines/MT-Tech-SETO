clear all; close all; clc
% generate sine wave to simulate solar generation cycle and could cover
% 4 minutes (240 seconds)
% 0-30      - no action
% 30-90     - ramp up
% 90-150    - hold peak
% 150-210   - ramp down
% 210-240   - no action

% cloud cover events
% 45-55 - 20% max gen
% 120-140 - 30% cover
% 180-190 - 15% cover

% running variables
ts = 0:1/60:240;            % 4 minutes of time
xs = zeros(length(ts),1);   % modulation signal

maxGen = 3;                 % peak of genration

% create modulation signal based on time (similar to PST simulation)
for n=1:length(ts)
    
    % Handle basic generation curve
    if ts(n) >= 30 && ts(n) < 90
        % up ramp
        % use first 1/4 cycle of a sin to ramp up
        xs(n) = maxGen*sin( (ts(n)-30) * 2*pi/240); 
    elseif ts(n) >= 150 && ts(n) < 210
        % down ramp
        % use second 1/4 cycle of a sin to ramp down
        xs(n) = maxGen*sin( (ts(n)-150) * 2*pi/240 + pi/2);   
    elseif ts(n) >= 90 && ts(n) < 150
        % held peak
        xs(n) = maxGen;
    else
        % no generation
        xs(n) = 0;
    end
    
    % Handle Cloud cover generation reductions
    if ts(n) >= 45 && ts(n) < 55
        xs(n) = maxGen*.2;
    elseif ts(n) >= 120 && ts(n) < 140
        xs(n) = xs(n)*0.7;
    elseif ts(n) >= 180 && ts(n) < 190
        xs(n) = xs(n)*0.85;
    end
    
end

% plot modulation signal
figure
plot(ts,xs)
clear ts xs w