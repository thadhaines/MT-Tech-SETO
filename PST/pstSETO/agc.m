function agc(k,flag)
%AGC performs actions related to automatic generation control.
% AGC performs actions related to automatic generation control.
%
% Syntax: agc(k,flag)
%
%   NOTES:  Idea is to piggy back AGC signals to the tg_sig variable..
%
%   Input:
%   k - data index
%   flag  - dictate what operation to perform
%       0 - initialize AGC
%       1 - calculate ACE and distribute to control generators
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/20/20    20:09   Thad Haines     Version 1

global g

if flag == 0
    % initialize
    
    % for each defined agc
    for n = 1:g.agc.n_agc
        % set initial ACE signal to zero
        g.agc.agc(n).aceSig = 0;
        
        % set icS intial condition
        g.area.area(g.agc.agc(n).area).icS = g.area.area(g.agc.agc(n).area).icS ...
            * g.area.area(g.agc.agc(n).area).icA(1);
        
        % calcualte each AGC B per user input
        if g.agc.agc(n).Btype == 1
            % calculate B as a percentage of maximum area capacity
            g.agc.agc(n).Bcalc = g.agc.agc(n).B* ...
                g.area.area(g.agc.agc(n).area).maxCapacity/100;
        else
            % use absolute Frequency bias
            g.agc.agc(n).Bcalc = g.agc.agc(n).B;
        end
        
        % set nextActionTime to start time
        g.agc.agc(n).nextActionTime = g.agc.agc(n).startTime;
    end
    
    % set inital states of generation control...
    
end

% No action if flag == 1

if flag == 2 && k ~=1 % skip first step
    
    % create index for most recent monitored values
    j = k-1;
    
    %for each area
    for n = 1:g.agc.n_agc
        % calcualte RACE using most recently monitored values
        delta_w = g.area.area(g.agc.agc(n).area).aveF(j) - g.sys.sys_freq(j);
        icError = g.area.area(g.agc.agc(n).area).icA(j) - ...
            g.area.area(g.agc.agc(n).area).icS(j);
        
        % handle frequency bias scaling
        fError = delta_w * g.agc.agc(n).Bcalc * 10 * (1+abs(delta_w)*g.agc.agc(n).kBv);
        
        g.agc.agc(n).race(k) = real(icError) + fError*(g.sys.Fbase/g.sys.basmva); % ensure PU value
        
        % check if due to distribute
        if g.sys.t(k) >= g.agc.agc(n).nextActionTime
            % update ACE signal via RACE
            % placeholder notification
            fprintf('*** t = %4.5f Distributing Area %d ACE...\n',g.sys.t(k), g.agc.agc(n).area);
            
            % conditional logic should go here
            
            g.agc.agc(n).aceSig = g.agc.agc(n).race(k);
            
            % increment nextActionTime
            g.agc.agc(n).nextActionTime = g.agc.agc(n).nextActionTime ...
                + g.agc.agc(n).actionTime;
        end
        
        % quick test to see if tg sig work - NOT final solution
        % essentially just a Proportional controller...
        g.agc.agc(n).ace2dist(k) = g.agc.agc(n).aceSig * g.agc.agc(n).Kp;
        
        % ensure agc signal always updated in tg model
        for gNdx = 1:g.agc.agc(n).n_ctrlGen
            g.tg.tg_sig(g.agc.agc(n).tgNdx(gNdx),k) = g.tg.tg_sig(g.agc.agc(n).tgNdx(gNdx),k) ...
                - g.agc.agc(n).ace2dist(k)*g.agc.agc(n).ctrlGen(gNdx).pF;
        end
        
    end
    
    
end
end% end function