function agc(k,flag)
%AGC performs actions related to automatic generation control.
% AGC performs actions related to automatic generation control.
%
% Syntax: agc(k,flag)
%
%   NOTES:  Idea is to piggy back AGC signals to the tg_sig variable.
%           No integrator windup is present in PI filter.
%           Coonditional logic not coded yet.
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
%   07/21/20    15:58   Thad Haines     Version 1

global g

if flag == 0
    % initialize
    
    % for each defined agc
    for n = 1:g.agc.n_agc
        % set initial ACE signal to zero
        g.agc.agc(n).aceSig = 0;
        g.agc.agc(n).iaceStartNdx = 1;
        g.agc.agc(n).iaceN = 0; % number of  iace vals in window intergrator
        
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
    
    % inital states of generation control set when zeros allocated in s_simu_batch
    
end

% No action if flag == 1

if flag == 2 && k ~=1 % skip first step
    
    % create index for most recent monitored values
    j = k-1;
    
    %for each area
    for n = 1:g.agc.n_agc
        % calcualte RACE using most recently monitored values via j index
        delta_w = g.area.area(g.agc.agc(n).area).aveF(j) - g.sys.sys_freq(j);
        icError = g.area.area(g.agc.agc(n).area).icA(j) - ...
            g.area.area(g.agc.agc(n).area).icS(j);
        
        % handle frequency bias scaling
        fError = delta_w * g.agc.agc(n).Bcalc * 10 * (1+abs(delta_w)*g.agc.agc(n).kBv);
        % calculate race using real power and PU fError
        g.agc.agc(n).race(k) = real(icError) + fError*(g.sys.Fbase/g.sys.basmva); % ensure PU value
        
        % set derivative of smoothed ace
        g.agc.agc(n).d_sace(k) = g.agc.agc(n).race(k);
        % handle output of pi filter
        g.agc.agc(n).sace(k) = g.agc.agc(n).Kp*(g.agc.agc(n).sace(j)*g.agc.agc(n).a ...
            + g.agc.agc(n).race(k)); 
        
        %TODO: add windup correction ?
        
        g.agc.agc(n).iace(k) = (g.agc.agc(n).iace(k) + g.agc.agc(n).iace(j))/2;
        g.agc.agc(n).iaceN = g.agc.agc(n).iaceN + 1;
        
        % check if due to distribute
        if g.sys.t(k) >= g.agc.agc(n).nextActionTime
            % update ACE signal via RACE
            % placeholder notification
            fprintf('*** t = %4.5f Distributing Area %d ACE...\n',g.sys.t(k), g.agc.agc(n).area);
            
            % TODO: create conditional logic here
            
            % aceSig is the output of the PI filter - updatd every actionTime
%             g.agc.agc(n).aceSig = g.agc.agc(n).aceSig + g.agc.agc(n).sace(k); % no IACE?

            % attempt at window intergration 'average window value'
            g.agc.agc(n).aceSig = g.agc.agc(n).aceSig + g.agc.agc(n).sace(k)...
                + g.agc.agc(n).Kiace * sum(g.agc.agc(n).iace(g.agc.agc(n).iaceStartNdx:k)) ...
                / g.agc.agc(n).iaceN;

            g.agc.agc(n).iaceStartNdx = k; % set next start index
            g.agc.agc(n).iaceN = 0; % reset counter
            
            % increment nextActionTime
            g.agc.agc(n).nextActionTime = g.agc.agc(n).nextActionTime ...
                + g.agc.agc(n).actionTime;
        end
        
        % Log ace2dist
        g.agc.agc(n).ace2dist(k) = g.agc.agc(n).aceSig * g.agc.agc(n).gain;
        
        % ensure agc signal updated in tg model for each agc ctrlGen
        for gNdx = 1:g.agc.agc(n).n_ctrlGen
            % send signal to input of lowpass filter
            g.agc.agc(n).ctrlGen(gNdx).input(k) =  g.agc.agc(n).ace2dist(k) ...
            * g.agc.agc(n).ctrlGen(gNdx).pF;
            % set derivative of filter
            g.agc.agc(n).ctrlGen(gNdx).dx(k) = (-g.agc.agc(n).ctrlGen(gNdx).x(k) ...
                +  g.agc.agc(n).ctrlGen(gNdx).input(k))/g.agc.agc(n).ctrlGen_con(gNdx,3);
            % send filtered signal to governor
            g.tg.tg_sig(g.agc.agc(n).tgNdx(gNdx),k) = g.tg.tg_sig(g.agc.agc(n).tgNdx(gNdx),k) ...
                - g.agc.agc(n).ctrlGen(gNdx).x(k);
        end
    end
end
end% end function