function agc(k,flag)
%AGC performs actions related to automatic generation control.
% AGC performs actions related to automatic generation control.
%
% Syntax: agc(k,flag)
%
%   NOTES:  Idea is to piggy back AGC signals to the tg_sig variable.
%           No integrator windup is present in PI filter.
%           IACE inclusion is pretty shoddy, but exists - use at own risk.
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
%   07/22/20    10:17   Thad Haines     Version 1.0.1 - Added conditional AGC
%   08/04/20    06:07   Thad Haines     Version 1.0.2 - adjusted for VTS

global g

if flag == 0
    % initialize
    
    % for each defined agc
    for n = 1:g.agc.n_agc
        % set initial ACE signal to zero for all time
        g.agc.agc(n).aceSig = zeros(1,size(g.agc.agc(n).ace2dist,2));
        % Create place holders for current and available capacity
        g.agc.agc(n).curGen = zeros(1,size(g.agc.agc(n).ace2dist,2)); % current AGC gen output
        
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
        
        % calculate total AGC capacity and percent/mw available
        for cg = 1:g.agc.agc(n).n_ctrlGen
            g.agc.agc(n).curGen(1) = g.agc.agc(n).curGen(1) ...
                + g.mac.pmech(g.agc.agc(n).macBusNdx(cg),1)* g.agc.agc(n).ctrlGen(cg).mBase;
        end
    end
    
    % inital states of generation control set when zeros allocated
    
end

% No action if flag == 1

if flag == 2 && k ~=1  % skip first step
    
    % create index for most recent monitored values
    j = k-1;
    
    %for each area
    for n = 1:g.agc.n_agc
        % calcualte RACE using most recently monitored values via j index
        delta_w = g.area.area(g.agc.agc(n).area).aveF(j) - g.sys.sys_freq(j);
        icError = g.area.area(g.agc.agc(n).area).icA(j) - ...
            g.area.area(g.agc.agc(n).area).icS(j);
        
        % handle frequency bias scaling
        fError = delta_w * g.agc.agc(n).Bcalc * 10 * (1+abs(delta_w)*g.agc.agc(n).Kbv);
        % calculate race using real power and PU fError
        g.agc.agc(n).race(k) = real(icError) + fError*(g.sys.Fbase/g.sys.basmva); % ensure PU value scaled correctly
        
        % set derivative of smoothed ace
        g.agc.agc(n).d_sace(k) = g.agc.agc(n).race(k);
        % handle output of pi filter
        g.agc.agc(n).sace(k) = g.agc.agc(n).Kp*(g.agc.agc(n).sace(j)*g.agc.agc(n).a ...
            + g.agc.agc(n).race(k));
        
        %TODO: add windup correction ?
        
        %         % handle integration of RACE -> IACE (not used...)
        %         g.agc.agc(n).iace(k) = g.agc.agc(n).iace(j) + (g.agc.agc(n).race(k) + g.agc.agc(n).race(j)) /2 ...
        %             * abs(g.sys.t(k)-g.sys.t(j));
        
        % check if due to distribute
        if (g.sys.t(k) >= g.agc.agc(n).nextActionTime)
            % update ACE signal via RACE
            % placeholder notification
            fprintf('*** t = %4.5f Distributing Area %d ACE...\n',g.sys.t(k), g.agc.agc(n).area);
            
            % Conditional logic
            if g.agc.agc(n).condAce == 1
                % check if current SACE is of the same sign as delta_w
                condOk = sign(g.agc.agc(n).sace(k)) == sign(delta_w);
            else
                condOk = 1;
            end
            %fprintf('DEBUG AGC %d, delta_w = %2.7f\t sign of RACE %d\t condOk = %d\n' , ...
            % n, delta_w, sign(g.agc.agc(n).sace(k)), condOk);
            
            % aceSig is signal sent to all generators (after a gain) % maybe not the most efficient method...
            g.agc.agc(n).aceSig(k:end) = g.agc.agc(n).aceSig(k) + ... % modified..
                g.agc.agc(n).sace(k)*condOk;
            
            % increment nextActionTime
            g.agc.agc(n).nextActionTime = g.agc.agc(n).nextActionTime ...
                + g.agc.agc(n).actionTime;
            
            clear condOk
        end
        
        % Gain aceSig and Log ace2dist
        g.agc.agc(n).ace2dist(k) = g.agc.agc(n).aceSig(k) * g.agc.agc(n).gain;
        
        % reset current generation value
        g.agc.agc(n).curGen(k) = 0;
        
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
                - g.agc.agc(n).ctrlGen(gNdx).x(k); % NOTE: negative sign is important here.
            % update curGen value
            g.agc.agc(n).curGen(k) = g.agc.agc(n).curGen(k) ...
                + g.mac.pmech(g.agc.agc(n).macBusNdx(gNdx),k)* g.agc.agc(n).ctrlGen(gNdx).mBase;
        end
    end
end
end% end function