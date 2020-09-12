function huensMethod(simTblock)
%HUENSMETHOD Simulates given time block using Huen's method
% HUENSMETHOD Simulates given time block using Huen's method
%
% Syntax: huensMethod(simTblock)
%
%   NOTES:
%
%   Input:
%   simTblock - Index of time block to simulate
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   09/12/20    11:31   Thad Haines     Version 1

%%
global g

fprintf(['*~~ Using Huen''s integration method ', ...
    'for time block %d\n*** t_block = [%7.4f, %7.4f ]\n'], ...
    simTblock, g.vts.fts{simTblock}(1), g.vts.fts{simTblock}(end))

% incorporate fixed time vector into system time vector
nSteps = length(g.vts.fts{simTblock});
g.sys.t(g.vts.dataN:g.vts.dataN+nSteps-1) = g.vts.fts{simTblock};

% account for pretictor last step time check
g.sys.t(g.vts.dataN+nSteps) = g.sys.t(g.vts.dataN+nSteps-1) ...
    + g.sys.sw_con(simTblock,7);

for cur_Step = 1:nSteps
    k = g.vts.dataN;
    j = k+1;
    
    % display k and t at every first, last, and 250th step
    if ( mod(k,250)==0 ) || cur_Step == 1 || cur_Step == nSteps
        fprintf('*** k = %5d, \tt(k) = %7.4f\n',k,g.sys.t(k))
    end
    
    %% Time step start
    initStep(k)
    
    %% Predictor Solution ======================================================
    networkSolutionVTS(k, g.sys.t(k))
    % most recent accepted network soln is at index k
    monitorSolution(k);
    dynamicSolution(k)
    dcSolution(k)
    % NOTE: g.k.h_sol updated in network solution i_simu call
    predictorIntegration(k, j, g.k.h_sol)
    
    %% Corrector Solution ======================================================
    networkSolutionVTS(j, g.sys.t(j))
    dynamicSolution(j)
    dcSolution(j)
    correctorIntegration(k, j, g.k.h_sol)
    
    %% Live plot call
    if g.sys.livePlotFlag
        livePlot(k)
    end
    
    g.vts.dataN = j;                        % inc data counter
    g.vts.tot_iter = g.vts.tot_iter  + 2;   % inc tot soln counter
    g.vts.slns(g.vts.dataN) = 2;            % track step solution
end
% Account for next time block using VTS
handleStDx(j, [], 3) % update g.vts.stVec to init. cond. of St
handleStDx(k, [], 1) % update g.vts.dxVec to init. cond. of dx

end % end huensMethod