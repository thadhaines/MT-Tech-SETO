function status = vtsOutputFcn(t,y,flag)
%VTSOUTPUTFCN performs associated flag actions with ODE solvers
% VTSOUTPUTFCN performs associated flag actions with ODE solvers.
%
% Syntax: status = vtsOutputFcn(t,y,flag)
%
%   NOTES: Handles 3 conditions: step (flag == empty), init, done
%   Uses globals for logging
%
%   Input:
%   t - simulation time
%   y - solution vector
%   flag - dictate function action
%
%   Output:
%   status - required for normal operation (return 1 to stop)
%
%   History:
%   Date        Time    Engineer        Description
%   07/27/20    11:19   Thad Haines     Version 1
%   08/03/20    10:24   Thad Haines     Version 1.1 - clean up actions

global g
persistent FirstRun
status = 0; % required for normal operation

%{
    MATLAB DOCUMENTATION:
The solver calls status = myOutputFcn(t,y,[]) after each integration
step for which output is requested. t contains points where output was
generated during the step, and y is the numerical solution at the points
in t. If t is a vector, then the ith column of y corresponds to the ith
element of t.

If length(tspan) > 2, then the output is produced at every point in tspan.
If length(tspan) = 2, then the output is produced according to the Refine option.
%}
if isempty(flag) % normal step completion
    
    % restore network to initial solution
    handleNetworkSln(g.vts.dataN ,2) % may cause issues with DC.
    % Save and restore first prediction of dx in a similar fashion as newtork soln?
    
    %% Line Monitoring and Area Calculations ===================================
    monitorSolution(g.vts.dataN);
    
    %% Live plot call
    if g.sys.livePlotFlag
        livePlot(g.vts.dataN)
    end
    
    % after each successful integration step by ODE solver:
    g.vts.dataN = g.vts.dataN+1;    % increment logged data index 'dataN'
    g.sys.t(g.vts.dataN) = t;       % log step time
    g.vts.stVec = y;                % update state vector
    handleStDx(g.vts.dataN, y, 2)   % place new solution results into associated globals
    
    % display dataN and t at every first, last, and 250th step
    if ( mod(g.vts.dataN,250)==0 ) || g.vts.iter > 100 || FirstRun
        fprintf('*** dataN: %5d, g.sys.t(dataN) = %10.6f\tperformed %d solutions\n', ...
            g.vts.dataN, t, g.vts.iter)
        FirstRun = 0;
    end
    
    g.vts.tot_iter = g.vts.tot_iter + g.vts.iter;  % update total iterations
    g.vts.slns(g.vts.dataN) = g.vts.iter;          % log solution step iterations
    g.vts.iter = 0;                                % reset iteration counter
    
    %{
    MATLAB DOCUMENTATION:
The solver calls myOutputFcn([tspan(1) tspan(end)],y0,'init') before
beginning the integration to allow the output function to initialize.
tspan and y0 are the input arguments to the ODE solver.
    %}
    
elseif flag(1) == 'i' % init solver for new time block
    FirstRun = 1;
    g.sys.t(g.vts.dataN) = t(1);    % log step time
    handleStDx(g.vts.dataN, y, 2)   % set initial state conditions
    % set initial dx conditions? - Shouldn't matter as input function calculates new derivatives for dataN...
    
    % debug display
    if g.sys.DEBUG
        disp('*** VTS Flag == init')
        fprintf('*** Data step: %d\t%3.5f\n', g.vts.dataN, t(1))
    end
    
    %{
    MATLAB DOCUMENTATION:
The solver calls myOutputFcn([],[],'done') once integration is
complete to allow the output function to perform cleanup tasks.
    %}
    
elseif flag(1) == 'd'
    % time block period complete
    
    % debug display
    if g.sys.DEBUG
        disp('*** VTS Flag == done')
        fprintf('*** Last complete data step: %d\t%8.10f\n', g.vts.dataN, g.sys.t(g.vts.dataN))
    end
end

end