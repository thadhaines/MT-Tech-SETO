function status = vtsOutputFcn(t,y,flag)
%VTSOUTPUTFCN performs associated flag actions with ODE solvers
% VTSOUTPUTFCN performs associated flag actions with ODE solvers.
% Handles 3 conditions: step (flag == empty), init, done
% uses globals for logging

global g % used for data logging

status = 0; % required for normal operation (return 1 to stop)

%{
    The solver calls status = myOutputFcn(t,y,[]) after each integration
step for which output is requested. t contains points where output was
generated during the step, and y is the numerical solution at the points
in t. If t is a vector, then the ith column of y corresponds to the ith
element of t.

If length(tspan) > 2, then the output is produced at every point in tspan.
If length(tspan) = 2, then the output is produced according to the Refine option.
%}
if isempty(flag) % normal step completion
%% =============================================================================
%% Line Monitoring and Area Calculations =======================================
    %% Line Monitoring
    if g.lmon.n_lmon~=0
        lmon(g.vts.dataN)
    end
    
    %% Average Frequency Calculation
    calcAveF(g.vts.dataN,1);
    
    %% Area Total Calcvulations
    if g.area.n_area ~= 0
        calcAreaVals(g.vts.dataN,1);
    end
    
    %% Live plot call
    if g.sys.livePlotFlag
        livePlot(g.vts.dataN)
    end
    
    % after each successful integration step by solver:
    g.vts.dataN = g.vts.dataN+1; % increment logged data index 'dataN'
    g.sys.t(g.vts.dataN) = t; % log step time
    g.vts.stVec = y; % update state vector
    % i.e. call handleStDx to place new solution results into associated globals
    handleStDx(g.vts.dataN, y, 2)
    
    fprintf('* Data step: %6d\tat time:\t%8.4f\tused %4d solutions...\n', g.vts.dataN, t, g.vts.iter)
    g.vts.tot_iter = g.vts.tot_iter  + g.vts.iter; % count total iterations...
    g.vts.iter = 0; % reset iteration (solution) counter
    
    %{
    The solver calls myOutputFcn([tspan(1) tspan(end)],y0,'init') before
beginning the integration to allow the output function to initialize.
tspan and y0 are the input arguments to the ODE solver.
    %}
elseif flag(1) == 'i' % init solver for time period t
    
    
    g.vts.dataN = g.vts.dataN+1; % increment logged data index 'dataN'
    g.sys.t(g.vts.dataN) = t(1); % log step time
    handleStDx(g.vts.dataN, y, 2) % log initial conditions
    
    % debug display
    disp('*** ')
    disp('Flag == init')
    fprintf('Data step: %d\t%3.5f\n', g.vts.dataN, t(1))

    %{
    The solver calls myOutputFcn([],[],'done') once integration is
complete to allow the output function to perform cleanup tasks.
    %}
elseif flag(1) == 'd' % time period complete
    % debug display
    disp('*** ')
    disp('Flag == done')
    fprintf('Last complete data step: %d\t%8.4f\n', g.vts.dataN, g.sys.t(g.vts.dataN))
    
    % NOTE: Called after last normal integration step
    %%%g.vts.dataN = g.vts.dataN+1; % increment logged data index 'dataN'
    % included here as init does not increment data index
    % NOTE: extra index must be removed after all simulated times are complete
end

end