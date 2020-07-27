function status = testOutF1(t,y,flag)
% testOutF1 to test output function associated with ODE solvers
% handles 3 conditions: step (flag == empty), init, done
% uses globals for logging

global yGlobal tGlobal dataN % used for data logging
global A B C D U % used for output Calculation

%{
myOutputFcn must return a status of 0 or 1.
If status = 1, then the solver halts integration.
You can use this mechanism, for instance, to implement a Stop button.
%}
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
    
    % after each successful integration step by solver:
    dataN = dataN+1; % increment logged data index 'dataN'
    tGlobal(dataN) = t; % log step time
    
    yGlobal(dataN) = C*y+D*U; % log result
    % i.e. call handleStDx to place new solution results into associated globals
    
    % debug display
    %if mod(dataN, 42) == 0 || dataN == 1
        disp('Flag is empty')
        dataN
        t
        y
    %end
    
    %{
    The solver calls myOutputFcn([tspan(1) tspan(end)],y0,'init') before
beginning the integration to allow the output function to initialize.
tspan and y0 are the input arguments to the ODE solver.
    %}
elseif flag(1) == 'i' % init solver for time period t
    
    tGlobal(dataN) = t(1); % log starting time
    yGlobal(dataN) = C*y+D*U; % log initial conditions
    
    % debug display
    disp('*** ')
    disp('Flag == init')
    t
    y
    dataN
    
    %{
    The solver calls myOutputFcn([],[],'done') once integration is
complete to allow the output function to perform cleanup tasks.
    %}
elseif flag(1) == 'd' % time period complete
    % debug display    
    disp('*** ')
    disp('Flag == done')
    t
    y
    dataN
    
    
    % NOTE: Called after last normal integration step
    dataN = dataN+1; % included here as init does not increment data index
    % NOTE: extra index must be removed after all simulated times are complete
end

end