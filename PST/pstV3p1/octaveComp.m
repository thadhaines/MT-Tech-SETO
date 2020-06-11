%% Octave specific work arounds for PST compatibility
% https://wiki.octave.org/Differences_between_Octave_and_Matlab

%% supress warning about | and & usage vs || and && 
% while this doesn't cause code to crash, it does cause many warning to be 
% output to the console which greately slows down simulations
warning('off', 'Octave:possible-matlab-short-circuit-operator'); 

%% Load built in package for state space functionality in linear analysis
pkg load control 

%% Always output MATLAB compatible data
save_default_options ('-mat7-binary')


%% Other Octave notes:
%{
    MATLAB simulation runs are generally faster than Octave by about 50%...
    
    Calls to save and load must include explicit file ending of .mat (or other).
    
    Octave does not run code sections as MATLAB does using ctrl+Enter,
    instead, highlight code to run and press F9.
    
    Octave does not handle legend fontsizes defined in legend calls.
    In General, MATLAB plots look nicer than Octave plots.
    
 %  Handled compatiblity issues:
    
    To enable MATLAB to read saved octave data, '-mat7-binary' must be used 
    as an option during save call e.g.:
      save('-mat7-binary','saveFile.mat', 'var1', 'var2', ... )
    Handled save_default_options code above.
      
    stepfun is not included in Octave -> similar funciton saved in root PST directory as stepfun
    
    
%}
