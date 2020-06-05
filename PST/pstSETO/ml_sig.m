function ml_sig(k)
% ML_SIG Defines modulation signal for lmod control
% Syntax: f = ml_sig(k)
%
%   History:
%   Date        Time    Engineer        Description
%   08/15/97    16:40   Graham Rogers   Version 1
%   06/05/20    09:02   Thad Haines     V2 - using global g, no t passed
%   in, no dummy varibale f passed out. Default behavior commented out.
%

%global g

%fprintf('%4.4f \t %d\n', t(k), k); % DEBUG

%if g.sys.t(k) > 1
    % load step
    %g.lmod.lmod_sig(1,k) = 0.25; % modify first load only
%end
end