function [dxVec] = vtsInputFcn(t, y)
%VTSINPUTFCN passed to ODE solver to perfrom required step operations
% VTSINPUTFCN passed to ODE solver to perfrom required step operations
%
% Syntax: vtsInputFcn(k)
%
%   NOTES: Updates g.vts.dxVec, and returns values
%
%   Input:
%   t - simulation time
%
%   Output:
%   dxVec - requried derivative vector for ODE solver
%
%   History:
%   Date        Time    Engineer        Description
%   07/27/20    11:19   Thad Haines     Version 1

global g

%% call handleStDx with flag==2 to update global states with newest passed in soln.
% write slnVec vector of values to associated states at index k
% i.e. update states at g.vts.dataN with newest solution
  handleStDx(g.vts.dataN, y, 2)

%% Line Monitoring and Area Calculations ==================================
monitorSolution(g.vts.dataN);

%% Start initStep action ==================================================
initStep(g.vts.dataN)

%% Start of Network Solution ==============================================
networkSolutionVTS(g.vts.dataN, t)

%% Start Dynamic Solution =================================================
dynamicSolution(g.vts.dataN )

%% Start of DC solution ===================================================
dcSolution(g.vts.dataN )

%% call handleStDx with flag==1 to update global dxVec
handleStDx(g.vts.dataN , [], 1) % update g.vts.dxVec (solution vector not needed)

dxVec = g.vts.dxVec; % return for ODE fcn requirements

g.vts.iter = g.vts.iter + 1;
end % end vtsInputFcn