function [ xdot ] = getXdot( t, x)
%getXdot return xdot from statespace for ODE45 use
%   t = filler variable
%   A = A matrix from system
%   x = initial state vector
%   B = B matrix from system
%   U = Input to system
global A B U
    xdot = A*x + B*U;
end

