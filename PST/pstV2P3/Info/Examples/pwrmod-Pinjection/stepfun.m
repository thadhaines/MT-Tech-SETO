function y = stepfun(t,t0)
%stepfun Unit step function.
%
%   stepfun(t,t0), where t is an increasing vector, t0 is a float or int value
%   returns a vector the same length as t with zeros where t < t0
%   and ones where t >= t0.

%   Based on MATLAB stepfun by:
%   J.N. Little 6-3-87
%   Copyright 1986-2002 The MathWorks, Inc. 

[m,n] = size(t);
y = zeros(m,n);
i = find(t>=t0); % return indices of t values larger than t0

if isempty(i)
  % t never greater than or equal to t0
    return
end

i = i(1); % first t>=t0
% create ones in output
if m == 1
  % row vector
    y(i:n) = ones(1,n-i+1);
else
  % column vector
    y(i:m) = ones(m-i+1,1);
end
