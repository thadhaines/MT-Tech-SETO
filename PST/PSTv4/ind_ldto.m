function ind_ldto(i,k)
%IND_LDTO is a template for Induction Motor Load Torque Calculation
% IND_LDTO is a template for Induction Motor Load Torque Calculation as
% a function of slip
%
% Syntax: ind_ldto(i,k)
%
%   NOTES:  Format for motor load data - mld_con
%           1 motor number
%           2 bus number
%           3 stiction load pu on motor base (f1)
%           4 stiction load coefficient (i1)
%           5 external load  pu on motor base(f2)
%           6 external load coefficient (i2)
%
%           load has the form
%           tload = f1*slip^i1 + f2*(1-slip)^i2
%
%   Input:
%   i - is the motor number
%          - 0 for vectorized computation
%   k - integer time (data index)
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   11/xx/95    xx:xx   Graham Rogers   Version 1.0
%   (c) copyright Joe Chow 1995-1996 All rights reserved
%   07/13/20    11:55   Thad Haines     Revised format of globals and internal function documentation

global g

if i == 0  % vector computation
    numot = length(g.ind.mld_con(:,2));
    g.ind.tload(:,k)=g.ind.mld_con(:,3).*(g.ind.slip(:,k).^g.ind.mld_con(:,4)) ...
        +g.ind.mld_con(:,5).*(ones(numot,1)-g.ind.slip(:,k)).^g.ind.mld_con(:,6);
else
    g.ind.tload(i,k)=g.ind.mld_con(i,3)*(g.ind.slip(i,k)^g.ind.mld_con(i,4)) ...
        + g.ind.mld_con(i,5)*(1-g.ind.slip(i,k))^g.ind.mld_con(i,6);
end

