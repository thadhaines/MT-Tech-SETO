function lmon(k)
%LMON Calculation and logging of line flows at index k
% LMON Calculation and logging of line flows at index k
%
% Syntax: LMON
%
%   NOTES:  Could/should be expanded to calculate line flows for area interchanges...
% 
%   Input: 
%   k - data index
%
%   Output: 
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/17/20    07:59   Thad Haines     Version 1

global g

% collect values for calculation
V1 = g.bus.bus_v(g.bus.bus_int(g.lmon.busFromTo(:,1)), k);
V2 = g.bus.bus_v(g.bus.bus_int(g.lmon.busFromTo(:,2)), k);
R = g.line.line(g.lmon.lmon_con,3); 
X = g.line.line(g.lmon.lmon_con,4); 
B = g.line.line(g.lmon.lmon_con,5);
tap = g.line.line(g.lmon.lmon_con,6); 
phi = g.line.line(g.lmon.lmon_con,7);
            
jay = sqrt(-1);
nline = size(V1,1);

for i = 1:nline
  if tap(i) == 0
    tap(i) = 1;
  end
end

% copied from line_pq
tps = tap.*exp(jay*phi*pi/180);
tpsi = diag(ones(nline,1)./tps);
tps = diag(tps);
z = R + jay*X;
y = diag(ones(nline,1)./z);
chg = diag(jay*B/2);

cur1 = tps*(y*(tpsi*V1-V2) + chg*V1); % iFrom
cur2 = y*(V2 - tpsi*V1) + chg*V2; % iTo
S1 = V1.*conj(cur1); % sFrom
S2 = V2.*conj(cur2); % sTo

% may not be ideal data structure for this...
for Lndx=1:g.lmon.n_lmon
    g.lmon.line(Lndx).iFrom(k) = cur1(Lndx);
    g.lmon.line(Lndx).iTo(k) = cur2(Lndx);
    g.lmon.line(Lndx).sFrom(k) = S1(Lndx);
    g.lmon.line(Lndx).sTo(k) = S2(Lndx);
end    

end

