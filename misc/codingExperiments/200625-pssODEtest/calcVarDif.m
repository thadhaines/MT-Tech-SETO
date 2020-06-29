function [ tF, dif ] = calcVarDif( t1,y1, t2, y2 )
%calcDif Calculates diffefenec between variable timestep data
%   t1 = time series info for data y1
%   y1 = data related to t1
%   t2 = time series infor for data y2
%   y2 = data related to t2

% Check input data to ensure row vectors
if size(t1,1)>size(t1,2)
    t1 = t1';
end
if size(t2,1)>size(t2,2)
    t2 = t2';
end
if size(y1,1)>size(y1,2)
    y1 = y1';
end
if size(y2,1)>size(y2,2)
    y2 = y2';
end

% find data with more data points (fixed) and identifiy new working variables
% Assumes variable time step data has less points
if length(t1) > length(t2)
    tV = t2;
    yV = y2;
    tF = t1;
    yF = y1;
else
    tV = t1;
    yV = y1;
    tF = t2;
    yF = y2;
end

dif = zeros(length(tF),1); % placeholder

vNdx = 1; % for iterating through variable time step data
for fNdx=1:length(tF)
    
    nextNdx = vNdx+1;
    if nextNdx > max(max(size(yV)))
        % avoid over adjusting
        nextNdx = max(max(size(yV)));
    end
    
    % Adjust variable index as required
     if tF(fNdx) >= tV(nextNdx)
        vNdx = nextNdx;
    end
    
    
    dif(fNdx) = yF(fNdx)-yV(vNdx);
end


end

