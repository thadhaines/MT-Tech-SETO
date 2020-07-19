function calcAreaVals(k,flag)
%CALCAREAVALS calculates total area generation, load, and interchange.
% CALCAREAVALS calculates total area generation, load, and interchange.
%
% Syntax: calcAreaVals(k,flag)
%
%   NOTES:  Sign convention is to use positive numbers that correspond
%           to category. For example, the convention is believed to be:
%           Positive generation is power injection,
%           positive load is real and reactive absorbtion
%           positive B shunt is capacitive VAR injection
%           positive G shunt is inductive VAR injection
%
%   Input:
%   k - data index
%   flag  - dictate what operation to perform
%       0 - Handle future intialization operations?
%       1 - calculate area totals
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/19/20    08:50   Thad Haines     Version 1

global g

if flag == 0
    % This is included, but maybe not useful.
end

% bus array col definitions

% col4 p_gen(pu)
% col5 q_gen(pu)
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)

if flag == 1
    for n=1:g.area.n_area
        
        % collect and sum generator electric p and q values
        pGen = sum(g.mac.pelect(g.area.area(n).macBusNdx,k));
        qGen = sum(g.mac.qelect(g.area.area(n).macBusNdx,k));
        g.area.area(n).totGen(k) = pGen + 1j*qGen;
        
        % assumed the bus array is updated with new vals....
        % However, that is not the case.
        % a Y-matrix is manipulated for load changes instead
        % could reflect lmod, or rlmod changes by indexing the respective
        % _con array and adding the output signal to the total sum...
        % but it seems like that may not be super useful as non-conforming
        % load changes wouldn't be accounted for.
        % The solution may be to alter/add code to nc_load that updates values?
        
       % pLoad = sum(g.bus.bus(g.area.area(n).loadBusNdx,6));
       % qLoad = sum(g.bus.bus(g.area.area(n).loadBusNdx,7));
        
       % shuntG = sum(g.bus.bus(g.area.area(n).loadBusNdx,8));
       % shuntB = sum(g.bus.bus(g.area.area(n).loadBusNdx,9));
        
       % g.area.area(n).totLoad(k) = pLoad + 1j*(qLoad + shuntG - shuntB);

        % calculate actual interchange        
        if g.area.area(n).n_export~=0
            % calculate exported power
            [sExport, ~] = line_pq2(g.area.area(n).exportLineNdx,k);
        else
            sExport = 0;
        end
        
        if g.area.area(n).n_import~=0
            % calculate imported power
            [~, sImport] = line_pq2(g.area.area(n).importLineNdx,k);
        else
            sImport = 0;
        end
        
        g.area.area(n).icA(k) = sum(sExport) + sum(sImport); % ( Positive number represents overgenration / export )
        
    end
end% end function