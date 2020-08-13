function agc_indx
%AGC_INDX creates indices and data structures for AGC
% AGC_INDX creates indices and data structures AGC
%
% Syntax: agc_indx
%
%   NOTES:
%
%   Input:
%   VOID
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/21/20    15:58   Thad Haines     Version 1.0.0

global g

if ~isempty(g.agc.agc)
    g.agc.n_agc = length(g.agc.agc);
    
    % for each area...
    for ndx = 1:g.agc.n_agc
        % Determine if agc area is valid
        if isempty( intersect(g.area.area_def, g.agc.agc(ndx).area))
            error(['Area ', int2str(g.agc.agc(ndx).area), ' not found'])
        end
        
        % count number of ctrl gens
        g.agc.agc(ndx).n_ctrlGen = size(g.agc.agc(ndx).ctrlGen_con,1);
        
        g.agc.agc(ndx).macBusNdx = [];
        g.agc.agc(ndx).tgNdx = [];
        g.agc.agc(ndx).ctrlGen =[];
        
        % for each areas controlled generator...
        for gNdx =1:g.agc.agc(ndx).n_ctrlGen
            % collect ctrl gen machine index
            g.agc.agc(ndx).macBusNdx = [ g.agc.agc(ndx).macBusNdx , ...
                find(g.mac.mac_con(:,2) == ...
                g.agc.agc(ndx).ctrlGen_con(gNdx,1))];
            % collect turbine governor index
            g.agc.agc(ndx).tgNdx = [ g.agc.agc(ndx).tgNdx, ...
                find(g.tg.tg_con(:,2) == ...
                g.agc.agc(ndx).macBusNdx(gNdx))];
            
            % put bus number into ctrlGen struct
            g.agc.agc(ndx).ctrlGen(gNdx).busNum = g.mac.mac_con( ...
                g.agc.agc(ndx).macBusNdx(gNdx), 2);
            % put participation factor into ctrlGen Struct
            g.agc.agc(ndx).ctrlGen(gNdx).pF = g.agc.agc(ndx).ctrlGen_con(gNdx, 2);
        end
    end
    
else
    g.agc.n_agc = 0;
end
end

