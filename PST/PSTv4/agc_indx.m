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
%   08/20/20    10:52   Thad Haines     Version 1.0.1 - added fields to calculate maxGen
%   08/25/20    12:06   Thad Haines     Version 1.0.2 - addition of warnings for empties and mis-matched machine arrays

%%
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
        g.agc.agc(ndx).maxGen = 0; % Maximum MW output of ctrl gens
        
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
            
            % store machine base
            g.agc.agc(ndx).ctrlGen(gNdx).mBase = g.mac.mac_con(g.agc.agc(ndx).macBusNdx(gNdx), 3);
            
            % calculate maximum MW output
            g.agc.agc(ndx).maxGen = g.agc.agc(ndx).maxGen + g.agc.agc(ndx).ctrlGen(gNdx).mBase;
            
            % notify of empties/mis-matches
            if isempty(g.agc.agc(ndx).macBusNdx)
                fprintf('*!* AGC %d has no indexed machines...\n', ndx)
            end
            if isempty(g.agc.agc(ndx).ctrlGen)
                fprintf('*!* AGC %d has no controlled machines...\n', ndx)
            end
            if size(g.agc.agc(ndx).macBusNdx) ~= size(g.agc.agc(ndx).tgNdx)
                fprintf('*!* AGC %d a size mis-match in number of machines and governors...\n', ndx)
            end
            
        end
    end
    
else
    g.agc.n_agc = 0;
end
end

