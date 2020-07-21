function area_indx
%AREA_INDX creates indices and data structures for areas monitoring
% AREA_INDX creates indices and data structures for areas monitoring
%
% Syntax: area_indx
%
%   NOTES:  Doesn't assume/require same order of area and bus array
%
%   Input:
%   VOID
%
%   Output:
%   VOID
%
%   History:
%   Date        Time    Engineer        Description
%   07/18/20    05:07   Thad Haines     Version 1.0.0
%   07/19/20    15:18   Thad Haines     Version 1.0.1 - minor additions to indices, area counts
%   07/20/20    19:42   Thad Haines     Version 1.0.2 - added max capacity to area init

global g

if ~isempty(g.area.area_def)
    areaNums = unique(g.area.area_def(:,2)); % collect unique area numbesr
    g.area.n_area = length(areaNums);
    macBus = g.mac.mac_con(:,2); % used for indexing machine_con
    
    % used for indexing the bus array
    genBus = g.bus.bus((g.bus.bus(:,10)==2) ,1);
    slackBus = g.bus.bus((g.bus.bus(:,10)==1) ,1);
    if ~isempty(slackBus)   
        genBus = [genBus; slackBus];
    end
    loadBus = g.bus.bus((g.bus.bus(:,10)==3),1);
    
    for areaN = 1:g.area.n_area
        g.area.area(areaN).number = areaNums(areaN);
        % collect area bus numbers
        g.area.area(areaN).areaBuses = g.area.area_def((g.area.area_def(:,2) == areaNums(areaN)),1);
        
        % doesn't assume/require same order of area and bus array
        
        % collect area machine bus numbers
        g.area.area(areaN).macBus = intersect(macBus, g.area.area(areaN).areaBuses);
        
        % collect area mac_con index numbers
        g.area.area(areaN).macBusNdx = [];
        for mNdx=1:length(g.area.area(areaN).macBus)
            % add index to global
            g.area.area(areaN).macBusNdx = [ g.area.area(areaN).macBusNdx ...
                , find(g.mac.mac_con(:,2) == g.area.area(areaN).macBus(mNdx))]; % for references to machine array values
        end
        % calculate area max capacity
        g.area.area(areaN).maxCapacity = sum(g.mac.mac_con(g.area.area(areaN).macBusNdx, 3));
        
        % collect area load bus numbers
        g.area.area(areaN).loadBus = intersect(loadBus, g.area.area(areaN).areaBuses);
        
        % collect area load bus index numbers
        g.area.area(areaN).loadBusNdx = []; % for references to bus array values
        for ndx=1:length(g.area.area(areaN).loadBus)
            % add index to global
            g.area.area(areaN).loadBusNdx = [ g.area.area(areaN).loadBusNdx ...
                , find(g.bus.bus(:,1) == g.area.area(areaN).loadBus(ndx))]; % for references to machine array values
        end
        
        % collect area generator bus numbers
        g.area.area(areaN).genBus = intersect(genBus, g.area.area(areaN).areaBuses);
        
        % collect area generator bus index numbers
        g.area.area(areaN).genBusNdx = []; % for references to bus array values
        for ndx=1:length(g.area.area(areaN).genBus)
            % add index to global
            g.area.area(areaN).genBusNdx = [ g.area.area(areaN).genBusNdx ...
                , find(g.bus.bus(:,1) == g.area.area(areaN).genBus(ndx))]; % for references to machine array values
        end
        
        % place holders for logged values
        g.area.area(areaN).totH = []; % to account for possible trips/changes in Inertia
        g.area.area(areaN).aveF = [];
        g.area.area(areaN).totGen = [];
        g.area.area(areaN).totLoad = [];
        g.area.area(areaN).icA = []; % Actual interchange
        g.area.area(areaN).icS = []; % Scheduled interchange (AGC?)
        
        % Create placeholders for line indices
        g.area.area(areaN).exportLineNdx = [];
        g.area.area(areaN).importLineNdx = [];
        
        % Find Interchange line indices
        for line=1:length(g.line.line)
            % From end in area
            if ~isempty(intersect(g.area.area(areaN).areaBuses, g.line.line(line,1)))
                % To End out of area
                if isempty(intersect(g.area.area(areaN).areaBuses, g.line.line(line,2)))
                    % add export line
                    g.area.area(areaN).exportLineNdx = ...
                        [g.area.area(areaN).exportLineNdx, line];
                end
                
            % To end in area
            elseif ~isempty(intersect(g.area.area(areaN).areaBuses, g.line.line(line,2)))
                % From end out of area
                if isempty(intersect(g.area.area(areaN).areaBuses, g.line.line(line,1)))
                    % add import line
                    g.area.area(areaN).importLineNdx = ...
                        [g.area.area(areaN).importLineNdx, line];
                end
            end
        end
        
        % count number of interchange lines
        g.area.area(areaN).n_export = length(g.area.area(areaN).exportLineNdx);
        g.area.area(areaN).n_import = length(g.area.area(areaN).importLineNdx);
    end
else
    g.area.n_area = 0;
end
end

