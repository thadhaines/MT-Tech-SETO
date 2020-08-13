function [s clearedVars] = cleanZeros( s , clearedVars, varargin)
%CLEANZEROS crawl input structure and remove any all zero entries
%   Returns cleaned structure and clearedVars cell of removed zeros
%
% Syntax: [s clearedVars] = cleanZeros( s , clearedVars, varargin)
%
%   NOTES:  Is a recursive alogrithm
%
%   Input:
%   s - data index
%   clearVars - choose between operations
%   varargin{1} - used to handle field names upon recursion
%
%   Output:
%   s - cleaned structure
%   clearedVars - cell array of cleared variable names
%
%   History:
%   Date        Time    Engineer        Description
%   07/27/20    10:18   Thad Haines     Version 1 with updated documentation

debugOutput = 0; % prints names of cleared items

%% clean struct
fNames = fieldnames(s); % returns field names as cell
% for each structured array
for sNdx = 1:size(s, 2)
    % for each field
    for fNdx = 1:size(fNames)
        % skip cells and function handles
        if iscell( s(sNdx).(fNames{fNdx}) ) || isa( s(sNdx).(fNames{fNdx}), 'function_handle')
            continue
            
        elseif ~isempty( s(sNdx).(fNames{fNdx}) )
            % check if struct, if so recursive call with additional name vield
            if isstruct( s(sNdx).(fNames{fNdx}) )
                if nargin <= 2
                    [s(sNdx).(fNames{fNdx}), clearedVars] = ...
                        cleanZeros( s(sNdx).(fNames{fNdx}) , clearedVars, fNames{fNdx});
                else
                    % keep appending field names
                    [s(sNdx).(fNames{fNdx}), clearedVars] = ...
                        cleanZeros( s(sNdx).(fNames{fNdx}) , clearedVars, [varargin{1},'.',fNames{fNdx}]);
                end
            else
                % check if zeros
                zeroTest = all( s(sNdx).(fNames{fNdx})==0 ) ; % check if all zeros
                if all(zeroTest)
                    s(sNdx).(fNames{fNdx}) =[] ; % clear variable
                    if nargin <= 2
                        if debugOutput
                            fprintf('cleared struct(%d).%s\n', sNdx, fNames{fNdx})
                        end
                        clearedVars{end+1} = [fNames{fNdx}]; % add name to cell for reference
                    else
                        if debugOutput
                            fprintf('cleared struct(%d).%s\n', sNdx, [varargin{1},'.',fNames{fNdx}])
                        end
                        clearedVars{end+1} = [varargin{1},'.',fNames{fNdx}]; % add name to cell for reference
                    end
                end
            end
        end
    end
end

