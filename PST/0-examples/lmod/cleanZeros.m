function [s clearedVars] = cleanZeros( s , clearedVars, varargin)
%CLEANZEROS crawl input structure and remove any all zero entries
%   Returns cleaned structure and clearedVars cell of removed zeros


%% clean struct
fNames = fieldnames(s); % returns field names as cell
% for each structured array
for sNdx = 1:size(s, 2)
    % for each field
    for fNdx = 1:size(fNames)
        if ~isempty( s(sNdx).(fNames{fNdx}) )
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
                        fprintf('cleared s(%d).%s\n', sNdx, fNames{fNdx})
                        clearedVars{end+1} = [fNames{fNdx}]; % add name to cell for reference
                    else
                        fprintf('cleared s(%d).%s\n', sNdx, [varargin{1},'.',fNames{fNdx}])
                        clearedVars{end+1} = [varargin{1},'.',fNames{fNdx}]; % add name to cell for reference
                    end
                end
            end
        end
    end
end

