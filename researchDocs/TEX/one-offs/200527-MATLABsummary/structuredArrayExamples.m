clear; clc; close all
%% Cell Creation
a{1,1} = rand(2,3)
a{2,1} = magic(4)
a{2,2} = 'This is a cell'

%% Structure Creation
s.a = a
s.t = 2
s.Name = 'First'

%% Addition to Structure
s(2).Name = 'Second'
% not each field is required to have the same type of data
s(4).a = 'This is a string' 
% structures do not have to be created sequentially
s(5).Name = 'Fifth' 
% definition of field via setfield function
setfield(s(5), 't', 56) 

%% Simple Dynamic Field Name Example
s(2).Name % display current data
s(2).('Name') % same as above, but using dynamic field format
f = 'Name'; % Using dynamic field names with variables is more useful
s(2).(f) %

%% Some structure manipulations
fNames = fieldnames(s); % returns field names as cell
% for each structured array
for sNdx = 1:size(s, 2)
    % Introduce each strucutred array
    fprintf('*** s.(%d) has Fields:\n', sNdx)
    % display all fields
    for fNdx = 1:size(fNames)
        fprintf('%s: \n', fNames{fNdx} ) % display field name
        if isempty( getfield(s(sNdx), fNames{fNdx}) )
            % clarify existance of empty fields
            disp('   Empty') 
        else
            % display data of each field using dynamic field names
            disp( s(sNdx).(fNames{fNdx}) ) 
            % same as above but using getfield function
            disp( getfield(s(sNdx), fNames{fNdx}) ) 
        end
    end
end