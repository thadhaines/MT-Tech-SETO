function [ outputVector ] = cleanVector( inputVector, unwanted )
%CLEANVECTOR remove unwanted values from inputVector
%   Detailed explanation goes here
%
outputVector = inputVector;

if ~isempty(unwanted)
    for n = 1:max(size(unwanted))
        outputVector = outputVector(outputVector~=unwanted(n));
    end
end

end

