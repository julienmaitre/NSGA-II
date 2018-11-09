%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Convert decimal to binary number                    %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 18th 2017                                                %
% Version : 1.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [populationInBinary] = convertDecimalToBinary(populationInDecimal, individualNumber, GAParameters, testFunctionParameters)

stepForDivision = (testFunctionParameters.upperLimit - testFunctionParameters.lowerLimit)/ ...
    (2^GAParameters.numberOfBits);
populationInBinary = cell(individualNumber, 1);

for i = 1:1:individualNumber
    
    populationInBinary{i,1} = de2bi(round((populationInDecimal(i,1:testFunctionParameters.dim) - testFunctionParameters.lowerLimit)./stepForDivision), GAParameters.numberOfBits);
    
end
