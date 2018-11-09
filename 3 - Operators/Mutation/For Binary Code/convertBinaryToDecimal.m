%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Convert binary to decimal number                    %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 18th 2017                                                %
% Version : 1                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [populationInDecimal] = convertBinaryToDecimal(populationInBinary, individualNumber, GAParameters, testFunctionParameters)

populationInDecimal = zeros(individualNumber, testFunctionParameters.dim);

for i = 1:1:individualNumber
    
    for j = 1:1:testFunctionParameters.dim
    
        [~,Column] = find(populationInBinary{i, 1}(j, 1:GAParameters.numberOfBits));
        Locus = Column-1;
        populationInDecimal(i,j) = (sum(power(2,Locus))/(2^GAParameters.numberOfBits))* ...
            (testFunctionParameters.upperLimit(1,j) - testFunctionParameters.lowerLimit(1,j)) + testFunctionParameters.lowerLimit(1,j);
        clearvars Column Locus
        
    end
end

