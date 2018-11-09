%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Uniform Crossover                            %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : Ocotober 20th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% References : Introduction to Evolutionary Algorithms
%              Xinjie Yu && Mitsuo Gen - Springer
%
%              Uniform Crossover in Genetic Algorithms
%              	Gilbert Syswerda	


function [offsprings] = uniformCrossover(parents, GAParameters, testFunctionParameters)

%% Phase 1 : Convert decimal to binary

parentsInBinary = convertDecimalToBinary(parents, 2, GAParameters, testFunctionParameters);

%% Phase 2 : Perform the uniform crossover

for i = 1:1:testFunctionParameters.dim
    
    for j = 1:1:GAParameters.numberOfBits

        if rand < 0.5
        
            offsprings{1, 1}(i, j) = parentsInBinary{1, 1}(i, j);
            offsprings{2, 1}(i, j) = parentsInBinary{2, 1}(i, j);
    
        else
            
            offsprings{1, 1}(i, j) = parentsInBinary{2, 1}(i, j);
            offsprings{2, 1}(i, j) = parentsInBinary{1, 1}(i, j);
        end
        
    end
    
end