%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Single-Point Crossover or Simple Crossover               %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : Ocotober 19th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% References : Introduction to Evolutionary Algorithms
%              Xinjie Yu && Mitsuo Gen - Springer
%
%              Genetic Algorithms in Search, Optimization & Machine Learning
%              Addison-Wesley - 1989

function [offsprings] = singlePointCrossover(parents, GAParameters, testFunctionParameters)

%% Phase 1 : Convert decimal to binary

parentsInBinary = convertDecimalToBinary(parents, 2, GAParameters, testFunctionParameters);

%% Phase 2 : Perform the single-point crossover

for i = 1:1:testFunctionParameters.dim
    
    % Phase 2.1 : Assign randomly a point between 1 and the (number of genes - 1)
    pointForCrossover = randi(GAParameters.numberOfBits - 1);
    
    % Phase 2.2 : Perform the crossover between parents
    offsprings{1, 1}(i, 1:GAParameters.numberOfBits) = [parentsInBinary{1, 1}(i, 1:pointForCrossover) parentsInBinary{2, 1}(i, pointForCrossover+1:GAParameters.numberOfBits)];
    offsprings{2, 1}(i, 1:GAParameters.numberOfBits) = [parentsInBinary{2, 1}(i, 1:pointForCrossover) parentsInBinary{1, 1}(i, pointForCrossover+1:GAParameters.numberOfBits)];
    
end

