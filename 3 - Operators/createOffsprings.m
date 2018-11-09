%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Create the Offspring Population                     %
%                   by using the Mating Pool Population                   %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 19th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer 

function [offspringPopulation] = createOffsprings(matingPool, currentGeneration, GAParameters, testFunctionParameters)

% Initialize the offspring population
offspringPopulation = zeros(GAParameters.offspringSize, testFunctionParameters.dim);

% Generate a random integer permutation
randomIntegerPermutation = randperm(GAParameters.matingPoolSize)';

% Create the offspring population from the mating pool
for i = 1:2:GAParameters.offspringSize
    
    parents(1, :) = matingPool(randomIntegerPermutation(i, 1), :);
    parents(2, :) = matingPool(randomIntegerPermutation(i+1, 1), :);

    % Test the probability if two parents cross over
    if rand < GAParameters.crossoverRate
        
        offsprings = choiceCrossoverOperator(parents, GAParameters, testFunctionParameters);
        
    else
        
        switch GAParameters.chromosomeRepresentation
            
            case 'Binary Code'
                
                parentsInBinary = convertDecimalToBinary(parents, 2, GAParameters, testFunctionParameters);
                offsprings = parentsInBinary;
                
            case 'Real Number'
                
                offsprings = parents;
                
        end
        
    end
    
    offsprings = choiceMutationOperator(offsprings, currentGeneration, GAParameters, testFunctionParameters);
    
    offspringPopulation(i:i+1, 1:testFunctionParameters.dim) = offsprings(1:2, 1:testFunctionParameters.dim);
    
    % Mettre une fonction de verification, si la solution est viable ou pas
    
end
