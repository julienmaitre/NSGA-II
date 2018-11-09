%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           Non Uniform Mutation                          %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 21th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer


function offsprings = nonUniformMutation(offsprings, currentGeneration, GAParameters, testFunctionParameters)

% Initialize and set variable
offspringNumber = size(offsprings, 1); % Compute the number of offsprings

Exponent = ((1 - currentGeneration)/GAParameters.maximumGeneration)^ ...
    GAParameters.controlParameterNUMutation; % Compute the exponent of the non uniform mutation

for i = 1:1:offspringNumber

    for j = 1:1:testFunctionParameters.dim
        
        if rand < GAParameters.mutationRate
        
            if rand >= 0.5
            
                Delta = (testFunctionParameters.upperLimit(1,j) - offsprings(i,j))*(1 - rand^Exponent);
                offsprings(i,j) = offsprings(i,j) + Delta;
            
            else
            
                Delta = (offsprings(i,j) - testFunctionParameters.lowerLimit(1,j))*(1 - rand^Exponent);
                offsprings(i,j) = offsprings(i,j) - Delta;
                
            end
            
        end
        
    end
    
end
            