%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Uniform Mutation                            %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 21th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer


function offsprings = uniformMutation(offsprings, GAParameters, testFunctionParameters)

% Initialize and set variable
offspringNumber = size(offsprings, 1); % Compute the number of offsprings

% For each offspring
for i = 1:1:offspringNumber
    
    % For each parameter of the solution (gene)
    for j = 1:1:testFunctionParameters.dim
        
        if rand < GAParameters.mutationRate

            offsprings(i,j) = testFunctionParameters.lowerLimit(1,j) + (testFunctionParameters.upperLimit(1,j) - testFunctionParameters.lowerLimit(1,j))*rand;
            
        end
        
    end
    
end

