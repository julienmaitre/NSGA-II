%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Boundary Mutation                            %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 21th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer


function offsprings = boundaryMutation(offsprings, GAParameters, testFunctionParameters)

% Initialize and set variable
offspringNumber = size(offsprings, 1); % Compute the number of offsprings

for i = 1:1:offspringNumber
    
    for j = 1:1:testFunctionParameters.dim
            
        if rand < GAParameters.mutationRate
            
            if rand <= 0.5
            
                offsprings(i,j) = testFunctionParameters.lowerLimit(1,j);
                
            else
                
                offsprings(i,j) = testFunctionParameters.upperLimit(1,j);
                
            end
            
        end
        
    end
    
end