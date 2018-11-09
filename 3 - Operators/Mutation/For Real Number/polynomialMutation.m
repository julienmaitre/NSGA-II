%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Boundary Mutation                            %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 24th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : A Fast Elitist Multiobjective Genetic Algorithm : NSGA II
%             Aravind Seshadri
%
%             Analysing Mutation Schemes for Real Parameter Genetic
%             Algorithms
%             Kalyanmoy Deb & Debayan Deb


function offsprings = polynomialMutation(offsprings, GAParameters, testFunctionParameters)

% Initialize and set variable
offspringNumber = size(offsprings, 1); % Compute the number of offsprings

for i = 1:1:offspringNumber
    
    for j = 1:1:testFunctionParameters.dim
            
        if rand < GAParameters.mutationRate
            
            u = rand;
    
            if u < 0.5
        
                Delta = ((2*u)^(1/(GAParameters.distributionIndexForPMX+1))) - 1;
                offsprings(i,j) = offsprings(i,j) + Delta*(testFunctionParameters.upperLimit(1,j) - testFunctionParameters.lowerLimit(1,j));
                
            else
        
                Delta = 1-(2*(1-u))^(1/(GAParameters.distributionIndexForPMX+1));
                offsprings(i,j) = offsprings(i,j) + Delta*(testFunctionParameters.upperLimit(1,j) - testFunctionParameters.lowerLimit(1,j));
                
            end

        end
        
    end
    
end