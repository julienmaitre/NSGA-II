%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        Simulated Binary Crossover                       %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 21th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference: Introduction to Evolutionary Algorithms
%            Xinji Yu & Mitsuo Gen - Springer
%
%            https://www.slideshare.net/paskorn/simulated-binary-crossover-presentation
%
%            Self-Adaptive Genetic Algorithms with Simulated Binary
%            Crossover
%            Kalyanmoy Deb & Hans-Georg Beyer

function offsprings = simulatedBinaryCrossover(parents, GAParameters, testFunctionParameters)

% Initialize 2 offspring solutions
offsprings = zeros(2, testFunctionParameters.dim);

% Perform the Simulated Binary Crossover 
for i = 1:1:testFunctionParameters.dim
    
    u = rand;
    
    if u <= 0.5
        
        Beta = (2*u)^(1/(GAParameters.controlParameterSBX + 1));
        
    else
        
        Beta = (1/(2*(1-u)))^(1/(GAParameters.controlParameterSBX + 1));
        
    end
    
    x = 0.5*(parents(1,i) + parents(2,i));
    
    if parents(2,i) > parents(1,i)
    
        offsprings(1,i) = x - 0.5*Beta*(parents(2,i) - parents(1,i));
        offsprings(2,i) = x + 0.5*Beta*(parents(2,i) - parents(1,i));
        
    else
    
        offsprings(1,i) = x - 0.5*(parents(1,i) - parents(2,i));
        offsprings(2,i) = x + 0.5*(parents(1,i) - parents(2,i));
        
    end

end
