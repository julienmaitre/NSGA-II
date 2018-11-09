%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%        Initialization of the Population for the first generation        %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 18th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference: Introduction to Evolutionary Algorithms
%            Xinji Yu & Mitsuo Gen - Springer


function [population] = initializePopulation(EAParameters, testFunctionParameters)

% Initialize the variable of the initial population.
population = zeros(EAParameters.populationSize, testFunctionParameters.dim);

% Initialize the population of the generation 0 with uniformly distributed 
% random numbers between the upper limit and lower limit of each chromosome 
% (or variable).
% If the upper limit and lower limit are vectors
for i = 1:1:EAParameters.populationSize
  
    population(i,:) = testFunctionParameters.lowerLimit + ...
        (testFunctionParameters.upperLimit - testFunctionParameters.lowerLimit).*rand(1, testFunctionParameters.dim);
    
end

