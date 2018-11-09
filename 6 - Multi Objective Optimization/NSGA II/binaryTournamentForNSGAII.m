%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      Binary Tournament For NSGA II                      %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 28th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan (2002). 
%             A Fast and Elitist Multiobjective Genetic Algorithm : NSGA II


function [matingPool, fitnessValues] = binaryTournamentForNSGAII(rank, crowdingDistance, population, fitnessValues, NSGAIIParameters, testFunctionParameters)

% Initialize the variables of the mating pool.
matingPool = zeros(NSGAIIParameters.matingPoolSize, testFunctionParameters.dim);
fitnessValues.matingPool = zeros(NSGAIIParameters.matingPoolSize, testFunctionParameters.objectiveNumber);

% Compute the size (the number of individuals) of the overall population.
[individualNumber,~] = size(population);

for i = 1:1:NSGAIIParameters.matingPoolSize
    
    % Generate 2 unique integers selected randomly (binary tournament) from 
    % 1 to the overall population size inclusive.
    indices = randperm(individualNumber, 2);
    
    % Determine the stronger individual according to their rank calculated 
    % from the NSGA II algorithm.
    if rank(indices(1,1)) < rank(indices(1,2))
        
        matingPool(i,:) = population(indices(1,1),:);
        fitnessValues.matingPool(i,:) = fitnessValues.population(indices(1,1),:);
        
    elseif rank(indices(1,1)) > rank(indices(1,2))
        
        matingPool(i,:) = population(indices(1,2),:);
        fitnessValues.matingPool(i,:) = fitnessValues.population(indices(1,2),:);
        
    else
        
        if crowdingDistance(indices(1,1)) < crowdingDistance(indices(1,2))
        
            matingPool(i,:) = population(indices(1,2),:);
            fitnessValues.matingPool(i,:) = fitnessValues.population(indices(1,2),:);
            
        elseif crowdingDistance(indices(1,1)) > crowdingDistance(indices(1,2))
            
            matingPool(i,:) = population(indices(1,1),:);
            fitnessValues.matingPool(i,:) = fitnessValues.population(indices(1,1),:);
        
        end
        
    end
     
end
