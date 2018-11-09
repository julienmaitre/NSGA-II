%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      Fill the Parent Population According to the NSGA II Algorithm      %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 30th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan (2002). 
%             A Fast and Elitist Multiobjective Genetic Algorithm : NSGA II


function [population, fitnessValues, newRank] = fillParentPopulation(Front, crowdingDistance, overallPopulation, fitnessValues, NSGAIIParameters, testFunctionParameters)

% Initialize variables
front = 1;
population = zeros(0,testFunctionParameters.dim);
objectiveValues = zeros(0,testFunctionParameters.objectiveNumber);
newRank = zeros(0,1);

% Fill the parent population of the next generation
while size(population,1) + length(Front(front).f) < NSGAIIParameters.populationSize
    
    population = [population ; overallPopulation(Front(front).f,:)];
    objectiveValues = [objectiveValues ; fitnessValues.overallPopulation(Front(front).f,:)];
    newRank = [newRank ; ones(1,length(Front(front).f))'*front];
    front = front +1;
    
end

currentNumberOfIndividuals = size(population,1);

if currentNumberOfIndividuals < NSGAIIParameters.populationSize
    
    [~,indices] = sort(crowdingDistance(Front(front).f),'descend');
    j = 0;
    
    for i = currentNumberOfIndividuals + 1:1:NSGAIIParameters.populationSize
        
        j = j+1;
        population(i,:) = overallPopulation(Front(front).f(indices(j)),:);
        objectiveValues(i,:) = fitnessValues.overallPopulation(Front(front).f(indices(j)),:);
        newRank (i,1) = front;
        
    end
    
end

fitnessValues.population = objectiveValues;

