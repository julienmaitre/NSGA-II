%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Crowding Distance Assignment for NSGA II                %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 28th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan (2002). 
%             A Fast and Elitist Multiobjective Genetic Algorithm : NSGA II

function [distance] = crowdingDistanceAssignment(Front, fitnessValues, testFunctionParameters)

%% Step 1 : Compute the crowding distance for each individual in each front

% Number of Fronts
frontsNumber = length(Front);

% Initialize the distance variable
distance = zeros(size(fitnessValues.overallPopulation,1), 1);

% For each Front i
for i = 1:1:frontsNumber-1

    % Number of individuals belonging to the Front i
    individualNumber = length(Front(i).f);
    
    % Initialize variables
    fmax = zeros(1, testFunctionParameters.objectiveNumber); % The maximum objective values in the Front i
    fmin = zeros(1, testFunctionParameters.objectiveNumber); % The minimum objective values in the Front i
    
    % Find the maximum and minimum values in the Front i for each objective
    for k = 1:1:testFunctionParameters.objectiveNumber
    
        fmax(1,k) = max(fitnessValues.overallPopulation(Front(i).f,k));
        fmin(1,k) = min(fitnessValues.overallPopulation(Front(i).f,k));
        
    end
    
    % Compute the crowding distance for each individual j in the Front i
    % according to each objective
    for k = 1:1:testFunctionParameters.objectiveNumber
        
        [~,indices] = sort(fitnessValues.overallPopulation(Front(i).f,k),'ascend');
        distance(Front(i).f(indices(1))) = Inf;
        distance(Front(i).f(indices(end))) = Inf;
        
        % Compute the crowding distance for each individual j in the Front 
        % i
        for j = 2:1:individualNumber-1
        
            distance(Front(i).f(indices(j))) = distance(Front(i).f(indices(j)))+...
                abs(fitnessValues.overallPopulation(Front(i).f(indices(j+1))) - fitnessValues.overallPopulation(Front(i).f(indices(j-1))))/(fmax(1,k)-fmin(1,k));
            
        end
        
    end
    
end
