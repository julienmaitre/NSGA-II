%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Create the Offspring Population                     %
%             by using the Mating Pool Population for NSGA II             %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 30th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : https://books.google.ca/books?id=uJ8NDgAAQBAJ&pg=PA70&lpg=PA70&dq=decision+variable+with+bounds+SBX&source=bl&ots=CS5_7VvDYN&sig=UWvdAykiERrj6-UwfayduqtQNaU&hl=fr&sa=X&ved=0ahUKEwje4IfM7ebXAhWh4IMKHcLEB80Q6AEISjAF#v=onepage&q=decision%20variable%20with%20bounds%20SBX&f=false


function [offspringPopulation] = createOffspringsForNSGAII(matingPool, currentGeneration, NSGAIIParameters, testFunctionParameters)

% Initialize variable
individualNumber = 0;

% Initialize the offspring population
Q = zeros(NSGAIIParameters.populationSize, testFunctionParameters.dim);

% Create the offspring population from the mating pool
while individualNumber < NSGAIIParameters.populationSize
  
    % Generate 2 random integers permutation
    indices = randperm(NSGAIIParameters.matingPoolSize, 2);
        
    % Get parent individuals according to the 2 random integers 
    parents(1,1:testFunctionParameters.dim) = matingPool(indices(1,1),:);
    parents(2,1:testFunctionParameters.dim) = matingPool(indices(1,2),:);  
    
    % Test the probability if two parents cross over
    if rand < NSGAIIParameters.crossoverRate
        
        offsprings = choiceCrossoverOperator(parents, NSGAIIParameters, testFunctionParameters);
        
    else
        
        switch NSGAIIParameters.chromosomeRepresentation
            
            case 'Binary Code'
                
                parentsInBinary = convertDecimalToBinary(parents, 2, NSGAIIParameters);
                offsprings = parentsInBinary;
                
            case 'Real Number'
                
                offsprings = parents;
                
        end
        
    end
    
    % Perform the mutation of the two new offspring
    offsprings = choiceMutationOperator(offsprings, currentGeneration, NSGAIIParameters, testFunctionParameters);
    
    
    for i = 1:1:testFunctionParameters.dim
        
        % Verify if the new solutions are within the boundaries of the search
        % space.   
        if offsprings(1,i) > testFunctionParameters.upperLimit(1,i)
        
            offsprings(1,i) = testFunctionParameters.upperLimit(1,i);
        
        elseif offsprings(1,i) < testFunctionParameters.lowerLimit(1,i)
        
            offsprings(1,i) = testFunctionParameters.lowerLimit(1,i);
        
        end
    
        if offsprings(2,i) > testFunctionParameters.upperLimit(1,i)
        
            offsprings(2,i) = testFunctionParameters.upperLimit(1,i);
        
        elseif offsprings(2,i) < testFunctionParameters.lowerLimit(1,i)
        
            offsprings(2,i) = testFunctionParameters.lowerLimit(1,i);
        
        end
        
    end
    
    % Assign the two new offspring mutated in the offspring population
    Q(individualNumber+1:individualNumber+2,1:testFunctionParameters.dim) = offsprings(1:2,1:testFunctionParameters.dim);
    individualNumber = individualNumber + 2;    

end

offspringPopulation = Q(1:NSGAIIParameters.populationSize,:);

