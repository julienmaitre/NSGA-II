%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   Fast Non Dominated Sort for NSGA II                   %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 27th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan (2002). 
%             A Fast and Elitist Multiobjective Genetic Algorithm : NSGA II


function [rank, Front, overallPopulation, fitnessValues] = fastNonDominatedSort_V2(population, offspringPopulation, fitnessValues, testFunctionParameters)

    %% Step 1 : Initialize

    % Initialize variables of Fronts
    front = 1; % Index of the current Front
    Front(front).f = [];

    % Combine parent and offspring population
    tmpOverallPopulation = [population ; offspringPopulation];
    tmpFitnessValues = [fitnessValues.population ; fitnessValues.offspringPopulation];

    % Delete doublons in the population 
    [tmpOverallPopulation, indices, ~] = unique(tmpOverallPopulation,'rows','stable');
    fitnessValues.tmpOverallPopulation = tmpFitnessValues(indices, 1:testFunctionParameters.objectiveNumber);
    
    % Delete individuals in the population with NaN objective values
    indices = find(ismember(isnan(fitnessValues.tmpOverallPopulation), zeros(1,testFunctionParameters.objectiveNumber), 'rows') == 1);
    overallPopulation = tmpOverallPopulation(indices,:);
    fitnessValues.overallPopulation = fitnessValues.tmpOverallPopulation(indices,:);
    
    % Compute the constrained values if the test problem is constrained
    switch testFunctionParameters.contrainedProblem
        
        case 'Contrained'

            % Compute the constraines
            constrainedValues = testFunctionParameters.cobj(overallPopulation, testFunctionParameters);
            
            % Delete individuals in the population with NaN constrained
            % values
            indices = find(ismember(isnan(constrainedValues), zeros(1,testFunctionParameters.constrainedNumber), 'rows') == 1);
            overallPopulation = overallPopulation(indices,:);
            fitnessValues.overallPopulation = fitnessValues.overallPopulation(indices,:);
            constrainedValues = constrainedValues(indices,:);
            
            % Compute the size of the overall Population (number of 
            % individuals)
            individualsNumber = size(overallPopulation, 1);
    
            % Initialize the rank variable 
            rank = zeros(individualsNumber, 1);
            
        otherwise
            
            constrainedValues = [];
            
            % Compute the size of the overall Population (number of 
            % individuals)
            individualsNumber = size(overallPopulation, 1);
    
            % Initialize the rank variable 
            rank = zeros(individualsNumber, 1);
            
    end

    %% Step 2 : Find the first front according to the Pareto dominance

    for i = 1:1:individualsNumber

        % Initialize variables
        individual(i).np = 0; % Number of individuals that dominate this individual
        individual(i).Sp = []; % Indices of individuals that this individual dominate

        for j = 1:1:individualsNumber

            if i ~= j

                % Initialize variables
                dominatingCounter = 0; % Counter of the number of objectives that are dominating
                dominatedCounter = 0; % Counter of the number of objectives that are dominated
                equalCounter = 0; % Counter of the number of objectives that are equal

                % Test the pareto dominance between an individual i and an 
                % individual j according to their number of objectives and the
                % number of contraines violated
                [individual] = paretoDominanceTest(dominatedCounter, dominatingCounter, equalCounter, i, j, fitnessValues, testFunctionParameters, individual, constrainedValues);

            end

        end

        % If the number of individuals that dominate the individual i is equal
        % to zero, i belongs to the first front
        if individual(i).np == 0

            rank(i) = 1;
            Front(front).f = [Front(front).f i];

        end

    end

    %% Step : 3 Find the subsequent fronts

    while ~isempty(Front(front).f) % While Front(front).f is not empty

        % Initialize a temporary variable used to store the members of the next 
        % front
        Q = [];

        % For each individual i of Front(front).f
        for i = 1:1:length(Front(front).f)

            % If Sp of i is not empty
            if ~isempty(individual(Front(front).f(i)).Sp)

                % For each individual j of Sp of i
                for j = 1:1:length(individual(Front(front).f(i)).Sp)

                    % nj = nj-1
                    individual(individual(Front(front).f(i)).Sp(j)).np = individual(individual(Front(front).f(i)).Sp(j)).np - 1;

                    % If nj of individual j is equal to zero
                    if individual(individual(Front(front).f(i)).Sp(j)).np == 0

                        % Individual j belongs to the next front
                        rank(individual(Front(front).f(i)).Sp(j)) = front+1;
                        Q = [Q individual(Front(front).f(i)).Sp(j)];

                    end

                end

            end

        end

        front = front+1;
        Front(front).f = Q;

    end
    
end
    
%% Pareto dominance test according the constraines and then to the objectives

function [individual] = paretoDominanceTest(dominatedCounter, dominatingCounter, equalCounter, i, j, fitnessValues, testFunctionParameters, individual, constrainedValues)

    switch testFunctionParameters.contrainedProblem
        
        case 'Not contrained'
            
            % Pareto dominance according only to the objectives
            [dominatedCounter, dominatingCounter, equalCounter] = paretoDominance(i, j, fitnessValues, testFunctionParameters, dominatedCounter, dominatingCounter, equalCounter);

            % Store the information about the Pareto dominance if i 
            % dominates j
            if dominatedCounter == 0 && equalCounter ~= testFunctionParameters.objectiveNumber
                    
                % Add j to the set of solutions dominated by i
                individual(i).Sp = [individual(i).Sp j];

            % If j dominates i    
            elseif dominatingCounter == 0 && equalCounter ~= testFunctionParameters.objectiveNumber

                % Increment the domination counter of i
                individual(i).np = individual(i).np + 1;

            end
                
        case 'Contrained'
            
            for k = 1:1:testFunctionParameters.constrainedNumber
                
                switch testFunctionParameters.constrainedType{1,k}
                    
                    case 'Superior or equal'

                        if constrainedValues(i,k)>=0
                        
                            constrainedValues(i,k) = 0;
                            
                        end
                        
                        if constrainedValues(j,:)>=0
                            
                            constrainedValues(j,k) = 0;
                            
                        end
                        
                    case 'Inferior or equal'
                        
                        if constrainedValues(i,k)<=0
                        
                            constrainedValues(i,k) = 0;
                            
                        end
                        
                        if constrainedValues(j,:)<=0
                            
                            constrainedValues(j,k) = 0;
                            
                        end
                        
                end
                
            end

            if sum(abs(constrainedValues(i,:))) < sum(abs(constrainedValues(j,:)))
                
                % Add j to the set of solutions dominated by i
                individual(i).Sp = [individual(i).Sp j];

            elseif sum(abs(constrainedValues(i,:))) > sum(abs(constrainedValues(j,:)))
                
                % Increment the domination counter of i
                individual(i).np = individual(i).np + 1;

            else
                
                [dominatedCounter, dominatingCounter, equalCounter] = paretoDominance(i, j, fitnessValues, testFunctionParameters, dominatedCounter, dominatingCounter, equalCounter);

                % Store the information about the Pareto dominance if i 
                % dominates j
                if dominatedCounter == 0 && equalCounter ~= testFunctionParameters.objectiveNumber
                    
                    % Add j to the set of solutions dominated by i
                    individual(i).Sp = [individual(i).Sp j];

                % If j dominates i    
                elseif dominatingCounter == 0 && equalCounter ~= testFunctionParameters.objectiveNumber

                    % Increment the domination counter of i
                    individual(i).np = individual(i).np + 1;

                end
                
            end
            
    end

end

%% Pareto dominance according only to the objectives

function [dominatedCounter, dominatingCounter, equalCounter] = paretoDominance(i, j, fitnessValues, testFunctionParameters, dominatedCounter, dominatingCounter, equalCounter)

        for k = 1:1:testFunctionParameters.objectiveNumber
        
            switch testFunctionParameters.optimizationType{1,k}
            
                case 'Maximization'
                
                    if fitnessValues.overallPopulation(i,k) > fitnessValues.overallPopulation(j,k)
                    
                        dominatingCounter = dominatingCounter+1;
                            
                    elseif fitnessValues.overallPopulation(i,k) < fitnessValues.overallPopulation(j,k)
                            
                        dominatedCounter = dominatedCounter+1;
                            
                    else
                            
                        equalCounter = equalCounter+1;
                            
                    end
                
                case 'Minimization'
                        
                    if fitnessValues.overallPopulation(i,k) < fitnessValues.overallPopulation(j,k)
                            
                        dominatingCounter = dominatingCounter+1;
                            
                    elseif fitnessValues.overallPopulation(i,k) > fitnessValues.overallPopulation(j,k)
                            
                        dominatedCounter = dominatedCounter+1;
                            
                    else
                            
                        equalCounter = equalCounter+1;
                            
                    end
                
            end
        
        end
    
end
