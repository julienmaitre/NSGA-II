%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      A Fast and Elitist Multiobjective Genetic Algorithm - NSGA II      %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 25th 2017                                               %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan (2002). 
%             A Fast and Elitist Multiobjective Genetic Algorithm : NSGA II


function [population, objectiveValues] = NSGAII(functionName, NSGAIIParameters, testFunctionParameters)

    %% Step 1 : Initialization

    % Step 1.1 : Assign and set the parameter(s) for NSGA II
    algorithmUsed = 'NSGA II';
    currentGeneration = 0;
    disp(['Generation ' num2str(currentGeneration) '/' num2str(NSGAIIParameters.maximumGeneration)])

    % Step 1.2 : Generate the initial population P0 (Parent Population)
    [population] = initializePopulation(NSGAIIParameters, testFunctionParameters);

    % Step 1.3 : Create the empty offspring population and the empty 
    % fitness values of the offspring population
    offspringPopulation = [];
    fitnessValues.offspringPopulation = [];

    %% Step 2 : Fitness, Rank, Set of Solution in the Fronts and Crowding Distance Assignment

    % Step 2.1 : Calculate fitness values of individuals in Population P0
    % according to the problem to optimize
    [fitnessValues.population] = testFunctionParameters.fobj(population, testFunctionParameters);

    % Step 2.2 : Calculate the rank and front of individuals according to 
    % the fast-non-dominated-sort of NSGA II
    [rank, Front, overallPopulation, fitnessValues] = fastNonDominatedSort_V2(population, offspringPopulation, fitnessValues, testFunctionParameters);

    % Step 2.3 : Calculate the crowding distance for each individual 
    % according to the NSGA II algorithm
    [crowdingDistance] = crowdingDistanceAssignment(Front, fitnessValues, testFunctionParameters);

    %% Step 3 : Save Results of the Generation 0

    % Step 3.1 : Create the path name of the results location
    dateTime = clock;
    dateTime = strcat(' (',num2str(dateTime(1)),num2str(dateTime(2)),num2str(dateTime(3)), ...
        '-',num2str(dateTime(4)),'h',num2str(dateTime(5)),'min',num2str(floor(dateTime(6))),'sec)');
    pathName_0 = strcat(NSGAIIParameters.folderName,'\',functionName,'\',algorithmUsed,'\',NSGAIIParameters.testIndex,' - ',dateTime);

    % Step 3.2 : Create the path name and the directory file of the file 
    % results of the generation 0
    dirName = strcat('Generation_',num2str(currentGeneration));
    pathName = strcat(pathName_0,'\',dirName);
    mkdir(pathName)
    fileName = strcat('resultsGeneration_',num2str(currentGeneration),'.mat');
    pathToSave = strcat(pathName,'\',fileName);

    % Step 3.3 : Get the necessary variables representing the results of
    % the generation 0
    population = overallPopulation;
    fitnessValues.population = fitnessValues.overallPopulation;
    objectiveValues = fitnessValues.population;

    % Step 3.4 : Save the results of the generation 0 in the directory 
    % created with pathToSave variable
    save(pathToSave ,'rank','Front','crowdingDistance','population','objectiveValues')

    %% Step 4 : Plot the Pareto Front of the Generation 0

    % The results of the Pareto Front of the generation 0 can only be 
    % displayed for test problem functions with 2 or 3 objectives
    displayResults(testFunctionParameters.ParetoOptimalFront, objectiveValues, rank, functionName, algorithmUsed, currentGeneration, NSGAIIParameters, testFunctionParameters)
    
    %% Step 5 : Termination Test

    while currentGeneration < NSGAIIParameters.maximumGeneration % Step 5.1 : Check if stopping criterion is satisfied

        %% Step 6 : Fill the Mating Pool 

        % Step 6.1 : Fill the mating pool according to the binary 
        % tournament based on the crowded-comparison operator
        [matingPool, fitnessValues] = binaryTournamentForNSGAII(rank, crowdingDistance, population, fitnessValues, NSGAIIParameters, testFunctionParameters);

        %% Step 7 : Variation

        % Step 7.1 : Apply recombination and mutation operators to the
        % mating pool and set Q(t+1) to the resulting population.
        [offspringPopulation] = createOffspringsForNSGAII(matingPool, currentGeneration, NSGAIIParameters, testFunctionParameters);

        %% Step 8 : Fitness, Rank, Set of Solution in the Fronts and Crowding Distance Assignment

        % Step 8.1 : Calculate fitness values of individuals in Population 
        % P_currentGeneration according to the optimizing problem
        [fitnessValues.offspringPopulation] = testFunctionParameters.fobj(offspringPopulation, testFunctionParameters);

        % Step 8.2 : Calculate the rank and front of individuals according  
        % to the fast-non-dominated-sort of NSGA II
        [rank, Front, overallPopulation, fitnessValues] = fastNonDominatedSort_V2(population, offspringPopulation, fitnessValues, testFunctionParameters);

        % Step 8.3 : Calculate the crowding distance for each individual 
        % according to the NSGA II algorithm
        [crowdingDistance] = crowdingDistanceAssignment(Front, fitnessValues, testFunctionParameters);

        % Step 8.4 : Fill the population (Parent Population) P_t+1 of the 
        % next generation
        [population, fitnessValues, newRank] = fillParentPopulation(Front, crowdingDistance, overallPopulation, fitnessValues, NSGAIIParameters, testFunctionParameters);

        % Step 8.5 : Increment the generation
        currentGeneration = currentGeneration+1;
        disp(['Generation ' num2str(currentGeneration) '/' num2str(NSGAIIParameters.maximumGeneration)])

        %% Step 9 : Save Results of the Generation currentGeneration

        % Step 9.1 : Create the path name and the directory file of the 
        % file results of the generation 0
        dirName = strcat('Generation_',num2str(currentGeneration));
        pathName = strcat(pathName_0,'\',dirName);
        mkdir(pathName)
        fileName = strcat('ResultsGeneration_',num2str(currentGeneration),'.mat');
        pathToSave = strcat(pathName,'\',fileName);

        % Step 9.2 : Get the necessary variables representing the results 
        % of the generation 0
        objectiveValues = fitnessValues.population;
        rank = newRank;

        % Step 9.2 : Save the results of the generation 0 in the directory 
        % created with pathToSave variable
        save(pathToSave ,'rank','Front','crowdingDistance','population','objectiveValues')

        %% Step 10 : Plot the Pareto Front of the Generation currentGeneration

        % The results of the Pareto Front of the generation 
        % currentGeneration can only be displayed for test problem 
        % functions with 2 or 3 objectives
        displayResults(testFunctionParameters.ParetoOptimalFront, objectiveValues, rank, functionName, algorithmUsed, currentGeneration, NSGAIIParameters, testFunctionParameters)

    end
    
end

function displayResults(optimalParetoFront, objectiveValues, rank, functionName, algorithmUsed, currentGeneration, NSGAIIParameters, testFunctionParameters)

    if isempty(optimalParetoFront) == 0

        if testFunctionParameters.objectiveNumber == 2

            indices = find(rank == 1);
            data = objectiveValues(indices,:);
            plot(data(:,1),data(:,2),'*r')
            hold on
            plot(optimalParetoFront(:,1),optimalParetoFront(:,2))
            hold off
            title(['Optimization of ', functionName, ' at generation ', num2str(currentGeneration), '/', num2str(NSGAIIParameters.maximumGeneration)]);
            legend(algorithmUsed, 'Optimal Pareto Front');
            xlabel('Objective 1');
            ylabel('Objective 2');
            drawnow

        elseif testFunctionParameters.objectiveNumber == 3

            indices = find(rank == 1);
            data = objectiveValues(indices,:);
            plot(data(:,1),data(:,2),data(:,3),'*r')
            grid on
            hold on
            plot3(optimalParetoFront(:,1),optimalParetoFront(:,2), optimalParetoFront(:,3), 'ob')
            hold off
            title(['Optimization of ', functionName, ' at generation ', num2str(currentGeneration), '/', num2str(NSGAIIParameters.maximumGeneration)]);
            legend(algorithmUsed, 'Optimal Pareto Front');
            xlabel('Objective 1');
            ylabel('Objective 2');
            zlabel('Objective 3');
            drawnow

        end

    else

        if testFunctionParameters.objectiveNumber == 2

            indices = find(rank == 1);
            data = objectiveValues(indices,:);
            plot(data(:,1),data(:,2),'*r')
            title(['Optimization of ', functionName, ' at generation ', num2str(currentGeneration), '/', num2str(NSGAIIParameters.maximumGeneration)]);
            legend(algorithmUsed);
            xlabel('Objective 1');
            ylabel('Objective 2');
            drawnow

        elseif testFunctionParameters.objectiveNumber == 3

            indices = find(rank == 1);
            data = objectiveValues(indices,:);
            plot3(data(:,1),data(:,2),data(:,3),'*r')
            grid on
            title(['Optimization of ', functionName, ' at generation ', num2str(currentGeneration), '/', num2str(NSGAIIParameters.maximumGeneration)]);
            legend(algorithmUsed);
            xlabel('Objective 1');
            ylabel('Objective 2');
            zlabel('Objective 3');
            drawnow

        end

    end
    
end

