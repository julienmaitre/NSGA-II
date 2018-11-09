%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%     Multi-Objective Test Problems to Compare Optimization Algorithms    %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : November 09th 2018                                               %
% Version : 2.2                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Information

% This function containts the entire information and implementations of the  
% benchmark functions from the literature 

% Inputs : 
%   functionName is the name of the function that the user has chosen to optimise
%   dim is the number of variables (dimension of the problem)
%   objectiveNumber is the number of functions (objectives) to optimize for one test problem
%   pointNumber is the number of points used for the representation of the Pareto-optimal front, surface or other

% Outputs :
%   fobj is the object that allows to evaluate the solutions according to the test problem chosen by the user 
%   lb is the lower bound value of variables 
%   up is the uppper bound value of variables
%   dim is the number of variables (dimension of the problem)
%   lowerLimit is the vector of lower bound values for each variable
%   upperLimit is the vector of upper bound values for each variable
%   objectiveNumber is the number of functions (objectives) to optimize for one test problem
%   optimizationType are the array of objectives (maximization or minimization) for each function to optimize
%   pointNumber is the number of points used for the representation of the Pareto-optimal front, surface or other
%   ParetoOptimalFront are the values of the Pareto-optimal front, surface or other according to the number of points 


function [testFunctionParameters] = multiObjectiveTestProblems(functionName, dim, objectiveNumber, pointNumber)

switch functionName
    
    %% Unconstrained Test Functions 
    
    case 'ZDT1'
        % Information about ZDT1
        testFunctionParameters.fobj = @ZDT1;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = 1;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Pareto-optimal front of ZDT1
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0,1,pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - sqrt(testFunctionParameters.ParetoOptimalFront(:,1));
        
    case 'ZDT2'
        % Information about ZDT2
        testFunctionParameters.fobj = @ZDT2;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = 1;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of ZDT2
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0,1,pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - testFunctionParameters.ParetoOptimalFront(:,1).^2;
        
    case 'ZDT3'
        % Information about ZDT3
        testFunctionParameters.fobj = @ZDT3;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = 1;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of ZDT3
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        pointNumberVector = linspace(1,pointNumber,6);
        newPointNumber = floor(pointNumberVector(1,2:end))-floor(pointNumberVector(1,1:end-1));
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = [linspace(0,0.0830015349,newPointNumber(1,1)) ...
            linspace(0.1822287280,0.2577623634,newPointNumber(1,2)) ...
            linspace(0.4093136748,0.4538821041,newPointNumber(1,3)) ...
            linspace(0.6183967944,0.6525117038,newPointNumber(1,4)) ...
            linspace(0.8233317983,0.8518328654,newPointNumber(1,5) + 1)];
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - sqrt(testFunctionParameters.ParetoOptimalFront(:,1)) - ...
            testFunctionParameters.ParetoOptimalFront(:,1).*sin(10*pi*testFunctionParameters.ParetoOptimalFront(:,1));
        
    case 'ZDT4'
        % Information about ZDT4
        testFunctionParameters.fobj = @ZDT4;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = -5;
        testFunctionParameters.ub = 5;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = [0 testFunctionParameters.lb*ones(1,dim - 1)];
        testFunctionParameters.upperLimit = [1 testFunctionParameters.ub*ones(1,dim - 1)];
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of ZDT4
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0,1,pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - sqrt(testFunctionParameters.ParetoOptimalFront(:,1));
        
    case 'ZDT6'
        % Information about ZDT6
        testFunctionParameters.fobj = @ZDT6;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = 1;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of ZDT6
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0.2807753191,1,pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - testFunctionParameters.ParetoOptimalFront(:,1).^2;
        
    case 'SCH'
        % Information about SCH
        testFunctionParameters.fobj = @SCH;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = -1000;
        testFunctionParameters.ub = 1000;
        testFunctionParameters.dim = 1;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of SCH
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0,4,pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = (sqrt(testFunctionParameters.ParetoOptimalFront(:,1)) - 2).^2;
        
    case 'FON'
        % Information about FON
        testFunctionParameters.fobj = @FON;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = -4;
        testFunctionParameters.ub = 4;
        testFunctionParameters.dim = 3;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of FON
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = zeros(pointNumber,testFunctionParameters.objectiveNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,1) = linspace(0,1-exp(-4),pointNumber);
        testFunctionParameters.ParetoOptimalFront(1:pointNumber,2) = 1 - exp(-(2 - sqrt(-log(1 - testFunctionParameters.ParetoOptimalFront(:,1)))).^2);
        
    case 'POL'
        % Information about POL
        testFunctionParameters.fobj = @POL;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = -pi;
        testFunctionParameters.ub = pi;
        testFunctionParameters.dim = 2;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of POL
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
        
    case 'KUR'
        % Information about KUR
        testFunctionParameters.fobj = @KUR;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = -5;
        testFunctionParameters.ub = 5;
        testFunctionParameters.dim = dim;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of KUR
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
        
    %% Constrained Test Functions 
    case 'CONSTR'
        % Information about CONSTR
        testFunctionParameters.fobj = @CONSTR;
        testFunctionParameters.cobj = @constrainedCONSTR;
        testFunctionParameters.contrainedProblem = 'Contrained';
        testFunctionParameters.lb = 'Different';
        testFunctionParameters.ub = 'Different';
        testFunctionParameters.dim = 2;
        testFunctionParameters.lowerLimit = [0.1 0];
        testFunctionParameters.upperLimit = [1 5];
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 2;
        testFunctionParameters.optimizationType = repmat({'Minimization'}, 1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = repmat({'Superior or equal'}, 1, testFunctionParameters.constrainedNumber);
        testFunctionParameters.constraints = [6 1];
        
        % Optimal Pareto front of CONSTR
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
        
    case 'SRN'
        % Information about SRN
        testFunctionParameters.fobj = @SRN;
        testFunctionParameters.cobj = @constrainedSRN;
        testFunctionParameters.contrainedProblem = 'Contrained';
        testFunctionParameters.lb = -20;
        testFunctionParameters.ub = 20;
        testFunctionParameters.dim = 2;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 2;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = repmat({'Inferior or equal'}, 1, testFunctionParameters.constrainedNumber);
        testFunctionParameters.constraints = [225 -10];
        
        % Optimal Pareto front of SRN
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
        
    case 'TNK'
        testFunctionParameters.fobj = @TNK;
        testFunctionParameters.cobj = @constrainedTNK;
        testFunctionParameters.contrainedProblem = 'Contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = pi;
        testFunctionParameters.dim = 2;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,testFunctionParameters.dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,testFunctionParameters.dim);
        testFunctionParameters.objectiveNumber = 2;
        testFunctionParameters.constrainedNumber = 2;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = repmat({'Inferior or equal'}, 1, testFunctionParameters.constrainedNumber);
        testFunctionParameters.constraints = [0 0.5];
        
        % Optimal Pareto front of TNK
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
        
    case 'WATER'
        testFunctionParameters.fobj = @WATER;
        testFunctionParameters.cobj = @constrainedWATER;
        testFunctionParameters.contrainedProblem = 'Contrained';
        testFunctionParameters.lb = 'Different';
        testFunctionParameters.ub = 'Different';
        testFunctionParameters.dim = 3;
        testFunctionParameters.lowerLimit = [0.01 0.01 0.01];
        testFunctionParameters.upperLimit = [0.45 0.10 0.10];
        testFunctionParameters.objectiveNumber = 5;
        testFunctionParameters.constrainedNumber = 7;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, 5);
        testFunctionParameters.constrainedType = repmat({'Inferior or equal'}, 1, testFunctionParameters.constrainedNumber);
        testFunctionParameters.constrainedType{1,2} = 'Inferior';
        testFunctionParameters.constraints = [1 1 50000 16000 10000 2000 550];
        
        % Optimal Pareto front of WATER
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
     
    %% Many-Objective and Unconstrained Test Problems
    case 'DTLZ1'
        % Information about DTLZ1
        testFunctionParameters.fobj = @DTLZ1;
        testFunctionParameters.cobj = 'None';
        testFunctionParameters.contrainedProblem = 'Not contrained';
        testFunctionParameters.lb = 0;
        testFunctionParameters.ub = 1;
        testFunctionParameters.dim = objectiveNumber + 5 - 1;
        testFunctionParameters.lowerLimit = testFunctionParameters.lb*ones(1,dim);
        testFunctionParameters.upperLimit = testFunctionParameters.ub*ones(1,dim);
        testFunctionParameters.objectiveNumber = objectiveNumber;
        testFunctionParameters.constrainedNumber = 0;
        testFunctionParameters.optimizationType = repmat({'Minimization'},1, testFunctionParameters.objectiveNumber);
        testFunctionParameters.constrainedType = 'None';
        testFunctionParameters.constraints = [];
        
        % Optimal Pareto front of DTLZ1
        testFunctionParameters.pointNumber = pointNumber;
        testFunctionParameters.ParetoOptimalFront = [];
                    
end

end

% ZDT1
function fitnessValues = ZDT1(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);

    h = sum(population(:,2:end),2);
    g = 1 + 9.*h./(size(population,2) - 1);

    fitnessValues(:,2) = g.*(1 - sqrt(population(:,1)./g));
    
end

% ZDT2
function [fitnessValues] = ZDT2(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);

    h = sum(population(:,2:end),2);
    g = 1 + 9.*h./(size(population,2) - 1);

    fitnessValues(:,2) = g.*(1 - ((population(:,1)./g).^2));

end

% ZDT3
function [fitnessValues] = ZDT3(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);

    h = sum(population(:,2:end),2);
    g = 1 + 9.*h./(size(population,2)-1);

    fitnessValues(:,2) = g.*(1 - sqrt(population(:,1)./g) - (population(:,1)./g).*sin(10*pi.*population(:,1)));

end

% ZDT4
function [fitnessValues] = ZDT4(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);

    h = 10*cos(4*pi.*population(:,2:end));
    h = (population(:,2:end).^2) - h;

    g = 1 + 10*(size(population,2) - 1) + sum(h,2);

    fitnessValues(:,2) = g.*(1 - sqrt(population(:,1)./g));
    
end

% ZDT6
function [fitnessValues] = ZDT6(population, testFunctionParameters)

    fitnessValues(:,1) = 1 - exp(-4.*population(:,1)).*sin(6*pi.*population(:,1)).^6;

    g = 1 + 9.*(sum(population(:,2:end),2)./(size(population,2) - 1)).^0.25;

    fitnessValues(:,2) = g.*(1 - ((fitnessValues(:,1)./g).^2));
    
end

% SCH
function [fitnessValues] = SCH(population, testFunctionParameters)

    fitnessValues(:,1) = population.^2;
    fitnessValues(:,2) = (population - 2).^2;
    
end

% FON
function [fitnessValues] = FON(population, testFunctionParameters)

    fitnessValues(:,1) = 1 - exp(-sum((population(:,1:end) - (1/sqrt(3))).^2,2));
    fitnessValues(:,2) = 1 - exp(-sum((population(:,1:end) + (1/sqrt(3))).^2,2));
    
end

% POL
function [fitnessValues] = POL(population, testFunctionParameters)

    A1 = 0.5*sin(1) - 2*cos(1) + sin(2) - 1.5*cos(2);
    A2 = 1.5*sin(1) - cos(1) + 2*sin(2) - 0.5*cos(2);
    B1 = 0.5*sin(population(:,1)) - 2*cos(population(:,1)) + sin(population(:,2)) - 1.5*cos(population(:,2));
    B2 = 1.5*sin(population(:,1)) - cos(population(:,1)) + 2*sin(population(:,2)) - 0.5*cos(population(:,2));
    
    fitnessValues(:,1) = 1 + (A1 - B1).^2 + (A2 - B2).^2;
    fitnessValues(:,2) = (population(:,1) + 3).^2 + (population(:,2) + 1).^2;
    
end

% KUR
function [fitnessValues] = KUR(population, testFunctionParameters)

    r = population.^2;
    h = -10*exp(-0.2.*sqrt(r(:,1:end-1) + r(:,2:end)));
    fitnessValues(:,1) = sum(h,2);

    fitnessValues(:,2) = sum((abs(population).^0.8) + (5.*sin(population.^3)),2);
    
end

% CONSTR
function [fitnessValues] = CONSTR(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);
    fitnessValues(:,2) = (1 + population(:,2))./population(:,1);
    
end

function [constrainedValues] = constrainedCONSTR(population, testFunctionParameters)

    constrainedValues(:,1) = population(:,2) + 9*population(:,1) - 6;
    constrainedValues(:,2) = -population(:,2) + 9*population(:,1) - 1;
    
end

% SRN
function [fitnessValues] = SRN(population, testFunctionParameters)
    
    fitnessValues(:,1) = (population(:,1) - 2).^2 + (population(:,2) - 1).^2 + 2;
    fitnessValues(:,2) = 9*population(:,1) - (population(:,2) - 1).^2;
    
end

function [constrainedValues] = constrainedSRN(population, testFunctionParameters)

    constrainedValues(:,1) = population(:,1).^2 + population(:,1).^2 - 225;
    constrainedValues(:,2) = population(:,1) - 3*population(:,2) + 10;
    
end

% TNK
function [fitnessValues] = TNK(population, testFunctionParameters)

    fitnessValues(:,1) = population(:,1);
    fitnessValues(:,2) = population(:,2);

end

function [constrainedValues] = constrainedTNK(population, testFunctionParameters)

    constrainedValues(:,1) = -population(:,1).^2 - population(:,2).^2 + 1 + ...
        0.1*cos(16*atan(population(:,1)./population(:,2)));
    constrainedValues(:,2) = (population(:,1) - 0.5).^2 + ...
        (population(:,2) - 0.5).^2 - 0.5;
    
end

% WATER
function [fitnessValues] = WATER(population, testFunctionParameters)
    
    fitnessValues(:,1) = 106780.37*(population(:,2) - population(:,3)) + 61704.67;
    fitnessValues(:,2) = 3000*population(:,1);
    fitnessValues(:,3) = ((305700)*2289*population(:,2))./((0.06*2289)^0.65);
    fitnessValues(:,4) = (250)*2289*exp(-39.75*population(:,2) + 9.9*population(:,3) + 2.74);
    fitnessValues(:,5) = 25*(1.39./(population(:,1).*population(:,2)) + 4940*population(:,3) - 80);
    
end

function [constrainedValues] = constrainedWATER(population, testFunctionParameters)

    constrainedValues(:,1) = 0.00139./(population(:,1).*population(:,2)) + ...
        4.94*population(:,3) - 0.08 - 1;
    constrainedValues(:,2) = 0.000306./(population(:,1).*population(:,2)) + ...
        1.082*population(:,3) - 0.0986 - 1;
    constrainedValues(:,3) = 12.307./(population(:,1).*population(:,2)) + ...
        49408.24*population(:,3) + 4051.02 - 50000;
    constrainedValues(:,4) = 2.098./(population(:,1).*population(:,2)) + ...
        8046.33*population(:,3) - 696.71 - 16000;
    constrainedValues(:,5) = 2.138./(population(:,1).*population(:,2)) + ...
        7883.39*population(:,3) - 705.04 - 10000;
    constrainedValues(:,6) = 0.417./(population(:,1).*population(:,2)) + ...
        1721.26*population(:,3) - 136.54 - 2000;
    constrainedValues(:,7) = 0.164./(population(:,1).*population(:,2)) + ...
        631.13*population(:,3) - 54.48 - 550;
    
end

% DTLZ1
function [fitnessValues] = DTLZ1(population, testFunctionParameters)

    k = 5; % as suggested by Deb
    
    objectiveNumber = testFunctionParameters.objectiveNumber;
    n = objectiveNumber + k - 1; % this is the number of variable by default
    
    xm = population(:,n - k + 1:end); % xm contains the last k variables
    g = 100*(k + sum((xm - 0.5).^2 - cos(20*pi*(xm - 0.5)),2));

    fitnessValues(:,1) = 1/2*prod(population(:,1:objectiveNumber - 1),2).*(1 + g);
    
    for i = 2:1:objectiveNumber-1
        
        fitnessValues(:,i) = (1/2)*prod(population(:,1:objectiveNumber - i),2).*(1 - population(:,objectiveNumber - i + 1)).*(1 + g);
        
    end
    
    fitnessValues(:,objectiveNumber) = (1/2)*(1 - population(:,1)).*(1 + g);
        
end


