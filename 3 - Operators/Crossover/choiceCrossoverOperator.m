%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    Choice of the Crossover Operator                     %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 19th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [offsprings] = choiceCrossoverOperator(parents, GAParameters, testFunctionParameters)

switch GAParameters.crossoverType
    
    case 'Single-Point Crossover' % Ok
       
        offsprings = singlePointCrossover(parents, GAParameters, testFunctionParameters);
        
    case 'Two-Point Crossover' % Ok
        
        offsprings = twoPointCrossover(parents, GAParameters, testFunctionParameters);
        
    case 'Uniform Crossover' % Ok
        
        offsprings = uniformCrossover(parents, GAParameters, testFunctionParameters);
        
    case 'Whole Arithmetic Crossover' % Ok
        
        offsprings = wholeArithmeticCrossover(parents, testFunctionParameters);
        
    case 'Local Arithmetic Crossover' % Ok
        
        offsprings = localArithmeticCrossover(parents, testFunctionParameters);
        
    case 'Blend Alpha Crossover' % Ok
        
        offsprings = blendAlphaCrossover(parents, GAParameters, testFunctionParameters);
    
    case 'Simulated Binary Crossover' % Ok
        
        offsprings = simulatedBinaryCrossover(parents, GAParameters, testFunctionParameters);
        
    case 'Unimodal Normal Distribution Crossover' % Not ok
        
        offsprings = unimodalNormalDistributionCrossover(parents, GAParameters);
   
end