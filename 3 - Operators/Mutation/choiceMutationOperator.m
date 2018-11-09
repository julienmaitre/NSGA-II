%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Choice of the Mutation Operator                     %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 19th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer

function [offsprings] = choiceMutationOperator(offsprings, currentGeneration, GAParameters, testFunctionParameters)

switch GAParameters.mutationType
    
    case 'Bit-Flip Mutation' % Ok
        
        offsprings = bitFlipMutation(offsprings, GAParameters, testFunctionParameters);
    
    case 'Uniform Mutation' % Ok
     
        offsprings = uniformMutation(offsprings, GAParameters, testFunctionParameters);
        
    case 'Boundary Mutation' % Ok
        
        offsprings = boundaryMutation(offsprings, GAParameters, testFunctionParameters);
        
    case 'Non Uniform Mutation' % Ok
        
        offsprings = nonUniformMutation(offsprings, currentGeneration, GAParameters, testFunctionParameters);
        
    case 'Normal Mutation' % Revoir
        
        offsprings = normalMutation(offsprings, GAParameters);
        
    case 'Polynomial Mutation' % Ok
        
        offsprings = polynomialMutation(offsprings, GAParameters, testFunctionParameters);
        
end