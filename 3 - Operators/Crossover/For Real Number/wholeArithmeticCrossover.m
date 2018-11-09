%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        Whole Arithmetic Crossover                       %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 20th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference: Introduction to Evolutionary Algorithms
%            Xinji Yu & Mitsuo Gen - Springer


function offsprings = wholeArithmeticCrossover(parents, testFunctionParameters)

alpha = rand;

offsprings(1, 1:testFunctionParameters.dim) = alpha*parents(1, 1:testFunctionParameters.dim) + (1 - alpha)*parents(2, 1:testFunctionParameters.dim);
offsprings(2, 1:testFunctionParameters.dim) = (1 - alpha)*parents(1, 1:testFunctionParameters.dim) + alpha*parents(2, 1:testFunctionParameters.dim);
