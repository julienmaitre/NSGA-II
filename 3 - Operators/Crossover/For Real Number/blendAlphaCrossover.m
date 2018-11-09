%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Blend Alpha Crossover                          %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 20th 2017                                                %
% Version : 2.0                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference: Introduction to Evolutionary Algorithms
%            Xinji Yu & Mitsuo Gen - Springer
%
%            Other reference : Eshelman and Schaffer [5] reported that
%            alpha=0.5 is a good choice - 1993


function offspring = blendAlphaCrossover(parents, GAParameters, testFunctionParameters)

offspring = zeros(1, testFunctionParameters.dim);

for i = 1:1:testFunctionParameters.dim
    
    if parents(1, i) < parents(2, i)
        
        A = parents(1,i) - GAParameters.BLXAlphaCoefficient*(parents(2, i) - parents(1, i));
        B = parents(2,i) + GAParameters.BLXAlphaCoefficient*(parents(2, i) - parents(1, i));
        
        offspring(1,i) = A + (B - A)*rand;
        
    else
        
        A = parents(2,i) - GAParameters.BLXAlphaCoefficient*(parents(1, i) - parents(2, i));
        B = parents(1,i) + GAParameters.BLXAlphaCoefficient*(parents(1, i) - parents(2, i));
        
        offspring(1,i) = A + (B - A)*rand;
        
    end
    
end

