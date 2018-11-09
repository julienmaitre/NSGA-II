%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            Bit-flip mutation                            %
%                                                                         %
% Author : Julien Maitre                                                  %
% Date : October 19th 2017                                                %
% Version : 1                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference : Introduction to Evolutionary Algorithms
%             Xinjie Yu && Mitsuo Gen - Springer


function [offsprings] = bitFlipMutation(offsprings, GAParameters, testFunctionParameters)

for i = 1:1:2
    
    for j = 1:1:testFunctionParameters.dim
        
        for k = 1:1:GAParameters.numberOfBits
            
            if rand < GAParameters.mutationRate

                offsprings{i,1}(j,k) = not(offsprings{i,1}(j,k));
            
            end

        end
        
    end
    
end

offsprings = convertBinaryToDecimal(offsprings, 2, GAParameters, testFunctionParameters);

