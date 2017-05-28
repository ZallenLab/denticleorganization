

function [fractionMatrix, sumofMTG] = GimmeFractions( matrixToGet, directionToSum, varargin)
% Sum up the elements in a matrix (either by column or by row) and get the
% fraction of the total that each element comprises. 
% GimmeFractions( matrixToGet, directionToSum, varargin)
%           varargin can be anything, use if you want to eliminate zero counts from the sum
% Make a fraction of total from a matrix given the matrix to convert and
% the max value to use as a denominator
% directionToSum variable: this is the value being used as the total to
% calculate the fraction each element is of the whole
% '1' = sum down the column; 
% '2' = sum across the row


narginchk(2,3)

    if nargin == 2, 
        sumofMTG = sum(matrixToGet,directionToSum);

    elseif nargin == 3, 

        if directionToSum == 1,   % sum of all rows in each of the columns
            for column = 1:size(matrixToGet,2), 
                temp = matrixToGet(matrixToGet(:,column)~=0);        % Get rid of any zeros (important for cases using histcounts)
                if isempty(temp),
                    temp = 0;
                end                    
                temp2 = sum(temp,1);
                sumofMTG(1,column) = temp2; %#ok<*AGROW>
            end
        
        elseif directionToSum == 2,  % sum of all columns over all rows
            for row = 1:size(matrixToGet,1),
                temp = matrixToGet(matrixToGet(row,:)~=0);
                temp2 = sum(temp,2);
                sumofMTG(row,1) = temp2;
            end
        end

    elseif nargin > 3
            error('Too many arguments!!!')

    elseif nargin < 2
            error('Not enough arguments!!!')

    end

    
fractionMatrix = zeros(size(matrixToGet));

     
    if directionToSum == 1,  % sum of all rows in each of the columns
        
        for column = 1:size(matrixToGet,2),
            for row = 1:size(matrixToGet,1)
                if matrixToGet(row,column) == 0, 
                    fractionMatrix(row,column) = 0;
                elseif matrixToGet(row,column) ~= 0,
                     fractionMatrix(row,column) = matrixToGet(row,column) / sumofMTG(1,column);
                end 
%                 fractionMatrix(row,column) = matrixToGet(row,column) / sumofMTG(1,column);
            end 
        end
        
    elseif directionToSum == 2,  % sum of all columns over all rows
        for row = 1:size(matrixToGet,1),

            for column = 1:size(matrixToGet,2)
                
                if matrixToGet(row,column) == 0,
                    fractionMatrix(row,column) = 0;
                elseif matrixToGet(row,column) ~= 0,
                    fractionMatrix(row,column) = matrixToGet(row,column) / sumofMTG(row,1);
                end
%                 fractionMatrix(row,column) = matrixToGet(row,column) / sumofMTG(row,1);
            end

        end
    end
    
end

