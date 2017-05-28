

function [matrixOUT] = ColumnShifter( matrixToGet, sortoncol, keep, shift)
% Column shifter - reorganize a pre-existing data matrix to get an XY type plot
    % matrixToGet, X coordinates to keep, block size to shift by, column to use for sorting (eg DIC column)
    % ( matrixToGet, sortoncol, keep, shift)


startrow = 1;  		% Row placer, start at row 1
rows = 0;  			% End row (based on size of currently bit off hunk of data) ((add on the total number of current data points))

% SHOULD INPUT SOME CHECKS SO THAT YOU KNOW YOU ARE EVENLY DIVIDING UP THE INPUT MATRIX AND YOU'RE NOT GETTTING ERRORS DUE TO REMAINDER DIVISION

block_start = keep + 1;
block_end = keep + shift;

matrixOUT = zeros(size(matrixToGet));
itemstosort = unique(matrixToGet(:,sortoncol));

for k = 1:length(itemstosort), 
	datain = matrixToGet(matrixToGet(:,sortoncol) == itemstosort(k),:);

	% Set-up to get appropriate staggering
    currentsize = size(datain,1);    % Number of rows that will be needed for this bit of data
    rows = rows + currentsize;    % Number of rows that will be needed for this bit of data.

    % Put the data into the output matrix
    matrixOUT(startrow:rows,(1:keep)) = datain(:,1:keep);
    matrixOUT(startrow:rows,(block_start:block_end)) = datain(:,(keep+1):end);

    % Set up for the next iteration
    startrow = rows + 1;            % move down one row to start placing the next chunk

    block_start = block_end + 1;
    block_end = block_end + shift;

end 
