

function [ averagesMatrix ] = GroupAvg(matrixToGet, basename, sortbycol, genotype)
% GroupAvg(matrixToGet, basename, sortbycol, genotype) takes the average of subsets of data that are specified in sortbycol and returns a file with average, stdev, and number of values
%   Detailed explanation goes here


% Set up new space to keep average, standard deviation, and number for each set
averagesMatrix = zeros(1,(size(matrixToGet,2)*3)+1);

% get the data and separate it into the correct subsets
types = unique(matrixToGet(:,sortbycol));

    for currtype = 1:(size(types,1)),
        dataset = matrixToGet(matrixToGet(:,sortbycol) == types(currtype,:),:);

        columnplace = 1;

        % Walk through the data subset to get average, standard deviation, and number for each column of data in the matrix
        % Do not include any rows where the value is zero (otherwise the average will get crazy off)
        for column = 1:size(dataset,2),

            tempAvgs(:,columnplace:columnplace+2) = ...
                [mean(matrixToGet(matrixToGet(:,column)~=0,column)),...
                std(matrixToGet(matrixToGet(:,column)~=0,column)),...
                size(matrixToGet(matrixToGet(:,column)~=0,column),1)];

            columnplace = columnplace +3;

        end

        averagesMatrix(currtype,:) = [currtype,tempAvgs];			% append the embryoID to the first column of the averages matrix

    end

    %name = sprintf('%1$s%2$s.%3$s', basename,genotype,'csv');
    name = SpiffyName('csv', basename, genotype);

    % if an averages file exists (from a previous round), append the new data to it
    if exist(name,'file') ~= 0,
        dlmwrite(name,averagesMatrix,'-append');

        % Otherwise if there is no file for this genotype to save to, create it and then add this data.
    elseif exist(name,'file') == 0,
        fileID = fopen(name,'w+');
        fprintf(fileID, '%1s,%2s,%3s,%4s\n', 'EmbryoID', 'Avg', 'StdDev', 'Number');
        dlmwrite(name,averagesMatrix,'delimiter',',','-append');
        fclose(fileID);

    end
end


