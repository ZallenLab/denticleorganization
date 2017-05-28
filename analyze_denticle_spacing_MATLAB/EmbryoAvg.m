

function [ averagesMatrix ] = EmbryoAvg(matrixToGet, basename, curID, genotype)
%EmbryoAvg(matrixToGet, basename, curID, genotype)
%   Averages based on embryoID to generate cleaner averaged plots


% Titlenames
AvgALL_titles = ['embryoID', reshape([arrayfun(@(x)sprintf('Avg_%i',x),1:size(matrixToGet,2),'uniformOutput',false);...
                                arrayfun(@(x)sprintf('StDev_%i',x),1:size(matrixToGet,2),'uniformOutput',false);...
                                arrayfun(@(x)sprintf('N_%i',x),1:size(matrixToGet,2),'uniformOutput',false)],1,[])];

AvgOnly_titles = ['embryoID', arrayfun(@(x)sprintf('Avg_%i',x),1:size(matrixToGet,2),'uniformOutput',false)];


% Set up new space to keep average, standard deviation, and number for each set
averagesMatrix = zeros(1,(size(matrixToGet,2)*3));
avgonlyMatrix = zeros(1,size(matrixToGet,2));

% Walk through the original matrix to get average, standard deviation, and number for each column of data in the original matrix
% Do not include any rows where the value is zero (otherwise the average will get crazy off)

columnplace = 1;
for column = 1:size(matrixToGet,2),
    
    averagesMatrix(:,columnplace:columnplace+2) = ...
                [mean(matrixToGet(matrixToGet(:,column)~=0,column)),...
                 std(matrixToGet(matrixToGet(:,column)~=0,column)),...
                 size(matrixToGet(matrixToGet(:,column)~=0,column),1)];   
             
    columnplace = columnplace +3;
end

for column = 1:size(matrixToGet,2),
     avgonlyMatrix(:,column) = mean(matrixToGet(matrixToGet(:,column)~=0,column));
end 

% append the embryoID to the first column of the averages matrix
averagesMatrix = [curID,averagesMatrix];
avgonlyMatrix = [curID,avgonlyMatrix];


%name = sprintf('%1$s%2$s.%3$s', basename,genotype,'csv');
averagesName = SpiffyName('csv', basename, genotype,'avg+std+N');

%     dlmwrite(SpiffyName('csv',basename,genotype,'avg+std+N','basic'),averagesMatrix,'-append')

    if exist(averagesName,'file') == 1,    % if an averages file exists (from a previous round), append the new data to it      
        dlmwrite(averagesName,averagesMatrix,'-append'); 
        
    elseif exist(averagesName,'file') == 0,	% Otherwise if there is no file for this genotype to save to, create it and then add this data. 
        % AddHeaders(averagesName,AvgALL_titles,averagesMatrix);
        dlmwrite(averagesName,averagesMatrix);
    end
    
    
avgonlyName = SpiffyName('csv', basename, genotype);
    
%     dlmwrite(SpiffyName('csv',basename,genotype,'basic'),avgonlyMatrix,'-append')
        %{
            0   name does not exist.
            1   name is a variable in the workspace.
            2  One of the following is true:    
                name exists on your MATLAB® search path as a file with extension .m.
                name is the name of an ordinary file on your MATLAB search path.
                name is the full pathname to any file.
        %}


    if exist(avgonlyName,'file') == 2,    % if an averages file exists (from a previous round), append the new data to it      
        dlmwrite(avgonlyName,avgonlyMatrix,'-append'); 

    elseif exist(avgonlyName,'file') == 0,          % Otherwise if there is no file for this genotype to save to, create it and then add this data. 
        % AddHeaders(avgonlyName,AvgOnly_titles,avgonlyMatrix);
        dlmwrite(avgonlyName,avgonlyMatrix);
    end 


end

