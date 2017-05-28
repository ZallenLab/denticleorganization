

function DenticleOrganization_ImportData(genotype, scale)



%%  GET DATA
for dentrow = 0:6,
    
    dirquerystr = sprintf('*%1$s*%2$c%3$i.%4$s',genotype,'r',dentrow,'txt');
    CellCounterData = dir([inputdirectory, filesep, dirquerystr]);
    
    
    % CellCounterData is a complex matrix with
    
    if dentrow == 6,        % Combine rows 6 and 7, since they're essentially the same (easier to save a row 7 than 2bxr6)
        dirquerystr2 = sprintf('*%1$s*%2$c%3$i.%4$s',genotype,'r',7,'txt');
        CellCounterData2 = dir([inputdirectory, filesep, dirquerystr2]);
        CellCounterData = cat(1,CellCounterData,CellCounterData2);
    end
    
    numfiles = length(CellCounterData);         % Get the number of .txt files found in the above
    
    % Calculate spacing data, cell size data
    % Get each data file individually in sequence, import the data portion
    % so you can run calculations on it
    for k = 1:numfiles,
        fileToRead = CellCounterData(k).name;
        dataImport = importdata(fileToRead);
        data = dataImport.data;
        
        ntotal = size(data,1);      % nt = the number of data points total (denticles and cell edges)
        dist = zeros;               %reset matrix (i hope)
        
        % Get name of file to add into matrices to allow for separation by embryo
        [filePath,fileName,fileExt] = fileparts(fileToRead);
        embryoID = 100*str2double(fileName(1:4))+str2double(fileName(6:7));
        % 1:4 is the embryoID#, 6:7 is the stage
        % Convert first 4 characters (ID number) from a string to a number so it can go into the established matrices without trouble
        % Changing from fileName(6:7) to just (7)
        % to eliminate a problem with scientific
        % notation in the csv file at the end...
        
        if isnan(embryoID) == 1,
            error('something is wrong with your input text file names')
        end
        
        beltID = str2double(fileName((length(fileName))-2));
        
        % SORT DATA
        % Sort positions (all objects) by order along row; points can be picked in any order in the CellCounter plugin
        data = sortrows(data,beltdir);
        
        % List of just denticles
        d_ind = find(data(:,1) == denticlemarker);
        num_dent = length(d_ind);
        denticles = data(d_ind,:);
        
        % List of just edges
        e_ind = find(data(:,1) == edgemarker);
        num_edges = length(e_ind);
        edges = data(e_ind,:);
        
        
        % Separate out on a cell basis
        ibccounter = 1;          % Counter for interbycell matrix
        aecounter = 1;
        
        bycell = zeros(1,20);       % Make empty matrix of a consistent column number
        intrabycell = zeros(1,7);
        cellnumber = 0;
        
        for i = 1:num_edges-1,
            cellnumber = cellnumber + 1;                
            a = 1;                          % Counter for inclDent matrix; resets with each new cell edge pair
            inclDent = zeros(20,5);         % Preallocate & reset for each iteration
            adjacentDent = zeros(2,5);      % Preallocate & reset for each iteration
            
            % Get edges of current cell
            ledge = edges(i,:);
            redge = edges(i+1,:);
            
            % Get X values from all the denticles in this row
            onlydentvals = denticles(:,beltdir);
            
            
            % Look at all the denticles in the dataset and determine which are contained within the cell boundaries
            for j = 1:size(onlydentvals,1),
                
                % Get X coordinate for the denticle
                dentinq = onlydentvals(j);
                
                % Compare the X coordinates of denticle 'j' to the left and
                % right edges. If the denticle is contained within a cell,
                % write all its associated data into the 'a'th row of
                % inclDent
                if (dentinq > ledge(:,beltdir)) && (dentinq < redge(:,beltdir)),
                    inclDent(a,1:5) = denticles(j,1:5);
                    a = a + 1;
                    
                    % If the denticle is not the first or the last in the dataset, get data on the denticles on either side
                    if (j ~= 1) && (j ~= size(onlydentvals,1)),
                        dentbefore = onlydentvals(j-1);
                        dentafter = onlydentvals(j+1);
                        
                        % If the neighboring denticle(s) are NOT in the same cell, write them into the appropriate matrix
                        if onlydentvals(j-1) < ledge(:,beltdir),
                            adjacentDent(1,1:5) = denticles(j-1,1:5);
                        end
                        
                        if onlydentvals(j+1) > redge(:,beltdir),
                            adjacentDent(2,1:5) = denticles(j+1,1:5);
                        end
                        
                        % Otherwise if the denticle is the first or the last in the set, write the adjacent coordinates as zeros
                    elseif j == 1,
                        adjacentDent(1,1:5) = [0 0 0 0 0];
                        
                    elseif j == size(onlydentvals,1),
                        adjacentDent(2,1:5) = [0 0 0 0 0];
                    end
                end
            end
            
            
            
            if max(max(inclDent)) == 0,
                cellnumber = cellnumber-1;
                %{
                    This gets the maximum value in the entire inclDent
                    matrix; if this is an empty matrix, the max will be 0.
                    Stop if there are no denticles assigned between the two
                    edges (indicates a break ex. at the midline; don't need
                    to save as multiple files any more)
                %}
            else
                inclDent = inclDent(any(inclDent,2),:);         % This eliminates zero rows in the 'inclDent' matrix so when it is counting the number of dent/cell you only count nonzero rows
                
                dentincell = size(inclDent,1);
                
                % Embryo identification data
                bycell(i,1) = embryoID;
                bycell(i,2) = dentrow;
                bycell(i,3) = beltID;
                
                % basic data about the cell
                bycell(i,4) = cellnumber;                                   % Cell number within the row
                bycell(i,5) = dentincell;                                   % number of denticles in the current cell
                bycell(i,6) = (norm((redge(3:4))-(ledge(3:4))))/scale;      % DV length
                
                % Denticle to edge distances
                bycell(i,7) = norm((inclDent(1,(3:4)))-(ledge(3:4)))/scale;
                bycell(i,8) = norm((redge(3:4))-(inclDent(dentincell,3:4)))/scale;
                
                % Adjacent/inter distance values
                if (adjacentDent(1,3) ~= 0) && (adjacentDent(1,4) ~= 0),        % If the X and Y coordinate values are non-zero
                    bycell(i,9) = norm((adjacentDent(1,3:4))-(inclDent(1,3:4)))/scale;      % Calculate inter distance using the denticle pair on the left
                else bycell(i,9) = 0;
                end
                
                if (adjacentDent(2,3) ~= 0) && (adjacentDent(2,4) ~= 0),
                    bycell(i,10) = norm((adjacentDent(2,3:4))-(inclDent(dentincell,3:4)))/scale; % Calculate inter distances using the denticle pair on the right
                else bycell(i,10) = 0;
                end
                
                
                % Same cell/intra distance values
                if dentincell > 1,
                    c = 11;
                    for j = 1:(dentincell - 1),
                        bycell(i,c) = norm(inclDent(j,3:4)-inclDent(j+1,3:4))/scale;
                        c = c+1;
                    end
                end
                
                
                % Write intra data into a separate matrix with one cell per
                % row, type 1
                if dentincell == 1,
                    intrabycell(ibccounter,1:7) = 0;
                elseif dentincell >1,
                    for j = 1:(dentincell-1),
                        intrabycell(ibccounter,1) = embryoID;
                        intrabycell(ibccounter,2) = dentrow;
                        intrabycell(ibccounter,3) = beltID;
                        intrabycell(ibccounter,4) = i; % Cell/iteration number
                        intrabycell(ibccounter,5) = dentincell;  % # of denticles in cell
                        intrabycell(ibccounter,6) = (norm((redge(3:4))-(ledge(3:4))))/scale; % DV length of cell
                        intrabycell(ibccounter,7) = norm(inclDent(j,3:4)-inclDent(j+1,3:4))/scale; % Intra distance
                        ibccounter = ibccounter+1;
                    end
                end
                
                % Write intra data into a separate matrix, type 2
                % (IntraInterSpread.m)
                if dentincell == 1,
                    aebydic(aecounter,1) = embryoID;
                    aebydic(aecounter,3) = dentrow;
                    aebydic(aecounter,4) = beltID;
                    aebydic(aecounter,5) = dentincell;
                    
                    aebydic(aecounter,7) = norm((inclDent(1,(3:4)))-(ledge(3:4)))/scale;
                aebydic(aecounter,8) = norm((redge(3:4))-(inclDent(1,3:4)))/scale;
                
                elseif dentincell >1,
                    for j = 1:dentincell,
                        aebydic(aecounter,1) = embryoID; %#ok<*SAGROW>
                        aebydic(aecounter,3)= dentrow;
                        aebydic(aecounter,4) = beltID;
                        aebydic(aecounter,5) = dentincell;
                        if j == 1,
                            aebydic(aecounter,7) = norm((inclDent(1,(3:4)))-(ledge(3:4)))/scale;
                            aebydic(aecounter,6) = norm(inclDent(j,3:4)-inclDent(j+1,3:4))/scale; % Intra distance
                        elseif (j ~= dentincell) && (j ~=1),
                            aebydic(aecounter,6) = norm(inclDent(j,3:4)-inclDent(j+1,3:4))/scale; % Intra distance
                        elseif j == dentincell
                            aebydic(aecounter,8) = norm((redge(3:4))-(inclDent(j,3:4)))/scale;
                        end

                        aecounter = aecounter+1;
                    end
                end
            
            aebydic = aebydic(any(aebydic,2),:);
            
            dlmwrite(SpiffyName('csv', 'IntraDV', genotype),aebydic,'-append');
            
        end
    end
    
    
    % OUTPUT FILES
    
    % save('Results.mat', 'bycell', 'intrabycell');
    
    cbcfilename = (SpiffyName('csv','cellbyCellStack', genotype));
    ibcfilename = (SpiffyName('csv','intrabyCellStack', genotype));
    % these cannot be changed to include the timerun in the filename or
    % everything else breaks!
    
    dlmwrite(cbcfilename,bycell,'-append');
    dlmwrite(ibcfilename,intrabycell,'-append');
    % Add the new data from this row to the end of the existing matrix
    
    
    
end
disp('done with row:')
disp(dentrow)
    
    
end

disp('Done with import and quantification')

%% Reformat data (column headers, reorder)

% Cellbycell
% embryoID,row,belt,cell,dentincell,Dvlen,dentEdgeL,dentEdgeR,AdjL,AdjR,Intra1-2,Intra2-3,Intra3-4,Intra4-5,Intra5-6,Intra6-7,Intra7-8,Intra8-9

cbctitles_raw = {'embryoID', 'row', 'belt', 'cell', 'dentincell', 'Dvlen', 'dentEdgeL', 'dentEdgeR', 'AdjL', ...
    'AdjR', 'Intra1-2', 'Intra2-3','Intra3-4','Intra4-5', 'Intra5-6', 'Intra6-7', 'Intra7-8', 'Intra8-9',...
    'Intra9-10', 'Intra10-11','Intra11-12','Intra12-13', 'Intra13-14', 'Intra14-15', 'Intra15-16', 'Intra16-17','Intra17-18', 'Intra18-19',...
    'Intra19-20','Intra20-21','Intra21-22','Intra22-23', 'Intra23-24', 'Intra24-25', 'Intra25-26', 'Intra26-27','Intra27-28', 'Intra28-29','Intra29-30',...
    '-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
cbctemp = dlmread(cbcfilename);
cbctitles = cbctitles_raw(1:(size(cbctemp,2)));
AddHeaders((SpiffyName('csv','CellbyCell', genotime)), cbctitles, cbctemp);

% Intrabycell
ibctitles = {'embryoID', 'row', 'belt', 'cell', 'dentincell', 'Dvlen', 'Intra','-','-','-','-'};
ibctemp = dlmread(ibcfilename);
AddHeaders(( SpiffyName('csv','IntrabyCell', genotime)), ibctitles, ibctemp);

% Rearrange data to look like it would across each cell
% Uses cbctemp data

maxintra = max(cbctemp(:,5)) -1;    % Max number of intra values that need to be accomodated is one fewer than the total number of denticles in the final output data

% Titles set-up and joining routine
cotitles_part1 = {'embryoID', 'row', 'belt', 'cell', 'dentincell', 'Dvlen', 'AdjL', 'dentEdgeL'};
cotitles_part2 = {'-', 'Intra1-2', 'Intra2-3','Intra3-4','Intra4-5', 'Intra5-6', 'Intra6-7', 'Intra7-8', 'Intra8-9',...
    'Intra9-10', 'Intra10-11','Intra11-12','Intra12-13', 'Intra13-14', 'Intra14-15', 'Intra15-16', 'Intra16-17','Intra17-18', 'Intra18-19',...
    'Intra19-20','Intra20-21','Intra21-22','Intra22-23', 'Intra23-24', 'Intra24-25', 'Intra25-26', 'Intra26-27','Intra27-28', 'Intra28-29','Intra29-30'};
% empty '-' at position 1 in cotitles2 is a contrivance so you don't have to subtract 1 from maxdent; the effect is the same as (1:(maxdent-1))
cotitles_part3 = {'dentEdgeR', 'AdjR'};
cotitles_part4 = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};

cotitles = [cotitles_part1,cotitles_part2(2:maxintra+1),cotitles_part3];

% Add pad headers ('-') if necessary - ie if there are blank rows
% at the end of the cbctemp matrix
if size(cotitles,2) < size(cbctemp,2)
    delta = size(cbctemp,2) - size(cotitles,2);
    cotitles = [cotitles,cotitles_part4(1:delta)];
end

cellorder = [cbctemp(:,1:6), cbctemp(:,9), cbctemp(:,7),cbctemp(:,11:(maxintra +10)), cbctemp(:,8), cbctemp(:,10)];
%  Intra values placement: 11(the start in cbctemp) : (maxdent -1 (one less intra value than there are denticles) +10 (number of positions before the intra goes in))
dlmwrite((SpiffyName('csv','cellOrderStack', genotype)),cellorder);
AddHeaders((SpiffyName('csv','CellOrder', genotime)), cotitles, cellorder)

% Rearrange Intra v dv
aetemp = dlmread(SpiffyName('csv','IntraDV', genotype));
aetemp = sortrows(aetemp,5);

dlmwrite(SpiffyName('csv','IntraDVsort', genotype),aetemp);


disp('Done with reformatting')

