 
% 

function DataDivider(inputdir,thingtodo, genotype, totalgenotypes)
% datadivider Subdivide into pooled or individual datasets
% 0 gives jsut statistics for pooled data
% 1 gives pooled data
% 2 gives pooled + individual embryos data
% 3 gives individual data only


%%
% Import data
% load('UserInputData.mat');
cellbyCellStackALL = dlmread([inputdir, filesep, SpiffyName('csv','cellbyCellStack', genotype)]);
intrabyCellStackALL = dlmread([inputdir, filesep, SpiffyName('csv','intrabyCellStack', genotype)]);
cellOrderStackALL = dlmread([inputdir, filesep, SpiffyName('csv','cellOrderStack', genotype)]);

global embryoIDs noEmbryos totalCells dentNumbers

embryoIDs = unique(cellbyCellStackALL(cellbyCellStackALL(:,1)~=0));
noEmbryos = length(embryoIDs);
totalCells = size(cellbyCellStackALL,1);
dentNumbers = unique(cellbyCellStackALL(cellbyCellStackALL(:,5)~=0,5));            % Cells in the data have how many denticles?
maxDent_all = max(dentNumbers);                   % Largest number of denticles/cell in the dataset


formatOut = 'yyyy.mm.dd_HHMM';
timerun = datestr(now,formatOut);

save('DataSet_ALL.mat','cellbyCellStackALL','intrabyCellStackALL','cellOrderStackALL','embryoIDs','noEmbryos','totalCells','dentNumbers','maxDent_all','thingtodo');

%size(temp(any(temp,2),:))
% Data set general info/statistics - Numbers, per row, per belt, dent no, etc
datasetStats = zeros(7,1);

    datasetStats(1,1) = noEmbryos;       % # of embryos analyzed
    datasetStats(2,1) = size(cellbyCellStackALL(any(cellbyCellStackALL,2),:),1);    % Cells analyzed
    datasetStats(3,1) = sum(cellbyCellStackALL(:,5),1);    % Total number of denticles
    datasetStats(4,1) = max(unique(cellbyCellStackALL(:,2)));      % # of rows analyzed
    datasetStats(5,1) = mean(cellbyCellStackALL(cellbyCellStackALL(:,6)~=0,6));   % mean_celllength mean
    datasetStats(6,1) = std(cellbyCellStackALL(cellbyCellStackALL(:,6)~=0,6));   % mean_celllength stdev
    datasetStats(7,1) = mean(intrabyCellStackALL(intrabyCellStackALL(:,7)~=0,7));  % mean_intraspacing mean
    datasetStats(8,1) = std(intrabyCellStackALL(intrabyCellStackALL(:,7)~=0,7));  % mean_intraspacing stdev
    datasetStats(9,1) = mean(cellbyCellStackALL(cellbyCellStackALL(:,5)~=0,5));   % mean_dentnumber mean 
    datasetStats(10,1) = std(cellbyCellStackALL(cellbyCellStackALL(:,5)~=0,5));   % mean_dentnumber stdev
    datasetStats(11,1) = min(cellbyCellStackALL(cellbyCellStackALL(:,5)~=0,6));   % min length
    datasetStats(12,1) = max(cellbyCellStackALL(cellbyCellStackALL(:,5)~=0,6));   % max length
    datasetStats(13,1) = min(intrabyCellStackALL(intrabyCellStackALL(:,7)~=0,7));   % min spacing
    datasetStats(14,1) = max(intrabyCellStackALL(intrabyCellStackALL(:,7)~=0,7));   % max spacing

    


for ds = 1:10,
    bythenos(1,ds) = length(cellbyCellStackALL(cellbyCellStackALL(:,5) == ds,5));  % cells with n denticles
    bythenos(2,ds) = mean(cellbyCellStackALL(cellbyCellStackALL(:,5) == ds,6)); % mean cell length for cells with ds denticles
    bythenos(3,ds) = std(cellbyCellStackALL(cellbyCellStackALL(:,5) == ds,6)); % mean cell length for cells with ds denticles
    bythenos(4,ds) = mean(intrabyCellStackALL(intrabyCellStackALL(:,5) == ds,7)); % mean intra spacing for cells with ds denticles
    bythenos(5,ds) = std(intrabyCellStackALL(intrabyCellStackALL(:,5) == ds,7)); % mean intra spacing for cells with ds denticles
end


Stats = [{'# of embryos            ',   num2str(datasetStats(1,1));...
        '# of cells analyzed       ',   num2str(datasetStats(2,1));...
        'Total # of denticles      ',   num2str(datasetStats(3,1));...
        '# of rows analyzed        ',   num2str(datasetStats(4,1));...
        'mean cell length          ',   strcat(sprintf('%0.2f',datasetStats(5,1)), ' ± ', sprintf('%0.2f',(datasetStats(6,1))));...    %strcat(num2str(datasetStats(5,1)),' ±  ',num2str(datasetStats(6,1)));...
        'min & max cell length     ',   strcat(sprintf('%0.2f',datasetStats(11,1)), ' ± ', sprintf('%0.2f',(datasetStats(12,1))));...     %strcat(num2str(datasetStats(11,1)),' to  ',num2str(datasetStats(12,1)));...
        'mean spacing distance     ',   strcat(sprintf('%0.2f',datasetStats(7,1)), ' ± ', sprintf('%0.2f',(datasetStats(8,1))));...     %strcat(num2str(datasetStats(7,1)),' ±  ',num2str(datasetStats(8,1)));...
        'min & max spacing distance',   strcat(sprintf('%0.2f',datasetStats(13,1)), ' ± ', sprintf('%0.2f',(datasetStats(14,1))));...     %strcat(num2str(datasetStats(13,1)),' to  ',num2str(datasetStats(14,1)));...
        'mean denticle number      ',   strcat(sprintf('%0.2f',datasetStats(9,1)), ' ± ', sprintf('%0.2f',(datasetStats(10,1))));...     %strcat(num2str(datasetStats(9,1)),' ±  ',num2str(datasetStats(10,1)));...
        'max denticles             ',   num2str(maxDent_all)}]; % 'embryoIDs',''}; embryoNameCell];

ByNumbers = [{'cells with [1:10] denticles',num2str(bythenos(1,:));...
              'mean cell length    ', num2str(bythenos(2,:)); '                   ± ',num2str(bythenos(3,:));...
              'mean separation dist', num2str(bythenos(4,:)); '                   ± ',num2str(bythenos(5,:))}];

fileID = fopen(SpiffyName('txt', 'Stats',genotype),'w');
formatSpec = '%s\t %s\n';
fprintf(fileID,formatSpec,['genotype is          ',genotype]);

formatSpec = '\n%s\t %s';
[nrows, ncols] = size(Stats);
for row = 1:nrows,
  fprintf(fileID,formatSpec,Stats{row,:});
end

% formatSpec = '\n\n%d %d %d %d %d %d %d %d %d %d\n';           
% fprintf(fileID,formatSpec,[cellswdent]);

formatSpec = '\n\n%s\t %s';
for row = 1:size(ByNumbers,1),
    fprintf(fileID,formatSpec,ByNumbers{row,:});
end

fclose(fileID);

% T = cell2table(Stats, 'VariableNames',{'Property', 'number'});
% writetable(T,SpiffyName('txt', 'Stats',genotype))
%

if thingtodo == 0, % just run for stats output
    % eID = 'all';
    % genoID = sprintf('%1s_%2s',genotype,eID); %#ok<*NASGU>
    % % genoID = sprintf('%1s_%2s_%3s',genotype,eID,timerun); %#ok<*NASGU>

    % numerical_eID = 0;

    % % Get all data
    % cellbyCellStack = cellbyCellStackALL;
    % intrabyCellStack = intrabyCellStackALL;
    % cellOrderStack = cellOrderStackALL;

    % % dentNumbers = unique(intrabyCellStack(intrabyCellStack(:,5)~=0,5));            % Cells in the data have how many denticles?
    % % maxDent_all = max(dentNumbers);                   % Largest number of denticles/cell in the dataset
    % maxDent = maxDent_all;

    % save('currentDataSet.mat','cellbyCellStack','intrabyCellStack','cellOrderStack','eID','numerical_eID','genoID','maxDent_all','maxDent','noEmbryos','embryoIDs','thingtodo')


    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    formatSpec = '\n\n %s';
    fprintf(fileID,formatSpec,'embryoIDs');
    embryoStats = embryoIDs;
    formatSpec = '\n %d';
    fprintf(fileID,formatSpec,embryoStats);
    fclose(fileID);



elseif thingtodo == 1,
    % Pooled dataset only

    % Identity data
    %     eID = sprintf('%1$s','all');
    %             genoID = sprintf('%1s_%2s',eID,genotime); %#ok<*NASGU>

    eID = 'all';
    genoID = sprintf('%1s_%2s',genotype,eID); %#ok<*NASGU>
    % genoID = sprintf('%1s_%2s_%3s',genotype,eID,timerun); %#ok<*NASGU>

    numerical_eID = 0;

    % Get all data
    cellbyCellStack = cellbyCellStackALL;
    intrabyCellStack = intrabyCellStackALL;
    cellOrderStack = cellOrderStackALL;

    % dentNumbers = unique(intrabyCellStack(intrabyCellStack(:,5)~=0,5));            % Cells in the data have how many denticles?
    % maxDent_all = max(dentNumbers);                   % Largest number of denticles/cell in the dataset
    maxDent = maxDent_all;

    % save('currentDataSet.mat','cellbyCellStack','intrabyCellStack','cellOrderStack','eID','numerical_eID','genoID','maxDent_all','maxDent','noEmbryos','embryoIDs','thingtodo')


    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    formatSpec = '\n\n %s';
    fprintf(fileID,formatSpec,'embryoIDs');
    embryoStats = embryoIDs;
    formatSpec = '\n\n %d';
    fprintf(fileID,formatSpec,embryoStats);
    fclose(fileID);


    %     slopes_wrow;
    DenticleCalculations(genotype, thingtodo, cellbyCellStack, intrabyCellStack, cellOrderStack, eID, numerical_eID, genoID, maxDent, maxDent_all, totalgenotypes);


    disp 'Done with calculations for pooled dataset'


elseif thingtodo == 2,
    % Pooled and individual
    load('DataSet_ALL.mat');

    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    formatSpec = '\n %s %s %s\n';
    embryoStats_labels = ['embryoID ', 'no.cells ','no.denticles ','[cellsbyCol]','[dentsbyCol]'];
    fprintf(fileID,formatSpec,embryoStats_labels);
    fclose(fileID);

    statssave = [0 0 0 0];


    for k = 0:noEmbryos,

        if k == 0,    % The 'all case'
            % Identity data
            eID = 'all';
            genoID = sprintf('%1s_%2s',genotype,eID); %#ok<*NASGU>
            % genoID = sprintf('%1s_%2s_%3s',genotype,eID,timerun); %#ok<*NASGU>
            numerical_eID = 0;
            disp('All embryos')

            % Get all data
            cellbyCellStack = cellbyCellStackALL;
            intrabyCellStack = intrabyCellStackALL;
            cellOrderStack = cellOrderStackALL;

            maxDent_all = max(unique(intrabyCellStack(intrabyCellStack(:,5)~=0,5)));                             % Largest number of denticles/cell in the dataset
            maxDent = maxDent_all;

            % save('currentDataSet.mat','cellbyCellStack','intrabyCellStack','cellOrderStack','eID','numerical_eID','genoID','maxDent_all','maxDent','thingtodo')

        elseif k ~= 0,    % The 'individual embryo' case
            % Identity data
            IDtoGet = embryoIDs(k);
            disp(IDtoGet)
            eID = sprintf('%1$i',IDtoGet);
            genoID = sprintf('%1s_%2s',genotype,eID); %#ok<*NASGU>
            % genoID = sprintf('%1s_%2s_%3s',genotype,eID,timerun);
            numerical_eID = IDtoGet;

            % Get subset of data
            cellbyCellStack = cellbyCellStackALL(cellbyCellStackALL(:,1) == IDtoGet,:);      % Get embryo cbc data
            intrabyCellStack = intrabyCellStackALL(intrabyCellStackALL(:,1) == IDtoGet,:);   % Get embryo ibc data
            cellOrderStack = cellOrderStackALL(cellOrderStackALL(:,1) == IDtoGet,:);   % Get embryo ibc data

            maxDent_all = max(unique(intrabyCellStackALL(intrabyCellStackALL(:,5)~=0,5)));                             % Largest number of denticles/cell in the dataset
            maxDent = max(unique(intrabyCellStack(intrabyCellStack(:,5)~=0,5)));                             % Largest number of denticles/cell in the data subset

            EmbryoAvg(cellOrderStack,'Avg_CellOrder',numerical_eID,genotype);

            % save('currentDataSet.mat','cellbyCellStack','intrabyCellStack','cellOrderStack','eID','numerical_eID','genoID','maxDent','maxDent_all','thingtodo')


            EID_cbc = 1;
            Row_cbc = 2;
            Belt_cbc = 3;
            Cell_cbc = 4;
            DIC_cbc = 5;
            Length_cbc = 6;


            statsoncolumns = zeros(1,14);
            for currentrow = 1:7,
                thisDataset = cellbyCellStack(cellbyCellStack(:,Row_cbc) == (currentrow-1),DIC_cbc);

                % # of cells (length)
                statsoncolumns(1,currentrow) = length(thisDataset);
                % # of denticles (sum of elements)
                statsoncolumns(1,(currentrow + 7)) = sum(thisDataset);

                celldata = statsoncolumns(:,1:7); dentdata = statsoncolumns(:,8:14);
                statsmaxmin = [min(celldata(any(celldata,1))),max(celldata(any(celldata,1))),min(dentdata(any(dentdata,1))),max(dentdata(any(dentdata,1)))];
            end

            % stats_embryos(k,1) = statsmaxmin(1); stats_embryos(k,2) = statsmaxmin(2);
            % stats_embryos(k,3) = statsmaxmin(3); stats_embryos(k,4) = statsmaxmin(4);
            % stats_embryos(k,5) = size(cellbyCellStack,1); stats_embryos(k,6) = sum(cellbyCellStack(:,5),1);
            % save('stats.mat','stats_embryos');


            fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
            embryoStats = [embryoIDs(k) size(cellbyCellStack,1) sum(cellbyCellStack(:,5),1) statsoncolumns statsmaxmin];
            formatSpec = '\n %d     %d %d,  %d %d %d %d %d %d %d,   %d %d %d %d %d %d %d,       %d : %d;  %d : %d';
            fprintf(fileID,formatSpec,embryoStats);
            fclose(fileID);


            statsTemp(k,:) = [size(cellbyCellStack,1), sum(cellbyCellStack(:,5),1)];
            save('statstemp.mat', 'statsTemp');

        end

        % Run the other three scripts to do actual calculations
        DenticleCalculations(genotype, thingtodo, cellbyCellStack, intrabyCellStack, cellOrderStack, eID, numerical_eID, genoID, maxDent, maxDent_all, totalgenotypes);

        % Get all the data again to set up for the next round
        load('DataSet_ALL.mat');

    end

    load('statstemp.mat')
    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    netstats = [min(statsTemp(:,1)), max(statsTemp(:,1)),min(statsTemp(:,2)), max(statsTemp(:,2))];
    formatSpec = '\n \n %d - %d,   %d - %d ';
    fprintf(fileID,formatSpec,netstats);


    clear statsTemp
    disp 'Done with calculations for individual datasets'


elseif thingtodo == 3,
    % Individual only

    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    formatSpec = '\n %s %s %s\n';
    embryoStats_labels = ['embryoID ', 'no.cells ','no.denticles '];
    fprintf(fileID,formatSpec,embryoStats_labels);
    fclose(fileID);


    load('DataSet_ALL.mat');
    for k = 1:noEmbryos,
        % Identity data
        IDtoGet = embryoIDs(k);
        disp(IDtoGet)
        eID = sprintf('%1$i',IDtoGet);
        genoID = sprintf('%1s_%2s',genotype,eID); %#ok<*NASGU>
        % genoID = sprintf('%1s_%2s_%3s',genotype,eID,timerun);
        numerical_eID = IDtoGet;

        % Get subset of data
        cellbyCellStack = cellbyCellStackALL(cellbyCellStackALL(:,1) == IDtoGet,:);      % Get embryo cbc data
        intrabyCellStack = intrabyCellStackALL(intrabyCellStackALL(:,1) == IDtoGet,:);   % Get embryo ibc data
        cellOrderStack = cellOrderStackALL(cellOrderStackALL(:,1) == IDtoGet,:);   % Get embryo ibc data

        maxDent_all = max(unique(intrabyCellStackALL(intrabyCellStackALL(:,5)~=0,5)));                             % Largest number of denticles/cell in the dataset
        maxDent = max(unique(intrabyCellStack(intrabyCellStack(:,5)~=0,5)));                             % Largest number of denticles/cell in the dataset

        EmbryoAvg(cellOrderStack,'Avg_CellOrder',numerical_eID,genotype);

        % save('currentDataSet.mat','cellbyCellStack','intrabyCellStack','cellOrderStack','eID','numerical_eID','genoID','maxDent','maxDent_all','thingtodo')


        EID_cbc = 1;
        Row_cbc = 2;
        Belt_cbc = 3;
        Cell_cbc = 4;
        DIC_cbc = 5;
        Length_cbc = 6;


        statsoncolumns = zeros(1,14);
        for currentrow = 1:7,
            thisDataset = cellbyCellStack(cellbyCellStack(:,Row_cbc) == (currentrow-1),DIC_cbc);

            % # of cells (length)
            statsoncolumns(1,currentrow) = length(thisDataset);
            % # of denticles (sum of elements)
            statsoncolumns(1,(currentrow + 7)) = sum(thisDataset);

            celldata = statsoncolumns(:,1:7); dentdata = statsoncolumns(:,8:14);
            statsmaxmin = [min(celldata(any(celldata,1))),max(celldata(any(celldata,1))),min(dentdata(any(dentdata,1))),max(dentdata(any(dentdata,1)))];
        end

        stats_cellmin(k,1) = statsmaxmin(1); stats_cellmax(k,2) = statsmaxmin(2); stats_dentmin(k,3) = statsmaxmin(3); stats_dentmax(k,4) = statsmaxmin(4);


        fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
        embryoStats = [embryoIDs(k) size(cellbyCellStack,1) sum(cellbyCellStack(:,5),1) statsoncolumns statsmaxmin];
        formatSpec = '\n %d     %d %d,  %d %d %d %d %d %d %d,   %d %d %d %d %d %d %d,       %d : %d;  %d : %d';
        fprintf(fileID,formatSpec,embryoStats);
        fclose(fileID);


        % Run the other three scripts to do actual calculations
        DenticleCalculations(genotype, thingtodo, cellbyCellStack, intrabyCellStack, cellOrderStack, eID, numerical_eID, genoID, maxDent, maxDent_all, totalgenotypes);

        % Get all the data again to set up for the next round
        load('DataSet_ALL.mat');

    end

    fileID = fopen(SpiffyName('txt', 'Stats',genotype),'a');
    netstats = [min(stats_cellmin(:,1)), max(stats_cellmin(:,2)),min(stats_cellmin(:,3)), max(stats_cellmin(:,4))];
    formatSpec = '\n %d - %d,   %d - %d ';
    fprintf(fileID,formatSpec,netstats);
    fclose(fileID);

    disp 'Done with calculations for individual datasets'


end 
