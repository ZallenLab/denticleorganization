


function DenticleCalculations(genotype, thingtodo, cellbyCellStack, intrabyCellStack, cellOrderStack, eID, numerical_eID, genoID, maxDent, maxDent_all, totalgenotypes)

%% Calculations to do to denticle organization data obtained from DenticleDist
%{
  Calculations based on row parameters
      - Data set statistics
      - Basic sort-data by row
      - DIC frequency
      - Cell size frequency
      - Relative postion across cell face
      - How close to centered are denticles?
      - Cell overlay
      - Cumulative distance to the next denticle (independent of cell boundaries)

  ------------------------
  Calculations based on cells only (or cell pairs)
      - Inter pair positioning (placemnt of outermost denticles in neighboring
      cells relative to the cell DV edge)
      - Cell order organization and misc calculations

  ------------------------
  Calculations based on denticle in cell numbers, independent of row 
      - Sort by DIC 
      - Combine data from cells with 5 or more denticles ('many' sort of category) 
      - DVlen vs Intra (sorted by DIC) 
      - DIC vs DVlen 
      - Cell size frequency


  ------------------------
  Input files: (genotype_*byCellStack.csv from DenticleDist_akls*.m)
      CellbyCell
          (1) embryoID    (2) row         (3 belt)        (4 cell)
          (5) dentincell (6) Dvlen       (7) dentEdgeL   (8) dentEdgeR   (9) AdjL
          (10) AdjR (11) Intra2     (12) Intra3     (13) Intra4     (14) Intra5
          (15) Intra6 (16) Intra7     (17) Intra8     (18) Intra9

      IntrabyCell
          (1) embryoID    (2) row         (3 belt)        (4 cell)
          (5) dentincell (6) Dvlen       (7) Intra

      CellOrder
          (1) embryoID    (2) row         (3 belt)        (4 cell)
          (5) dentincell (6) Dvlen       (7) AdjL        (8) dentEdgeL   (9) Intra2
          (10) Intra3 (11) Intra4     (12) Intra5     (13 Intra6)     (14) Intra7
          (15) Intra8 (16) Intra9     (17) dentEdgeR  (18 AdjR)
%}



% load('UserInputData.mat');
% load('currentDataSet.mat');
% load('DataSet_ALL.mat');


%% SET PARAMETERS AND TITLES 
stdbin_DVlen = 4;
standardBin = 1:5:40;

EID = 1;
ROW = 2;
BELT = 3;
CELLNUMBER = 4;
DICNUMBER = 5;
DVLENGTH = 6;



INTRA = 7;  % intra by cell
INTRA1 = 11; % cell by cell

% cell order stack
adjL_co = 7; 
deL_co = 8; 
intra2_co = 9;

% Cell by cell
DE_L = 7;
DE_R = 8;

adjL_cbc = 9; 
adjL_cbc = 10;





%% ------------------
% Sort by Row - output file with all data (DIC, DVlength, dent-edge, inter, intra) sorted into groups by cell row (column)

maxNeeded = 2*(size(cellbyCellStack,1));   % Find max number of data rows in the entire dataset; this *should* be the maximum value, ever.

blocksize = 8;  % Number of columns that will be needed per data 'block' - ie. Intra, DVlen + 1 for a blank separator
byrow = zeros(maxNeeded,5*blocksize);
columnplace = 0;

for k = 0:max(unique(cellbyCellStack(:,2))), % cell column counter
    data = cellbyCellStack(cellbyCellStack(:,2) == k,:);       % logical indexing to get data from embryo k, row m only
    dataIntra = intrabyCellStack(intrabyCellStack(:,2) == k,:);
    
    h = size(data,1);   
    % Becuase the number of rows will change for each dataset (cell #/column is never the same, you need to specify the length of the block of data you are putting into the whole matrix (which already has the maximum size needed))
    
    byrow(1:h,(1 +columnplace)) = data(:,5);                                      % DIC
    byrow(1:h,(1 +blocksize +columnplace)) = data(:,6);                           % DVlength
    byrow(1:2*h,(1 + 2*blocksize +columnplace)) = [data(:,7);data(:,8)];          % d-e
    byrow(1:2*h,(1 +3*blocksize +columnplace)) = [data(:,9);data(:,10)];          % Inter/Next cell data
    
    h = size(dataIntra,1);  
    % Same as above, but the size of the data for intra varies with the DIC # (where the other parameters always have a constant number of data points)
    
    byrow(1:h,(1 +4*blocksize +columnplace)) = dataIntra(:,7);                    % Intra/same cell data
    
    columnplace = columnplace + 1;    % Each row that gets copied, shift over to a new column to get in the right place for the next data set 
    
end

% % For datasets that include 'individual embryo' separations, take the average for each embryo 
% if thingtodo >= 2,
%     EmbryoAvg(byrow,'Averages_Set_sortbyRow',numerical_eID, genotype);
% end

% Save it all with column headers
minrow = 0;
maxrow = max(6,max(unique(cellbyCellStack(:,2))));

titles_byrow = [arrayfun(@(x)sprintf('DIC row %i',x),minrow:6,'uniformOutput',false),'-',arrayfun(@(x)sprintf('DVlength row %i',x),minrow:6,'uniformOutput',false),'-',...
               arrayfun(@(x)sprintf('dentEdge row %i',x),minrow:6,'uniformOutput',false),'-',arrayfun(@(x)sprintf('Inter/NextCell row %i',x),minrow:6,'uniformOutput',false),'-',...
               arrayfun(@(x)sprintf('Intra/SameCell row %i',x),minrow:6,'uniformOutput',false),'-'];

AddHeaders(SpiffyName('csv','SortbyRow', genoID),titles_byrow,byrow);





%% ------------------------------------
% Intra and Intercellular Spacing lists
intraall = cellbyCellStack(:,11:end);  % Take all the intra values
intraall = reshape(intraall,[],1);      % Reorganize into a single column
intraall = intraall(intraall ~= 0);     % Remove the zeros

interall = cellbyCellStack(:,9:10);     % Take all the inter values
interall = reshape(interall,[],1);      % Reorganize into a single column
interall = interall(interall ~= 0);     % Remove all the zeros

a = size(intraall,1);                   % Size data to make a combined matrix
b = size(interall,1);

ealists = zeros(max(a,b),2);            % New matrix with 2 cols, the max needed # of rows
ealists(1:a,1) = intraall;
ealists(1:b,2) = interall;

titles_intrainterlist = {'Intracellular Spacing', 'Intercellular Spacing'};
dlmwrite(SpiffyName('csv','IntraInterStack', genotype), ealists);
AddHeaders(SpiffyName('csv','IntraInterList', genoID), titles_intrainterlist, ealists);







%% -----------------
% Sort all the data available by # of denticles in the cell

% Sort by DIC (to maximum number)
columnplace = 0;  % Column shifter
blocksizeDIC = maxDent_all +1;  % Number of columns that will be needed per data 'block' - ie. Intra, DVlen for denticles 1:maxDent + 1 for a blank separator

byDIC = zeros(1000,5*blocksizeDIC);    % Preallocation

    for k = 1:maxDent_all,

        dicall_data = cellbyCellStack(cellbyCellStack(:,5) == k,:);
        dicall_data2 = intrabyCellStack(intrabyCellStack(:,5) == k,:);

        h = size(dicall_data,1);

        byDIC(1:h,(1 +columnplace)) = dicall_data(:,5);     %DIC
        byDIC(1:h,(1 +blocksizeDIC +columnplace)) = dicall_data(:,6);        % DVlen
        byDIC(1:2*h,(1 +2*blocksizeDIC +columnplace)) = [dicall_data(:,7);dicall_data(:,8)];  % d-e
        byDIC(1:2*h,(1 +3*blocksizeDIC +columnplace)) = [dicall_data(:,9);dicall_data(:,10)]; % Next/Inter

        h = size(dicall_data2,1);

        byDIC(1:h,(1 +4*blocksizeDIC +columnplace)) = dicall_data2(:,7);     % Same/Intra

        columnplace = columnplace + 1;

    end

    if thingtodo >= 2,  % if you're doing pool + individual or individual only cases
      EmbryoAvg(byDIC, 'Average_Set_sortbyDIC', numerical_eID, genotype);
    end 
    
    % Save
    titles_byDIC = [arrayfun(@(x)sprintf('DIC %i',x),1:maxDent_all,'uniformOutput',false),'-',arrayfun(@(x)sprintf('DVlength %i',x),1:maxDent_all,'uniformOutput',false),'-',...
                     arrayfun(@(x)sprintf('dentEdge %i',x),1:maxDent_all,'uniformOutput',false),'-',arrayfun(@(x)sprintf('Inter/NextCell %i',x),1:maxDent_all,'uniformOutput',false),'-',...
                     arrayfun(@(x)sprintf('SameCell/Intra %i',x),1:maxDent_all,'uniformOutput',false),'-'];   

    AddHeaders(SpiffyName('csv', 'SortbyDIC', genoID), titles_byDIC, byDIC)
    





%% ----------------------------
% Relative position in the cell
    % Calculates the postion (as a percentage of the DV cell length) of each
    % denticle in the cell. Uses DIC value to sort the data into blocks.

sumofcell = cumsum(cellOrderStack(:,8:size(cellOrderStack,2) -1),2);        
    % Sum across the row (this one is weird for some reason) - get 'positional' measurement
    % This sum is from dent-edge LEFT, through all the intra measurements
    % (col 9 to the 3rd from last column) up to dent-edge RIGHT (2nd to
    % last column). 

relativeposition = zeros(size(sumofcell));
    % Calculate the fraction of the total (column 10) that the previous m
    % columns equal
    for row = 1:size(sumofcell,1), 
       relativeposition(row,:) = 100*(sumofcell(row,:)/sumofcell(row,(size(sumofcell,2))));       % fraction taken based on the final column of the sumofcell matrix, which should be the total
    end
    
    relativeposition = [cellOrderStack(:,5), relativeposition];     % Add DIC value to the data table


relativepositionzeroed = zeros(size(relativeposition));
    % Zero out excess intra columns
    for row = 1:size(relativeposition,1),
        relativepositionzeroed(row,:) = relativeposition(row,:);
        
        if relativeposition(row,1) == maxDent_all, 
            continue
            
        elseif relativeposition(row,1) < maxDent_all,
            relativepositionzeroed(row,(relativeposition(row,1) +2):(size(relativeposition,2) -1)) = 0;

        elseif relativeposition(row,1) > maxDent_all,
            error('Im borked')
        end
    end
    
    
% Individual data points for each denticle in each cell
totalrpSet = zeros(size(relativeposition,1),0);   % Zero columns, the max number of rows you could possibly need

for dents = 1:maxDent,
    relpossubset = relativeposition(relativeposition(:,1) == dents,:);  % only data for x denticles
    
    % if there aren't x denticles, then make a zeros array to pad the final output
    if isempty(relpossubset) == 1, 
        relpossubset = zeros(size(relativeposition,1),dents+1);
    end         
    
    % Get rid of the 'dic' data that is in column 1
    relposSet = relpossubset(:,2:dents+1);  
    
    % add rows of zeros to the end if needed so this set can mesh with the total set
    if size(totalrpSet,1) > size(relposSet,1),
        relposSet = padarray(relposSet,(size(totalrpSet,1) - size(relposSet,1)),0,'post');
    end
    
    % add this new data to the growing total dataset
    totalrpSet = [totalrpSet relposSet]; %#ok<AGROW>
end


%Save
titles_relposition_pts = '1_d1';
  for dicval = 2:maxDent_all,
      str = strcat(num2str(dicval),'_',sprintf('d%i',1));
      for j = 2:min(dicval,maxDent_all),
          str = [str,{strcat(num2str(dicval),'_', sprintf('d%i',j))}]; %#ok<AGROW>
      end
      titles_relposition_pts = [titles_relposition_pts,str]; %#ok<AGROW>
  end

AddHeaders(SpiffyName('csv','RelativePosition', genoID), titles_relposition_pts, totalrpSet);






%% ------------------------------------------    
% DV length vs Intra, spacing separated by DIC 
spacingvsLength_byDIC = [0,0,0,0; intrabyCellStack(:,1), intrabyCellStack(:,5:7)];
titles_intraDVbyDIC = ['embryoID', 'DIC' ,'DVlength',arrayfun(@(x)sprintf('Intra%i_DIC',x),1:maxDent,'uniformOutput',false)];
AddHeaders(SpiffyName('csv','spacingvsLength_byDIC',genoID),titles_intraDVbyDIC,ColumnShifter(spacingvsLength_byDIC,2,3,1));
dlmwrite(SpiffyName('csv','spacingvsLength_noheaders',genoID),ColumnShifter(spacingvsLength_byDIC,2,3,1));



%% --------------------------------------------
% DVlength vs dent-edge for 1DIC cells
dvlength_edge = [cellbyCellStack(cellbyCellStack(:,DICNUMBER)== 1, DVLENGTH), cellbyCellStack(cellbyCellStack(:,DICNUMBER)== 1,DE_L); cellbyCellStack(cellbyCellStack(:,DICNUMBER)== 1,DVLENGTH), cellbyCellStack(cellbyCellStack(:,DICNUMBER)== 1,DE_R)];
titles_DV_edge = {'DVlength','dent-edge'};

AddHeaders(SpiffyName('csv','DVlength_dentedge_1dic',genoID), titles_DV_edge, dvlength_edge)





%% -----------------------------------------------
% Organize for cell length versus denticle number

titles_cellbyDIC = ['denticles/cell','-',arrayfun(@(x)sprintf('Cell length %i',x),1:maxDent,'uniformOutput',false)];

lengthanddic = [cellbyCellStack(:,5), cellbyCellStack(:,6)];
AddHeaders(SpiffyName('csv','cellLengths_by_DIC',genoID), titles_cellbyDIC, ColumnShifter(lengthanddic,1,1,(length(maxDent))));



%% -----------------------------------------------
%         predicted_spacing = cellbyCellStack(cell_k,DVlength)/(cellbyCellStack(cell_k,DIC) + 1);

if thingtodo < 3 && numerical_eID == 0,

   for typeofanimal = 1:2,  % typeofanimals are embryos or larvae 
      basisequations = {'embryoEquations','larvalEquations'};

      % Root mean squared error - WT equations
      RMS_calculated = zeros(size(cellbyCellStack,1),5);
      % RMS_calculated positions
      % titles_rmse = {'eID','DIC','DVLength','predicted distance','rmse'};

      tempdata = cellbyCellStack(cellbyCellStack(:,DICNUMBER) > 1,:);

      for cell_k = 1:size(tempdata,1),

         if typeofanimal == 1,     % ie embryos
            % WT EMBYRO SPACING EQUATIONS
            if tempdata(cell_k,DICNUMBER) == 2,
               predicted_spacing = 0.34*tempdata(cell_k,DVLENGTH) + 0.77;
            elseif tempdata(cell_k,DICNUMBER) == 3,
               predicted_spacing = 0.24*tempdata(cell_k,DVLENGTH) + 0.64;
            elseif tempdata(cell_k,DICNUMBER) == 4,
               predicted_spacing = 0.22*tempdata(cell_k,DVLENGTH) + 0.14;
            elseif tempdata(cell_k,DICNUMBER) == 5,
               predicted_spacing = 0.16*tempdata(cell_k,DVLENGTH) + 0.16;
            elseif tempdata(cell_k,DICNUMBER) >= 6,
               predicted_spacing = tempdata(cell_k,DVLENGTH)/(tempdata(cell_k,DICNUMBER) + 1);
            end   
            
            label = 'embryoEquations';

       
         else %if typeofanimal == 2,   % ie larvae
            % LARVAL SPACING EQUATIONS
            if tempdata(cell_k,DICNUMBER) == 2,
               predicted_spacing = 0.37*tempdata(cell_k,DVLENGTH) + 0.63;
            elseif tempdata(cell_k,DICNUMBER) == 3,
               predicted_spacing = 0.27*tempdata(cell_k,DVLENGTH) + 0.23;
            elseif tempdata(cell_k,DICNUMBER) == 4,
               predicted_spacing = 0.22*tempdata(cell_k,DVLENGTH) + 0.1;
            elseif tempdata(cell_k,DICNUMBER) == 5,
               predicted_spacing = 0.18*tempdata(cell_k,DVLENGTH) + 0.07;
            elseif tempdata(cell_k,DICNUMBER) == 6,
               predicted_spacing = 0.14*tempdata(cell_k,DVLENGTH) + 0.07;
            elseif tempdata(cell_k,DICNUMBER) == 7,
               predicted_spacing = 0.13*tempdata(cell_k,DVLENGTH) + 0.18;
            elseif tempdata(cell_k,DICNUMBER) == 8,
               predicted_spacing = 0.11*tempdata(cell_k,DVLENGTH) + 0.25;
            elseif tempdata(cell_k,DICNUMBER) == 9,
               predicted_spacing = 0.09*tempdata(cell_k,DVLENGTH) + 0.42;
            elseif tempdata(cell_k,DICNUMBER) == 10,
               predicted_spacing = 0.078*tempdata(cell_k,DVLENGTH) + 0.59;
            elseif tempdata(cell_k,DICNUMBER) >= 11,
               predicted_spacing = tempdata(cell_k,DVLENGTH)/(tempdata(cell_k,DICNUMBER) + 1);
            end   
            label = 'larvalEquations';
         end 
         

         observed_spacing = tempdata(cell_k,INTRA1:(INTRA1+tempdata(cell_k,DICNUMBER)-2));
         
         for intravals = 1:length(observed_spacing),
             mean_square(intravals) = (observed_spacing(intravals) - predicted_spacing)^2;
         end
         
         RMS_calculated(cell_k,1) = tempdata(cell_k,EID);
         RMS_calculated(cell_k,2) = tempdata(cell_k,DICNUMBER);
         RMS_calculated(cell_k,3) = tempdata(cell_k,DVLENGTH);
         RMS_calculated(cell_k,4) = predicted_spacing;
         RMS_calculated(cell_k,5) = sqrt(mean(mean_square));

      end 
      
        dlmwrite(SpiffyName('csv', label, genoID,'RMSE_self'),RMS_calculated);
        dlmwrite(SpiffyName('csv', label, genoID,'RMSE_sorted_self'),ColumnShifter(RMS_calculated,2,4,1));

        RMS_GT = [totalgenotypes*ones(size(RMS_calculated,1),1), RMS_calculated];
        for dentnumber= 2:6, 
            dentStr = num2str(dentnumber);
            temp = RMS_GT(RMS_GT(:,3)==dentnumber,:);
            
            dlmwrite(SpiffyName('csv', label, dentStr,'RMSE'),temp,'-append');
        end

      
   end 

end



% % Sort to give data on an embryo-by-embryo basis
% % titles_rmse_byEmbryo = {'eID','DIC','rmse'};

% columnplace = 0;  % Column shifter
% columnplace2 = 0;

% blocksizeRMS = noEmbryos +1;  % Number of columns that will be needed per data 'block' - ie. Intra, DVlen for denticles 1:maxDent + 1 for a blank separator

% RMS_subset = [RMS_calculated(:,1:2) RMS_calculated(:,5)];

% for dentno = 2:maxDent,
%     rms_sorttemp1 = RMS_subset(RMS_subset(:,2) == dentno,:);
    
%     for embryonumber = 1:noEmbryos,
        
%         rms_sorttemp = rms_sorttemp1(rms_sorttemp1(:,1) == embryoIDs(embryonumber),:);
        
%         h = size(rms_sorttemp,1);
        
%         RMS_byEmbryo(1:h,(1 +columnplace)) = rms_sorttemp(:,1);     % eID
%         RMS_byEmbryo(1:h,(1 +blocksizeRMS +columnplace)) = rms_sorttemp(:,2);        % dic
%         RMS_byEmbryo(1:h,(1 +2*blocksizeRMS +columnplace)) = rms_sorttemp(:,3);     % rms
        
%         columnplace = columnplace + 1;
        
%         RMS_byEmbryo_rmsonly(1:h,(1 +columnplace2)) = rms_sorttemp(:,3);     % rms
%         columnplace2 = columnplace2 + 1;
        
%     end
    
%     columnplace = size(RMS_byEmbryo,2) + 1;
%     columnplace2 = columnplace2 + 1;
    
% end

% dlmwrite(SpiffyName('csv','RMSE_byEmbryo_wtequations',genoID), RMS_byEmbryo_rmsonly);

% end



%% -------------
% Root mean square error
% using intra by cell 

% if thingtodo < 3 && numerical_eID == 0, 
    
%     tempdata = intrabyCellStack(intrabyCellStack(:,DICNUMBER) >1,:);
    
%     for cell_k = 1:size(tempdata,1), 
%         cell_actualspacing = tempdata(cell_k,INTRA);
        
%         % Predictions WT EMBYRO + LARVA EQUATIONS   
%         if tempdata(cell_k,DICNUMBER) == 2,
%             predicted_spacing = 0.35*tempdata(cell_k,DVLENGTH) + 0.71;
%         elseif tempdata(cell_k,DICNUMBER) == 3,
%             predicted_spacing = 0.26*tempdata(cell_k,DVLENGTH) + 0.50;
%         elseif tempdata(cell_k,DICNUMBER) ==  4,
%             predicted_spacing = 0.21*tempdata(cell_k,DVLENGTH) + 0.29;
%         elseif tempdata(cell_k,DICNUMBER) ==  5,
%             predicted_spacing = 0.18*tempdata(cell_k,DVLENGTH) + 0.09;
%         elseif tempdata(cell_k,DICNUMBER) >=  6,
%             predicted_spacing = tempdata(cell_k,DVLENGTH)/(tempdata(cell_k,DICNUMBER) + 1);
%         end
        
%         % Predictions - wt embryo equations
% %         if tempdata(cell_k,DICNUMBER) == 2,
% %             predicted_spacing = 0.34*tempdata(cell_k,DVLENGTH) + 0.77;
% %         elseif tempdata(cell_k,DICNUMBER) == 3,
% %             predicted_spacing = 0.24*tempdata(cell_k,DVLENGTH) + 0.64;
% %         elseif tempdata(cell_k,DICNUMBER) == 4,
% %             predicted_spacing = 0.22*tempdata(cell_k,DVLENGTH) + 0.14;
% %         elseif tempdata(cell_k,DICNUMBER) == 5,
% %             predicted_spacing = 0.22*tempdata(cell_k,DVLENGTH) + 0.14;
% %         elseif tempdata(cell_k,DICNUMBER) >= 6,
% %             predicted_spacing = tempdata(cell_k,DVLENGTH)/(tempdata(cell_k,DICNUMBER) + 1);
% %         end
% % 
       
        
%         rmsError(cell_k,:) = [tempdata(cell_k,DICNUMBER), tempdata(cell_k,DVLENGTH), tempdata(cell_k,INTRA) - predicted_spacing];
        
%     end 
    
    
%     for j = 1:maxDent,
%         temp = rmsError(rmsError(:,1) == j,:);
%         rootmeansqerr(j) = sqrt(mean(temp(:,3).^2));
%     end 
    
%    dlmwrite(SpiffyName('csv','rmse_individual','gtrows'),rootmeansqerr,'-append')
   

    
% end 





%% ---------------------
% Noise-in-data 
  % Based off of equations in Swain et al 2002, PNAS  and Elowitz et al, 2002, Science
  % .* operator gives element-wise multiplication of two matrices
  % .^ operator gives element-wise raised to the nth power of the matrix

for dentnumber = 2:maxDent, 

  separations = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber,INTRA);
  celllengths = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber, DVLENGTH);

  if length(separations) < 2, 
    continue

  else
    n_int_square = ((mean((separations - celllengths).^2))) / (2*mean(separations)*mean(celllengths)); 
    n_int = sqrt(n_int_square);

    n_ext_square = (mean(separations.*celllengths) - (mean(separations)*mean(celllengths))) / (mean(separations)*mean(celllengths));
    n_ext = sqrt(n_ext_square);

    n_tot_square = (mean((separations.*separations) + (celllengths.*celllengths)) - (2*mean(separations)*mean(celllengths))) / ((2*mean(separations)*mean(celllengths)));
    n_tot = sqrt(n_tot_square);


    dentStr = num2str(dentnumber);
    dlmwrite(SpiffyName('csv',dentStr,'Noise'),[0,numerical_eID, n_tot, n_int, n_ext],'-append');
  %   dlmwrite(SpiffyName('csv',dentStr,'Noise_CHECK'),[numerical_eID, n_tot_square, (n_int_square + n_ext_square)],'-append');
  end 

end 





%% ---------------------
% Correlation coefficient calculations 

% ADD: iff statement to do only when runing individual embryo groups

for dentnumber = 2:maxDent, 
  separations = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber,INTRA);
  celllengths = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber, DVLENGTH);

  if length(separations) < dentnumber, 
    continue

  else
    [r,p] = corrcoef(separations,celllengths);
    r2 = r.^2;

    dentStr = num2str(dentnumber);
    dlmwrite(SpiffyName('csv',dentStr,'CorrCoef'),[totalgenotypes,numerical_eID, r(1,2), r2(1,2), p(1,2)],'-append');
  end 


end



%% -------- 
% Calculate X error to plot both X and Y error on separation v. length plots

% [totalgenotypes numerical_eID std()]
if numerical_eID ~= 0,
    for dentnumber = 2:maxDent, 
      dentStr = num2str(dentnumber);

      separations = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber,INTRA);
      celllengths = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber, DVLENGTH);

      sd_err = [mean(celllengths) std(celllengths) mean(separations) std(separations) length(separations)];

      dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSD'),[totalgenotypes,numerical_eID, dentnumber, sd_err],'-append');
      dlmwrite(SpiffyName('csv','all', 'spacingvsLength_XerrorSD'),[totalgenotypes,numerical_eID, dentnumber, sd_err],'-append');

      sem_err = [mean(celllengths) std(celllengths)/sqrt(length(celllengths)) mean(separations) std(separations) length(separations)];

      dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSEM'),[totalgenotypes,numerical_eID, dentnumber, sem_err],'-append');
      dlmwrite(SpiffyName('csv','all', 'spacingvsLength_XerrorSEM'),[totalgenotypes,numerical_eID, dentnumber, sem_err],'-append');

      %   dlmwrite(SpiffyName('csv',dentStr,'Noise_CHECK'),[numerical_eID, n_tot_square, (n_int_square + n_ext_square)],'-append');

    end 


elseif numerical_eID == 0, 
    for dentnumber = 2:maxDent, 
        dentStr = num2str(dentnumber);

      separations = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber,INTRA);
      celllengths = intrabyCellStack(intrabyCellStack(:,DICNUMBER) == dentnumber, DVLENGTH);

      sd_err = [mean(celllengths) std(celllengths) mean(separations) std(separations) length(separations)];

      dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSD'),[totalgenotypes,numerical_eID, dentnumber, sd_err],'-append');
      dlmwrite(SpiffyName('csv','all', 'spacingvsLength_all_XerrorSD'),[totalgenotypes,numerical_eID, dentnumber, sd_err],'-append');

      sem_err = [mean(celllengths) std(celllengths)/sqrt(length(celllengths)) mean(separations) std(separations) length(separations)];

      dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSEM'),[totalgenotypes, dentnumber, sem_err],'-append');
      dlmwrite(SpiffyName('csv','all', 'spacingvsLength_all_XerrorSEM'),[totalgenotypes, dentnumber, sem_err],'-append');

      %   dlmwrite(SpiffyName('csv',dentStr,'Noise_CHECK'),[numerical_eID, n_tot_square, (n_int_square + n_ext_square)],'-append');

    end 
end 




%% ---------
% Edge portion 

if numerical_eID == 0, 
  edgetotal = cellbyCellStack(:,DE_L) + cellbyCellStack(:,DE_R);
  edgeratio = edgetotal./cellbyCellStack(:,DVLENGTH);
  edgedata = [cellbyCellStack(:,DICNUMBER), edgetotal, edgeratio];


  for dentnumber = 2:maxDent, 
    dentStr = num2str(dentnumber);

    dlmwrite(SpiffyName('csv', dentStr, 'edgefraction_all'),[totalgenotypes*(ones(size(edgedata(edgedata(:,1) == dentnumber,:),1),1)), edgedata(edgedata(:,1) == dentnumber,:)],'-append');
  end
end



disp('finished. next?')





