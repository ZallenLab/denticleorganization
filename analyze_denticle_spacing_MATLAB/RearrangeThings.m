

function RearrangeThings(genotypes)



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
deL_cbc = 7; 
deL_cbc = 8; 
adjL_cbc = 9; 
adjL_cbc = 10;



% Make denticles/cell, cell size conglomerates

dentpercell = zeros(5000, length(genotypes));
cellsize = zeros(5000, length(genotypes));
intra2DIC = zeros(5000, length(genotypes));


for k = 1:length(genotypes),
	genotype = genotypes{k};
	
	temp = dlmread(sprintf('%1$s_%2$s.%3$s',genotype,'cellOrderStack','csv'));
	temp2 = dlmread(sprintf('%1$s_%2$s.%3$s',genotype,'intrabyCellStack','csv'));	


	dentpercell(1:length(temp),k) = temp(:,DICNUMBER);
	cellsize(1:length(temp),k) = temp(:,DVLENGTH);
	intra2DIC(1:length(temp2(temp2(:,DICNUMBER)==2)),k) = temp2(temp2(:,DICNUMBER)==2,INTRA);

end 


dlmwrite(SpiffyName('csv','allgenotypes','DenticlesperCell'),dentpercell);
dlmwrite(SpiffyName('csv','allgenotypes','CellSize'),cellsize);
dlmwrite(SpiffyName('csv','allgenotypes','Intra_2DIC'),intra2DIC);











% Rearrange Noise_n files
for dentno = 2:6, 
	dentStr = num2str(dentno);
	temp = dlmread(SpiffyName('csv',dentStr,'Noise'));

	% columns by genotype
	dlmwrite(SpiffyName('csv',dentStr,'Noise_bygenotype'), ColumnShifter(temp,1,2,3));
	

	% (x,y) plots
	temp2 = [temp(:,1), temp(:,4:5)];
	dlmwrite(SpiffyName('csv',dentStr,'Noise_ext'), ColumnShifter(temp2,1,2,1));

	temp2pad = [temp(:,1),temp(:,5),temp(:,4)];
	dlmwrite(SpiffyName('csv',dentStr,'Noise_int'), ColumnShifter(temp2pad,1,2,1));
    
    % x, shifted plots
    temp4 = [temp(:,1),temp(:,3)];
	dlmwrite(SpiffyName('csv',dentStr,'Noise_tot'), ColumnShifter(temp4,1,1,1));
	

end 







clear temp temp2 temp3 temp4

% Rearrange CorrCoef files
% [gt,numerical_eID, r, r2, p]

for dentno = 2:6, 
	dentStr = num2str(dentno);
	temp = dlmread(SpiffyName('csv',dentStr,'CorrCoef'));

	% columns by genotype
	dlmwrite(SpiffyName('csv',dentStr,'CorrCoef_bygenotype'), ColumnShifter(temp,1,2,3));
	

	% (x,y) plot
	temp2 = [temp(:,1), temp(:,3:4)];
	dlmwrite(SpiffyName('csv',dentStr,'CorrCoef_r2'), ColumnShifter(temp2,1,2,1));

	temp2pad = [temp(:,1),temp(:,4),temp(:,3)];
	dlmwrite(SpiffyName('csv',dentStr,'CorrCoef_r'), ColumnShifter(temp2pad,1,2,1));
	

end 






clear temp temp2 temp3 temp4

% Rearrange separationVlength_Xerror files
% [gt,numerical_eID, dentnumber, sd_err(0 0 0 0 0 0 )]

% temp = dlmread(SpiffyName('csv','all','spacingvsLength_XerrorSD'));
% temp = temp(temp(:,3) < 7,:);
% 
% % columns by denticle number
% temp2 = ColumnShifter(temp,3,5,3);
% dlmwrite(SpiffyName('csv','all','spacingvsLength_Xerror_byDIC'), temp2);
% 
% % denticle number then genotype
% temp3 = ColumnShifter(temp2,1,5,15);
% dlmwrite(SpiffyName('csv','DIC_GT','separationVlength_Xerror'), temp3);


testdir = dir(sprintf('%1$s*.%2$s','spacingvsLength_XerrorSD','csv'));
if isempty(testdir) == 1,
	0

else
	for dentno = 2:6, 
		dentStr = num2str(dentno);

		temp = dlmread(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSD'));

		% columns by genotype
		dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSD_byDIC'), ColumnShifter(temp,1,5,3));
	end 



	for dentno = 2:6, 
		dentStr = num2str(dentno);
		temp2 = dlmread(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSEM'));

		% columns by genotype
		dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_XerrorSEM_byDIC'), ColumnShifter(temp2,1,5,3));
	end 




	 for dentno = 2:6, 
		dentStr = num2str(dentno);
		temp2pad = dlmread(SpiffyName('csv','RMSE_embryoEquations',dentStr));

		% columns by genotype
		dlmwrite(SpiffyName('csv',dentStr,'RMSE_embryoEquations_bygenotype'), ColumnShifter(temp2pad,1,5,1));
     end 
    
     
	 for dentno = 2:6, 
		dentStr = num2str(dentno);
		temp2pad = dlmread(SpiffyName('csv','RMSE_larvalEquations',dentStr));

		% columns by genotype
		dlmwrite(SpiffyName('csv',dentStr,'RMSE_larvalEquations_bygenotype'), ColumnShifter(temp2pad,1,5,1));
	end 


end 
 
 clear temp temp2 temp3 temp4



% 
% testdir = dir(sprintf('%1$s*.%2$s','spacingvsLength_all_XerrorSD','csv'));
% if isempty(testdir) == 1,
%     continue
% else
% 	for dentno = 2:6, 
% 		dentStr = num2str(dentno);
% 
% 		temp = dlmread(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSD'));
% 
% 		% columns by genotype
% 		dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSD_byDIC'), ColumnShifter(temp,1,5,3));
% 	end 
% 
% 
% 
% 	for dentno = 2:6, 
% 		dentStr = num2str(dentno);
% 		temp2 = dlmread(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSEM'));
% 
% 		% columns by genotype
% 		dlmwrite(SpiffyName('csv',dentStr,'spacingvsLength_all_XerrorSEM_byDIC'), ColumnShifter(temp2,1,5,3));
% 	end 
% 
% 
% 
% 
% 	 for dentno = 2:6, 
% 		dentStr = num2str(dentno);
% 		temp3 = dlmread(SpiffyName('csv',dentStr,'RMSE'));
% 
% 		% columns by genotype
% 		dlmwrite(SpiffyName('csv',dentStr,'RMSE_bygenotype'), ColumnShifter(temp3,1,5,1));
% 	end 
% 
%  
% end 
%  
 clear temp temp2 temp3 temp4

% Rearrange Edge files
% [gt,dic, tot, %]

for dentno = 2:6, 
	dentStr = num2str(dentno);
	temp = dlmread(SpiffyName('csv', dentStr, 'edgefraction_all'));

	% columns by genotype
	dlmwrite(SpiffyName('csv',dentStr,'edgefraction_bygenotype'), ColumnShifter(temp,1,3,1));
	
end






clear temp temp2 temp3 temp4

% Combine larvae sets

filenames = {'all_spacingvsLength_noheaders', 'cellbyCellStack', 'intrabyCellStack','all_RMSE_sorted_self_embryoEquations','all_RMSE_sorted_self_larvalEquations'};  

for rnd = 1:length(filenames),
	temp = dlmread(SpiffyName('csv',filenames{rnd},'L1'));
	temp2 = dlmread(SpiffyName('csv',filenames{rnd},'1stInstar'));

	padding_temp2 = size(temp,2) - size(temp2,2);

	temp2pad = [temp2,zeros(size(temp2,1),padding_temp2)];

	temp4 = [temp2pad;temp];
	dlmwrite(SpiffyName('csv',filenames{rnd},'Larvae'),temp4)
end 




 clear temp temp2 temp3 temp4

% Combine embryo and larvae sets

filenames = {'all_spacingvsLength_noheaders', 'cellbyCellStack', 'intrabyCellStack'};  

for rnd = 1:length(filenames),
	temp = dlmread(SpiffyName('csv',filenames{rnd},'L1'));
	temp2 = dlmread(SpiffyName('csv',filenames{rnd},'1stInstar'));
    temp3 = dlmread(SpiffyName('csv',filenames{rnd},'yw'));

	padding_temp2 = size(temp,2) - size(temp2,2);
	temp2pad = [temp2,zeros(size(temp2,1),padding_temp2)];
    
    padding_temp3 = size(temp,2) - size(temp3,2);
	temp3pad = [temp3,zeros(size(temp3,1),padding_temp3)];
    
	temp4 = [temp3pad;temp2pad;temp];

	dlmwrite(SpiffyName('csv',filenames{rnd},'Embryos_with_Larvae'),temp4)
end 




 



