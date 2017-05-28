STD = 1.75;

larva = dlmread('1stInstar_all_spacingvsLength_noheaders.csv');
larva2 = dlmread('L1_all_spacingvsLength_noheaders.csv');
larvanet = [larva(:,1:8);larva2(:,1:8)];

embryos = dlmread('yw_all_spacingvsLength_noheaders.csv');


for number = 0:3,
    % first, sort out what data to use
    dentcol = 5+number;
    celllengthcol = 3;
    
    embryodata = [embryos(embryos(:,dentcol)~=0,celllengthcol), embryos(embryos(:,dentcol)~=0,dentcol)];
    larvadata = [larvanet(larvanet(:,dentcol)~=0,celllengthcol), larvanet(larvanet(:,dentcol)~=0,dentcol)];
    
    egroup = ones(size(embryodata,1),1);
    lgroup = 2*ones(size(larvadata,1),1);

    
    
    data =[embryodata;larvadata];
    group = [egroup; lgroup];
    
    types = {'embryos','larvae'};
    
    % scatterplot (needs statistics toolbox)
    gscatter(data(:,1),data(:,2),group,'kr','..',15,15);
    hold on

    for groupnumber=1:2
        %# indices of points in this group
        index = ( group == groupnumber );
        
        %# subtract mean
        Mu = mean(data(index,:) );
        X0 = bsxfun(@minus, data(index,:), Mu);
        % bsxfun applies element-wise operation (in this case minus) to the
        % specified arrays
        
        %    STD = 1.75;                     %# 2 standard deviations
        conf = 2*normcdf(STD)-1;     %# covers around 95% of population
        % printout conf -> approximate % of population covered
        scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
        
        Cov = cov(X0) * scale;
        [V, D] = eig(Cov);
        
        %     [D, order] = sort(diag(D), 'descend');
        %     D = diag(D);
        %     V = V(:, order);
        %
        t = linspace(0,2*pi,100);
        circle = [cos(t) ; sin(t)];        %# unit circle
        VV = V*sqrt(D);               %# scale eigenvectors
        ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space
        
        %# plot cov and major/minor axes
        plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');
        
        % save into a csv file
        dlmwrite(SpiffyName('csv',num2str(number + 2),'ellipsecoordinates',types{groupnumber}),ellipsetodraw);


    end

end 

% 
% 
% 
%     
% l2 = [larvanet(larvanet(:,5)~=0,3), larvanet(larvanet(:,5)~=0,5)];
% l3 = [larvanet(larvanet(:,6)~=0,3), larvanet(larvanet(:,6)~=0,6)];
% l4 = [larvanet(larvanet(:,7)~=0,3), larvanet(larvanet(:,7)~=0,7)];
% 
% e2 = [embryos(embryos(:,5)~=0,3), embryos(embryos(:,5)~=0,5)];
% e3 = [embryos(embryos(:,6)~=0,3), embryos(embryos(:,6)~=0,6)];
% e4 = [embryos(embryos(:,7)~=0,3), embryos(embryos(:,7)~=0,7)];
% 
% lgroup = 2*ones(size(l2,1),1);
% egroup = ones(size(e2,1),1);
% 
% data =[e2;l2];
% group = [egroup; lgroup];
% 
% % scatterplot (needs statistics toolbox)
% gscatter(data(:,1),data(:,2),group,'kr','..',15,15);
% hold on
% 
% 
% 
% for groupnumber=1:2
%     %# indices of points in this group
%     index = ( group == groupnumber );
% 
%     %# subtract mean
%     Mu = mean(data(index,:) );
%     X0 = bsxfun(@minus, data(index,:), Mu);
%     % bsxfun applies element-wise operation (in this case minus) to the
%     % specified arrays
%     
% %    STD = 1.75;                     %# 2 standard deviations
%     conf = 2*normcdf(STD)-1     %# covers around 95% of population
%         % printout conf -> approximate % of population covered
%     scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
% 
%     Cov = cov(X0) * scale;
%     [V, D] = eig(Cov)
% 
% %     [D, order] = sort(diag(D), 'descend');
% %     D = diag(D);
% %     V = V(:, order);
% % 
%     t = linspace(0,2*pi,100);
%     circle = [cos(t) ; sin(t)];        %# unit circle
%     VV = V*sqrt(D);               %# scale eigenvectors
%     ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space
% 
%     %# plot cov and major/minor axes
%     plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');
%     %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
%     %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
% end
% 
% 
% 
% % STD = 2;                     %# 2 standard deviations
% % conf = 2*normcdf(STD)-1;     %# covers around 95% of population
% % scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
% % 
% % Cov = cov(X0) * scale;
% % [V, D] = eig(Cov);
% 
% 
% %   for groupnumber=1:2
% %     %# indices of points in this group
% %     index = ( group2 == groupnumber );
% % 
% %     %# subtract mean
% %     % Shift the data to be centered around zero (much easier to calculate
% %     % stuff from
% %     Mu = mean(data2(index,:) );
% %     dataround0 = bsxfun(@minus, data2(index,:), Mu);
% %     % bsxfun applies element-wise operation (in this case minus) to the
% %     % specified arrays
% %     
% %     %# eigen decomposition [sorted by eigen values]
% %     % This gives you the scaling factors needed
% %     [V, D] = eig( dataround0'*dataround0 ./ (sum(index)-1) );     %#' cov(X0)
% %     [D, order] = sort(diag(D), 'descend');
% %     D = diag(D);
% %     V = V(:, order);
% % 
% %     
% %     t = linspace(0,2*pi,100);
% %     circle = [cos(t) ; sin(t)];        %# unit circle
% %     VV = V*sqrt(D);               %# scale eigenvectors
% %     ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space
% % 
% %     %# plot cov and major/minor axes
% %     plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');
% % 
% %     %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
% %     %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
% % end
% 














% figure

% larva = dlmread('1stInstar_all_spacingvsLength_noheaders.csv');
% larva2 = dlmread('L1_all_spacingvsLength_noheaders.csv');
% larvanet = [larva(:,1:8);larva2(:,1:8)];

% l2 = [larvanet(larvanet(:,5)~=0,3), larvanet(larvanet(:,5)~=0,5)];
% l3 = [larvanet(larvanet(:,6)~=0,3), larvanet(larvanet(:,6)~=0,6)];
% l4 = [larvanet(larvanet(:,7)~=0,3), larvanet(larvanet(:,7)~=0,7)];

% lgroup2 = 2*ones(size(l2,1),1);


% embryos = dlmread('yw_all_spacingvsLength_noheaders.csv');

% e2 = [embryos(embryos(:,5)~=0,3), embryos(embryos(:,5)~=0,5)];
% e3 = [embryos(embryos(:,6)~=0,3), embryos(embryos(:,6)~=0,6)];
% e4 = [embryos(embryos(:,7)~=0,3), embryos(embryos(:,7)~=0,7)];

% egroup2 = ones(size(e2,1),1);

% data2 =[e2;l2];
% group2 = [egroup2; lgroup2];

% % scatterplot (needs statistics toolbox)
% gscatter(data2(:,1),data2(:,2),group2,'kr','..',15,15);
% hold on



% for groupnumber=1:2
%     %# indices of points in this group
%     index = ( group2 == groupnumber );

%     %# subtract mean
%     Mu = mean(data2(index,:) );
%     X0 = bsxfun(@minus, data2(index,:), Mu);
%     % bsxfun applies element-wise operation (in this case minus) to the
%     % specified arrays
    
%     STD = 1.75;                     %# 2 standard deviations
%     conf = 2*normcdf(STD)-1     %# covers around 95% of population
%         % printout conf -> approximate % of population covered
%     scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

%     Cov = cov(X0) * scale;
%     [V, D] = eig(Cov)

% %     [D, order] = sort(diag(D), 'descend');
% %     D = diag(D);
% %     V = V(:, order);
% % 
%     t = linspace(0,2*pi,100);
%     circle = [cos(t) ; sin(t)];        %# unit circle
%     VV = V*sqrt(D);               %# scale eigenvectors
%     ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space

%     %# plot cov and major/minor axes
%     plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');
%     %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
%     %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
% end



% % STD = 2;                     %# 2 standard deviations
% % conf = 2*normcdf(STD)-1;     %# covers around 95% of population
% % scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
% % 
% % Cov = cov(X0) * scale;
% % [V, D] = eig(Cov);


% %   for groupnumber=1:2
% %     %# indices of points in this group
% %     index = ( group2 == groupnumber );
% % 
% %     %# subtract mean
% %     % Shift the data to be centered around zero (much easier to calculate
% %     % stuff from
% %     Mu = mean(data2(index,:) );
% %     dataround0 = bsxfun(@minus, data2(index,:), Mu);
% %     % bsxfun applies element-wise operation (in this case minus) to the
% %     % specified arrays
% %     
% %     %# eigen decomposition [sorted by eigen values]
% %     % This gives you the scaling factors needed
% %     [V, D] = eig( dataround0'*dataround0 ./ (sum(index)-1) );     %#' cov(X0)
% %     [D, order] = sort(diag(D), 'descend');
% %     D = diag(D);
% %     V = V(:, order);
% % 
% %     
% %     t = linspace(0,2*pi,100);
% %     circle = [cos(t) ; sin(t)];        %# unit circle
% %     VV = V*sqrt(D);               %# scale eigenvectors
% %     ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space
% % 
% %     %# plot cov and major/minor axes
% %     plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');
% % 
% %     %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
% %     %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
% % end

















% %%% NO CONTROL OF ANYTHING

% larva = dlmread('1stInstar_all_spacingvsLength_noheaders.csv');
% larva2 = dlmread('L1_all_spacingvsLength_noheaders.csv');
% larvanet = [larva(:,1:8);larva2(:,1:8)];

% l2 = [larvanet(larvanet(:,5)~=0,3), larvanet(larvanet(:,5)~=0,5)];
% l3 = [larvanet(larvanet(:,6)~=0,3), larvanet(larvanet(:,6)~=0,6)];
% l4 = [larvanet(larvanet(:,7)~=0,3), larvanet(larvanet(:,7)~=0,7)];

% lgroup2 = 2*ones(size(l2,1),1);


% embryos = dlmread('yw_all_spacingvsLength_noheaders.csv');

% e2 = [embryos(embryos(:,5)~=0,3), embryos(embryos(:,5)~=0,5)];
% e3 = [embryos(embryos(:,6)~=0,3), embryos(embryos(:,6)~=0,6)];
% e4 = [embryos(embryos(:,7)~=0,3), embryos(embryos(:,7)~=0,7)];

% egroup2 = ones(size(e2,1),1);

% data2 =[e2;l2];
% group2 = [egroup2; lgroup2];

% % scatterplot (needs statistics toolbox)
% gscatter(data2(:,1),data2(:,2),group2,'kr','..',15,15);
% hold on


%   for groupnumber=1:2
%     %# indices of points in this group
%     index = ( group2 == groupnumber );

%     %# subtract mean
%     % Shift the data to be centered around zero (much easier to calculate
%     % stuff from
%     Mu = mean(data2(index,:) );
%     dataround0 = bsxfun(@minus, data2(index,:), Mu);
%     % bsxfun applies element-wise operation (in this case minus) to the
%     % specified arrays
    
%     %# eigen decomposition [sorted by eigen values]
%     % This gives you the scaling factors needed
%     [V, D] = eig( dataround0'*dataround0 ./ (sum(index)-1) );     %#' cov(X0)
%     [D, order] = sort(diag(D), 'descend');
%     D = diag(D);
%     V = V(:, order);

%     t = linspace(0,2*pi,100);
%     circle = [cos(t) ; sin(t)];        %# unit circle
%     VV = V*sqrt(D);               %# scale eigenvectors
%     ellipsetodraw = bsxfun(@plus, VV*circle, Mu'); %#' project circle back to orig space

%     %# plot cov and major/minor axes
%     plot(ellipsetodraw(1,:), ellipsetodraw(2,:), 'LineWidth',4,'Color','b');

%     %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
%     %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
% end


