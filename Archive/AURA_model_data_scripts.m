%% AURA model-data comparison with NO2

% Load directories and other constants

fig_directory = '/Users/yshiga/Documents/Research/AURA/Figures/';

model_dir = '/Users/yshiga/Documents/Research/AURA/Data/Model_output/';


% Load model output (omi&gmi data)

load([ model_dir 'aura_omi_L3_surf_01_v1.mat'])

lat_lon_grid = meshgrid(lons_01_NH,lats_01_NH);


% Load lat lon for EPA sites and

lat_lon_EPA_site = readtable([ model_dir 'lat_lon_EPA_site.csv']);

ll_epa = table2array(lat_lon_EPA_site(:,2:3));  % store lat lon in double array

% Search by lat lon for index in model for each EPA site 

for i = 1 : length(ll_epa)
    
    ind_ll(i,1) = find(min(abs(ll_epa(i,1)-lats_01))==(abs(ll_epa(i,1)-lats_01))); % find nearest latitude to EPA latitude
    
    ind_ll(i,2) = find(min(abs(ll_epa(i,2)-lons_01))==(abs(ll_epa(i,2)-lons_01))); % find nearest longitude to EPA longitude
    
end

[latlon, ~, inds] = unique([ind_ll(:,1), ind_ll(:,2)], 'rows', 'stable'); % Identify unique latlon location and their index for averaging later

%% Find model output for corresponding lat lon points
% summer mean

data_model_summer = mean(cf_aug,cf_july,cf_june);

for i = 1 : length(ind_ll)
    
    data_EPA(i,:) = data_NO2_01(ind_ll(i,1),ind_ll(i,2),:);

end

% Load EPA NO2 summer average data from 2005-2020

summer_daily_mean_all = readtable([ model_dir 'summer_daily_mean_all.csv']);

% summer_month_mean_all = readtable([model_dir 'summer_month_mean_all.csv']);


% create matrix where rows are site and columns years from 2005 to 2020

list_of_sites = unique(lat_lon_EPA_site.stn); %create list of unique station numbers

all_years = sscanf(sprintf(' %s',summer_daily_mean_all.year{:}),'%f',[1,Inf]);

%% Loop over sites to reorganize EPA data into matrix

% Daily first

% create place holder
summer_daily_mean_EPA = nan(length(list_of_sites),length(2005:2020));

year_list = 2005:2020;

for i = 1 : length(list_of_sites)
    
    site_ind = strcmp(summer_daily_mean_all.stn, lat_lon_EPA_site.stn{i});
    
    years_temp = all_years(site_ind);
    
    data_temp = summer_daily_mean_all.seas_mean(site_ind);
    
    year_index = find(ismember(year_list,years_temp));

    for j = 1 : length(years_temp)
        
       year_index_temp = year_index(j);
        
        summer_daily_mean_EPA(i,year_index_temp) = data_temp(j);
        
    end
    
end
        
summer_month_mean_EPA = nan(length(list_of_sites),length(2005:2020));

year_list = 2005:2020;

for i = 1 : length(list_of_sites)
    
    site_ind = strcmp(summer_month_mean_all.stn, lat_lon_EPA_site.stn{i});
    
    years_temp = all_years(site_ind);
    
    data_temp = summer_month_mean_all.seas_mean(site_ind);
    
    year_index = find(ismember(year_list,years_temp));

    for j = 1 : length(years_temp)
        
       year_index_temp = year_index(j);
        
        summer_month_mean_EPA(i,year_index_temp) = data_temp(j);
        
    end
    
end
  


%% stats
[RHO2,P] = corrcoef(summer_daily_mean_EPA(:),data_EPA(:),'rows','pairwise');

CorrC(1) = RHO2(1,2);

[RHO2,P] = corrcoef(summer_month_mean_EPA(:),data_EPA(:),'rows','pairwise');

CorrC(2) = RHO2(1,2);

Bias(1) = nanmean(summer_daily_mean_EPA(:)-data_EPA(:));
Bias(2) = nanmean(summer_month_mean_EPA(:)-data_EPA(:));


RMSE(1) = sqrt(nanmean((summer_daily_mean_EPA(:)-data_EPA(:)).^2));
RMSE(2) = sqrt(nanmean((summer_month_mean_EPA(:)-data_EPA(:)).^2));


%% figures correlation
figure;

subplot(1,2,1)

% plot(summer_daily_mean_EPA(:),data_EPA(:),'o')
dscatter(summer_daily_mean_EPA(isfinite(summer_daily_mean_EPA(:))),data_EPA(isfinite(summer_daily_mean_EPA(:))))
title({'Model-data comparison for', '2005-2020 summer averages', 'using Daily data'})
str=sprintf('r = %1.4f\nBias = %e\nRMSE = %e',CorrC(1),Bias(1),RMSE(1) );

T = text(min(get(gca, 'xlim'))+2, max(get(gca, 'ylim')), str); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(1,2,2)

% plot(summer_month_mean_EPA(:),data_EPA(:),'o')
dscatter(summer_month_mean_EPA(isfinite(summer_month_mean_EPA(:))),data_EPA(isfinite(summer_month_mean_EPA(:))))
title({'Model-data comparison for', '2005-2020 summer averages', 'using Monthly averages'})
str=sprintf('r = %1.4f\nBias = %e\nRMSE = %e',CorrC(2),Bias(2),RMSE(2) );

T = text(min(get(gca, 'xlim'))+2, max(get(gca, 'ylim')), str); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


set(gcf, 'PaperPosition', [0 0 12 5]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [12 5 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_model_data_summer.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_model_data_summer.png'],'-dpng') % save as png


%% mappping

figure;
load coastlines

subplot(1,2,1)
axesm('MapProjection','miller','MapLatLimit',[15 55],'MapLonLimit',[-130 -60])
framem on
gridm on;
axis off;
mlabel; 
plabel;
t=pcolorm(lats_01,lons_01,mean(data_NO2_01,3));
caxis([0 10*10^15])
hold on
plotm(coastlat, coastlon)
colorbar

subplot(1,2,2)
axesm('MapProjection','miller','MapLatLimit',[15 55],'MapLonLimit',[-130 -60])
framem on
gridm on;
axis off;
mlabel; 
plabel;scatterm(ll_epa(:,1),ll_epa(:,2),nanmean(summer_daily_mean_EPA,2)*5,nanmean(summer_daily_mean_EPA,2),'filled')
hold on
plotm(coastlat, coastlon)
colorbar

set(gcf, 'PaperPosition', [0 0 10 6]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 6 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_model_data_summer_spatial_averages.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_model_data_summer_spatial_averages.png'],'-dpng') % save as png




%% HCHO

% load model output
load('/Users/yshiga/Documents/Research/AURA/Data/Model_output/aura_omi_L3_no2_01_v2.mat')

% load cleaned EPA data
summer_daily_mean_all_HCHO = readtable('/Users/yshiga/Documents/Research/AURA/Data/Model_output/summer_daily_mean_all_HCHO.csv');

summer_month_mean_all_HCHO = readtable('/Users/yshiga/Documents/Research/AURA/Data/Model_output/summer_month_mean_all_HCHO.csv');

lat_lon_EPA_site_HCHO = readtable('/Users/yshiga/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_HCHO.csv');


ll_epa = table2array(lat_lon_EPA_site_HCHO(:,2:3));

for i = 1 : length(ll_epa)
    
    ind_ll(i,1) = find(min(abs(ll_epa(i,1)-lats_01))==(abs(ll_epa(i,1)-lats_01)));
    
    ind_ll(i,2) = find(min(abs(ll_epa(i,2)-lons_01))==(abs(ll_epa(i,2)-lons_01)));
    
end

% model data for specific grid

for i = 1 : length(ind_ll)
    
    data_EPA(i,:) = data_NO2_01(ind_ll(i,1),ind_ll(i,2),:);

end



% create matrix where rows are site and columns years from 2005 to 2020

list_of_sites = unique(lat_lon_EPA_site_HCHO.stn); 

all_years = sscanf(sprintf(' %s',summer_daily_mean_all_HCHO.year{:}),'%f',[1,Inf]);

% Loop over site
% Daily first

% create place holder
summer_daily_mean_EPA = nan(length(list_of_sites),length(2005:2020));

year_list = 2005:2020;

for i = 1 : length(list_of_sites)
    
    site_ind = strcmp(summer_daily_mean_all_HCHO.stn, lat_lon_EPA_site_HCHO.stn{i});
    
    years_temp = all_years(site_ind);
    
    data_temp = summer_daily_mean_all_HCHO.seas_mean(site_ind);
    
    year_index = find(ismember(year_list,years_temp));

    for j = 1 : length(years_temp)
        
       year_index_temp = year_index(j);
        
        summer_daily_mean_EPA(i,year_index_temp) = data_temp(j);
        
    end
    
end
        
summer_month_mean_EPA = nan(length(list_of_sites),length(2005:2020));

year_list = 2005:2020;

for i = 1 : length(list_of_sites)
    
    site_ind = strcmp(summer_month_mean_all_HCHO.stn, lat_lon_EPA_site_HCHO.stn{i});
    
    years_temp = all_years(site_ind);
    
    data_temp = summer_month_mean_all_HCHO.seas_mean(site_ind);
    
    year_index = find(ismember(year_list,years_temp));

    for j = 1 : length(years_temp)
        
       year_index_temp = year_index(j);
        
        summer_month_mean_EPA(i,year_index_temp) = data_temp(j);
        
    end
    
end
 
%% stats
[RHO2,P] = corrcoef(summer_daily_mean_EPA(:),data_EPA(:),'rows','pairwise');

CorrC(1) = RHO2(1,2);

[RHO2,P] = corrcoef(summer_month_mean_EPA(:),data_EPA(:),'rows','pairwise');

CorrC(2) = RHO2(1,2);

Bias(1) = nanmean(summer_daily_mean_EPA(:)-data_EPA(:));
Bias(2) = nanmean(summer_month_mean_EPA(:)-data_EPA(:));


RMSE(1) = sqrt(nanmean((summer_daily_mean_EPA(:)-data_EPA(:)).^2));
RMSE(2) = sqrt(nanmean((summer_month_mean_EPA(:)-data_EPA(:)).^2));


%% figures correlation
figure;

subplot(1,2,1)

% plot(summer_daily_mean_EPA(:),data_EPA(:),'o')
dscatter(summer_daily_mean_EPA(isfinite(summer_daily_mean_EPA(:))),data_EPA(isfinite(summer_daily_mean_EPA(:))))
title({'HCHO model-data comparison for', '2005-2020 summer averages', 'using Daily data'})
str=sprintf('r = %1.4f\nBias = %e\nRMSE = %e',CorrC(1),Bias(1),RMSE(1) );

T = text(min(get(gca, 'xlim'))+2, max(get(gca, 'ylim')), str); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(1,2,2)

% plot(summer_month_mean_EPA(:),data_EPA(:),'o')
dscatter(summer_month_mean_EPA(isfinite(summer_month_mean_EPA(:))),data_EPA(isfinite(summer_month_mean_EPA(:))))
title({'HCHO model-data comparison for', '2005-2020 summer averages', 'using Monthly averages'})
str=sprintf('r = %1.4f\nBias = %e\nRMSE = %e',CorrC(2),Bias(2),RMSE(2) );

T = text(min(get(gca, 'xlim'))+2, max(get(gca, 'ylim')), str); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


set(gcf, 'PaperPosition', [0 0 12 5]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [12 5 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_model_data_summer.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_model_data_summer.png'],'-dpng') % save as png


%% mappping

figure;
load coastlines

subplot(1,2,1)
axesm('MapProjection','miller','MapLatLimit',[15 55],'MapLonLimit',[-130 -60])
framem on
gridm on;
axis off;
mlabel; 
plabel;
% t=pcolorm(lats_01,lons_01,mean(data_NO2_01,3));
% caxis([0 10*10^15])
hold on
plotm(coastlat, coastlon)
colorbar

subplot(1,2,2)
axesm('MapProjection','miller','MapLatLimit',[15 55],'MapLonLimit',[-130 -60])
framem on
gridm on;
axis off;
mlabel; 
plabel;scatterm(ll_epa_filter(:,1),ll_epa_filter(:,2),nanmean(summer_daily_mean_EPA_filter,2)*6,nanmean(summer_daily_mean_EPA_filter,2),'filled')
hold on
plotm(coastlat, coastlon)
colorbar

set(gcf, 'PaperPosition', [0 0 10 6]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 6 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_model_data_summer_spatial_averages_filter.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_model_data_summer_spatial_averages_filter.png'],'-dpng') % save as png
