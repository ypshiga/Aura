clear
clc
close all

%% PLOTTING scripts for AURA
% Loads in-situ EPA data from R code and OMI/GMI model output to produce plots
% of NO2 and HCHO over the US
% 

%% Constants
% load ('C:\Users\sphilip\Desktop\MATLAB\ams3\aura_omi_L3_surf_01_v1.mat');

disp('Starting model-data-comparison code for AURA NO2 and HCHO... Loading constants')

model_dir = '/Users/yshiga/Documents/Research/AURA/Data/Model_output/';

fig_directory = '/Users/yshiga/Documents/Research/AURA/Figures/';

tracerN = {'HCHO','NO2'}; % names of tracers

% tracer = 2; %

ssnN = {'2005-2019', '2005-2009','2010-2014','2015-2019'}; % year ranges

stats = {'r','slope','Mean Bias','RMSE'};  % stat labels

var_names = {'gmi_hcho_surf_01','gmi_hcho_vc_01','gmi_no2_surf_01','gmi_no2_vc_01','omi_hcho_surf_01','omi_hcho_vc_01','omi_no2_surf_01','omi_no2_vc_01'};

%% Load model output (omi&gmi data)

load([ model_dir 'aura_omi_L3_surf_01_v1.mat'])

% restrict to North America

lat1=20;  % lower lat bound

lat2=55; % upper lat bound

lon1=-140;  % lower lon bound

lon2=-40;  % upper lon bound

lats = lats_01_NH(1, find(lats_01_NH<=lat1, 1, 'last' ) : find(lats_01_NH<=lat2, 1, 'last' )); % subset lat

lons = lons_01_NH(1, find(lons_01_NH<=lon1, 1, 'last' ) : find(lons_01_NH<=lon2, 1, 'last' )); % subset lon

lat_lon_grid = meshgrid(lons_01_NH,lats_01_NH); % create grid of lats and lons for in-situ search

%% Load in-situ locations

% Load lat lon for EPA sites

lat_lon_EPA_site = readtable([ model_dir 'lat_lon_EPA_site_hourly_NO2_v5_corrected.csv']);

ll_epa = table2array(lat_lon_EPA_site(:,2:3));  % store lat lon in double array

ll_epa(:,4) = str2double(lat_lon_EPA_site.stn);

% Search by lat lon for index in model for each EPA site

for i = 1 : length(ll_epa)
    
    temp_lat_find=find(min(abs(ll_epa(i,1)-lats_01_NH))==(abs(ll_epa(i,1)-lats_01_NH))); % find nearest latitude to EPA latitude
    
    ind_ll(i,1) = temp_lat_find(1);% pick first if equidistant
    
    temp_lon_find=find(min(abs(ll_epa(i,2)-lons_01_NH))==(abs(ll_epa(i,2)-lons_01_NH))); % find nearest longitude to EPA latitude
    
    ind_ll(i,2) = temp_lon_find(1); % pick first if equidistant
    
end

[latlon, ~, inds] = unique([ind_ll(:,1), ind_ll(:,2)], 'rows', 'stable'); % Identify unique latlon location and their index for averaging later

%% Find model output for corresponding lat lon points

for j = 1 : 8
    
    for i = 1 : length(ind_ll) %loop over indices of data points
        
        MODEL_DATA(i,:,j) = eval([var_names{j} '(ind_ll(i,1),ind_ll(i,2),:)']);
        
    end
    
end
% Load EPA NO2 summer average data from 2005-2020

summer_daily_mean_all = readtable([ model_dir 'summer_hourly_mean_NO2_v5_corrected.csv']);

% summer_month_mean_all = readtable([model_dir 'summer_month_mean_all.csv']);


% create matrix where rows are site and columns years from 2005 to 2020

list_of_sites = lat_lon_EPA_site.stn; %create list of unique station numbers

all_years = sscanf(sprintf(' %s',summer_daily_mean_all.year{:}),'%f',[1,Inf]); % get list of years (need to convert string to number)

%% Loop over sites to reorganize EPA data into matrix

% Daily first

% create place holder for observation data
OBS = nan(length(list_of_sites),length(2005:2019));

year_list = 2005:2019;

for i = 1 : length(list_of_sites) % loop over sites
    
    site_ind = strcmp(summer_daily_mean_all.stn, lat_lon_EPA_site.stn{i}); % identify index of data for given site by station id number
    
    years_temp = all_years(site_ind); % use index to pull out years when data is available
    
%          data_temp = summer_daily_mean_all.seas_mean(site_ind); % also use index to pull out data
    
    data_temp = summer_daily_mean_all.seas_mean_c(site_ind); % also use index to pull out data
    
    
    year_index = find(ismember(year_list,years_temp)); % use find to create index of years when data is available
    
    for j = 1 : length(years_temp) % loop over number (length) of years available
        
        year_index_temp = year_index(j); % for given year use year index to obtain correct year placement
        
        OBS(i,year_index_temp) = data_temp(j); % place data in matrix with row equal to site and col equal to year
        
    end
    
end

%% Clean data : 1. remove 2020; 2. remove sites without "all" years of data; 3. average sites within same model grid

% OBS(:,end)=[]; % remove 2020
%
% miss_data_ind = sum(isnan(OBS),2)>0; % find sites with missing data between 2005:2019
%
% OBS(miss_data_ind,:)=[]; % remove sites with missing data between 2005:2019 (~750 site to ~200 sites)
%
% MODEL_DATA(miss_data_ind,:,:) = []; % remove sites from model data

OBS_LL = ll_epa; % create new lat lon variable

% OBS_LL(miss_data_ind,:) = [];% remove sites from lat lon

% inds_overlap = inds; % create new var for overlap search
%
% inds_overlap(miss_data_ind) = []; % % remove missing data
%
ind_ll_2 = ind_ll;

% ind_ll_2(miss_data_ind,:)=[];

[latlon2, ~, inds2] = unique([ind_ll_2(:,1), ind_ll_2(:,2)], 'rows', 'stable'); % Identify unique latlon location and their index for averaging later

[~, ind_temp] = unique(inds2, 'rows'); % find uniqe values

% duplicate indices
duplicate_ind = setdiff(1:size(inds2, 1), ind_temp); % find indices of duplicate values

% duplicate values
duplicate_value = inds2(duplicate_ind); % find values (index) of overlapping values

% average sites within same grid cell

ind_remove_later = []; % placeholder for indices that will be removed later

for i = 1 : length(duplicate_value) % loop through indices of overlapping values
    
    over_ind = inds2==duplicate_value(i);
    
    obs_temp = OBS(over_ind,:); % grab values for sites in same mode grid cell
    
    ll_temp = OBS_LL(over_ind,:); % grab lat lon
    
    ind_remove_later = [ind_remove_later;find(over_ind)]; % store these indices - will be removed later
    
    OBS_add(i,:) = nanmean(obs_temp); % average obs values
    
    OBS_LL_add(i,:) = nanmean(ll_temp); % average lat lon
    
    temp_model = MODEL_DATA(over_ind,:,:);
    
    MODEL_add(i,:,:) = nanmean(temp_model);
    
end

OBS(ind_remove_later(:),:) = []; % remove overlapping sites from obs

OBS_LL(ind_remove_later(:),:) = [];% remove overlapping sites from lat lon


OBS =[OBS; OBS_add]; % add back in averaged values

OBS_LL =[ OBS_LL ;OBS_LL_add]; % add back in averaged lat lon

MODEL_DATA(ind_remove_later(:),:,:)=[];% remove overlapping sites from MODEL

MODEL_DATA = [MODEL_DATA;MODEL_add]; %  add back in "averaged" MODEL


% double check nan values for MODEL

miss_data_model = sum(sum(isnan(MODEL_DATA),2),3)>0; % find sites with missing data between 2005:2019

MODEL_DATA(miss_data_model,:,:) = [];

OBS_LL(miss_data_model,:) = [];

OBS(miss_data_model,:) = [];

clear OBS_add OBS_LL_add temp_model obs_temp ll_temp MODEL_add

%% HCHO

% load cleaned EPA data
summer_daily_mean_all_HCHO = readtable('/Users/yshiga/Documents/Research/AURA/Data/Model_output/summer_daily_mean_all_HCHO.csv');

lat_lon_EPA_site_HCHO = readtable('/Users/yshiga/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_HCHO.csv');

ll_epa = table2array(lat_lon_EPA_site_HCHO(:,2:3));

ll_epa(:,4) = str2double(lat_lon_EPA_site_HCHO.stn);


clear ind_ll

for i = 1 : length(ll_epa)
    
    ind_ll(i,1) = find(min(abs(ll_epa(i,1)-lats_01_NH))==(abs(ll_epa(i,1)-lats_01_NH))); % find nearest latitude to EPA latitude
    
    ind_ll(i,2) = find(min(abs(ll_epa(i,2)-lons_01_NH))==(abs(ll_epa(i,2)-lons_01_NH))); % find nearest longitude to EPA longitude
    
end

[latlon, ~, inds] = unique([ind_ll(:,1), ind_ll(:,2)], 'rows', 'stable'); % Identify unique latlon location and their index for averaging later

%% Find model output for corresponding lat lon points

for j = 1 : 8
    
    for i = 1 : length(ind_ll) %loop over indices of data points
        
        MODEL_DATA_HCHO(i,:,j) = eval([var_names{j} '(ind_ll(i,1),ind_ll(i,2),:)']);
        
    end
    
end

%% Loop over sites to reorganize EPA data into matrix

% create matrix where rows are site and columns years from 2005 to 2020

list_of_sites = unique(lat_lon_EPA_site_HCHO.stn);

all_years = sscanf(sprintf(' %s',summer_daily_mean_all_HCHO.year{:}),'%f',[1,Inf]);

% Loop over site
% Daily first

% create place holder
OBS_HCHO = nan(length(list_of_sites),length(2005:2020));

year_list = 2005:2020;

for i = 1 : length(list_of_sites)
    
    site_ind = strcmp(summer_daily_mean_all_HCHO.stn, lat_lon_EPA_site_HCHO.stn{i});
    
    years_temp = all_years(site_ind);
    
    data_temp = summer_daily_mean_all_HCHO.seas_mean(site_ind);
    
    year_index = find(ismember(year_list,years_temp));
    
    for j = 1 : length(years_temp)
        
        year_index_temp = year_index(j);
        
        OBS_HCHO(i,year_index_temp) = data_temp(j);
        
    end
    
end

%% Clean data : 1. remove 2020; 2. remove sites without "all" years of data; 3. average sites within same model grid

OBS_HCHO(:,end)=[]; % remove 2020

% miss_data_ind = sum(isnan(OBS_HCHO),2)>0; % find sites with missing data between 2005:2019

% remove outlier

outlier_hi_vals = sum(OBS_HCHO>15,2)>0; % pick arbitrary "high" long term mean value 15 ppb

% miss_data_ind = miss_data_ind | outlier_hi_vals;

miss_data_ind = outlier_hi_vals;

OBS_HCHO(miss_data_ind,:)=[]; % remove sites with missing data between 2005:2019 (~350 site to ~50 sites)

MODEL_DATA_HCHO(miss_data_ind,:,:) = []; % remove sites from model data

OBS_LL_HCHO = ll_epa; % create new lat lon variable

OBS_LL_HCHO(miss_data_ind,:) = [];% remove sites from lat lon

% inds_overlap = inds; % create new var for overlap search
%
% inds_overlap(miss_data_ind) = []; % % remove missing data
%
ind_ll_2 = ind_ll;

ind_ll_2(miss_data_ind,:)=[];

[latlon2, ~, inds2] = unique([ind_ll_2(:,1), ind_ll_2(:,2)], 'rows', 'stable'); % Identify unique latlon location and their index for averaging later

[~, ind_temp] = unique(inds2, 'rows'); % find uniqe values

% duplicate indices
duplicate_ind = setdiff(1:size(inds2, 1), ind_temp); % find indices of duplicate values

% duplicate values
duplicate_value = inds2(duplicate_ind); % find values (index) of overlapping values

% average sites within same grid cell

clear ind_remove_later_HCHO OBS_add_HCHO OBS_LL_add_HCHO

ind_remove_later = []; % placeholder for indices that will be removed later

for i = 1 : length(duplicate_value) % loop through indices of overlapping values
    
    over_ind = inds2==duplicate_value(i);
    
    obs_temp = OBS_HCHO(over_ind,:); % grab values for sites in same mode grid cell
    
    ll_temp = OBS_LL_HCHO(over_ind,:); % grab lat lon
    
    ind_remove_later = [ind_remove_later;find(over_ind)]; % store these indices - will be removed later
    
    OBS_add_HCHO(i,:) = nanmean(obs_temp); % average obs values
    
    OBS_LL_add_HCHO(i,:) = nanmean(ll_temp); % average lat lon
    
    temp_model = MODEL_DATA_HCHO(over_ind,:,:);
    
    MODEL_add(i,:,:) = nanmean(temp_model);
    
end

OBS_HCHO(ind_remove_later(:),:) = []; % remove overlapping sites from obs

OBS_LL_HCHO(ind_remove_later(:),:) = [];% remove overlapping sites from lat lon


OBS_HCHO =[OBS_HCHO ; OBS_add_HCHO]; % add back in averaged values

OBS_LL_HCHO =[ OBS_LL_HCHO ;OBS_LL_add_HCHO]; % add back in averaged lat lon

MODEL_DATA_HCHO(ind_remove_later(:),:,:)=[];% remove overlapping sites from MODEL

MODEL_DATA_HCHO = [MODEL_DATA_HCHO;MODEL_add]; %  add back in "averaged" MODEL


% double check nan values for MODEL

miss_data_model = sum(sum(isnan(MODEL_DATA_HCHO),2),3)>0; % find sites with missing data between 2005:2019

MODEL_DATA_HCHO(miss_data_model,:,:) = [];

OBS_LL_HCHO(miss_data_model,:) = [];

OBS_HCHO(miss_data_model,:) = [];


% double check nan values for MODEL

miss_data_model = sum(sum(isnan(MODEL_DATA_HCHO),2),3)>0; % find sites with missing data between 2005:2019

MODEL_DATA_HCHO(miss_data_model,:,:) = [];

OBS_LL_HCHO(miss_data_model,:) = [];

OBS_HCHO(miss_data_model,:) = [];


clear OBS_add_HCHO OBS_LL_add_HCHO temp_model obs_temp ll_temp

% double check nan values for MODEL

miss_data_model = sum(sum(isnan(MODEL_DATA_HCHO),2),3)>0; % find sites with missing data between 2005:2019

MODEL_DATA_HCHO(miss_data_model,:,:) = [];

OBS_LL_HCHO(miss_data_model,:) = [];

OBS_HCHO(miss_data_model,:) = [];

%%  NO2 city processing 
% NYC and other city lat lon here

city_names ={'New York','Los Angeles','Chicago','DC','Pittsburg','Atlanta','Houston','USA'};

city_lat_lon = [40.7128, 74.0060;34.0522, 118.2437;41.8781, 87.6298;38.9072, 77.0369;40.4406, 79.9959;33.7490, 84.3880;29.7604, 95.3698];


% Search by lat lon for index in model for each EPA site

clear ind_ll_city

for i = 1 : size(city_lat_lon,1)
    
    temp_lat_find=abs(OBS_LL(:,1)-city_lat_lon(i,1))<.25; % find sites within 0.5 latitude
        
    temp_lon_find=abs(OBS_LL(:,2)-(-city_lat_lon(i,2)))<.25; % find sites within 0.5 latitude
    
    ind_ll_city{i} = find(temp_lat_find & temp_lon_find); % pick first if equidistant
    
end

% Stats for timeseries mean and std

for i = 1 : length(ind_ll_city)+1
    
    if i==length(ind_ll_city)+1
        
        city_stats(i,1,:) = nanmean(OBS(:,:));
        
        city_stats(i,2,:) = nanstd(OBS(:,:));
        
        city_mod_stats(i,1,:,:) = squeeze(nanmean(MODEL_DATA(:,:,:)));
        
        city_mod_stats(i,2,:,:) = squeeze(nanstd(MODEL_DATA(:,:,:)));
        
        n_obs(i,:) = sum(~isnan(OBS(:,:)),1);

        
    else
        
        clear vals_temp
        
        vals_temp = ind_ll_city{i}(:);
        
        for j = 1 : length(vals_temp)
            
            setting_type{i,j} = unique(summer_daily_mean_all.Location_Setting(OBS_LL(ind_ll_city{i}(j),4)==str2double(summer_daily_mean_all.stn)));
            if ~isempty(setting_type{i,j})
                
                is_rural(i,j) = strcmp(setting_type{i,j}{:},'RURAL');
                
            else
                
                is_rural(i,j)=0;
                
            end
        end
        ind_ll_city{i} = ind_ll_city{i}(~is_rural(i,1:length(ind_ll_city{i})));
        % Time series stats for all and URBAN only
        n_obs(i,:) = sum(~isnan(OBS(ind_ll_city{i},:)),1);

        city_stats(i,1,:) = nanmean(OBS(ind_ll_city{i},:),1);
        
        city_stats(i,2,:) = nanstd(OBS(ind_ll_city{i},:),0,1);
        
        city_mod_stats(i,1,:,:) = squeeze(nanmean(MODEL_DATA(ind_ll_city{i},:,:),1));
        
        city_mod_stats(i,2,:,:) = squeeze(nanstd(MODEL_DATA(ind_ll_city{i},:,:),1));
        
    end
end

%% PLots for CITY and NA NO2 timeseries column
figure;

mod = 8;

for i = 1 : size(city_stats,1)
    
    subaxis(2,8,i,'SpacingVert',0.05,'SpacingHoriz',0.02)
    
    errorbar(squeeze( city_mod_stats(i,1,:,mod)),squeeze( city_mod_stats(i,2,:,mod)))
    
    xticks(0:5:15)
    
    xlim([0 16])
    
    xticklabels({'2004','2009','2014','2019'})
    
    grid on
    
    title(city_names{i})
    
    if i==1
        
        ylabel('OMI column NO_2 [DU]')
        
    end
    
    subaxis(2,8,i+8,'SpacingVert',0.05,'SpacingHoriz',0.02)
    
    errorbar(squeeze( city_stats(i,1,:)),squeeze( city_stats(i,2,:)),'color',[0.9100    0.4100    0.1700])
    
    if i==1
        
        ylabel('In-situ (corrected) NO_2 [ppb]')
        
    end
    
    xticks(0:5:15)
    
    xlim([0 16])
    
    xticklabels({'2004','2009','2014','2019'})
    %    if i==2
    %         legend ('In-situ','OMI Column')
    %    end
    
    grid on
    
end

set(gcf, 'PaperPosition', [0 0 18 6]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [18 6 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_time_series_city_col.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_time_series_city_col.png'],'-dpng') % save as png

%% PLots for CITY and NA NO2 timeseries column NORMALIZED
figure;
mod = 8;
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    errorbar(normalize(squeeze( city_stats(i,1,:))),normalize(squeeze( city_stats(i,2,:))))
    hold on
    errorbar(normalize(squeeze( city_mod_stats(i,1,:,mod))),normalize(squeeze( city_mod_stats(i,2,:,mod))))
   if i==2
        legend ('In-situ','OMI Surface')
   end
    if i==1
       ylabel('Standard deviations')
    end
       ylim([-4 4])

   xticks(0:5:15)
    xlim([0 16])
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
        grid on

end



set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_time_series_city_norm_col.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_time_series_city_norm_col.png'],'-dpng') % save as png

%% PLots for CITY and NA NO2 timeseries SURFACE GMI
figure;
mod = 3;
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    errorbar(normalize(squeeze( city_stats(i,1,:))),normalize(squeeze( city_stats(i,2,:))))
    hold on
    errorbar(normalize(squeeze( city_mod_stats(i,1,:,mod))),normalize(squeeze( city_mod_stats(i,2,:,mod))))
   if i==2
        legend ('In-situ','GMI Surface')
   end
    if i==1
       ylabel('Standard deviations')
    end
       ylim([-4 4])

   xticks(0:5:15)
    xlim([0 16])
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
        grid on

end



set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_time_series_city_norm_surf_GMI.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_time_series_city_norm_surf_GMI.png'],'-dpng') % save as png
        
%% HCHO city processing 
% Search by lat lon for index in model for each EPA site
clear ind_ll_city city_mod_stats_hcho
clear city_stats_hcho setting_type

for i = 1 : size(city_lat_lon,1)
    
    temp_lat_find=abs(OBS_LL_HCHO(:,1)-city_lat_lon(i,1))<.25; % find sites within 0.5 latitude
        
    temp_lon_find=abs(OBS_LL_HCHO(:,2)-(-city_lat_lon(i,2)))<.25; % find sites within 0.5 latitude
    
    ind_ll_city{i} = find(temp_lat_find & temp_lon_find); % pick first if equidistant
    
end

% Stats for timeseries mean and std

for i = 1 : length(ind_ll_city)+1
    
     if i==length(ind_ll_city)+1
         
        city_stats_hcho(i,1,:) = nanmean(OBS_HCHO(:,:));
        
        city_stats_hcho(i,2,:) = nanstd(OBS_HCHO(:,:));
        
        city_mod_stats_hcho(i,1,:,:) = squeeze(nanmean(MODEL_DATA_HCHO(:,:,:)));
        
        city_mod_stats_hcho(i,2,:,:) = squeeze(nanstd(MODEL_DATA_HCHO(:,:,:)));
        
        n_obs_hcho(i,:) = sum(~isnan(OBS_HCHO(:,:)),1);

        
     else
    
        clear vals_temp
    
        vals_temp = ind_ll_city{i}(:);
    
        for j = 1 : length(vals_temp)
    
            setting_type_hcho{i,j} = unique(summer_daily_mean_all_HCHO.Location_Setting(OBS_LL_HCHO(ind_ll_city{i}(j),4)==str2double(summer_daily_mean_all_HCHO.stn)));
            if ~isempty(setting_type_hcho{i,j})

                is_rural(i,j) = strcmp(setting_type_hcho{i,j}{:},'RURAL');
                
            else
                
                is_rural(i,j)=0;
                
            end
        end
        
    ind_ll_city{i} = ind_ll_city{i}(~is_rural(i,1:length(ind_ll_city{i})));
    
    % Time series stats for all and URBAN only 
    
    n_obs_hcho(i,:) = sum(~isnan(OBS_HCHO(ind_ll_city{i},:)),1);
    
    city_stats_hcho(i,1,:) = nanmean(OBS_HCHO(ind_ll_city{i},:),1);
    
    city_stats_hcho(i,2,:) = nanstd(OBS_HCHO(ind_ll_city{i},:),0,1);

    city_mod_stats_hcho(i,1,:,:) = squeeze(nanmean(MODEL_DATA_HCHO(ind_ll_city{i},:,:),1));
    
    city_mod_stats_hcho(i,2,:,:) = squeeze(nanstd(MODEL_DATA_HCHO(ind_ll_city{i},:,:),0,1));
    
     end
end
%% PLots HCHO SURF
figure;
mod = 5;
for i = 1 : size(city_stats_hcho,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    errorbar(squeeze( city_stats_hcho(i,1,:)),squeeze( city_stats_hcho(i,2,:)))
    hold on
    errorbar(squeeze( city_mod_stats_hcho(i,1,:,mod)),squeeze( city_mod_stats_hcho(i,2,:,mod)))
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI Surface')
   end
    if i==1
       ylabel('HCHO [ppb]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_time_series_city_surf.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_time_series_city_surf.png'],'-dpng') % save as png

%% PLots HCHO SURF NORMALIZED
    
figure;
mod = 6;
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    errorbar(normalize(squeeze( city_stats_hcho(i,1,:))),normalize(squeeze( city_stats_hcho(i,2,:))))
    hold on
    errorbar(normalize(squeeze( city_mod_stats_hcho(i,1,:,mod))),normalize(squeeze( city_mod_stats_hcho(i,2,:,mod))))
   if i==2
        legend ('In-situ','OMI Surface')
   end
    if i==1
       ylabel('Standard deviations')
    end
       ylim([-4 4])

   xticks(0:5:15)
    xlim([0 16])
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
        grid on

end



set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_time_series_city_norm_surf.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_time_series_city_norm_surf.png'],'-dpng') % save as png

%% PLots HCHO Column
figure;
mod = 6;
for i = 1 : size(city_stats_hcho,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    errorbar(squeeze( city_stats_hcho(i,1,:)),squeeze( city_stats_hcho(i,2,:)))
    hold on
    errorbar(squeeze( city_mod_stats_hcho(i,1,:,mod)),squeeze( city_mod_stats_hcho(i,2,:,mod)))
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI Column')
   end
    if i==1
       ylabel('HCHO [ppb]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_time_series_city_col.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_time_series_city_col.png'],'-dpng') % save as png

%% PLots HCHO Column Normalized
     
figure;
mod = 6;
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    errorbar(normalize(squeeze( city_stats_hcho(i,1,:))),normalize(squeeze( city_stats_hcho(i,2,:))))
    hold on
    errorbar(normalize(squeeze( city_mod_stats_hcho(i,1,:,mod))),normalize(squeeze( city_mod_stats_hcho(i,2,:,mod))))
   if i==2
        legend ('In-situ','OMI Column')
   end
    if i==1
       ylabel('Standard deviations')
    end
       ylim([-4 4])

   xticks(0:5:15)
    xlim([0 16])
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
        grid on

end



set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_time_series_city_norm_col.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_time_series_city_norm_col.png'],'-dpng') % save as png

%% PLots Ratio Surf FNRs
figure;
mod = 7;
modhcho = 5;

for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    errorbar(vals_z ,errs__z)
    hold on
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    errs__z_mod = vals_z.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    errorbar(vals_z_mod,errs__z_mod)
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI surface','location','northwest')
    end
    if i==1
        ylabel('Ratio HCHO/NO_2 [unitless]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'RATIO_time_series_city_surf.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'RATIO_time_series_city_surf.png'],'-dpng') % save as png

%% Plots for RATIO Column FNRs
figure;
mod = 8;
modhcho = 6;

for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    errorbar(vals_z ,errs__z)
    hold on
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    errs__z_mod = vals_z.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    errorbar(vals_z_mod,errs__z_mod)
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI column','location','northwest')
    end
    if i==1
        ylabel('Ratio HCHO/NO_2 [unitless]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'RATIO_time_series_city_col.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'RATIO_time_series_city_col.png'],'-dpng') % save as png

%% PLots Ratio Surf NORMALZIED
figure;
mod = 7;
modhcho = 5;

for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    errorbar(normalize(vals_z) ,normalize(errs__z))
    hold on
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    errs__z_mod = vals_z.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod))
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI surface','location','northwest')
    end
    if i==3
        ylabel('Standard deviation Ratio HCHO/NO_2 [unitless]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'RATIO_time_series_city_surf_norm.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'RATIO_time_series_city_surf_norm.png'],'-dpng') % save as png

%% RATIO Column NORMALIZED
figure;
mod = 8;
modhcho = 6;

for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    errorbar(normalize(vals_z) ,normalize(errs__z))
    hold on
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    errs__z_mod = vals_z.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod))
    xticks(0:5:15)
    xlim([0 16])
    if i==2
        legend ('In-situ','OMI column','location','northwest')
    end
    if i==3
        ylabel('Standard deviation Ratio HCHO/NO_2 [unitless]')
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'RATIO_time_series_city_col_norm.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'RATIO_time_series_city_col_norm.png'],'-dpng') % save as png

%% COMBO plot zero mean unit std (for paper)
figure;

mod = 8; % column

modhcho = 6; % column

mod_s = 7;% surf

modhcho_s = 5; % surf

newcolors1 = [   222,235,247
158,202,225
49,130,189]./256;

newcolors2 = [   254,230,206
253,174,107
230,85,13]./256;
         
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    clear vals_hs vals_ns errs_hs errs_ns

    vals_hs = squeeze(city_stats_hcho(i,1,:));
    
    vals_ns = squeeze(city_stats(i,1,:));
    
    errs_hs = squeeze(city_stats_hcho(i,2,:));
    
    errs_ns = squeeze(city_stats(i,2,:));
    
    errorbar(normalize(vals_ns) ,normalize(errs_ns),'-o','color',newcolors1(3,:),'linewidth',1.5,'markerfacecolor',newcolors1(3,:))
    
    hold on
    
    errorbar(normalize(vals_hs) ,normalize(errs_hs),'-o','color',newcolors2(3,:),'linewidth',1.5,'markerfacecolor',newcolors2(3,:))
    
                clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

                vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod_s));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho_s));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod_s));
    
    errorbar(normalize(vals_ns_mod),normalize(errs_ns_mod),'-o','color',newcolors1(2,:))
    
    errorbar(normalize(vals_hs_mod),normalize(errs_hs_mod),'-o','color',newcolors2(2,:))
    
               clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

               vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod));
    
    errorbar(normalize(vals_ns_mod),normalize(errs_ns_mod),'--s','color',newcolors1(1,:),'linewidth',1)
    
    errorbar(normalize(vals_hs_mod),normalize(errs_hs_mod),'--s','color',newcolors2(1,:),'linewidth',1)
       
    % Surface EPA formaldheyde to no2 ratio FNR
    
    clear vals_z errs__z
    
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    
    errorbar(normalize(vals_z) ,normalize(errs__z),'--ok','linewidth',1.5,'markerfacecolor','k')
  
    % Surface OMI formaldheyde to no2 ratio FNR
    clear vals_z_mod errs__z_mod

    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s)./city_mod_stats(i,1,:,mod_s));
    
    errs__z_mod = vals_z_mod.*sqrt((squeeze(city_mod_stats(i,2,:,mod_s))./squeeze(city_mod_stats(i,1,:,mod_s))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho_s))./squeeze(city_mod_stats_hcho(i,1,:,modhcho_s))).^2);
    
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod),'--ok','color',[139,69,19]./256,'linewidth',1.5)
    
    % Column OMI formaldheyde to no2 ratio FNR
    clear vals_z_mod errs__z_mod

    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    
    errs__z_mod = vals_z_mod.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod),'--s','color',[128,128,128]./256,'linewidth',1.5)
   
    ylim([-4 4])
    xticks(0:5:15)
    xlim([0 16])
    text(8,-3.4,{['Obs per year (min,max) = NO_2 (' num2str(min(n_obs(i,:))) ',' num2str(max(n_obs(i,:))) '),  HCHO (' num2str(min(n_obs_hcho(i,:)))  ',' num2str(max(n_obs_hcho(i,:))) ')']},'HorizontalAlignment', 'center')

    if i==8
        leg_pos = legend ('In-situ NO_2','In-situ HCHO','OMI NO_2 Surface','OMI HCHO Surface','OMI NO_2 Column','OMI HCHO Column','AQS surface FNRs','OMI surface FNRs','OMI column FNRs','location','best','NumColumns',5);
       
        leg_pos.BoxFace.ColorType='truecoloralpha';
        
        leg_pos.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        rect = get(leg_pos, 'Position');

        set(leg_pos, 'Position', rect+[-.21 -.175 0 0])

    end
    if i==3
        label_h = ylabel('Standard deviation HCHO, NO_2, and FNRs');
        label_h.Position(2) = -5; % change vertical position of ylabel

    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'Combo_time_series_city_norm_paper.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'Combo_time_series_city_norm_paper.png'],'-dpng') % save as png

%% COMBO plot norm mean (for paper)
figure;

mod = 8; % column

modhcho = 6; % column

mod_s = 7;% surf

modhcho_s = 5; % surf

newcolors1 = [   222,235,247
158,202,225
49,130,189]./256;

newcolors2 = [   254,230,206
253,174,107
230,85,13]./256;
         
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    clear vals_hs vals_ns errs_hs errs_ns

    
    vals_hs = squeeze(city_stats_hcho(i,1,:));
    
    vals_ns = squeeze(city_stats(i,1,:));
    
    errs_hs = squeeze(city_stats_hcho(i,2,:));
    
    errs_ns = squeeze(city_stats(i,2,:));
    
    
    errorbar(vals_ns./nanmean(vals_ns) ,(errs_ns)./nanmean(vals_ns),'-o','color',newcolors1(3,:),'linewidth',1.5,'markerfacecolor',newcolors1(3,:))
    
    hold on
    
    
    errorbar(vals_hs./nanmean(vals_hs) ,(errs_hs)./nanmean(vals_hs),'-o','color',newcolors2(3,:),'linewidth',1.5,'markerfacecolor',newcolors2(3,:))
    
        clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

    vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod_s));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho_s));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod_s));
    
    
    errorbar(vals_ns_mod./nanmean(vals_ns_mod),(errs_ns_mod)./nanmean(vals_ns_mod),'-o','color',newcolors1(2,:))
    
    
    errorbar(vals_hs_mod./nanmean(vals_hs_mod),(errs_hs_mod)./nanmean(vals_hs_mod),'-o','color',newcolors2(2,:))
    
            clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

    vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod));
    
    
    errorbar(vals_ns_mod./nanmean(vals_ns_mod),(errs_ns_mod)./nanmean(vals_ns_mod),'-s','color',newcolors1(1,:),'linewidth',1)
    
    
    errorbar(vals_hs_mod./nanmean(vals_hs_mod),(errs_hs_mod)./nanmean(vals_hs_mod),'-s','color',newcolors2(1,:),'linewidth',1)
       
    % Surface EPA formaldheyde to no2 ratio FNR
    
    clear vals_z errs__z


    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    
    errs__z = vals_z./nanmean(vals_z).*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    
   
    errorbar(vals_z./nanmean(vals_z) ,(errs__z)./nanmean(vals_z),'--ok','linewidth',1.5,'markerfacecolor','k')
  
   
    % Surface OMI formaldheyde to no2 ratio FNR
    clear vals_z_mod errs__z_mod
    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s)./city_mod_stats(i,1,:,mod_s));
    
    errs__z_mod = vals_z_mod./nanmean(vals_z_mod).*sqrt((squeeze(city_mod_stats(i,2,:,mod_s))./squeeze(city_mod_stats(i,1,:,mod_s))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho_s))./squeeze(city_mod_stats_hcho(i,1,:,modhcho_s))).^2);
    
    
    errorbar(vals_z_mod./nanmean(vals_z_mod),(errs__z_mod)./nanmean(vals_z_mod),'--ok','color',[139,69,19]./256,'linewidth',1.5)
   
   
    % Column OMI formaldheyde to no2 ratio FNR
    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    
    errs__z_mod = vals_z_mod./nanmean(vals_z_mod).*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    
   
    errorbar(vals_z_mod./nanmean(vals_z_mod),(errs__z_mod)./nanmean(vals_z_mod),'--s','color',[128,128,128]./256,'linewidth',1.5)
   
    ylim([0 3])
    xticks(0:5:15)
    xlim([0 16])
    
        text(8,2.7,{['Obs per year (min,max) = NO_2 (' num2str(min(n_obs(i,:))) ',' num2str(max(n_obs(i,:))) '),  HCHO (' num2str(min(n_obs_hcho(i,:)))  ',' num2str(max(n_obs_hcho(i,:))) ')']},'HorizontalAlignment', 'center')

    if i==8
        leg_pos = legend ('In-situ NO_2','In-situ HCHO','OMI NO_2 Surface','OMI HCHO Surface','OMI NO_2 Column','OMI HCHO Column','AQS surface FNRs','OMI surface FNRs','OMI column FNRs','location','best','NumColumns',5);
       
        leg_pos.BoxFace.ColorType='truecoloralpha';
        
        leg_pos.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        rect = get(leg_pos, 'Position');

        set(leg_pos, 'Position', rect+[-.21 -.18 0 0])

    end
    if i==3
        label_h = ylabel('Normalized mean HCHO, NO_2, and FNRs (unitless)');
%         label_h.Position(1) = ; % change horizontal position of ylabel
        label_h.Position(2) = -1; % change vertical position of ylabel
    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'Combo_time_series_city_norm_mean_paper.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'Combo_time_series_city_norm_mean_paper.png'],'-dpng') % save as png

%% COMBO plot zero mean unit std (for paper) -- USING GMI
figure;

mod = 4; % column

modhcho = 2; % column

mod_s = 3;% surf

modhcho_s = 1; % surf

newcolors1 = [   222,235,247
158,202,225
49,130,189]./256;

newcolors2 = [   254,230,206
253,174,107
230,85,13]./256;
         
for i = 1 : size(city_stats,1)
    
    
    subaxis(4,2,i,'SpacingVert',0.05,'SpacingHoriz',0.05)
    
    clear vals_hs vals_ns errs_hs errs_ns

    vals_hs = squeeze(city_stats_hcho(i,1,:));
    
    vals_ns = squeeze(city_stats(i,1,:));
    
    errs_hs = squeeze(city_stats_hcho(i,2,:));
    
    errs_ns = squeeze(city_stats(i,2,:));
    
    errorbar(normalize(vals_ns) ,normalize(errs_ns),'-o','color',newcolors1(3,:),'linewidth',1.5,'markerfacecolor',newcolors1(3,:))
    
    hold on
    
    errorbar(normalize(vals_hs) ,normalize(errs_hs),'-o','color',newcolors2(3,:),'linewidth',1.5,'markerfacecolor',newcolors2(3,:))
    
                clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

                vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod_s));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho_s));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod_s));
    
    errorbar(normalize(vals_ns_mod),normalize(errs_ns_mod),'-o','color',newcolors1(2,:))
    
    errorbar(normalize(vals_hs_mod),normalize(errs_hs_mod),'-o','color',newcolors2(2,:))
    
               clear vals_hs_mod vals_ns_mod errs_hs_mod errs_ns_mod

               vals_hs_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho));
    
    vals_ns_mod = squeeze(city_mod_stats(i,1,:,mod));
    
    errs_hs_mod = squeeze(city_mod_stats_hcho(i,2,:,modhcho));
    
    errs_ns_mod = squeeze(city_mod_stats(i,2,:,mod));
    
    errorbar(normalize(vals_ns_mod),normalize(errs_ns_mod),'--s','color',newcolors1(1,:),'linewidth',1)
    
    errorbar(normalize(vals_hs_mod),normalize(errs_hs_mod),'--s','color',newcolors2(1,:),'linewidth',1)
       
    % Surface EPA formaldheyde to no2 ratio FNR
    
    clear vals_z errs__z
    
    vals_z = squeeze(city_stats_hcho(i,1,:)./city_stats(i,1,:));
    
    errs__z = vals_z.*sqrt((squeeze(city_stats(i,2,:))./squeeze(city_stats(i,1,:))).^2 + (squeeze(city_stats_hcho(i,2,:))./squeeze(city_stats_hcho(i,1,:))).^2);
    
    errorbar(normalize(vals_z) ,normalize(errs__z),'--ok','linewidth',1.5,'markerfacecolor','k')
  
    % Surface OMI formaldheyde to no2 ratio FNR
    clear vals_z_mod errs__z_mod

    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho_s)./city_mod_stats(i,1,:,mod_s));
    
    errs__z_mod = vals_z_mod.*sqrt((squeeze(city_mod_stats(i,2,:,mod_s))./squeeze(city_mod_stats(i,1,:,mod_s))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho_s))./squeeze(city_mod_stats_hcho(i,1,:,modhcho_s))).^2);
    
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod),'--ok','color',[139,69,19]./256,'linewidth',1.5)
    
    % Column OMI formaldheyde to no2 ratio FNR
    clear vals_z_mod errs__z_mod

    
    vals_z_mod = squeeze(city_mod_stats_hcho(i,1,:,modhcho)./city_mod_stats(i,1,:,mod));
    
    errs__z_mod = vals_z_mod.*sqrt((squeeze(city_mod_stats(i,2,:,mod))./squeeze(city_mod_stats(i,1,:,mod))).^2 + (squeeze(city_mod_stats_hcho(i,2,:,modhcho))./squeeze(city_mod_stats_hcho(i,1,:,modhcho))).^2);
    
    errorbar(normalize(vals_z_mod),normalize(errs__z_mod),'--s','color',[128,128,128]./256,'linewidth',1.5)
   
    ylim([-4 4])
    xticks(0:5:15)
    xlim([0 16])
    text(8,-3.4,{['Obs per year (min,max) = NO_2 (' num2str(min(n_obs(i,:))) ',' num2str(max(n_obs(i,:))) '),  HCHO (' num2str(min(n_obs_hcho(i,:)))  ',' num2str(max(n_obs_hcho(i,:))) ')']},'HorizontalAlignment', 'center')

    if i==8
        leg_pos = legend ('In-situ NO_2','In-situ HCHO','GMI NO_2 Surface','GMI HCHO Surface','GMI NO_2 Column','GMI HCHO Column','AQS surface FNRs','GMI surface FNRs','GMI column FNRs','location','best','NumColumns',5);
       
        leg_pos.BoxFace.ColorType='truecoloralpha';
        
        leg_pos.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        rect = get(leg_pos, 'Position');

        set(leg_pos, 'Position', rect+[-.21 -.175 0 0])

    end
    if i==3
        label_h = ylabel('Standard deviation HCHO, NO_2, and FNRs');
        label_h.Position(2) = -5; % change vertical position of ylabel

    end
    xticklabels({'2004','2009','2014','2019'})
    title(city_names{i})
    grid on
end


set(gcf, 'PaperPosition', [0 0 10 10]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [10 10 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'Combo_time_series_city_norm_paper_GMI.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'Combo_time_series_city_norm_paper_GMI.png'],'-dpng') % save as png

%% plot barplot
%% HCHO % Search by lat lon for index in model for each EPA site


city_names ={'New York','Los Angeles','Chicago','DC','Pittsburg','Atlanta','Houston','Phillidelphia','USA'};
city_lat_lon = [40.7128, 74.0060;34.0522, 118.2437;41.8781, 87.6298;38.9072, 77.0369;40.4406, 79.9959;33.7490, 84.3880;29.7604, 95.3698;39.9526, 75.1652];

clear ind_ll_city city_mod_stats_hcho
clear city_stats_hcho setting_type is_rural

for i = 1 : size(city_lat_lon,1)
    
    temp_lat_find=abs(OBS_LL_HCHO(:,1)-city_lat_lon(i,1))<.25; % find sites within 0.5 latitude
        
    temp_lon_find=abs(OBS_LL_HCHO(:,2)-(-city_lat_lon(i,2)))<.25; % find sites within 0.5 latitude
    
    ind_ll_city{i} = find(temp_lat_find & temp_lon_find); % pick first if equidistant
    
end

% Stats for timeseries mean and std

for i = 1 : length(ind_ll_city)+1
    
     if i==length(ind_ll_city)+1
         
        city_stats_hcho(i,1,:) = nanmean(OBS_HCHO(:,:));
        
        city_stats_hcho(i,2,:) = nanstd(OBS_HCHO(:,:));
        
        city_mod_stats_hcho(i,1,:,:) = squeeze(nanmean(MODEL_DATA_HCHO(:,:,:)));
        
        city_mod_stats_hcho(i,2,:,:) = squeeze(nanstd(MODEL_DATA_HCHO(:,:,:)));
        
     else
    
        clear vals_temp
    
        vals_temp = ind_ll_city{i}(:);
    
        for j = 1 : length(vals_temp)
    
            setting_type_hcho{i,j} = unique(summer_daily_mean_all_HCHO.Location_Setting(OBS_LL_HCHO(ind_ll_city{i}(j),4)==str2double(summer_daily_mean_all_HCHO.stn)));
            if ~isempty(setting_type_hcho{i,j})

                is_rural(i,j) = strcmp(setting_type_hcho{i,j}{:},'RURAL');
                
            else
                
                is_rural(i,j)=0;
                
            end
        end
    ind_ll_city{i} = ind_ll_city{i}(~is_rural(i,1:length(ind_ll_city{i})));
    % Time series stats for all and URBAN only 
    n_obs_hcho(i,:) = sum(~isnan(OBS_HCHO(ind_ll_city{i},:)));
    
    city_stats_hcho(i,1,:) = nanmean(OBS_HCHO(ind_ll_city{i},:),1);
    
    city_stats_hcho(i,2,:) = nanstd(OBS_HCHO(ind_ll_city{i},:),0,1);

    city_mod_stats_hcho(i,1,:,:) = squeeze(nanmean(MODEL_DATA_HCHO(ind_ll_city{i},:,:),1));
    
    city_mod_stats_hcho(i,2,:,:) = squeeze(nanstd(MODEL_DATA_HCHO(ind_ll_city{i},:,:),0,1));
    
     end
end

data_zhu = [ 3.0387596899224825;NaN;NaN;NaN;NaN;9.279069767441863; 4.775193798449612; 3.74418604651163; NaN];
order = [4,1,2,3];
data_vals_zhu = [9.279069767441863
4.775193798449612
3.74418604651163
 6.077519379844963
 3.0387596899224825];

data_err_upper_zhu=[11.88372093023256
 6.891472868217055
 5.046511627906978
 7.596899224806203
 4.558139534883722];

zhu_std = data_err_upper_zhu-data_vals_zhu;

zhu_std = [zhu_std(4); NaN;NaN;NaN;NaN; zhu_std(1:3); NaN];

data=nanmean(city_stats_hcho(:,1,:),3);
err=sqrt(nansum((city_stats_hcho(:,2,:).^2),3))./squeeze(sum(~isnan(city_stats_hcho(:,1,:)),3));

data_mod=nanmean(city_mod_stats_hcho(:,1,:,modhcho_s),3);
err_mod=sqrt(nansum(squeeze((city_mod_stats_hcho(:,2,:,modhcho_s).^2)),2))./squeeze(sum(~isnan(city_mod_stats_hcho(:,1,:,modhcho_s)),3));
figure
% bar(1:length(data_zhu),data_zhu)   
bar([data_zhu,data,data_mod])                

hold on
% bar(1:length(data),data)                

hold on
% bar(1:length(data_mod),data_mod)                
% 
er = errorbar(1:length(data),data,err,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
% 
er = errorbar((1:length(data_mod))+.22,data_mod,err_mod,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er = errorbar((1:length(data_mod))-.22,data_zhu,zhu_std,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
 xticklabels(city_names)
 set(gca,'XTickLabelRotation',45)

grid on
ylabel('HCHO [ppb]')
title('EPA in-situ HCHO summer (JJA) average 2005-2019') 
hold off
legend('Zhu et al. 2017','AQS','OMI surface')
% ylim([0 14])

set(gcf, 'PaperPosition', [0 0 6 4]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [6 4 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'HCHO_bar.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'HCHO_bar.png'],'-dpng') % save as png


%% Loop through model output types in "var_names"

% check 15 average criteria
p=1;

for k= 4:-1:1
    
    if k == 1
        
        whichYrs = 1:15;
        
    elseif k == 2
        
        whichYrs = 1:5;
        
    elseif k == 3
        
        whichYrs = 6:10;
        
    elseif k == 4
        
        whichYrs = 11:15;
        
    end
    
    if k>1
        
        ind_3yr_HCHO = sum(isnan(OBS_HCHO(:,whichYrs)),2)<=2;
        
        ind_for_all_HCHO(:,p) = ind_3yr_HCHO;
        
        ind_3yr = sum(isnan(OBS(:,whichYrs)),2)<=2;
        
        ind_for_all(:,p) = ind_3yr;
        
        p=p+1;
        
    else
        
        ind_15yr_HCHO = sum(ind_for_all_HCHO,2)==3;
        
        ind_15yr = sum(ind_for_all,2)==3;
        
    end
end

var_names = {'gmi_hcho_surf_01','gmi_hcho_vc_01','gmi_no2_surf_01','gmi_no2_vc_01','omi_hcho_surf_01','omi_hcho_vc_01','omi_no2_surf_01','omi_no2_vc_01'};

c_val_max = [5,5,5,5]; % set max values for color bars (first four are GMI then OMI)

c_val_max = repmat(c_val_max,1,2); % repeat max vals

unit_name_in_situ = 'ppb'; % units for in situ obs

for i = 1 : 8
    
    mod_name = var_names{i}(1:3); % pul model name "var_names"
    
    if strcmp(var_names{i}(5),'h') % identify name of measurement (val_name) and obs type (d_type) along with model unit (unit_name)
        
        val_name = 'HCHO';
        
        ind_all = ind_15yr_HCHO;
        
    else
        
        val_name = 'NO2';
        
        ind_all = ind_15yr;
        
    end
    
    if strcmp(var_names{i}(end-3),'c')
        
        d_type = 'Column';
        
        unit_name = 'DU';
        
    else
        
        d_type = 'Surface';
        
        unit_name = 'ppb';
        
    end
    
    % ASSIGN MODEL OUTPUT
    
    MOD_TEMP = eval(var_names{i}); % eval current model variable
    
    % subset to North America
    
    MDL = MOD_TEMP(find(lats_01_NH<=lat1, 1, 'last' ) : find(lats_01_NH<=lat2, 1, 'last' ), find(lons_01_NH<=lon1, 1, 'last' ) : find(lons_01_NH<=lon2, 1, 'last' ), :);
    
    f1 = figure; % create figure
    
    p=1; % counter for 5yr average index
    
    clear ind_for_all
    
    for k = 4:-1:1
        
        clear obs_plot obs_plot_ll model_plot
        
        if k == 1
            
            whichYrs = 1:15;
            
        elseif k == 2
            
            whichYrs = 1:5;
            
        elseif k == 3
            
            whichYrs = 6:10;
            
        elseif k == 4
            
            whichYrs = 11:15;
            
        end
        
        clear ind_3yr
        
        if strcmp(val_name,'HCHO')
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists)
            
            if k>1
                
                ind_3yr = sum(isnan(OBS_HCHO(:,whichYrs)),2)<=2;
                
                ind_for_all(:,p) = ind_3yr;
                
                p=p+1;
                
            else
             
                ind_3yr = ind_all;
                
            end
            
            temp_obs_plot = OBS_HCHO(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA_HCHO(ind_3yr,whichYrs,i); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL_HCHO(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA_HCHO;
            
            obs_ll_all = OBS_LL_HCHO;
            
            obs_all = OBS_HCHO;
            
            c_in_situ_max = 8;
            
        else
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists
            
            if k>1
                
                ind_3yr = sum(isnan(OBS(:,whichYrs)),2)<=2;
                
                ind_for_all(:,p) = ind_3yr;
                
                p=p+1;
                
            else
                
                ind_3yr = ind_all;
                
            end
            
            temp_obs_plot = OBS(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA(ind_3yr,whichYrs,i); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA;
            
            obs_all = OBS;
            
            obs_ll_all = OBS_LL;
            
            c_in_situ_max = 25;
            
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Plot Model OUTPUT
        %%%%%%%%%%%%%%%%%%%
        
        kk=k;
        
        subplot(5,4,kk); % place in subplot
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold'); % create projection
        
        setm(gca,'MLabelParallel','south'); % label parallels
        
        surfm(lats-0.1/2, lons-0.1/2, nanmean(MDL(:,:,whichYrs),3)); hold on % plot "surface" of model ouput
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );   % set color bar limit
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
        
        h(kk).Label.String = unit_name; % label color bar
        
        colormap('parula'); % set color map
        
        load coast; % load coasts
        
        plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
        
        setm(gca,'FlineWidth',1); % set line width
        
        title(sprintf('%s \n %s \n %s %s', ssnN{k}, val_name,upper(mod_name),d_type)); % title
        
        %%%%%%%%%%%%%%%%%%%
        % Plot Model Output (only at in-situ locations)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+4;
        
        subplot(5,4,kk); % place in subplot
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);% create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');  % create projection
        
        setm(gca,'MLabelParallel','south');  % label parallels
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),15,model_plot,'filled') % plot "scatter" of model ouput points
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' ); % set color bar limit
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
        
        % label color bar
        
        colormap('parula'); % set color map
        
        load coast; % load coasts
        
        plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
        
        setm(gca,'FlineWidth',1); % set line width
        
        title(sprintf('%s %s',upper(mod_name),val_name)); % title
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Plot in-situ data here
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        kk=k+8;
        
        subplot(5,4,kk);
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
        
        setm(gca,'MLabelParallel','south');
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),15,obs_plot,'filled')
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
        
        h(kk).Label.String = unit_name_in_situ;
        
        colormap('parula');
        
        load coast; plotm(lat,long,'k-','LineWIdth',1.5); setm(gca,'FlineWidth',1);
        
        title(sprintf('%s %s','EPA in-situ ',val_name));
        
        %%%%%%%%%%%%%%%%%%%
        % Plot scatter plot with best fit line comparing model/data average
        % values over time (essentially statistics on spatial relationship/correlation)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+12;
        
        subplot(5,4,kk);
        
        [slope, offset, R, sigmaSlope, sigmaInt] = HirshGilroy(obs_plot, model_plot);
        
        MeanBias = sum(model_plot - obs_plot) / length(obs_plot);
        
        StdDev   = std(model_plot - obs_plot, 0, 1);
        
        RMSE = sqrt(sum(((model_plot - obs_plot).^2)) / length(obs_plot));
        
        xmax = max(obs_plot(:));
        
        xmin = min(obs_plot(:));
        
        ymax = max(model_plot(:));
        
        ymin = min(model_plot(:));
        
        if offset < 0
            
            pm = '-';
            
        else
            
            pm = '+';
            
        end
        
        yhat0=slope*[0 xmax*1.2]+offset;                % evaluate over same range from origin
        
        v1=scatter(obs_plot, model_plot,20,['b' 'o'],'filled');  hold on % scatter numbers
        
        v1.MarkerFaceAlpha = .25; % make transparent
        
        v2=plot([0 xmax*1.2],yhat0,['r' '-'],'LineWidth',1.2) ;hold on % best fit line
        
        v3=plot([0 ylim*2], [0 ylim*2], '--k','LineWidth',1.2);% 1:1 line
        
        title(sprintf('%s average \n %s %s (%s)  vs in situ(%s)', ssnN{k},upper(mod_name),d_type,val_name, val_name),'FontWeight', 'bold');
        
        xlim([xmin xmax]); % fixing the x-y axis limits within a range
        
        ylim([ymin ymax]);
        
        h = xlabel('in situ (ppb)','FontWeight', 'bold'); %pos = get(h,'pos');
        
        h = ylabel([upper(mod_name) ' (ppb)'],'FontWeight', 'bold'); %pos = get(h,'pos');
        
        ltext = sprintf('  y = %.2f x %s %.1f  \n  r = %.2f  \n  N = %2.0f \n Mean Bias  = %2.0f (ppb) \n Std. Dev. = %2.0f (ppb) \n RMSE = %2.0f ', slope, pm, abs(offset), R, length(obs_plot),  MeanBias, StdDev, RMSE);
        
        hold on
        
        lgnd=legend(v2,ltext);
        
        hold on
        
        lgnd.BoxFace.ColorType='truecoloralpha';
        
        lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        %%%%%%%%%%%%%%%%%%%
        % Spatial map with scatter plot showing different statistics at
        % each site for "15 year" mean
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+16;
        
        subplot(5,4,kk);
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
        
        setm(gca,'MLabelParallel','south');
        
        % calculate 15-yr R, RMSE etc for each in situ site + add scatterm
        
        model_temp_1 = model_all(ind_all,whichYrs,i); % grab model data for given data/model type
        
        obs_15_yr = obs_all(ind_all,whichYrs);
        
        clear stats_15_yr
        
        for j = 1:size(obs_15_yr,1) % loop through all sites
            
            obs_temp= obs_15_yr(j,:); % get site data
            
            model_temp_2 = model_temp_1(j,:); % get model data
            
            [slope_temp, ~, r_temp, ~, ~] = HirshGilroy(obs_temp, model_temp_2);
            
            stats_15_yr(1,j) =  r_temp; % corr coef
            
            stats_15_yr(2,j) =  slope_temp; % corr coef
            
            stats_15_yr(3,j) = mean(obs_temp-model_temp_2); % mean Bias
            
            stats_15_yr(4,j) = sqrt(mean((obs_temp-model_temp_2).^2)); % rmse
            
        end
                
        scatterm(obs_ll_all(ind_all,1),obs_ll_all(ind_all,2),15,stats_15_yr(k,:),'filled')
        
        hold on
        
        load coast;
        
        plotm(lat,long,'k-','LineWIdth',1.5);
        
        setm(gca,'FlineWidth',1);
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
        
        title(sprintf('Indiviual site statistics: \n %s for 15 yrs ', stats{k}));
        
    end
    
    box on;
    
    set(gcf,'color',[1 1 1]);
    
    set(gcf, 'InvertHardCopy', 'off');
    
    f1.WindowState = 'maximized';
    
    print([fig_directory 'Fig_Scatter_ave_hourly_2d_4_' var_names{i}],'-dpng') % save as png
    
    print(gcf,'-dtiff','-r300', [fig_directory 'Fig_Scatter_ave_hourly_2d_4_' var_names{i}]);
    
    clf
    
    close all
    
end

close all

%% plot for PAPER/report NO2 and HCHO comparisons to EPA data (USING OMI)

% check 15 average criteria
p=1;

for k= 4:-1:1
    
    if k == 1
        
        whichYrs = 1:15;
        
    elseif k == 2
        
        whichYrs = 1:5;
        
    elseif k == 3
        
        whichYrs = 6:10;
        
    elseif k == 4
        
        whichYrs = 11:15;
        
    end
    
    if k>1
        
        ind_3yr_HCHO = sum(isnan(OBS_HCHO(:,whichYrs)),2)<=2;
        
        ind_for_all_HCHO(:,p) = ind_3yr_HCHO;
        
        ind_3yr = sum(isnan(OBS(:,whichYrs)),2)<=2;
        
        ind_for_all(:,p) = ind_3yr;
        
        p=p+1;
        
    else
        
        ind_15yr_HCHO = sum(ind_for_all_HCHO,2)==3;
        
        ind_15yr = sum(ind_for_all,2)==3;
        
    end
end

var_names = {'omi_hcho_surf_01','omi_no2_surf_01'};

c_val_max = [4,1.2,2,2]; % set max values for color bars (first four are GMI then OMI)

c_val_max = repmat(c_val_max,1,2); % repeat max vals


f1 = figure; % create figure

mod_ind =[5,7];

for i = 1 : 2
    
    
    mod_name = var_names{i}(1:3); % pul model name "var_names"
    disp(mod_name)
    disp(var_names{i}(5))
    if strcmp(var_names{i}(5),'h') % identify name of measurement (val_name) and obs type (d_type) along with model unit (unit_name)
        
        val_name = 'HCHO';
        
        ind_all = ind_15yr_HCHO;
        
        unit_name_in_situ = 'HCHO (ppb)'; % units for in situ obs

        
    else
        
        val_name = 'NO_2';
        
        ind_all = ind_15yr;
        
        unit_name_in_situ = 'NO_2 (ppb)'; % units for in situ obs

    end
    
    if strcmp(var_names{i}(end-3),'c')
        
        d_type = 'Column';
        
        unit_name = 'DU';
        
    else
        
        d_type = 'Surface';
        
        unit_name = 'ppb';
        
    end
    disp(val_name)
    
    % ASSIGN MODEL OUTPUT
    
    MOD_TEMP = eval(var_names{i}); % eval current model variable
    
    % subset to North America
    
    MDL = MOD_TEMP(find(lats_01_NH<=lat1, 1, 'last' ) : find(lats_01_NH<=lat2, 1, 'last' ), find(lons_01_NH<=lon1, 1, 'last' ) : find(lons_01_NH<=lon2, 1, 'last' ), :);
        
    p=1; % counter for 5yr average index
    
    clear ind_for_all
    
    for k = 1
        
        clear obs_plot obs_plot_ll model_plot
        
        whichYrs = 1:15;
        
        if strcmp(val_name,'HCHO')
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists)
            
            ind_3yr = ind_all;
            
            temp_obs_plot = OBS_HCHO(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA_HCHO(ind_3yr,whichYrs,mod_ind(i)); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL_HCHO(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA_HCHO;
            
            obs_ll_all = OBS_LL_HCHO;
            
            obs_all = OBS_HCHO;
            
            c_in_situ_max = 8;
            
        else
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists
            
            ind_3yr = ind_all;
            
            temp_obs_plot = OBS(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA(ind_3yr,whichYrs,mod_ind(i)); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA;
            
            obs_all = OBS;
            
            obs_ll_all = OBS_LL;
            
            c_in_situ_max = 25;
            
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Plot Model OUTPUT
        %%%%%%%%%%%%%%%%%%%
        
%         kk=k;
%         
%         subplot(3,2,kk); % place in subplot
%         
%         hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
%         
%         setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold'); % create projection
%         
%         setm(gca,'MLabelParallel','south'); % label parallels
%         
%         surfm(lats-0.1/2, lons-0.1/2, nanmean(MDL(:,:,whichYrs),3)); hold on % plot "surface" of model ouput
%         
%         set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );   % set color bar limit
%         
%         h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
%         
%         h(kk).Label.String = unit_name; % label color bar
%         
%         colormap('parula'); % set color map
%         
%         load coast; % load coasts
%         
%         plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
%         
%         setm(gca,'FlineWidth',1); % set line width
%         
%         title(sprintf('%s \n %s \n %s %s', ssnN{k}, val_name,upper(mod_name),d_type)); % title
%         
        %%%%%%%%%%%%%%%%%%%
        % Plot Model Output (only at in-situ locations)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+(i-1);
        
        h1 = subplot(3,2,kk); % place in subplot
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);% create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');  % create projection
        
        setm(gca,'MLabelParallel','south');  % label parallels
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),45,model_plot,'filled') % plot "scatter" of model ouput points
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 15, 'FontWeight','bold' ); % set color bar limit
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
        
%         if i==1
%             
%             h(kk).Position = h(kk).Position + [-1.44e-1 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         else
%             
%             h(kk).Position = h(kk).Position + [2.5e-2 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         end
        
        h(kk).Label.String = unit_name_in_situ;

%         label color bar
        
        colormap('jet'); % set color map
        
        load coast; % load coasts
        
        plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
        
        setm(gca,'FlineWidth',1); % set line width
        
        originalSize1 = get(gca, 'Position');
        
        title(sprintf('%s %s surface concentration',upper(mod_name),val_name),'FontSize', 15); % title
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Plot in-situ data here
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        kk=k+2+(i-1);
        
        h2 = subplot(3,2,kk);
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
        
        setm(gca,'MLabelParallel','south');
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),45,obs_plot,'filled')
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
        
%         if i==1
%             
%             h(kk).Position = h(kk).Position + [-1.44e-1 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         else
%             
%             h(kk).Position = h(kk).Position + [2.5e-2 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         end
        
        h(kk).Label.String = unit_name_in_situ;
        
        colormap('jet');
        
        load coast; 
        plotm(lat,long,'k-','LineWIdth',1.5); 
        setm(gca,'FlineWidth',1);
        
        title(sprintf('%s %s %s','EPA-AQS in-situ',val_name, 'surface concentration'),'FontSize', 15, 'FontWeight','bold');
        
        %%%%%%%%%%%%%%%%%%%
        % Plot scatter plot with best fit line comparing model/data average
        % values over time (essentially statistics on spatial relationship/correlation)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+4 +(i-1);
        
        h3 = subplot(3,2,kk);
        
        [slope, offset, R, sigmaSlope, sigmaInt] = HirshGilroy(obs_plot, model_plot);
        
        MeanBias = sum(model_plot - obs_plot) / length(obs_plot);
        
        StdDev   = std(model_plot - obs_plot, 0, 1);
        
        RMSE = sqrt(sum(((model_plot - obs_plot).^2)) / length(obs_plot));
        
        xmax = max(obs_plot(:));
        
        xmin = min(obs_plot(:));
        
        ymax = max(model_plot(:));
        
        ymin = min(model_plot(:));
        
        if offset < 0
            
            pm = '-';
            
        else
            
            pm = '+';
            
        end
        
        yhat0=slope*[0 xmax*1.5]+offset;                % evaluate over same range from origin
        
        v1=scatter(obs_plot, model_plot,35,['b' 'o'],'filled');  hold on % scatter numbers
        
        v1.MarkerFaceAlpha = .45; % make transparent
        
        v2=plot([0 xmax*1.5],yhat0,['r' '-'],'LineWidth',1.2) ;hold on % best fit line
        
        v3=plot([0 ylim*2], [0 ylim*2], '--k','LineWidth',1.2);% 1:1 line
        
        title(sprintf('%s versus EPA-AQS %s surface concentration',upper(mod_name), val_name),'FontWeight', 'bold','FontSize', 15);
        
        h = xlabel('EPA-AQS in situ (ppb)','FontWeight', 'bold','FontSize', 15); %pos = get(h,'pos');
        
        if i >0
            
            ylim([0 4]);
            yticklabels ({''})
            yyaxis right
            
            xlim([0 10]); % fixing the x-y axis limits within a range
            
            ylim([0 4]);
            
            h = ylabel([upper(mod_name) ' (ppb)'],'FontWeight', 'bold','color','k','FontSize', 15); %pos = get(h,'pos');
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';
        else
             
            xlim([0 10]); % fixing the x-y axis limits within a range
            
            ylim([0 4]);
            
            h = ylabel([upper(mod_name) ' (ppb)'],'FontWeight', 'bold','color','k','FontSize', 15); %pos = get(h,'pos');
        end
        
        
%         ltext = sprintf('  y = %.2f x %s %.1f  \n  r = %.2f  \n  N = %2.0f \n Mean Bias  = %2.0f (ppb) \n Std. Dev. = %2.0f (ppb) \n RMSE = %2.0f ', slope, pm, abs(offset), R, length(obs_plot),  MeanBias, StdDev, RMSE);
        
        ltext = sprintf('slope = %.2f, r = %.2f, N = %2.0f (%s)', slope, R, length(obs_plot), ssnN{1});

        hold on
        
        lgnd=legend(v2,ltext,'location','southeast');
        
        hold on
        
        lgnd.BoxFace.ColorType='truecoloralpha';
        
        lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        originalSize2 = get(gca, 'Position');

        set(h3, 'Position', [originalSize2(1)+0.7e-2 originalSize2(4)-10e-2 originalSize1(3).*.93 originalSize1(4)*.95]); % Can also use gca instead of h1 if h1 is still active.

        if i==1
            p1 = get(h1, 'Position');
            p2 = get(h2, 'Position');
            p3 = get(h3, 'Position');
            
            set(h1, 'Position', p1 + [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h2, 'Position', p2+ [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h3, 'Position', p3+ [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
        else
            p1 = get(h1, 'Position');
            p2 = get(h2, 'Position');
            p3 = get(h3, 'Position');
            
            set(h1, 'Position', p1 - [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h2, 'Position', p2- [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h3, 'Position', p3- [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
        end
        box on
        
        %%%%%%%%%%%%%%%%%%%
        % Spatial map with scatter plot showing different statistics at
        % each site for "15 year" mean
        %%%%%%%%%%%%%%%%%%%
        
%         kk=k+16;
%         
%         subplot(5,4,kk);
%         
%         hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);
%         
%         setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
%         
%         setm(gca,'MLabelParallel','south');
%         
%         calculate 15-yr R, RMSE etc for each in situ site + add scatterm
%         
%         model_temp_1 = model_all(ind_all,whichYrs,i); % grab model data for given data/model type
%         
%         obs_15_yr = obs_all(ind_all,whichYrs);
%         
%         clear stats_15_yr
%         
%         for j = 1:size(obs_15_yr,1) % loop through all sites
%             
%             obs_temp= obs_15_yr(j,:); % get site data
%             
%             model_temp_2 = model_temp_1(j,:); % get model data
%             
%             [slope_temp, ~, r_temp, ~, ~] = HirshGilroy(obs_temp, model_temp_2);
%             
%             stats_15_yr(1,j) =  r_temp; % corr coef
%             
%             stats_15_yr(2,j) =  slope_temp; % corr coef
%             
%             stats_15_yr(3,j) = mean(obs_temp-model_temp_2); % mean Bias
%             
%             stats_15_yr(4,j) = sqrt(mean((obs_temp-model_temp_2).^2)); % rmse
%             
%         end
%                 
%         scatterm(obs_ll_all(ind_all,1),obs_ll_all(ind_all,2),15,stats_15_yr(k,:),'filled')
%         
%         hold on
%         
%         load coast;
%         
%         plotm(lat,long,'k-','LineWIdth',1.5);
%         
%         setm(gca,'FlineWidth',1);
%         
%         h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
%         
%         title(sprintf('Indiviual site statistics: \n %s for 15 yrs ', stats{k}));
%         
    end
    
   
    
end

box on;

set(gcf,'color',[1 1 1]);

set(gcf, 'InvertHardCopy', 'off');

f1.WindowState = 'maximized';
pause(1)
% 
% set(gcf, 'PaperPosition', [0 0 12 5]); %Position plot at left hand corner with width 10 and height 7.
% 
% set(gcf, 'PaperSize', [12 5 ]); %Set the paper to have width 10 and height 7.


print([fig_directory 'Fig_report_1b_orig'],'-dpng') % save as png

print(gcf,'-dtiff','-r300', [fig_directory 'Fig_report_1b_orig']);

%% plot for PAPER/report OMI NO2 and HCHO comparisons to EPA data (USING GMI)

% check 15 average criteria
p=1;

for k= 4:-1:1
    
    if k == 1
        
        whichYrs = 1:15;
        
    elseif k == 2
        
        whichYrs = 1:5;
        
    elseif k == 3
        
        whichYrs = 6:10;
        
    elseif k == 4
        
        whichYrs = 11:15;
        
    end
    
    if k>1
        
        ind_3yr_HCHO = sum(isnan(OBS_HCHO(:,whichYrs)),2)<=2;
        
        ind_for_all_HCHO(:,p) = ind_3yr_HCHO;
        
        ind_3yr = sum(isnan(OBS(:,whichYrs)),2)<=2;
        
        ind_for_all(:,p) = ind_3yr;
        
        p=p+1;
        
    else
        
        ind_15yr_HCHO = sum(ind_for_all_HCHO,2)==3;
        
        ind_15yr = sum(ind_for_all,2)==3;
        
    end
end

var_names = {'gmi_hcho_surf_01','gmi_no2_surf_01'};

c_val_max = [4,1.2,2,2]; % set max values for color bars (first four are GMI then OMI)

c_val_max = repmat(c_val_max,1,2); % repeat max vals


f1 = figure; % create figure

mod_ind =[1,3];

for i = 1 : 2
    
    
    mod_name = var_names{i}(1:3); % pul model name "var_names"
    disp(mod_name)
    disp(var_names{i}(5))
    if strcmp(var_names{i}(5),'h') % identify name of measurement (val_name) and obs type (d_type) along with model unit (unit_name)
        
        val_name = 'HCHO';
        
        ind_all = ind_15yr_HCHO;
        
        unit_name_in_situ = 'HCHO (ppb)'; % units for in situ obs

        
    else
        
        val_name = 'NO_2';
        
        ind_all = ind_15yr;
        
        unit_name_in_situ = 'NO_2 (ppb)'; % units for in situ obs

    end
    
    if strcmp(var_names{i}(end-3),'c')
        
        d_type = 'Column';
        
        unit_name = 'DU';
        
    else
        
        d_type = 'Surface';
        
        unit_name = 'ppb';
        
    end
    disp(val_name)
    
    % ASSIGN MODEL OUTPUT
    
    MOD_TEMP = eval(var_names{i}); % eval current model variable
    
    % subset to North America
    
    MDL = MOD_TEMP(find(lats_01_NH<=lat1, 1, 'last' ) : find(lats_01_NH<=lat2, 1, 'last' ), find(lons_01_NH<=lon1, 1, 'last' ) : find(lons_01_NH<=lon2, 1, 'last' ), :);
        
    p=1; % counter for 5yr average index
    
    clear ind_for_all
    
    for k = 1
        
        clear obs_plot obs_plot_ll model_plot
        
        whichYrs = 1:15;
        
        if strcmp(val_name,'HCHO')
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists)
            
            ind_3yr = ind_all;
            
            temp_obs_plot = OBS_HCHO(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA_HCHO(ind_3yr,whichYrs,mod_ind(i)); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL_HCHO(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA_HCHO;
            
            obs_ll_all = OBS_LL_HCHO;
            
            obs_all = OBS_HCHO;
            
            c_in_situ_max = 8;
            
        else
            
            % Filter for at least 3 of 5 years
            % Create separate case for 15 years - mean of 5 years (only if
            % mean for each 5 year exists
            
            ind_3yr = ind_all;
            
            temp_obs_plot = OBS(ind_3yr,whichYrs); % obs data to average
            
            temp_mod_plot = MODEL_DATA(ind_3yr,whichYrs,mod_ind(i)); %model data to average
            
            % set model years to NaN where obs data in NaN
            
            temp_mod_plot(isnan(temp_obs_plot))=nan;
            
            obs_plot = nanmean(temp_obs_plot,2); % mean
            
            obs_plot_ll = OBS_LL(ind_3yr,:);
            
            model_plot = nanmean(temp_mod_plot,2);
            
            model_all = MODEL_DATA;
            
            obs_all = OBS;
            
            obs_ll_all = OBS_LL;
            
            c_in_situ_max = 25;
            
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Plot Model OUTPUT
        %%%%%%%%%%%%%%%%%%%
        
%         kk=k;
%         
%         subplot(3,2,kk); % place in subplot
%         
%         hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
%         
%         setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold'); % create projection
%         
%         setm(gca,'MLabelParallel','south'); % label parallels
%         
%         surfm(lats-0.1/2, lons-0.1/2, nanmean(MDL(:,:,whichYrs),3)); hold on % plot "surface" of model ouput
%         
%         set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );   % set color bar limit
%         
%         h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
%         
%         h(kk).Label.String = unit_name; % label color bar
%         
%         colormap('parula'); % set color map
%         
%         load coast; % load coasts
%         
%         plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
%         
%         setm(gca,'FlineWidth',1); % set line width
%         
%         title(sprintf('%s \n %s \n %s %s', ssnN{k}, val_name,upper(mod_name),d_type)); % title
%         
        %%%%%%%%%%%%%%%%%%%
        % Plot Model Output (only at in-situ locations)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+(i-1);
        
        h1 = subplot(3,2,kk); % place in subplot
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);% create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');  % create projection
        
        setm(gca,'MLabelParallel','south');  % label parallels
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),45,model_plot,'filled') % plot "scatter" of model ouput points
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 15, 'FontWeight','bold' ); % set color bar limit
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold'); % add color bar
        
%         if i==1
%             
%             h(kk).Position = h(kk).Position + [-1.44e-1 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         else
%             
%             h(kk).Position = h(kk).Position + [2.5e-2 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         end
        
        h(kk).Label.String = unit_name_in_situ;

%         label color bar
        
        colormap('jet'); % set color map
        
        load coast; % load coasts
        
        plotm(lat,long,'k-','LineWIdth',1.5); % plot coustlines
        
        setm(gca,'FlineWidth',1); % set line width
        
        originalSize1 = get(gca, 'Position');
        
        title(sprintf('%s %s surface concentration',upper(mod_name),val_name),'FontSize', 15); % title
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Plot in-situ data here
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        kk=k+2+(i-1);
        
        h2 = subplot(3,2,kk);
        
        hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]); % create North America map
        
        setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
        
        setm(gca,'MLabelParallel','south');
        
        scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),45,obs_plot,'filled')
        
        set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );
        
        h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
        
%         if i==1
%             
%             h(kk).Position = h(kk).Position + [-1.44e-1 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         else
%             
%             h(kk).Position = h(kk).Position + [2.5e-2 -1e-2 -h(kk).Position(3)*.5 0];
%             
%         end
        
        h(kk).Label.String = unit_name_in_situ;
        
        colormap('jet');
        
        load coast; 
        plotm(lat,long,'k-','LineWIdth',1.5); 
        setm(gca,'FlineWidth',1);
        
        title(sprintf('%s %s %s','EPA-AQS in-situ',val_name, 'surface concentration'),'FontSize', 15, 'FontWeight','bold');
        
        %%%%%%%%%%%%%%%%%%%
        % Plot scatter plot with best fit line comparing model/data average
        % values over time (essentially statistics on spatial relationship/correlation)
        %%%%%%%%%%%%%%%%%%%
        
        kk=k+4 +(i-1);
        
        h3 = subplot(3,2,kk);
        
        [slope, offset, R, sigmaSlope, sigmaInt] = HirshGilroy(obs_plot, model_plot);
        
        MeanBias = sum(model_plot - obs_plot) / length(obs_plot);
        
        StdDev   = std(model_plot - obs_plot, 0, 1);
        
        RMSE = sqrt(sum(((model_plot - obs_plot).^2)) / length(obs_plot));
        
        xmax = max(obs_plot(:));
        
        xmin = min(obs_plot(:));
        
        ymax = max(model_plot(:));
        
        ymin = min(model_plot(:));
        
        if offset < 0
            
            pm = '-';
            
        else
            
            pm = '+';
            
        end
        
        yhat0=slope*[0 xmax*1.5]+offset;                % evaluate over same range from origin
        
        v1=scatter(obs_plot, model_plot,35,['b' 'o'],'filled');  hold on % scatter numbers
        
        v1.MarkerFaceAlpha = .45; % make transparent
        
        v2=plot([0 xmax*1.5],yhat0,['r' '-'],'LineWidth',1.2) ;hold on % best fit line
        
        v3=plot([0 ylim*2], [0 ylim*2], '--k','LineWidth',1.2);% 1:1 line
        
        title(sprintf('%s versus EPA-AQS %s surface concentration',upper(mod_name), val_name),'FontWeight', 'bold','FontSize', 15);
        
        h = xlabel('EPA-AQS in situ (ppb)','FontWeight', 'bold','FontSize', 15); %pos = get(h,'pos');
        
        if i >0
            
            ylim([0 4]);
            yticklabels ({''})
            yyaxis right
            
            xlim([0 10]); % fixing the x-y axis limits within a range
            
            ylim([0 4]);
            
            h = ylabel([upper(mod_name) ' (ppb)'],'FontWeight', 'bold','color','k','FontSize', 15); %pos = get(h,'pos');
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';
        else
             
            xlim([0 10]); % fixing the x-y axis limits within a range
            
            ylim([0 4]);
            
            h = ylabel([upper(mod_name) ' (ppb)'],'FontWeight', 'bold','color','k','FontSize', 15); %pos = get(h,'pos');
        end
        
        
%         ltext = sprintf('  y = %.2f x %s %.1f  \n  r = %.2f  \n  N = %2.0f \n Mean Bias  = %2.0f (ppb) \n Std. Dev. = %2.0f (ppb) \n RMSE = %2.0f ', slope, pm, abs(offset), R, length(obs_plot),  MeanBias, StdDev, RMSE);
        
        ltext = sprintf('slope = %.2f, r = %.2f, N = %2.0f (%s)', slope, R, length(obs_plot), ssnN{1});

        hold on
        
        lgnd=legend(v2,ltext,'location','southeast');
        
        hold on
        
        lgnd.BoxFace.ColorType='truecoloralpha';
        
        lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
        
        originalSize2 = get(gca, 'Position');

        set(h3, 'Position', [originalSize2(1)+0.7e-2 originalSize2(4)-10e-2 originalSize1(3).*.93 originalSize1(4)*.95]); % Can also use gca instead of h1 if h1 is still active.

        if i==1
            p1 = get(h1, 'Position');
            p2 = get(h2, 'Position');
            p3 = get(h3, 'Position');
            
            set(h1, 'Position', p1 + [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h2, 'Position', p2+ [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h3, 'Position', p3+ [.105 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
        else
            p1 = get(h1, 'Position');
            p2 = get(h2, 'Position');
            p3 = get(h3, 'Position');
            
            set(h1, 'Position', p1 - [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h2, 'Position', p2- [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
            set(h3, 'Position', p3- [.005 0 0 0]); % Can also use gca instead of h1 if h1 is still active.
        end
        box on
        
        %%%%%%%%%%%%%%%%%%%
        % Spatial map with scatter plot showing different statistics at
        % each site for "15 year" mean
        %%%%%%%%%%%%%%%%%%%
        
%         kk=k+16;
%         
%         subplot(5,4,kk);
%         
%         hh(kk) = worldmap([lats(1) lats(end)],[lons(1) lons(end)]);
%         
%         setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold');
%         
%         setm(gca,'MLabelParallel','south');
%         
%         calculate 15-yr R, RMSE etc for each in situ site + add scatterm
%         
%         model_temp_1 = model_all(ind_all,whichYrs,i); % grab model data for given data/model type
%         
%         obs_15_yr = obs_all(ind_all,whichYrs);
%         
%         clear stats_15_yr
%         
%         for j = 1:size(obs_15_yr,1) % loop through all sites
%             
%             obs_temp= obs_15_yr(j,:); % get site data
%             
%             model_temp_2 = model_temp_1(j,:); % get model data
%             
%             [slope_temp, ~, r_temp, ~, ~] = HirshGilroy(obs_temp, model_temp_2);
%             
%             stats_15_yr(1,j) =  r_temp; % corr coef
%             
%             stats_15_yr(2,j) =  slope_temp; % corr coef
%             
%             stats_15_yr(3,j) = mean(obs_temp-model_temp_2); % mean Bias
%             
%             stats_15_yr(4,j) = sqrt(mean((obs_temp-model_temp_2).^2)); % rmse
%             
%         end
%                 
%         scatterm(obs_ll_all(ind_all,1),obs_ll_all(ind_all,2),15,stats_15_yr(k,:),'filled')
%         
%         hold on
%         
%         load coast;
%         
%         plotm(lat,long,'k-','LineWIdth',1.5);
%         
%         setm(gca,'FlineWidth',1);
%         
%         h(kk) = colorbar('vert','FontSize', 15, 'FontWeight','bold');
%         
%         title(sprintf('Indiviual site statistics: \n %s for 15 yrs ', stats{k}));
%         
    end
    
   
    
end

box on;

set(gcf,'color',[1 1 1]);

set(gcf, 'InvertHardCopy', 'off');

f1.WindowState = 'maximized';
pause(1)
% 
% set(gcf, 'PaperPosition', [0 0 12 5]); %Position plot at left hand corner with width 10 and height 7.
% 
% set(gcf, 'PaperSize', [12 5 ]); %Set the paper to have width 10 and height 7.


print([fig_directory 'Fig_report_1b_GMI'],'-dpng') % save as png

print(gcf,'-dtiff','-r300', [fig_directory 'Fig_report_1b_GMI']);


    
%% histograms

figure;
histogram(nanmean(OBS(:,1:5),2),40);
hold on
histogram(nanmean(MODEL_DATA(:,1:5,3),2),20);
histogram(nanmean(MODEL_DATA(:,1:5,7),2),20);
xline(mean(nanmean(OBS(:,1:5))),'--','linewidth',2)
xline(mean(nanmean(MODEL_DATA(:,1:5,3),2)),'--','linewidth',2)
xline(mean(nanmean(MODEL_DATA(:,1:5,7),2)),'--','linewidth',2)
title('2005-2009 average JJA [ppb]');
legend('In-situ','GMI Surf','OMI surf')

set(gcf, 'PaperPosition', [0 0 12 5]); %Position plot at left hand corner with width 10 and height 7.

set(gcf, 'PaperSize', [12 5 ]); %Set the paper to have width 10 and height 7.

print([ fig_directory 'NO2_2005_2009_hist_summer_mod_data_d_4.pdf'],'-dpdf') % save as pdf

print([ fig_directory 'NO2_2005_2009_hist_summer_mod_data_d_4.png'],'-dpng') % save as png

figure;
subplot(2,1,1);
histogram(nanmean(OBS_HCHO(:,1:5),2));
title('2005-2009 average JJA EPA in-situ HCHO [ppb]');
subplot(2,1,2);
histogram(nanmean(MODEL_DATA_HCHO(:,1:5,1),2));
title('2005-2009 average JJA GMI Surface HCHO (Sampled at in-situ points) [ppb]');



% %% plots for 2005 - 2009 for GMI surface
%
% figure;
% whichYrs=1:5;
% subplot(2,1,1)
% worldmap([lats(1) lats(end)],[lons(1) lons(end)]);
%
% setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold'); setm(gca,'MLabelParallel','south');
%
% surfm(lats-0.1/2, lons-0.1/2, nanmean(MDL(:,:,whichYrs),3)); hold on
%
% % set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );   % for NO2 surface (ppb)
%
% colorbar('vert','FontSize', 15, 'FontWeight','bold');
%
% colormap(viridis);
%
% load coast;
%
% plotm(lat,long,'k-','LineWIdth',1.5); setm(gca,'FlineWidth',1);
%
% title(sprintf('NO2 GMI Surf'));
%
% subplot(2,1,2)
%
% worldmap([lats(1) lats(end)],[lons(1) lons(end)]);
%
% setm(gca,'MapProjection','miller','parallellabel','off','meridianlabel','off','FontSize', 6, 'FontWeight','bold'); setm(gca,'MLabelParallel','south');
%
% % put scatterm here!
%
% scatterm(obs_plot_ll(:,1),obs_plot_ll(:,2),50,nanmean(obs_plot(:,whichYrs),2),'filled')
%
% % set(hh(kk),'clim', [0 c_val_max(i)], 'FontSize', 10, 'FontWeight','bold' );
%
%
% load coast;
%
% plotm(lat,long,'k-','LineWIdth',1.5); setm(gca,'FlineWidth',1);
% colorbar('vert','FontSize', 15, 'FontWeight','bold');
%
%
% title(sprintf('EPA in-situ NO2'));