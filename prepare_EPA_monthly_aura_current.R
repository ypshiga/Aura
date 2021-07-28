#!/usr/bin/env Rscript

# This script uses a template script from Allan Just to average daily NO2 measurements from EPA data to be used for model data copmarison for AURA project (PI Sajeeev Philip)

# Downloaded daily EPA data from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Daily
# I simply manually clicked on links for 2005-2020 (accessed on 1-25-21)
# For NO2, Ozone and HAPS (for HCHO)


# Libraries ---------------------------------------------------------------

suppressMessages(library(data.table))
suppressMessages(library(fst))
library(sf)
library(nngeo)
library(ggplot2)
library(dplyr)
library(lubridate)
library(zoo)
library(seas)
library(R.matlab)
library(sp)
library(rgeos)

fig_dir = '/Users/yshiga/Documents/Research/AURA/Figures/'

data_dir = "~/Documents/Research/AURA/Data/Daily/NO2"

# Open Measurement Daily Data ---------------------------------------------

# Daily
storewd <- getwd()
# change working directory to where the data live
filepath <- data_dir # "data/aqs_raw/epa_daily"
setwd(filepath)
#list of data in folder
zipfiles <- list.files(pattern = "daily.*\\.zip")
# Load all daily data in folder in list
dailydt <- rbindlist(lapply(zipfiles, function(x) fread(cmd=paste0("unzip -cq ", x))))

# restore the working directory
setwd(storewd)
#clean up
rm(zipfiles, filepath)

# how big is the data
paste("The loaded daily PM measurement data table has dimensions", paste(dim(dailydt), collapse = ",")) 
# pryr::object_size(dailydt) #  For NO2 = 425 MB

# clean the day
dailydt[, day := as.Date(`Date Local`)]

# 
# # look at 2018 all sites
# dailydt <- dailydt[!(`State Name` %in% states_non_contig) & day > as.Date("2017-12-31") & day < as.Date("2019-01-01"), ]
# paste("PM Daily measurements table subset to 2018 has dimensions", paste(dim(dailydt), collapse = ",")) # 1243560,30

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Variable construction. 
# Create unique station identifier
# Pull Land Use and Location Setting info from aqs_sites.cvs (from EPA) ####
# Downloaded aqs_sites.zip from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Meta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Daily
# create a unique station identifier with uniform length
dailydt[, stn := paste0(
  stringr::str_pad(`State Code`, 2, "left", pad = "0"), # two digit for state
  stringr::str_pad(`County Code`, 3, "left", pad = "0"),# three digits for county
  stringr::str_pad(`Site Num`, 4, "left", pad = "0"))]  # four digits site

# rename the main parameter
# setnames(dailydt, "Arithmetic Mean", "NO2")

# Open site description file
sites <- fread("~/Documents/Research/AURA/Data/Daily/aqs_sites.csv")

# In site table, create a unique station identifier with uniform length
sites[, stn := paste0(
  stringr::str_pad(`State Code`, 2, "left", pad = "0"), # two digit for state
  stringr::str_pad(`County Code`, 3, "left", pad = "0"),# three digits for county
  stringr::str_pad(`Site Number`, 4, "left", pad = "0"))]  # four digits site

# Join Land Use and Location Setting fields to the daily measurements
setkey(dailydt, stn)
setkey(sites, stn)
dailydt[sites, c("Land_Use", "Location_Setting") := .(`Land Use`, `Location Setting`)]
rm(sites)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fixing Inconsistent Datums. Hourly, then Daily      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# general function for converting CRS
toEPSG <- function(x,y, from_epsg, to_epsg){ 
  fromPoint = st_sfc(st_point(c(x, y)), crs = from_epsg)
  toPoint = st_transform(fromPoint, crs = to_epsg)
  return(toPoint[[1]])}

#
# Transform coordinates to WGS84 where the Datum is known
uniquelonglatdatum <- unique(dailydt[, .(Longitude, Latitude, Datum)])
uniquelonglatdatum[Datum == "NAD27", datum_epsg := 4267]
uniquelonglatdatum[Datum == "NAD83", datum_epsg := 4269]
uniquelonglatdatum[Datum == "WGS84", datum_epsg := 4326]
uniquelonglatdatum[Datum %in% c("NAD27", "NAD83"),
                   c("long_wgs84", "lat_wgs84"):= as.data.table(t(mapply(Longitude, Latitude,
                                                                         from_epsg = datum_epsg, to_epsg = 4326, FUN = toEPSG)))]
# Merge these unique transformed coords back to measurement table
setkey(dailydt, Longitude, Latitude, Datum)
setkey(uniquelonglatdatum, Longitude, Latitude, Datum)
dailydt[uniquelonglatdatum, c("long_wgs84", "lat_wgs84") := list(long_wgs84, lat_wgs84)]
# # if WGS84 was original Datum, copy the coords to the _wgs84 columns
dailydt[Datum == "WGS84", c("long_wgs84", "lat_wgs84") := list(Longitude, Latitude)]
rm(uniquelonglatdatum)

# Number of stations with UNKNOWN Datums (9)
dailydt[Datum == "UNKNOWN", uniqueN(stn)]
# How many of these stations sometimes did have a known datum? (0)
dailydt[stn %in% dailydt[Datum == "UNKNOWN", unique(stn)] & Datum != "UNKNOWN" ]
# Zero "unknown" for NO2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Monthly averages by unique station number      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Create df for seasonal calculations from daily data
seasonaldt <- dailydt[,1:35]

# create "month" columns in list for monthly average - a placeholder setting date first day of given month
dailydt[, month := floor_date(as_date(dailydt$`Date Local`), "month")]

# create new variable by averaging all daily obs in a month for a unique station
dailydt[,monthly_mean := mean(`Arithmetic Mean`, na.rm=T), .(stn, month)]

# remove duplicates and create new data frame with mothly average data
monthlydt <- dailydt[!duplicated(dailydt[,31:37])]

# creat a year variable in seasonal data frame
seasonaldt$year <- format(seasonaldt$day,"%Y")

# create a season variable in seasonal data frame
seasonaldt$seas <- mkseas(x = seasonaldt$day,width='DJF')

# seasonal mean using daily data - take average by season and year and station
seasonaldt[,seas_mean := mean(`Arithmetic Mean`, na.rm=T), .(stn, seas, year)]

#remove duplicates so only seasonal average remain
seasonaldt<- seasonaldt[!duplicated(seasonaldt[,31:38])]

# seasonal mean using monthly averages copy monthly data frame
seasonaldt_month <-monthlydt

# creat a year variable for seasonal data frame
seasonaldt_month$year <- format(seasonaldt_month$day,"%Y")

# create a season variable for seasonal data frame
seasonaldt_month$seas <- mkseas(x = seasonaldt_month$day,width='DJF')

# seasonal mean using monthly mean -  averaging by season and year and station
seasonaldt_month[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]

# remove monthly mean variable and month variable
seasonaldt_month[,36:37]<-NULL

# remove duplicates so only seasonal averages remain
seasonaldt_month_unique <- seasonaldt_month[!duplicated(seasonaldt_month[,31:38])]

# remove data frame with duplicates
rm(seasonaldt_month)

summer_seasonldt <-seasonaldt[which(seasonaldt$seas == 'JJA'),]

summer_daily_mean_all <- summer_seasonldt[,c("stn","year","seas_mean")]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SAVE OUTPUT FOR PLOTTING

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(summer_daily_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_daily_mean_all.csv')


summer_seasonl_mon <-seasonaldt_month_unique[which(seasonaldt_month_unique$seas == 'JJA'),]

summer_month_mean_all <- summer_seasonl_mon[,c("stn","year","seas_mean")]

write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_month_mean_all.csv')

#create data frame of the lat lon for each unique site (stn)
lat_lon_EPA_site <- unique(seasonaldt[,.(Latitude,Longitude,stn)])

write.csv(lat_lon_EPA_site,'~/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site.csv')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial Plots of Data                                     ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot to check differences between seasonal mean using daily data vs monthly average
ggplot() +
  geom_point(aes(x=seasonaldt$seas_mean,y=seasonaldt_month_unique$seas_mean),alpha = .6) + 
  labs(x='Daily data',y='Monthly Averages') +
  ggtitle("Seasonal Means of NO2 using Daily data vs Monthly averages")
ggsave(paste0(fig_dir,"Daily_vs_Monthly_seasonal_averages.png"), width = 6, height = 3.5)

ggplot() +
  geom_point(aes(x=seasonaldt$seas_mean[seasonaldt_month_unique$seas=='JJA'],y=seasonaldt_month_unique$seas_mean[seasonaldt_month_unique$seas=='JJA']),alpha = .6) + 
  labs(x='Daily data',y='Monthly Averages') +
  ggtitle("Summer Means of NO2 using Daily data vs Monthly averages")
ggsave(paste0(fig_dir,"Daily_vs_Monthly_seasonal_JJA_averages.png"), width = 6, height = 3.5)


# plot number of sites with monthly average data per month over time
# use number/count of "month" to indicate unique sites 
ggplot(monthlydt) + aes(x=month) + geom_bar() + labs(x='Time [month]',y='Number of EPA sites with NO2 obs')
ggsave(paste0(fig_dir,"EPA_NO2_n_monthly_obs.png"), width = 6, height = 3.5)

#plot mean, min, max of monthly average data per month across all sites over time
monthlydt %>% group_by(month) %>%
  summarise(min = min(`Arithmetic Mean`, na.rm = TRUE),
            max = max(`Arithmetic Mean`, na.rm = TRUE),
            avg = mean(`Arithmetic Mean`,na.rm = TRUE)) %>%
  gather(metric, value, -month) %>%
  ggplot(.,aes(x = month, y = value, 
               group = metric, color = metric)) + labs(x='Time [month]',y='Monthly average NO2 [ppb]') +
  geom_line()
ggsave(paste0(fig_dir,"Ave_min_max_NO2_monthly.png"), width = 6, height = 3.5)

# Load OMI Model data 
#OMI_NO2 <- readMat("/Users/yshiga/Documents/Research/AURA/Data/Model_output/aura_omi_L3_no2_01_v2.mat")
OMI_NO2_EPA <- readMat("/Users/yshiga/Documents/Research/AURA/Data/Model_output/data_EPA.mat")

# Will need to find closest OMI Lat Lon for each EPA data site

# First organize OMI data (remove NaNs, gather lat lon info)
# find index of rows that are finite (not NaNs)
rows_to_keep_index <- !rowSums(!is.finite(OMI_NO2$data.NO2.01))

# expand lat lon to create grid (similar to meshgrid from Matlab)
omi_grid_lat_lon <- expand.grid(lon = unique(OMI_NO2$lons.01),lat = unique(OMI_NO2$lats.01))

# remove lat lon corresponding to data NaNs
omi_grid_lat_lon_finite = omi_grid_lat_lon[rows_to_keep_index,]

# search nearest grid cell to each EPA site
lat_lon_EPA_site$ind_omi <- 0

for (i in 1:nrow(lat_lon_EPA_site)){
lat_inds<-which(abs((omi_grid_lat_lon_finite$lat-lat_lon_EPA_site$Latitude[i])) == min(abs((omi_grid_lat_lon_finite$lat-lat_lon_EPA_site$Latitude[i]))))

lon_inds<-which(min(abs(omi_grid_lat_lon_finite$lon-lat_lon_EPA_site$Longitude[i])) == abs(omi_grid_lat_lon_finite$lon-lat_lon_EPA_site$Longitude[i]))
ind_pick <- intersect(lat_inds,lon_inds)

lat_lon_EPA_site$ind_omi[i] <- ind_pick
}

head(lat_lon_EPA_site)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## HCHO #######

# Libraries ---------------------------------------------------------------

suppressMessages(library(data.table))
suppressMessages(library(fst))
library(sf)
library(nngeo)
library(ggplot2)
library(dplyr)
library(lubridate)
library(zoo)
library(seas)
library(R.matlab)
library(sp)
library(rgeos)


data_dir = "~/Documents/Research/AURA/Data/Daily/HAPs"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Daily
storewd <- getwd()
# change working directory to where the data live
filepath <- data_dir # "data/aqs_raw/epa_daily"
setwd(filepath)
#list of data in folder
zipfiles <- list.files(pattern = "daily.*\\.zip")
# Load all daily data in folder in list
dailydt <- rbindlist(lapply(zipfiles, function(x) fread(cmd=paste0("unzip -cq ", x))))

# restore the working directory
setwd(storewd)
#clean up
rm(zipfiles, filepath)

# how big is the data
paste("The loaded daily HAPs measurement data table has dimensions", paste(dim(dailydt), collapse = ",")) 
# pryr::object_size(dailydt) #  For HAPs = 1.23 GB

dailydt <- dailydt[(`Parameter Name` %in%  'Formaldehyde'),]
paste("Formaldehyde measurements subset data table has dimensions", paste(dim(dailydt), collapse = ",")) 
pryr::object_size(dailydt) #  For HAPs (HCHO only) = 26.3  MB


# clean the day
dailydt[, day := as.Date(`Date Local`)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Variable construction. 
# Create unique station identifier
# Pull Land Use and Location Setting info from aqs_sites.cvs (from EPA) ####
# Downloaded aqs_sites.zip from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Meta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Daily
# create a unique station identifier with uniform length
dailydt[, stn := paste0(
  stringr::str_pad(`State Code`, 2, "left", pad = "0"), # two digit for state
  stringr::str_pad(`County Code`, 3, "left", pad = "0"),# three digits for county
  stringr::str_pad(`Site Num`, 4, "left", pad = "0"))]  # four digits site

# Open site description file
sites <- fread("~/Documents/Research/AURA/Data/Daily/aqs_sites.csv")

# In site table, create a unique station identifier with uniform length
sites[, stn := paste0(
  stringr::str_pad(`State Code`, 2, "left", pad = "0"), # two digit for state
  stringr::str_pad(`County Code`, 3, "left", pad = "0"),# three digits for county
  stringr::str_pad(`Site Number`, 4, "left", pad = "0"))]  # four digits site

# Join Land Use and Location Setting fields to the daily measurements
setkey(dailydt, stn)
setkey(sites, stn)
dailydt[sites, c("Land_Use", "Location_Setting") := .(`Land Use`, `Location Setting`)]
rm(sites)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fixing Inconsistent Datums. Hourly, then Daily      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# general function for converting CRS
toEPSG <- function(x,y, from_epsg, to_epsg){ 
  fromPoint = st_sfc(st_point(c(x, y)), crs = from_epsg)
  toPoint = st_transform(fromPoint, crs = to_epsg)
  return(toPoint[[1]])}

#
# Transform coordinates to WGS84 where the Datum is known
uniquelonglatdatum <- unique(dailydt[, .(Longitude, Latitude, Datum)])
uniquelonglatdatum[Datum == "NAD27", datum_epsg := 4267]
uniquelonglatdatum[Datum == "NAD83", datum_epsg := 4269]
uniquelonglatdatum[Datum == "WGS84", datum_epsg := 4326]
uniquelonglatdatum[Datum %in% c("NAD27", "NAD83"),
                   c("long_wgs84", "lat_wgs84"):= as.data.table(t(mapply(Longitude, Latitude,
                                                                         from_epsg = datum_epsg, to_epsg = 4326, FUN = toEPSG)))]
# Merge these unique transformed coords back to measurement table
setkey(dailydt, Longitude, Latitude, Datum)
setkey(uniquelonglatdatum, Longitude, Latitude, Datum)
dailydt[uniquelonglatdatum, c("long_wgs84", "lat_wgs84") := list(long_wgs84, lat_wgs84)]
# # if WGS84 was original Datum, copy the coords to the _wgs84 columns
dailydt[Datum == "WGS84", c("long_wgs84", "lat_wgs84") := list(Longitude, Latitude)]
rm(uniquelonglatdatum)

# Number of stations with UNKNOWN Datums (9)
dailydt[Datum == "UNKNOWN", uniqueN(stn)]
# How many of these stations sometimes did have a known datum? (0)
dailydt[stn %in% dailydt[Datum == "UNKNOWN", unique(stn)] & Datum != "UNKNOWN" ]
# Zero "unknown" for NO2

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Monthly averages by unique station number      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Create df for seasonal calculations from daily data
seasonaldt <- dailydt[,1:35]

# create "month" columns in list for monthly average - a placeholder setting date first day of given month
dailydt[, month := floor_date(as_date(dailydt$`Date Local`), "month")]

# create new variable by averaging all daily obs in a month for a unique station
dailydt[,monthly_mean := mean(`Arithmetic Mean`, na.rm=T), .(stn, month)]

# remove duplicates and create new data frame with mothly average data
monthlydt <- dailydt[!duplicated(dailydt[,31:37])]

# creat a year variable in seasonal data frame
seasonaldt$year <- format(seasonaldt$day,"%Y")

# create a season variable in seasonal data frame
seasonaldt$seas <- mkseas(x = seasonaldt$day,width='DJF')

# seasonal mean using daily data - take average by season and year and station
seasonaldt[,seas_mean := mean(`Arithmetic Mean`, na.rm=T), .(stn, seas, year)]

#remove duplicates so only seasonal average remain
seasonaldt<- seasonaldt[!duplicated(seasonaldt[,31:38])]

# seasonal mean using monthly averages copy monthly data frame
seasonaldt_month <-monthlydt

# creat a year variable for seasonal data frame
seasonaldt_month$year <- format(seasonaldt_month$day,"%Y")

# create a season variable for seasonal data frame
seasonaldt_month$seas <- mkseas(x = seasonaldt_month$day,width='DJF')

# seasonal mean using monthly mean -  averaging by season and year and station
seasonaldt_month[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]

# remove monthly mean variable and month variable
seasonaldt_month[,36:37]<-NULL

# remove duplicates so only seasonal averages remain
seasonaldt_month_unique <- seasonaldt_month[!duplicated(seasonaldt_month[,31:38])]

# remove data frame with duplicates
rm(seasonaldt_month)

summer_seasonldt <-seasonaldt[which(seasonaldt$seas == 'JJA'),]

summer_daily_mean_all <- summer_seasonldt[,c("stn","year","seas_mean","County Name","City Name","Location_Setting")]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## SAVE OUTPUT FOR PLOTTING

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(summer_daily_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_daily_mean_all_HCHO.csv')


summer_seasonl_mon <-seasonaldt_month_unique[which(seasonaldt_month_unique$seas == 'JJA'),]
summer_month_mean_all <- summer_seasonl_mon[,c("stn","year","seas_mean")]

write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_month_mean_all_HCHO.csv')

#create data frame of the lat lon for each unique site (stn)
lat_lon_EPA_site <- unique(seasonaldt[,.(Latitude,Longitude,stn)])


write.csv(lat_lon_EPA_site,'~/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_HCHO.csv')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial Plots of Data                                     ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot to check differences between seasonal mean using daily data vs monthly average
ggplot() +
  geom_point(aes(x=seasonaldt$seas_mean,y=seasonaldt_month_unique$seas_mean),alpha = .6) + 
  labs(x='Daily data',y='Monthly Averages') +
  ggtitle("Seasonal Means of HCHO using Daily data vs Monthly averages")
ggsave(paste0(fig_dir,"Daily_vs_Monthly_seasonal_averages_HCHO.png"), width = 6, height = 3.5)

ggplot() +
  geom_point(aes(x=seasonaldt$seas_mean[seasonaldt_month_unique$seas=='JJA'],y=seasonaldt_month_unique$seas_mean[seasonaldt_month_unique$seas=='JJA']),alpha = .6) + 
  labs(x='Daily data',y='Monthly Averages') +
  ggtitle("Summer Means of HCHO using Daily data vs Monthly averages")
ggsave(paste0(fig_dir,"Daily_vs_Monthly_seasonal_JJA_averages_HCHO.png"), width = 6, height = 3.5)


# plot number of sites with monthly average data per month over time
# use number/count of "month" to indicate unique sites 
ggplot(monthlydt) + aes(x=month) + geom_bar() + labs(x='Time [month]',y='Number of EPA sites with HCHO obs')
ggsave(paste0(fig_dir,"EPA_HCHO_n_monthly_obs.png"), width = 6, height = 3.5)

#plot mean, min, max of monthly average data per month across all sites over time
monthlydt %>% group_by(month) %>%
  summarise(min = min(`Arithmetic Mean`, na.rm = TRUE),
            max = max(`Arithmetic Mean`, na.rm = TRUE),
            avg = mean(`Arithmetic Mean`,na.rm = TRUE)) %>%
  gather(metric, value, -month) %>%
  ggplot(.,aes(x = month, y = value, 
               group = metric, color = metric)) + labs(x='Time [month]',y='Monthly average HCHO [ppb]') +
  geom_line()
ggsave(paste0(fig_dir,"Ave_min_max_HCHO_monthly.png"), width = 6, height = 3.5)


## Pick out cities
# Atlanta, Houston, Philadelphia, Nashville, NYC
# none for Atlanta or Nashville - need to look at county level maybe

sum(summer_seasonl_mon$`City Name`=="Houston")
# 15 in houston

sum(summer_seasonl_mon$`City Name`=="Philadelphia")
# 69 in Philadelphia

sum(summer_seasonl_mon$`City Name`=="New York")
# 81 in New York

