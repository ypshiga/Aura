#!/usr/bin/env Rscript

# This script uses a template script from Allan Just to average daily NO2 measurements from EPA data to be used for model data copmarison for AURA project (PI Sajeeev Philip)

# Downloaded daily EPA data from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Daily
# I simply manually clicked on links for 2005-2020 (accessed on 1-25-21)
# For NO2, Ozone and HAPS (for HCHO)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries                                         ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages(library(data.table))
suppressMessages(library(fst))
library(sf)
library(nngeo)
library(ggplot2)
library(dplyr)
library(lubridate)
fig_dir = '/Users/yshiga/Documents/Research/AURA/Figures/'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Open Measurement Daily Data.                        ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Daily
storewd <- getwd()
# change working directory to where the data live
filepath <- "~/Documents/Research/AURA/Data/Daily/NO2" # "data/aqs_raw/epa_daily"
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

# Assume UNKNOWN datums are WGS84. Assuming being the error is low if datum is incorrect
# dailydt[Datum == "UNKNOWN", c("long_wgs84", "lat_wgs84") := list(Longitude, Latitude)]

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Monthly averages by unique station number      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# create "month" columns in list for monthly average - a placeholder setting date first day of given month
dailydt[, month := floor_date(as_date(dailydt$`Date Local`), "month")]

# tried other ways to average but didn't use these
# df_month <- data.frame(Date=floor_date(as_date(dailydt$`Date Local`), "month"))
# monthlydt <- aggregate(dailydt$`Arithmetic Mean`,list(format(dailydt$day, "%Y-%m"),dailydt$`Site Num`),mean)

# create new variable by averaging all daily obs in a month for a unique station
dailydt[,monthly_mean := mean(`Arithmetic Mean`, na.rm=T), .(stn, month)]

# remove duplicates and create new data frame with mothly average data
monthlydt <- dailydt[!duplicated(dailydt[,31:37])]

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output                                            ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


output_file = "/Users/yshiga/Documents/Research/MtSinai/Data/monitors/epa_daily/dailybest_2018.fst"
write_fst(dailybest, output_file, compress = 100)
paste("Wrote output file:", output_file)
paste0("Final table had dimensions ", paste0(dim(dailybest), collapse = ","),
       ", included dates from ", min(dailybest$day), " to ", max(dailybest$day), 
       ", and included ", dailybest[,uniqueN(stn)], " unique stations.")




output_file = "/Users/yshiga/Documents/Research/MtSinai/Data/monitors/epa_daily/hourlydt_2018.fst"
write_fst(hourlydt, output_file, compress = 100)
paste("Wrote output file:", output_file)
paste0("Final table had dimensions ", paste0(dim(hourlydt), collapse = ","),
       ", included dates from ", min(hourlydt$day), " to ", max(hourlydt$day), 
       ", and included ", hourlydt[,uniqueN(stn)], " unique stations.")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filtering - no zero values - only 24 hour
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add number of obs per station
hourlydt[, N_stn := .N, by = stn] # n obs per station
dailybest[, N_stn := .N, by = stn] # n obs per station

# Mean and STD Total num. of non-NA obs. per station
hourlydt[, mean_pm25 := mean(pm25, na.rm =T), by = stn]
hourlydt[, std_pm25 := sd(pm25, na.rm =T), by = stn]
hourlydt[, mean_daily_std_pm25 := mean(daily_std, na.rm =T), by = stn]
hourlydt[, max_daily_std_pm25 := max(daily_std, na.rm =T), by = stn]
hourlydt[, mean_daily_range_pm25 := mean(daily_range, na.rm =T), by = stn]
hourlydt[, max_daily_range_pm25 := max(daily_range, na.rm =T), by = stn]
hourlydt[, mean_max_diff_avg := mean(daily_max-daily_mean, na.rm =T), by = stn]
hourlydt[, mean_max_diff_max := max(daily_max-daily_mean, na.rm =T), by = stn]

dailybest[,mean_pm25_pos := mean(pm25[pm25>=0], na.rm=T), by = stn]
dailybest[,std_pm25_pos := sd(pm25[pm25>=0], na.rm=T), by = stn]

# without zeros (positive or pos)
hourlydt[,mean_pm25_pos := mean(pm25[pm25>=0], na.rm=T), by = stn]
hourlydt[,std_pm25_pos := sd(pm25[pm25>=0], na.rm=T), by = stn]

hourlydt[, mean_daily_std_pm25_pos := mean(daily_std[pm25>=0], na.rm =T), by = stn]
hourlydt[, max_daily_std_pm25_pos := max(daily_std[pm25>=0], na.rm =T), by = stn]
hourlydt[, mean_daily_range_pm25_pos := mean(daily_range[pm25>=0], na.rm =T), by = stn]
hourlydt[, max_daily_range_pm25_pos := max(daily_range[pm25>=0], na.rm =T), by = stn]

hourlydt[, mean_max_diff_avg_pos := mean(daily_max[pm25>=0]-daily_mean[pm25>=0], na.rm =T), by = stn]
hourlydt[, mean_max_diff_max_pos := max(daily_max[pm25>=0]-daily_mean[pm25>=0], na.rm =T), by = stn]

dailybest[, mean_pm25 := mean(pm25, na.rm =T), by = stn]
dailybest[, std_pm25 := sd(pm25, na.rm =T), by = stn]

# # only use 24 hour measures
# 
# dailybest[, mean_pm25_24 := mean(pm25[`Sample Duration` == "24 HOUR"], na.rm =T), by = stn]
# dailybest[, std_pm25_24 := sd(pm25[`Sample Duration` == "24 HOUR"], na.rm =T), by = stn]
# 
# hourlydt[, mean_pm25_24 := mean(pm25[day %in% dailybest$day[dailybest$`Sample Duration` == "24 HOUR"]], na.rm =T), by = stn]
# hourlydt[, std_pm25_24 := sd(pm25[day %in% dailybest$day[dailybest$`Sample Duration` == "24 HOUR"]], na.rm =T), by = stn]
# 
# # only use 24 hour measures and positive
# 
# dailybest[, mean_pm25_24_pos := mean(pm25[`Sample Duration` == "24 HOUR" & pm25>=0], na.rm =T), by = stn]
# dailybest[, std_pm25_24_pos := sd(pm25[`Sample Duration` == "24 HOUR" & pm25>=0], na.rm =T), by = stn]
# 
# 
# hourlydt[, mean_pm25_24_pos := mean(pm25[day %in% dailybest$day[dailybest$`Sample Duration` == "24 HOUR"] & pm25>=0], na.rm =T), by = stn]
# hourlydt[, std_pm25_24_pos := sd(pm25[day %in% dailybest$day[dailybest$`Sample Duration` == "24 HOUR"] & pm25>=0], na.rm =T), by = stn]


# unique by station id
# filter by positive values and 24 hour sample duration
DS_hourly <- unique(hourlydt[,.(Latitude, Longitude, mean_pm25, std_pm25, mean_pm25_pos, std_pm25_pos,mean_daily_std_pm25_pos,mean_daily_range_pm25_pos,max_daily_std_pm25_pos,max_daily_range_pm25_pos , mean_max_diff_avg_pos,mean_max_diff_max_pos,stn , N_stn,Location_Setting, `State Name`)])
DS_daily <- unique(dailybest[,.(Latitude, Longitude, mean_pm25, std_pm25, mean_pm25_pos, std_pm25_pos,  stn,N_stn,Location_Setting,`State Name`)])
#DS_daily <- DS_daily[DS_daily$stn %in% DS_hourly$stn  & !is.nan(DS_daily$mean_pm25_24_pos)]
#DS_hourly <- DS_hourly[DS_hourly$stn %in% DS_daily$stn[!is.nan(DS_daily$mean_pm25_24_pos)]]
DS_daily <- DS_daily[DS_daily$stn %in% DS_hourly$stn]
DS_hourly <- DS_hourly[DS_hourly$stn %in% DS_daily$stn]


ggplot() + 
  geom_point(aes(x=DS_hourly$mean_pm25_pos,y=DS_daily$mean_pm25_pos)) + 
  labs(x='Annual mean using hourly data',y='Annual mean using daily data')
ggsave(paste0(fig_dir,"Compare_annual_mean_hourly_vs_daily.png"), width = 6, height = 3.5)

ggplot() + 
  geom_point(aes(x=DS_hourly$mean_pm25_pos,y=DS_hourly$mean_daily_std_pm25_pos)) + 
  labs(x='Annual mean using hourly data',y='Average daily std using hourly data')
ggsave(paste0(fig_dir,"Compare_annual_mean_hourly_vs_std.png"), width = 6, height = 3.5)

ggplot() + 
  geom_point(aes(x=DS_hourly$mean_pm25_pos,y=DS_hourly$mean_daily_range_pm25_pos)) + 
  labs(x='Annual mean using hourly data',y='Average daily range using hourly data')
ggsave(paste0(fig_dir,"Compare_annual_mean_hourly_vs_range.png"), width = 6, height = 3.5)

ggplot() + 
  geom_point(aes(x=DS_daily$std_pm25_pos,y=DS_hourly$mean_daily_std_pm25_pos)) + 
  labs(x='Annual std using daily data',y='Avg daily std using hourly data')
ggsave(paste0(fig_dir,"Compare_annual_std_daily_vs_avg_daily_std_hourly.png"), width = 6, height = 3.5)

ggplot() + 
  geom_point(aes(x=DS_hourly$std_pm25_pos,y=DS_hourly$max_daily_std_pm25_pos)) + 
  labs(x='Annual std using hourly data',y='Max daily std using hourly data')
ggsave(paste0(fig_dir,"Compare_annual_std_hourly_vs_max_daily_std_hourly.png"), width = 6, height = 3.5)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Stats for plotting                                ####
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# hourlydt[, N_stn := .N, by = stn] # n obs per station
# dailybest[, N_stn := .N, by = stn] # n obs per station
# 
# # Mean and STD Total num. of non-NA obs. per station
# hourlydt[, mean_pm25 := mean(pm25, na.rm =T), by = stn]
# hourlydt[, std_pm25 := sd(pm25, na.rm =T), by = stn]
# 
# dailybest[, mean_pm25 := mean(pm25, na.rm =T), by = stn]
# dailybest[, std_pm25 := sd(pm25, na.rm =T), by = stn]
# 
# 
# # unique by station id:
# DS_hourly <- unique(hourlydt[,.(Latitude, Longitude, mean_pm25, std_pm25, stn , N_stn,Location_Setting, `State Name`)])
# DS_daily <- unique(dailybest[,.(Latitude, Longitude, mean_pm25, std_pm25, stn,N_stn,Location_Setting,`State Name`)])
# DS_daily <- DS_daily[DS_daily$stn %in% DS_hourly$stn]
# DS_hourly <- DS_hourly[DS_hourly$stn %in% DS_daily$stn]

# plot with log10 scale
my_breaks_mean <- round(10^c(.4,.6,.8,1.0,1.2,1.4))

cat("Mean PM2.5 by site:\n")
ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  scale_colour_viridis_c() +
#  scale_colour_viridis_c(name = "mean_pm25", trans = "log10", breaks = my_breaks_mean, labels = my_breaks_mean) +
  labs(x = "Longitude", y = "Latitude", color= "Mean PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data mean PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_mean_hourly_2018_filter.png"), width = 6, height = 3.5)


cat("Mean PM2.5 by site:\n")
ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_daily[,.(mean_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 4) + 
#  scale_colour_viridis_c(name = "mean_pm25", trans = "log10", breaks = my_breaks_mean, labels = my_breaks_mean) +
  scale_colour_viridis_c() +
  
  labs(x = "Longitude", y = "Latitude", color= "Mean PM2.5", 
       size = "Number of\nObservations") +
  theme_minimal() +
  ggtitle("Daily data mean PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
  #theme(legend.position ='right') 

ggsave(paste0(fig_dir,"PM25_mean_daily_2018_filter.png"), width = 6, height = 3.5)

# plot with log10 scale
my_breaks_std <- round(10^c(.2,.4,.6,.8,1.0,1.2,1.4,1.6))

cat("STD PM2.5 by site:\n")
ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(std_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = std_pm25_pos, size = N_stn),
             alpha = 1) +
  #scale_colour_viridis_c(option = "magma", name = "std_pm25", trans = "log10", breaks = my_breaks_std, labels = my_breaks_std) +
  scale_colour_viridis_c() +
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  labs(x = "Longitude", y = "Latitude", color= "STD PM2.5", 
       size = "Number of\nObservations") +
  theme_minimal() +
  ggtitle("Hourly data standard deviation PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_STD_hourly_2018_filter.png"), width = 6, height = 3.5)


cat("STD PM2.5 by site:\n")
ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_daily[,.(std_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = std_pm25_pos, size = N_stn),
             alpha = 1) + 
 # scale_colour_viridis_c(option = "magma", name = "std_pm25", trans = "log10", breaks = my_breaks_std, labels = my_breaks_std) +
  scale_colour_viridis_c() +
# coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 4) + 
  labs(x = "Longitude", y = "Latitude", color= "STD PM2.5", 
       size = "Number of\nObservations") +
  theme_minimal() +
  ggtitle("Daily data standard deviation PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_STD_daily_2018_filter.png"), width = 6, height = 3.5)

ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_daily_std_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_daily_std_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Avg daily std PM2.5", trans = "log10") +
  
  #  scale_colour_viridis_c(name = "mean_pm25", trans = "log10", breaks = my_breaks_mean, labels = my_breaks_mean) +
  labs(x = "Longitude", y = "Latitude", color= "Avg daily std PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data mean PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_avg_std_hourly_2018_filter.png"), width = 6, height = 3.5)

ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_daily_range_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_daily_range_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Avg daily range PM2.5", trans = "log10") +
  labs(x = "Longitude", y = "Latitude", color= "Avg daily std PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data mean PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_avg_range_hourly_2018_filter.png"), width = 6, height = 3.5)


ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(max_daily_range_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = max_daily_range_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Max daily range PM2.5", trans = "log10") +
  labs(x = "Longitude", y = "Latitude", color= "Max daily range PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data mean PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_max_range_hourly_2018_filter.png"), width = 6, height = 3.5)

ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_daily_std_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_daily_std_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Avg daily std PM2.5", trans = "log10") +
  labs(x = "Longitude", y = "Latitude", color= "Avg daily std PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data avg daily std PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_avg_std_hourly_2018_filter.png"), width = 6, height = 3.5)


ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(max_daily_range_pm25_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = max_daily_range_pm25_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Max daily range PM2.5", trans = "log10") +
  labs(x = "Longitude", y = "Latitude", color= "Max daily range PM2.5", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Hourly data: max daily range PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_max_range_hourly_2018_filter.png"), width = 6, height = 3.5)


ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_max_diff_avg_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_max_diff_avg_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = expression(paste(Delta,
    "PM2.5")), trans = "log10") +
  #scale_colour_viridis_c(option="magma",name = expression(paste(Delta, "PM2.5"))) +
  labs(x = "Longitude", y = "Latitude", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Avg difference between mean and max PM2.5") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_avg_mean_max_diff_hourly_2018_filter_log_scale.png"), width = 6, height = 3.5)


ggplot() +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               fill = 'NA', color = "grey") + 
  geom_point(data = na.omit(unique(DS_hourly[,.(mean_max_diff_max_pos, Latitude, Longitude, stn, N_stn)])),
             aes(x=Longitude, y=Latitude, color = mean_max_diff_max_pos, size = N_stn),
             alpha = 1) + 
  # coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = T) + 
  scale_size_area(max_size = 5) + 
  #scale_colour_viridis_c() +
  scale_colour_viridis_c(option="magma",name = "Difference PM2.5", trans = "log10") +
  labs(x = "Longitude", y = "Latitude", 
       size = "Number of\nObservations") +
  # scale_size("RMSE") + 
  theme_minimal() +
  ggtitle("Difference between mean and max PM2.5 by station") +
  guides(color = guide_colourbar(order=1),
         size = guide_legend(order=2))
ggsave(paste0(fig_dir,"PM25_max_mean_max_diff_hourly_2018_filter.png"), width = 6, height = 3.5)

