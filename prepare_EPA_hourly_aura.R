  #!/usr/bin/env Rscript
  
  # This script uses a template script from Allan Just to average daily NO2 measurements from EPA data to be used for model data copmarison for AURA project (PI Sajeeev Philip)
  
  # Downloaded hourly EPA data from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Raw
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
  library(tidyr)
  
  fig_dir = '/Users/yshiga/Documents/Research/AURA/Figures/'
  
  data_dir = "~/Documents/Research/AURA/Data/Hourly/"
  
  # Open Measurement Daily Data ---------------------------------------------
  
  # Hourly
  storewd <- getwd()
  
  
  # change working directory to where the data live
  filepath <- data_dir 
  
  setwd(filepath)
  
  #list of data in folder
  zipfiles <- list.files(pattern = "hourly.*\\.zip")
  
  # look at all sites in lower 48
  states_non_contig <- c("Alaska", "Hawaii","Puerto Rico")
  
  # only mid-afternoon data
  time_range <- c("13:00", "14:00","15:00")
  
  # Load all daily data in folder in list
  datalist = list()
  # pre-process due to memory limit
  # loop over 2005:2019 data files - only save contiguous US sites and 1pm 2pm 3pm data
  for (i in 1:15) {
    # ... make some data
    hourlydt <- rbindlist(lapply(zipfiles[i], function(x) fread(cmd=paste0("unzip -cq ", x))))
    
    # only look at mid-afternoon data
    hourlydt <- hourlydt[!(`State Name` %in% states_non_contig) & (`Time Local` %in% time_range), ]
    datalist[[i]] <- hourlydt # add it to your list
  }
  hourlydt = do.call(rbind, datalist)
  
rm(datalist)
paste("PM Hourly 1pm-3pm measurements table subset has dimensions", paste(dim(hourlydt), collapse = ",")) # 1243560,30
pryr::object_size(hourlydt) #  For hourly for 2005 NO2 = 563 MB; For hourly 2005:2009 NO2 = 2.72 GB;For hourly 2005:2014 NO2 = 5.52 GB

# restore the working directory
setwd(storewd)
#clean up
rm(zipfiles, filepath)

# how big is the data
paste("The loaded daily PM measurement data table has dimensions", paste(dim(hourlydt), collapse = ",")) 
pryr::object_size(hourlydt) #  For hourly for 2005 NO2 = 563 MB; For hourly 2005:2009 NO2 = 2.72 GB;For hourly 2005:2014 NO2 = 5.52 GB

# plot number of data
# ggplot(hourlydt) + aes(x=`Date Local`) + geom_bar()  + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x='Time',y='Number of hourly EPA data per day')
# ggsave(paste0(fig_dir,"EPA_NO2_n_hourly_obs_2005_2019.png"), width = 6, height = 3.5)

# clean the day
hourlydt[, day := as.Date(`Date Local`)]

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

# Hourly
# create a unique station identifier with uniform length
hourlydt[, stn := paste0(
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
setkey(hourlydt, stn)
setkey(sites, stn)
hourlydt[sites, c("Land_Use", "Location_Setting") := .(`Land Use`, `Location Setting`)]
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
uniquelonglatdatum <- unique(hourlydt[, .(Longitude, Latitude, Datum)])
uniquelonglatdatum[Datum == "NAD27", datum_epsg := 4267]
uniquelonglatdatum[Datum == "NAD83", datum_epsg := 4269]
uniquelonglatdatum[Datum == "WGS84", datum_epsg := 4326]
uniquelonglatdatum[Datum %in% c("NAD27", "NAD83"),
                   c("long_wgs84", "lat_wgs84"):= as.data.table(t(mapply(Longitude, Latitude,
                                                                         from_epsg = datum_epsg, to_epsg = 4326, FUN = toEPSG)))]
# Merge these unique transformed coords back to measurement table
setkey(hourlydt, Longitude, Latitude, Datum)
setkey(uniquelonglatdatum, Longitude, Latitude, Datum)
hourlydt[uniquelonglatdatum, c("long_wgs84", "lat_wgs84") := list(long_wgs84, lat_wgs84)]
# # if WGS84 was original Datum, copy the coords to the _wgs84 columns
hourlydt[Datum == "WGS84", c("long_wgs84", "lat_wgs84") := list(Longitude, Latitude)]
rm(uniquelonglatdatum)

# Number of stations with UNKNOWN Datums (9)
hourlydt[Datum == "UNKNOWN", uniqueN(stn)]
# How many of these stations sometimes did have a known datum? (0)
hourlydt[stn %in% hourlydt[Datum == "UNKNOWN", unique(stn)] & Datum != "UNKNOWN" ]
# Zero "unknown" for NO2

# Assume UNKNOWN datums are WGS84. Assuming being the error is low if datum is incorrect
hourlydt[Datum == "UNKNOWN", c("long_wgs84", "lat_wgs84") := list(Longitude, Latitude)]
#

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Start filtering code     ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filter Step 2.
# Remove all "Qualifier" data
100*length(hourlydt$`Units of Measure`[hourlydt$Qualifier!=""])/length(hourlydt$POC)
# removes ~1.8% of data
hourlydt <-hourlydt[hourlydt$Qualifier=="",]

# Filter Step 3.
# only use FRM 
100*length(hourlydt$`Units of Measure`[hourlydt$`Method Type`!="FRM"])/length(hourlydt$POC)
# removes ~3.6% of remaining data
hourlydt <-hourlydt[hourlydt$`Method Type`=="FRM"]

# Filter Step 4.
# remove multiple POC > use lowest POC
# setkey to order POC then unique only keeps lowest POC
setkey(hourlydt, POC) 
temp_length<-length(hourlydt$POC)
hourlydt<-unique(setDT(hourlydt), by = 6:13) # removes ~500 data points (less than .01 % of data)
100*(temp_length-length(hourlydt$POC))/temp_length

# Filter Step 4.a
# remove negative
100*length(hourlydt$`Units of Measure`[hourlydt$`Sample Measurement`<0])/length(hourlydt$POC)
# removes ~0.8% of data
hourlydt <-hourlydt[hourlydt$`Sample Measurement`>=0,]
# remove nans
sum(is.nan(hourlydt$`Sample Measurement`))
# no NANs

# Filter Step 5.
# remove >1.5*IQR from 76% (per month x site x year)
# create "month" columns in list for monthly average - a placeholder setting date first day of given month
hourlydt[, month := floor_date(as_date(hourlydt$`Date Local`), "month")]
temp_length<-length(hourlydt$POC)
hourlydt_1p5iqr_filt <-hourlydt %>% group_by(stn,month) %>% filter(`Sample Measurement` < (quantile(`Sample Measurement`, probs=c(.75), na.rm = FALSE)+1.5*IQR(`Sample Measurement`, na.rm = TRUE)))
100*(temp_length-length(hourlydt_1p5iqr_filt$POC))/temp_length
# removes ~7.5% of data

rm(hourlydt) # clear old variable

# create data.table
hourlydt_1p5iqr_filt <- data.table(hourlydt_1p5iqr_filt)

# Processing Step 3.
# create new variable by averaging all hourly obs in a day for a unique station
hourlydt_1p5iqr_filt[,daily_mean := mean(`Sample Measurement`, na.rm=T), .(stn, day)]

# Processing Step 4.

# First count days of data per month
hourlydt_1p5iqr_filt[,mon_day_count := length(unique(`day`)), .(stn, month)]

temp_length<-length(hourlydt_1p5iqr_filt$POC) # count data length

# Filter for at least 5 days of data per month
hourlydt_1p5iqr_filt <-hourlydt_1p5iqr_filt[hourlydt_1p5iqr_filt$`mon_day_count`>=5,]
100*(temp_length-length(hourlydt_1p5iqr_filt$POC))/temp_length # calc data change
# removes ~0.05% of data

# create new variable `monthly_mean` by averaging all daily data in a month for a unique station
hourlydt_1p5iqr_filt[,monthly_mean := mean(`daily_mean`, na.rm=T), .(stn, month)]

# Processing Step 5.

# remove duplicates and create new data frame with monthly average data
monthlydt_1p5iqr__filt <- hourlydt_1p5iqr_filt[!duplicated(hourlydt_1p5iqr_filt[,26:31])]

rm(hourlydt_1p5iqr_filt) # clear old variables

# *****POTENTIAL FILTER ON MONTHLY DATA****
#monthlydt_1p5iqr__filt <-monthlydt_1p5iqr__filt %>% group_by(month) %>% filter(`monthly_mean` < (quantile(`monthly_mean`, probs=c(.75), na.rm = FALSE)+1.5*IQR(`monthly_mean`, na.rm = TRUE)))
# removes ~5.7% data

# monthlydt_3sd_filt <- hourlydt_3sd_filt[!duplicated(hourlydt_3sd_filt[,26:32])]
monthlydt_1p5iqr__filt<-data.table(monthlydt_1p5iqr__filt)
# Add correction HERE
# load data
correction <- readMat('/Users/yshiga/Documents/Research/AURA/Data/Model_output/aura_no2_insitu_correction_v1.mat')

# lat lon temp variable for each unique stn
lat_lon_EPA_site_temp <- monthlydt_1p5iqr__filt[!duplicated(monthlydt_1p5iqr__filt[,c("stn")])]

# create new variable of only lat lon and stn
ll_epa <- lat_lon_EPA_site_temp[,c("Latitude","Longitude","stn")]

rm(lat_lon_EPA_site_temp) # remove temp variable

# Search for nearest grid cell per site
# loop over lat lon and identify lat and lon index for model output
for (i in 1:length(ll_epa$Latitude)){
  min_temp_lat<-(abs(ll_epa$Latitude[i]-correction[["lats.gmi"]]))
  min_temp_lon<-(abs(ll_epa$Longitude[i]-correction[["lons.gmi"]]))
  
  ll_epa$lat_ind[i] <- which(min_temp_lat == min(min_temp_lat))
  ll_epa$lon_ind[i] <- which(min_temp_lon == min(min_temp_lon))
}

# add year variable
monthlydt_1p5iqr__filt$year <- format(monthlydt_1p5iqr__filt$day, "%Y")
# monthlydt_3sd_filt$year <- format(monthlydt_3sd_filt$day, "%Y")

# add season variable
monthlydt_1p5iqr__filt$seas <- mkseas(x = monthlydt_1p5iqr__filt$day, width = "DJF")
# monthlydt_3sd_filt$seas <- mkseas(x = monthlydt_3sd_filt$day, width = "DJF")

# remove non JJA
monthlydt_1p5iqr__filt_JJA <- monthlydt_1p5iqr__filt[seas=='JJA']

rm(monthlydt_1p5iqr__filt) # clear old variable


monthlydt_1p5iqr__filt_JJA$CF <- numeric()
# loop over each site
for (i in 1:length(ll_epa$stn)){
  # loop over all data points per site
  for (j in 1:sum(monthlydt_1p5iqr__filt_JJA$stn==ll_epa$stn[i])){
    temp<-monthlydt_1p5iqr__filt_JJA[monthlydt_1p5iqr__filt_JJA$stn==ll_epa$stn[i],]
    yr_ind <- year(temp$day[j])-2004
    mon_ind <- tolower(month.abb[month(temp$day[j])])
    st_id = ll_epa$stn[i]
    master_ind = monthlydt_1p5iqr__filt_JJA$stn==st_id & monthlydt_1p5iqr__filt_JJA$day==temp$day[j]
    mon_name_cf <- paste("cf.",mon_ind, sep="")
    cf_temp <- correction[[mon_name_cf]][ll_epa$lat_ind[i],ll_epa$lon_ind[i],yr_ind]
    monthlydt_1p5iqr__filt_JJA$CF[master_ind] <- cf_temp
  }
}

monthlydt_1p5iqr__filt_JJA[, monthly_mean_c := CF * monthly_mean]

#monthlydt_1p5iqr__filt_JJA$month_name <- format(monthlydt_1p5iqr__filt_JJA$day, "%m")

## PLOT HISTOGRAMS MONTHLY AVERAGE FOR JJA CORRECTED AND ORIG
mn_val = mean(monthlydt_1p5iqr__filt_JJA$monthly_mean_c)
md_val = median(monthlydt_1p5iqr__filt_JJA$monthly_mean_c)

ggplot(monthlydt_1p5iqr__filt_JJA, aes(x=monthly_mean_c)) +  geom_histogram(binwidth = 1) +
  geom_vline(aes(xintercept=mean(monthly_mean_c),color=sprintf("Mean = %.2f", mn_val)), linetype="dashed", size=1,show.legend=TRUE) +
  geom_vline(aes(xintercept=median(monthly_mean_c),color=sprintf("Median = %.2f", md_val)), linetype="dashed", size=1,show.legend=TRUE) +
  theme(legend.title=element_blank()) +
  xlim(c(0, 40))
ggsave(paste0(fig_dir,"EPA_NO2_hist_monthyl_corrected_obs_2005_2019_filter_.png"), width = 7, height = 3.5)

mn_val = mean(monthlydt_1p5iqr__filt_JJA$monthly_mean)
md_val = median(monthlydt_1p5iqr__filt_JJA$monthly_mean)

ggplot(monthlydt_1p5iqr__filt_JJA, aes(x=monthly_mean)) +  geom_histogram(binwidth = 1) +
  geom_vline(aes(xintercept=mean(monthly_mean),color=sprintf("Mean = %.2f", mn_val)), linetype="dashed", size=1,show.legend=TRUE) +
  geom_vline(aes(xintercept=median(monthly_mean),color=sprintf("Median = %.2f", md_val)), linetype="dashed", size=1,show.legend=TRUE) +
  theme(legend.title=element_blank())+
  xlim(c(0, 40))
ggsave(paste0(fig_dir,"EPA_NO2_hist_monthyl_corig_obs_2005_2019_filter_d.png"), width = 7, height = 3.5)

##
## Calculate Seasonal Averages
##

# seasonal mean using averages  - take average by season and year and station method 
monthlydt_1p5iqr__filt_JJA[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]
monthlydt_1p5iqr__filt_JJA[,seas_mean_c := mean(`monthly_mean_c`, na.rm=T), .(stn, seas, year)]

# remove duplicates so only seasonal average remain
seasonal_1p5iqr__filt_JJA<- monthlydt_1p5iqr__filt_JJA[!duplicated(monthlydt_1p5iqr__filt_JJA[,c("stn","year","seas")])]

# Trim variables
summer_month_mean_all <- seasonal_1p5iqr__filt_JJA[,c("stn","year","Latitude","Longitude","seas_mean","seas_mean_c","County Name","Location_Setting")]

# write summer mean to csv
write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_hourly_mean_v4_corrected.csv')

# write to lat lon stn data to CSV
lat_lon_EPA_site <- unique(summer_month_mean_all[,c("Latitude","Longitude","stn")])

lat_lon_EPA_site_temp <- summer_month_mean_all[!duplicated(summer_month_mean_all[,c("stn")])]

lat_lon_EPA_site_2 <- lat_lon_EPA_site_temp[,c("Latitude","Longitude","stn","County Name","Location_Setting")]

write.csv(lat_lon_EPA_site_2,'~/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_hourly_v4_corrected.csv')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# BUT FIRST...
# intermediate plots

# # histogram plot
# ggplot(hourlydt, aes(x=hourlydt$`Sample Measurement`)) +  geom_histogram(binwidth = 1) +
#   xlim(c(0, 40))
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_filter_1234.png"), width = 6, height = 3.5)

#time series
# hourlydt %>% group_by(day) %>%
#   summarise(min = min(`Sample Measurement`, na.rm = TRUE),
#             median = median(`Sample Measurement`, na.rm = TRUE),
#             avg = mean(`Sample Measurement`,na.rm = TRUE)) %>%
#   gather(metric, value, -day) %>%
#   ggplot(.,aes(x = day, y = value, 
#                group = metric, color = metric)) + labs(x='Time [day]',y='Mid-afternoon NO2 [ppb]') +
#   geom_line()
# ggsave(paste0(fig_dir,"Ave_min_med_NO2_hourly_1234.png"), width = 6, height = 3.5)
# 
# 
# hourlydt[hourlydt$`State Name`=="New York" & hourlydt$`County Name`=="New York" ,] %>% group_by(day) %>%
#   summarise(min = min(`Sample Measurement`, na.rm = TRUE),
#             median = median(`Sample Measurement`, na.rm = TRUE),
#             avg = mean(`Sample Measurement`,na.rm = TRUE)) %>%
#   gather(metric, value, -day) %>%
#   ggplot(.,aes(x = day, y = value, 
#                group = metric, color = metric)) + labs(x='Time [day]',y='Mid-afternoon NO2 [ppb]') +
#   geom_line()
# ggsave(paste0(fig_dir,"Ave_min_med_NO2_hourly_1234_NYC.png"), width = 6, height = 3.5)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Monthly averages by unique station number      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Create df for seasonal calculations from hourly data
# seasonaldt <- hourlydt[,1:30]

# tried other ways to average but didn't use these
# df_month <- data.frame(Date=floor_date(as_date(dailydt$`Date Local`), "month"))
# monthlydt <- aggregate(dailydt$`Arithmetic Mean`,list(format(dailydt$day, "%Y-%m"),dailydt$`Site Num`),mean)


# create new variable by averaging all hourly obs in a month for a unique station
#hourlydt[,monthly_mean := mean(`Sample Measurement`, na.rm=T), .(stn, month)]


# remove duplicates and create new data frame with mothly average data
#monthlydt <- hourlydt[!duplicated(hourlydt[,26:32])]


# add year
#monthlydt$year <- format(monthlydt$day, "%Y")

# add season
#monthlydt$seas <- mkseas(x = monthlydt$day, width = "DJF")

# calculate mean per season within each year
#seasonaldt <- aggregate(`Sample Measurement` ~ seas + year + stn, data = monthlydt, mean)

# AND AGAIN....
# intermediate plot histogram of summer averages

# ggplot(seasonaldt[seasonaldt$seas=="JJA",], aes(x=seasonaldt$`Sample Measurement`[seasonaldt$seas=="JJA"])) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(seasonaldt$`Sample Measurement`[seasonaldt$seas=="JJA"]),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(seasonaldt$`Sample Measurement`[seasonaldt$seas=="JJA"]),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Summer (JJA) averages NO2 in-situ (no std filter)')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_nofilter_1234_JJA.png"), width = 6, height = 3.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# filter per site per month using
#hourlydt_3sd_filt <-hourlydt %>% group_by(stn,month) %>% filter(`Sample Measurement` < (mean(`Sample Measurement`, na.rm = TRUE) + 3*sd(`Sample Measurement`, na.rm = TRUE)))
#hourlydt_1p5iqr_filt <-hourlydt %>% group_by(stn,month) %>% filter(`Sample Measurement` < (quantile(`Sample Measurement`, probs=c(.75), na.rm = FALSE)+1.5*IQR(`Sample Measurement`, na.rm = TRUE)))
#hourlydt_2sd_filt <-hourlydt %>% group_by(stn,month) %>% filter(`Sample Measurement` < (mean(`Sample Measurement`, na.rm = TRUE) + 2*sd(`Sample Measurement`, na.rm = TRUE)))

# ggplot(hourlydt, aes(x=hourlydt$`Sample Measurement`)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(hourlydt$`Sample Measurement`),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(hourlydt$`Sample Measurement`),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Hourly values (all months) NO 3sigma filter NO2 in-situ')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_no_filter_1234_ALL.png"), width = 6, height = 3.5)
# 
# ggplot(hourlydt_3sd_filt, aes(x=hourlydt_3sd_filt$`Sample Measurement`)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(hourlydt_3sd_filt$`Sample Measurement`),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(hourlydt_3sd_filt$`Sample Measurement`),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Hourly values (all months) with 3sigma filter NO2 in-situ')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_3_sig_filter_1234_ALL.png"), width = 6, height = 3.5)
# 
# 
# ggplot(hourlydt_1p5iqr_filt, aes(x=hourlydt_1p5iqr_filt$`Sample Measurement`)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(hourlydt_1p5iqr_filt$`Sample Measurement`),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(hourlydt_1p5iqr_filt$`Sample Measurement`),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Hourly values (all months) 75pct+3IQR filter NO2 in-situ')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_1p5iqr_filter_1234_ALL.png"), width = 6, height = 3.5)

# 
# 
# 
# monthlydt_1p5iqr__filt %>% 
#   filter(monthlydt_1p5iqr__filt$seas=='JJA') %>%
#   ggplot(.,aes(x = month, y = monthly_mean)) + labs(x='Time [month]',y='Mid-afternoon NO2 [ppb]') +
#   geom_point(alpha = 0.1 )
# ggsave(paste0(fig_dir,"ALL_NO2_monthly_IQR_JJA.png"), width = 6, height = 3.5)
# 
# ggplot() + 
#   geom_point(data=monthlydt_1p5iqr__filt[monthlydt_1p5iqr__filt$seas=='JJA',], aes(x=month, y=monthly_mean), color='blue',alpha = 0.1 ) + 
#   geom_point(data=monthlydt_1p5iqr__filt_2[monthlydt_1p5iqr__filt_2$seas=='JJA',], aes(x=month, y=monthly_mean), color='red',alpha = 0.1 ) 
# ggsave(paste0(fig_dir,"ALL_NO2_monthly_IQR_JJA_combo.png"), width = 6, height = 3.5)
# 
# 
# 
# 
# ggplot(monthlydt_1p5iqr__filt,aes(x = month, y = monthly_mean)) + labs(x='Time [month]',y='Mid-afternoon NO2 [ppb]') +
#   geom_point(alpha = 0.1 )
# 
# ggsave(paste0(fig_dir,"ALL_NO2_monthly_IQR.png"), width = 6, height = 3.5)
# 
# # remove mon
# 
# monthlydt_1p5iqr__filt %>% group_by(month) %>%
#   filter(any(monthlydt_1p5iqr__filt$seas=='JJA')) %>%
#   summarise(min = min(`monthly_mean`, na.rm = TRUE),
#             median = median(`monthly_mean`, na.rm = TRUE),
#             avg = mean(`monthly_mean`,na.rm = TRUE)) %>%
#   gather(metric, value, -month) %>%
#   ggplot(.,aes(x = month, y = value,
#                group = metric, color = metric)) + labs(x='Time [month]',y='Mid-afternoon NO2 [ppb]') +
#   geom_line()
# ggsave(paste0(fig_dir,"Ave_min_med_NO2_monthly_IQR_JJA.png"), width = 6, height = 3.5)
# 
# # remove monthly_mean  and year column
# #hourlydt_1p5iqr_filt = subset(hourlydt_1p5iqr_filt, select = -c(monthly_mean,year) )
# 
# #hourlydt_3sd_filt<- data.table(hourlydt_3sd_filt)
# 
# 
# # create new variable by averaging all hourly obs in a day for a unique station
# # hourlydt_1p5iqr_filt[,daily_noon_mean := mean(`Sample Measurement`, na.rm=T), .(stn, day)]
# 
# # create new variable by averaging all hourly obs in a month for a unique station
# #hourlydt_1p5iqr_filt[,monthly_mean := mean(`Sample Measurement`, na.rm=T), .(stn, month)]
# #hourlydt_3sd_filt[,monthly_mean := mean(`Sample Measurement`, na.rm=T), .(stn, month)]
# 
# 
# 
# # create new variable by averaging all daily data in a month for a unique station
# #hourlydt_1p5iqr_filt[,monthly_mean := mean(`daily_noon_mean`, na.rm=T), .(stn, month)]
# #hourlydt_3sd_filt[,monthly_mean := mean(`Sample Measurement`, na.rm=T), .(stn, month)]
# 
# # create new variable by counting number of data monthly average for a unique station
# #hourlydt_1p5iqr_filt[,monthly_count := length(`Sample Measurement`), .(stn, month)]
# #hourlydt_3sd_filt[,monthly_count := length(`Sample Measurement`), .(stn, month)]
# 
# # create new variable by calculating the percentage of data in a monthly average for a unique station
# # hourlydt_1p5iqr_filt[,monthly_pct := (hourlydt_1p5iqr_filt$monthly_count)/(3*as.numeric(as.Date(as.yearmon(hourlydt_1p5iqr_filt$month), frac = 1) - as.Date(as.yearmon(hourlydt_1p5iqr_filt$month)) 
# #                                                                                         + 1))]
# # hourlydt_3sd_filt[,monthly_pct := (hourlydt_3sd_filt$monthly_count)/(3*as.numeric(as.Date(as.yearmon(hourlydt_3sd_filt$month), frac = 1) - as.Date(as.yearmon(hourlydt_3sd_filt$month)) + 1))]
# 
# # remove monthly averages with less than 5% of hourly data
# # hourlydt_1p5iqr_filt <- hourlydt_1p5iqr_filt[hourlydt_1p5iqr_filt$monthly_pct>=.05,]
# # hourlydt_3sd_filt <- hourlydt_3sd_filt[hourlydt_3sd_filt$monthly_pct>=.05,]
# 
# # remove duplicates and create new data frame with monthly average data
# # monthlydt_1p5iqr__filt <- hourlydt_1p5iqr_filt[!duplicated(hourlydt_1p5iqr_filt[,26:31])]
# # monthlydt_3sd_filt <- hourlydt_3sd_filt[!duplicated(hourlydt_3sd_filt[,26:32])]
# 
# # Add correction HERE
# #correction <- readMat('/Users/yshiga/Documents/Research/AURA/Data/Model_output/aura_no2_insitu_correction_v1.mat')
# 
# #lat_lon_EPA_site_temp <- monthlydt_1p5iqr__filt[!duplicated(monthlydt_1p5iqr__filt[,c("stn")])]
# 
# #ll_epa <- lat_lon_EPA_site_temp[,c("Latitude","Longitude","stn")]
# 
# #rm(lat_lon_EPA_site_temp)
# # search for nearest grid cell per site
# 
# # for i = 1 : length(ll_epa)
# # 
# # ind_ll(i,1) = find(min(abs(ll_epa(i,1)-lats_01_NH))==(abs(ll_epa(i,1)-lats_01_NH))); % find nearest latitude to EPA latitude
# # 
# # ind_ll(i,2) = find(min(abs(ll_epa(i,2)-lons_01_NH))==(abs(ll_epa(i,2)-lons_01_NH))); % find nearest longitude to EPA longitude
# # 
# # end
# 
# # loop over lat lon and identify lat and lon index for model output
# 
# for (i in 1:length(ll_epa$Latitude)){
# min_temp_lat<-(abs(ll_epa$Latitude[i]-correction[["lats.gmi"]]))
# min_temp_lon<-(abs(ll_epa$Longitude[i]-correction[["lons.gmi"]]))
# 
# ll_epa$lat_ind[i] <- which(min_temp_lat == min(min_temp_lat))
# ll_epa$lon_ind[i] <- which(min_temp_lon == min(min_temp_lon))
# 
# }
# 
# # need to save correction factor for each site/month/year
# #
# 
# # add year
# monthlydt_1p5iqr__filt$year <- format(monthlydt_1p5iqr__filt$day, "%Y")
# # monthlydt_3sd_filt$year <- format(monthlydt_3sd_filt$day, "%Y")
# 
# # add season
# monthlydt_1p5iqr__filt$seas <- mkseas(x = monthlydt_1p5iqr__filt$day, width = "DJF")
# # monthlydt_3sd_filt$seas <- mkseas(x = monthlydt_3sd_filt$day, width = "DJF")
# 
# # calculate mean per season within each year
# seasonaldt_1p5iqr__filt <- aggregate(monthly_mean~ seas + year + stn, data = monthlydt_1p5iqr__filt, mean)
# # seasonaldt_3sd_filt <- aggregate(monthly_mean~ seas + year + stn, data = monthlydt_3sd_filt, mean)
# 
# # seasonal mean using averages  - take average by season and year and station method 2
# monthlydt_1p5iqr__filt[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]
# monthlydt_3sd_filt[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]
# 
# #remove duplicates so only seasonal average remain
# monthlydt_1p5iqr__filt_uniq<- monthlydt_1p5iqr__filt[!duplicated(monthlydt_1p5iqr__filt[,c("stn","year","seas")])]
# monthlydt_3sd_filt_uniq<- monthlydt_3sd_filt[!duplicated(monthlydt_3sd_filt[,c("stn","year","seas")])]
# 
# # only summer
# summer_seasonaldt_1p5iqr__filt_uniq <-monthlydt_1p5iqr__filt_uniq[which(monthlydt_1p5iqr__filt_uniq$seas == 'JJA'),]
# summer_seasonaldt_3sd_filt_uniq <-monthlydt_3sd_filt_uniq[which(monthlydt_3sd_filt_uniq$seas == 'JJA'),]
# 
# # only summer
# summer_seasonaldt_1p5iqr__filt <-seasonaldt_1p5iqr__filt[which(seasonaldt_1p5iqr__filt$seas == 'JJA'),]
# summer_seasonaldt_3sd_filt_uniq <-seasonaldt_3sd_filt[which(seasonaldt_3sd_filt$seas == 'JJA'),]
# 
# summer_month_mean_all <- summer_seasonaldt_1p5iqr__filt_uniq[,c("stn","year","seas_mean")]
# 
# # write summer mean to csv
# write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_hourly_mean_v2.csv')
# 
# # write to lat lon stn data to CSV
# lat_lon_EPA_site <- unique(summer_seasonaldt_1p5iqr__filt_uniq[,c("Latitude","Longitude","stn")])
# 
# lat_lon_EPA_site_temp <- summer_seasonaldt_1p5iqr__filt_uniq[!duplicated(summer_seasonaldt_1p5iqr__filt_uniq[,c("stn")])]
# 
# lat_lon_EPA_site_2 <- lat_lon_EPA_site_temp[,c("Latitude","Longitude","stn")]
# 
# write.csv(lat_lon_EPA_site_2,'~/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_hourly_v3.csv')
# 
# # sigma data
# summer_month_mean_all <- summer_seasonaldt_3sd_filt_uniq
# write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_hourly_mean_v2_sigma.csv')
# 
# # 
# # intermediate plots
# ggplot(seasonaldt_3sd_filt[seasonaldt_3sd_filt$seas=="JJA",], aes(x=seasonaldt_3sd_filt$monthly_mean[seasonaldt_3sd_filt$seas=="JJA"])) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(seasonaldt_3sd_filt$monthly_mean[seasonaldt_3sd_filt$seas=="JJA"]),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(seasonaldt_3sd_filt$monthly_mean[seasonaldt_3sd_filt$seas=="JJA"]),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Summer (JJA) averages NO2 in-situ (3 std filter)')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_5p_3sd_filter_1234_JJA.png"), width = 6, height = 3.5)
# # 
# 
# # intermediate plots
# ggplot(monthlydt_1p5iqr__filt[monthlydt_1p5iqr__filt$seas=='JJA',], aes(x=monthlydt_1p5iqr__filt$monthly_mean[monthlydt_1p5iqr__filt$seas=='JJA'])) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(monthlydt_1p5iqr__filt$monthly_mean[monthlydt_1p5iqr__filt$seas=='JJA']),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(monthlydt_1p5iqr__filt$monthly_mean[monthlydt_1p5iqr__filt$seas=='JJA']),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Monthly averages (JJA) NO2 in-situ (1.5*IQR std filter)')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_1p5iqr_monthly_JJA.png"), width = 6, height = 3.5)
# 
# 
# # intermediate plots
# ggplot(monthlydt_1p5iqr__filt_uniq[monthlydt_1p5iqr__filt_uniq$seas=="JJA",], aes(x=monthlydt_1p5iqr__filt_uniq$seas_mean[seasonaldt_1p5iqr__filt$seas=="JJA"])) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(aes(xintercept=mean(monthlydt_1p5iqr__filt_uniq$seas_mean[monthlydt_1p5iqr__filt_uniq$seas=="JJA"]),color="mean"), linetype="dashed", size=1,show.legend=TRUE) +
#   geom_vline(aes(xintercept=median(monthlydt_1p5iqr__filt_uniq$seas_mean[monthlydt_1p5iqr__filt_uniq$seas=="JJA"]),color="median"), linetype="dashed", size=1,show.legend=TRUE) +
#   xlim(c(0, 40)) +
#   theme(legend.box = 'horizontal') +
#   guides(fill = guide_legend(override.aes = list(linetype = 0))) +
#   scale_color_manual(name = "Statistics", values = c(median = "blue", mean = "red")) +
#   xlab('Summer (JJA) averages NO2 in-situ (3 std filter)')
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_1p5iqr_JJA.png"), width = 6, height = 3.5)
# 
# # 
# ggplot(seasonaldt_3sd_filt[seasonaldt_3sd_filt$seas=="JJA",], aes(x=seasonaldt_3sd_filt$monthly_mean[seasonaldt_3sd_filt$seas=="JJA"])) +
#   geom_histogram(binwidth = 1) +
# ggplot(monthlydt_1p5iqr__filt_uniq[monthlydt_1p5iqr__filt_uniq$seas=="JJA",], aes(x=monthlydt_1p5iqr__filt_uniq$seas_mean[seasonaldt_1p5iqr__filt$seas=="JJA"])) +
#   geom_histogram(binwidth = 1) 
# 
# library(plyr)
# 
# DF <- rbind(data.frame(fill="3sig", obs=seasonaldt_3sd_filt$monthly_mean[seasonaldt_3sd_filt$seas=="JJA"]),
#             data.frame(fill="1p5iqr", obs=monthlydt_1p5iqr__filt_uniq$seas_mean[monthlydt_1p5iqr__filt_uniq$seas=="JJA"]))
# 
# mu <- ddply(DF, "fill", summarise, grp.mean=mean(obs))
# 
# ggplot(DF, aes(x=obs, color=fill)) +
#   geom_histogram(aes(y=..density..), position="identity", alpha=0.25)+
#   geom_density(alpha=0.36)+
#   geom_vline(data=mu, aes(xintercept=grp.mean, color=fill),
#              linetype="dashed")+
#   scale_color_manual(values=c( "#E69F00", "#56B4E9","#999999"))+
#   scale_fill_manual(values=c( "#E69F00", "#56B4E9","#999999"))+
#   labs(title="Filter Comparison plot",x="NO2(ppb)", y = "Density")+
#   theme_classic()
# 
# ggsave(paste0(fig_dir,"EPA_NO2_hist_hourly_obs_2005_2019_JJA_compare_filt.png"), width = 6, height = 3.5)
