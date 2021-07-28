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
  
  # data_dir = "~/Documents/Research/AURA/Data/Hourly/"
  
  
  data_dir = "~/Documents/Research/AURA/Data/Hourly/HAPs/"
  
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
paste("HCHO Hourly 1pm-3pm measurements table subset has dimensions", paste(dim(hourlydt), collapse = ",")) # 454568,24
pryr::object_size(hourlydt) #  For hourly for HAPs = 76.7 MB 2005:2015

# restore the working directory
setwd(storewd)
#clean up
rm(zipfiles, filepath, states_non_contig,time_range )

# plot number of data
# ggplot(hourlydt) + aes(x=`Date Local`) + geom_bar()  + scale_x_discrete(labels = NULL, breaks = NULL) + labs(x='Time',y='Number of hourly EPA data per day')
# ggsave(paste0(fig_dir,"EPA_NO2_n_hourly_obs_2005_2019.png"), width = 6, height = 3.5)

hourlydt <- hourlydt[(`Parameter Name` %in%  'Formaldehyde'),]
paste("Formaldehyde measurements subset data table has dimensions", paste(dim(hourlydt), collapse = ",")) 

# clean the day
hourlydt[, day := as.Date(`Date Local`)]

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
sum(is.na(hourlydt$Qualifier))
sum(hourlydt$Qualifier[!is.na(hourlydt$Qualifier)] =="")
# removes ~0% of data

# Filter Step 3.
# only use FRM **** all data are non-FRM **** DID NOT APPLY 

# Filter Step 4.
# remove multiple POC > use lowest POC
# setkey to order POC then unique only keeps lowest POC
setkey(hourlydt, POC) 
temp_length<-length(hourlydt$POC)
hourlydt<-unique(setDT(hourlydt), by = 6:13) # removes 0 data points 
100*(temp_length-length(hourlydt$POC))/temp_length

# Filter Step 4.a
# remove negative
100*length(hourlydt$`Units of Measure`[hourlydt$`Sample Measurement`<0])/length(hourlydt$POC)
# removes ~0 of data
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
# removes ~3.8% of data

rm(hourlydt,temp_length) # clear old variables

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
# removes ~2.5% of data

# create new variable `monthly_mean` by averaging all daily data in a month for a unique station
hourlydt_1p5iqr_filt[,monthly_mean := mean(`daily_mean`, na.rm=T), .(stn, month)]

# Processing Step 5.

# remove duplicates and create new data frame with monthly average data
monthlydt_1p5iqr__filt <- hourlydt_1p5iqr_filt[!duplicated(hourlydt_1p5iqr_filt[,26:31])]

rm(hourlydt_1p5iqr_filt, temp_length) # clear old variables

# *****POTENTIAL FILTER ON MONTHLY DATA****
#monthlydt_1p5iqr__filt <-monthlydt_1p5iqr__filt %>% group_by(month) %>% filter(`monthly_mean` < (quantile(`monthly_mean`, probs=c(.75), na.rm = FALSE)+1.5*IQR(`monthly_mean`, na.rm = TRUE)))
# removes ~5.7% data

monthlydt_1p5iqr__filt<-data.table(monthlydt_1p5iqr__filt)

# add year variable
monthlydt_1p5iqr__filt$year <- format(monthlydt_1p5iqr__filt$day, "%Y")
# monthlydt_3sd_filt$year <- format(monthlydt_3sd_filt$day, "%Y")

# add season variable
monthlydt_1p5iqr__filt$seas <- mkseas(x = monthlydt_1p5iqr__filt$day, width = "DJF")
# monthlydt_3sd_filt$seas <- mkseas(x = monthlydt_3sd_filt$day, width = "DJF")

# remove non JJA
monthlydt_1p5iqr__filt_JJA <- monthlydt_1p5iqr__filt[seas=='JJA']

rm(monthlydt_1p5iqr__filt) # clear old variable

#monthlydt_1p5iqr__filt_JJA$month_name <- format(monthlydt_1p5iqr__filt_JJA$day, "%m")

## PLOT HISTOGRAMS MONTHLY AVERAGE FOR JJA 

mn_val = mean(monthlydt_1p5iqr__filt_JJA$monthly_mean)
md_val = median(monthlydt_1p5iqr__filt_JJA$monthly_mean)

ggplot(monthlydt_1p5iqr__filt_JJA, aes(x=monthly_mean)) +  geom_histogram(binwidth = 1) +
  geom_vline(aes(xintercept=mean(monthly_mean),color=sprintf("Mean = %.2f", mn_val)), linetype="dashed", size=1,show.legend=TRUE) +
  geom_vline(aes(xintercept=median(monthly_mean),color=sprintf("Median = %.2f", md_val)), linetype="dashed", size=1,show.legend=TRUE) +
  theme(legend.title=element_blank())
ggsave(paste0(fig_dir,"EPA_HCHOhist_monthyl_corig_obs_2005_2019_filter_d.png"), width = 7, height = 3.5)

##
## Calculate Seasonal Averages
##

# seasonal mean using averages  - take average by season and year and station method 
monthlydt_1p5iqr__filt_JJA[,seas_mean := mean(`monthly_mean`, na.rm=T), .(stn, seas, year)]

# remove duplicates so only seasonal average remain
seasonal_1p5iqr__filt_JJA<- monthlydt_1p5iqr__filt_JJA[!duplicated(monthlydt_1p5iqr__filt_JJA[,c("stn","year","seas")])]

# Trim variables
summer_month_mean_all <- seasonal_1p5iqr__filt_JJA[,c("stn","year","Latitude","Longitude","seas_mean","County Name","Location_Setting")]

# write summer mean to csv
write.csv(summer_month_mean_all,'~/Documents/Research/AURA/Data/Model_output/summer_hourly_mean_v4_HCHO.csv')

# write to lat lon stn data to CSV
lat_lon_EPA_site <- unique(summer_month_mean_all[,c("Latitude","Longitude","stn")])

lat_lon_EPA_site_temp <- summer_month_mean_all[!duplicated(summer_month_mean_all[,c("stn")])]

lat_lon_EPA_site_2 <- lat_lon_EPA_site_temp[,c("Latitude","Longitude","stn")]

write.csv(lat_lon_EPA_site_2,'~/Documents/Research/AURA/Data/Model_output/lat_lon_EPA_site_hourly_v4_HCHO.csv')

