

### Load libraries ----

library(here)
library(tidyverse)
library(dplyr)
library(caret)
library(visreg)  
library(mgcv)    
library(kernlab) 
library(glmnet) 
library(car)
library(gridExtra)
library(FSSgam) # https://github.com/beckyfisher/FSSgam
library(stringr)


# Clear environment ----
rm(list=ls())

### Set directories ----
m.dir <- here()
d.dir <- here("data")
o.dir <- here("outputs")


# 1. Load southeast sites ----
sites <- read.csv(paste(d.dir, "giant_kelp_site_region_metadata.csv", sep = '/')) %>%
  mutate_at(vars(site_campus_unique_ID, region), list(as.factor)) %>%
  dplyr::filter(region == "southeast") %>%
  glimpse()


# 2. Load data ----
df <- read.csv(paste(d.dir, "Kelp_Predictors_All_CA.csv", sep ='/')) %>% 
  mutate_at(vars(campus, survey_year, site_campus_unique_ID, zone, transect, region), list(as.factor)) %>%
  dplyr::filter(site_campus_unique_ID %in% sites$site_campus_unique_ID) %>%
  droplevels() %>%
  glimpse() 


# 3. Load info on years RCCA ----
years <- read.csv(paste(d.dir, "No_survey_years_per_site.csv", sep ='/')) %>%
  mutate_at(vars(site_campus_unique_ID, region), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(preMHW > 2) %>%
  droplevels() %>%
  glimpse() 


# get the sites for Central Coast model --
# sites with 3 or more survey years previous to MHW
df.region <- df %>% 
  #dplyr::select(-c(latitude, longitude)) %>%
  #right_join(scsites, by = c('site_campus_unique_ID')) %>%
  dplyr::filter(site_campus_unique_ID %in% years$site_campus_unique_ID) %>%
  glimpse() # 


# remove years without wave data --

df.region <- df.region %>%
  dplyr::filter(survey_year != "1999",
                survey_year != "2000",
                survey_year != "2001",
                survey_year != "2002",
                survey_year != "2003") %>%
  droplevels() %>%
  glimpse()


# 4. Choose variables and transform needed ----
names(df.region)

dat1 <- df.region %>%
  mutate_at(vars(campus, survey_year, site_campus_unique_ID, zone, transect, region), list(as.factor)) %>% 
  dplyr::select(
    # Factors 
    latitude, longitude,
    site_campus_unique_ID, survey_year, transect, zone, region,
    den_MACSTIPES, # stipes of Macrocystis pyrifera
    # Temperature
    'Min_Monthly_Temp',
    'Mean_Monthly_Upwelling_Temp',
    'Max_Monthly_Anomaly_Upwelling_Temp', # Log
    'Days_18C', # Log
    'Days_19C', # Log
    'Days_20C', # Log
    'Days_21C', # Log
    'Degree_Days_23C', # log
    # Nitrate
    'Mean_Monthly_Summer_Nitrate',
    'Max_Monthly_Anomaly_Summer_Nitrate', # log
    'Min_Monthly_Anomaly_Upwelling_Nitrate',
    'Days_2N', # Log
    'Days_3N', # Log
    'Days_4N', # Log
    # Bio
    'den_STRPURAD', # log
    'log_prev_year_spores',
    # Waves
    'wh_mean', # log
    'wh_99prc', # Log
    'UBR_Max', # Log
    # Substrate,
    'depth_mean',
    'mean_vrm', # log
    # NPP
    'Mean_Monthly_NPP',
    'Mean_Monthly_NPP_Upwelling',
    'Min_Monthly_NPP_Upwelling', # Log
    'Max_Monthly_NPP_Upwelling'
    ) %>% #glimpse()
  # Temperature transformations
  mutate(log_Max_Monthly_Anomaly_Upwelling_Temp = log(Max_Monthly_Anomaly_Upwelling_Temp + 1),
         log_Days_18C = log(Days_18C + 1),
         log_Days_19C = log(Days_19C + 1),
         log_Days_20C = log(Days_20C + 1),
         log_Days_21C = log(Days_21C + 1),
         log_Degree_Days_23C = log(Degree_Days_23C + 1)
         ) %>%
  dplyr::select(-c(Days_18C, Days_19C, Days_20C, Days_21C, Degree_Days_23C, Max_Monthly_Anomaly_Upwelling_Temp)) %>%
  # Nitrate transformations
  mutate(log_Days_2N = log(Days_2N + 1),
         log_Days_3N = log(Days_3N + 1),
         log_Days_4N = log(Days_4N + 1),
         log_Max_Monthly_Anomaly_Summer_Nitrate = log(Max_Monthly_Anomaly_Summer_Nitrate + 1)) %>%
  dplyr::select(-c(Days_2N, Days_3N, Days_4N, Max_Monthly_Anomaly_Summer_Nitrate)) %>%
  
  # Bio transformations
  mutate(log_den_MACSTIPES = log(den_MACSTIPES + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_MACSTIPES,
                   den_STRPURAD)) %>%
  # Wave transformation
  mutate(log_wh_mean = log(wh_mean + 1),
         log_wh_99prc = log(wh_99prc + 1),
         log_UBR_Max = log(UBR_Max + 1)) %>%
  dplyr::select(-c(wh_mean, wh_99prc, UBR_Max)) %>%
  # Substrate transformations
  mutate(log_mean_vrm = log(mean_vrm + 1)) %>%
  dplyr::select(-c(mean_vrm)) %>%
  # NPP transformations
  mutate(log_Min_Monthly_NPP_Upwelling = log(Min_Monthly_NPP_Upwelling + 1)) %>%
  dplyr::select(-c(Min_Monthly_NPP_Upwelling)) %>%
  drop_na() %>%
  glimpse() 




# 5. Divide data into train and test ----

inTraining <- createDataPartition(dat1$log_den_MACSTIPES, p = 0.75, list = FALSE)
train.gam <- dat1[ inTraining,]
test.gam  <- dat1[-inTraining,]


# 6. Set parameters to save outputs ----

species <- 'giant_kelp'
region <- 'southeast'
file_name <- paste(paste(species, region, sep = '_'), sep ='/')
file_name



# 7. Define predictor variables ----

names(train.gam)

pred.vars <- names(train.gam)[8:33]

# 8. Define Null model ----

model.null <- gam(log_den_MACSTIPES ~ 
                  s(site_campus_unique_ID, zone, bs = 're') +
                  s(survey_year, bs = 're') ,
                data = train.gam, 
                family = tw(),
                method = "REML") 

# 9. Define model set up ----

model.set <- generate.model.set(use.dat = train.gam,
                                test.fit = model.null,
                                pred.vars.cont = pred.vars,
                                max.predictors = 7,
                                cov.cutoff = 0.65, # cut off for correlations
                                k = 4,
                                null.terms = "s(site_campus_unique_ID, zone, bs = 're') +
                                s(survey_year, bs = 're')")


# 10. Run the full subset model selection ---- 

out.list <- fit.model.set(model.set,
                          max.models= 500,
                          parallel=T)



# 11. Model fits and importance ----
out.all=list()
var.imp=list()


# 12. Put results in a table ----
mod.table <- out.list$mod.data.out
mod.table <- mod.table[order(mod.table$AICc),]
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
out.i <- mod.table[which(mod.table$delta.AICc<=3),]
nrow(out.i)

out.all <- c(out.all,list(out.i)) 
var.imp <- c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))

names(out.all) <- 'Macro.density'
names(var.imp) <- 'Macro.density'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)


# 13. Save model fits and importance ----
write.csv(mod.table, file = paste(o.dir, paste(file_name, "all.mod.fits.csv", sep = '_'), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(file_name,  "best_models.csv", sep = '_'), sep="/"))
write.csv(all.var.imp, file=paste(o.dir, paste(file_name,  "all.var.imp.csv", sep = '_'), sep="/"))


