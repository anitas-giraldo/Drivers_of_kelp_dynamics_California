

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


# 1. Load north-central sites ----
sites <- read.csv(paste(d.dir, "bull_kelp_site_region_metadata.csv", sep = '/')) %>%
  mutate_at(vars(site_campus_unique_ID, region), list(as.factor)) %>%
  dplyr::filter(region == "north-central") %>%
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
    den_NERLUE, # stipes of Bull kelp
    # Temperature
    'Mean_Monthly_Temp',
    'Mean_Monthly_Summer_Temp',
    'Mean_Monthly_Upwelling_Temp',
    'Max_Monthly_Anomaly_Upwelling_Temp', 
    'MHW_Upwelling_Days', # log
    # Nitrate
    'Mean_Monthly_Nitrate',
    'Min_Monthly_Nitrate',
    'Max_Monthly_Nitrate',
    'Mean_Monthly_Upwelling_Nitrate',
    'Mean_Monthly_Summer_Nitrate',
    'Max_Monthly_Anomaly_Nitrate', 
    'Days_8N', # Log
    # Bio
    'den_STRPURAD', # log
    'log_prev_year_spores',
    # Waves
    'wh_mean', 
    'wh_max', 
    'wh_95prc', # log
    'UBR_Mean', # Log
    'UBR_Max', # Log
    # Substrate,
    'depth_mean',
    'mean_vrm', # log
    # NPP
    'Mean_Monthly_NPP',
    'Min_Monthly_NPP',
    'Max_Monthly_NPP_Upwelling'
    ) %>% #glimpse()
  # Temperature transformations
  mutate(log_MHW_Upwelling_Days = log(MHW_Upwelling_Days + 1)
         ) %>%
  dplyr::select(-c(MHW_Upwelling_Days)) %>%
  # Nitrate transformations
  mutate(log_Days_8N = log(Days_8N + 1)) %>%
  dplyr::select(-c(Days_8N)) %>%
  # Bio transformations
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1)) %>%
  dplyr::select(-c(den_NERLUE,
                   den_STRPURAD)) %>%
  # Wave transformation
  mutate(log_wh_95prc = log(wh_95prc + 1),
         log_UBR_Mean = log(UBR_Mean + 1),
         log_UBR_Max = log(UBR_Max + 1)) %>%
  dplyr::select(-c(wh_95prc, UBR_Mean, UBR_Max)) %>%
  # Substrate transformations
  mutate(log_mean_vrm = log(mean_vrm + 1)) %>%
  dplyr::select(-c(mean_vrm)) %>%
  drop_na() %>%
  glimpse() 




# 5. Divide data into train and test ----

inTraining <- createDataPartition(dat1$log_den_NERLUE, p = 0.75, list = FALSE)
train.gam <- dat1[ inTraining,]
test.gam  <- dat1[-inTraining,]


# 6. Set parameters to save outputs ----

species <- 'bull_kelp'
region <- 'north-central'
file_name <- paste(paste(species, region, sep = '_'), sep ='/')
file_name



# 7. Define predictor variables ----

names(train.gam)

pred.vars <- names(train.gam)[8:32]

# 8. Define Null model ----

model.null <- gam(log_den_NERLUE ~ 
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

names(out.all) <- 'Nereo.density'
names(var.imp) <- 'Nereo.density'
all.mod.fits <- do.call("rbind",out.all)
all.var.imp <- do.call("rbind",var.imp)


# 13. Save model fits and importance ----
write.csv(mod.table, file = paste(o.dir, paste(file_name, "all.mod.fits.csv", sep = '_'), sep ='/'))
write.csv(out.i, file=paste(o.dir, paste(file_name,  "best_models.csv", sep = '_'), sep="/"))
write.csv(all.var.imp, file=paste(o.dir, paste(file_name,  "all.var.imp.csv", sep = '_'), sep="/"))


