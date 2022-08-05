# Version in github

suppressPackageStartupMessages({
  library("mgcv")
  library("dplyr")
})

# Retrieving script input
args = commandArgs(TRUE)
s_input = args[1]
bootstrapResults <- data.frame(s_input)
bootstrapResults$bootSea <- rep(NA,length(1))
bootstrapResults$bootInt <- rep(NA,length(1))

# Initialization - Importing data
data <- read.csv("/DATA_LOCATION/data_filename.csv") # Read in data

calfin <- data %>% select(,c("year","julian","region","old_strata", 
                             "lat","lon","bathymetry","calfin_100m3"))
input <- calfin
input$calfin_100m3 <- log((input$calfin_100m3)/100 + 1)


  # Running all three models together
  knots <- list(julian = c(0.5, 366.5))
  Xvar <- input %>% filter(old_strata == s_input)
  Xvar <- na.omit(Xvar)
  
  calfinGam_unres <- bam(calfin_100m3 ~
                           s(julian, k = 12, bs = 'cc') + 
                           s(year, k = 10) +
                           s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                           s(bathymetry,k=10) +
                           ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                           ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                              m = list(c(1, 0.5), NA), k = c(25, 12)) + 
                           ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                              m = list(c(1, 0.5), NA), k = c(25, 15)), 
                         data = Xvar, method = 'fREML', knots = knots, 
                         discrete = TRUE)
  RSS1 <- summary(calfinGam_unres)$scale
  
  calfinGam_resSea <- bam(calfin_100m3 ~
                            s(julian, k = 12, bs = 'cc') + 
                            s(year, k = 10) +
                            s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                            s(bathymetry,k=10) +
                            ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                            ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                               m = list(c(1, 0.5), NA), k = c(25, 15)), 
                          data = Xvar, method = 'fREML', knots = knots, 
                          discrete = TRUE)
  RSSo_sea <- summary(calfinGam_resSea)$scale
  
  calfinGam_resInt <- bam(calfin_100m3 ~
                            s(julian, k = 12, bs = 'cc') + 
                            s(year, k = 10) +
                            s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                            s(bathymetry,k=10) +
                            ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                            ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                               m = list(c(1, 0.5), NA), k = c(25, 12)),
                          data = Xvar, method = 'fREML', knots = knots, 
                          discrete = TRUE)
  RSSo_int <- summary(calfinGam_resInt)$scale
  
  n <- 1
  LRS_sea <- n*log(RSSo_sea/RSS1)
  LRS_int <- n*log(RSSo_int/RSS1)
  
  
# Data Prep for Bootstrapping
  rep <- 1:100
  alpha <- 0.05
  tol <- end(rep)[1]*alpha
  
  res_unres <- residuals(calfinGam_unres, type = "response")
  
  fitted_Sea <- predict(calfinGam_resSea,type="response")
  res_resSea <- residuals(calfinGam_resSea, type = "response")
  
  fitted_Int <- predict(calfinGam_resInt,type="response")
  res_resInt <- residuals(calfinGam_resInt, type = "response")
  
  newData_sea <- Xvar
  newData_int <- Xvar
  LRSs_seaBoot <- numeric(length(rep))
  LRSs_intBoot <- numeric(length(rep))
  count_sea = 0
  count_int = 0
  
  for (r in rep) {
    
    # Bootstrapping - Seasonal
    if (count_sea <= 10) {
      boot_Sea = sample(res_resSea,length(res_resSea), replace = TRUE)
      newData_sea$calfin_100m3 = boot_Sea + fitted_Sea
      
      calfinGam_unresBoot <- bam(calfin_100m3 ~
                                   s(julian, k = 12, bs = 'cc') + 
                                   s(year, k = 10) +
                                   s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                   s(bathymetry,k=10) +
                                   ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                                   ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                                      m = list(c(1, 0.5), NA), k = c(25, 12)) + 
                                   ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                                      m = list(c(1, 0.5), NA), k = c(25, 15)), 
                                 data = newData_sea, method = 'fREML', knots = knots, 
                                 discrete = TRUE)
      RSS1Boot_sea <- summary(calfinGam_unresBoot)$scale
      
      calfinGam_resSeaBoot <- bam(calfin_100m3 ~
                                    s(julian, k = 12, bs = 'cc') + 
                                    s(year, k = 10) +
                                    s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                    s(bathymetry,k=10) +
                                    ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                                    ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                                       m = list(c(1, 0.5), NA), k = c(25, 15)), 
                                  data = newData_sea, method = 'fREML', knots = knots, 
                                  discrete = TRUE)
      RSSo_seaBoot <- summary(calfinGam_resSeaBoot)$scale
      
      LRS_seaBoot <- n*log(RSSo_seaBoot/RSS1Boot_sea)
      if (LRS_seaBoot > LRS_sea) {
        count_sea = count_sea+1
      }
    }
    
    # Bootstrapping - Interannual
    if (count_int <= 10) {
      boot_Int = sample(res_resInt,length(res_resInt), replace = TRUE)
      newData_int$calfin_100m3 = boot_Int + fitted_Int
      
      calfinGam_unresBoot <- bam(calfin_100m3 ~
                                   s(julian, k = 12, bs = 'cc') + 
                                   s(year, k = 10) +
                                   s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                   s(bathymetry,k=10) +
                                   ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                                   ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                                      m = list(c(1, 0.5), NA), k = c(25, 12)) + 
                                   ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                                      m = list(c(1, 0.5), NA), k = c(25, 15)), 
                                 data = newData_int, method = 'fREML', knots = knots, 
                                 discrete = TRUE)
      RSS1Boot_int <- summary(calfinGam_unresBoot)$scale
      
      calfinGam_resIntBoot <- bam(calfin_100m3 ~
                                    s(julian, k = 12, bs = 'cc') + 
                                    s(year, k = 10) +
                                    s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                    s(bathymetry,k=10) +
                                    ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                                    ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                                       m = list(c(1, 0.5), NA), k = c(25, 12)),
                                  data = newData_int, method = 'fREML', knots = knots, 
                                  discrete = TRUE)
      RSSo_intBoot <- summary(calfinGam_resIntBoot)$scale
      
      LRS_intBoot <- n*log(RSSo_intBoot/RSS1Boot_int)
      LRSs_intBoot[r] = LRS_intBoot
      
      if (LRS_intBoot > LRS_int) {
        count_int = count_int+1
      }
    }
    
  }
  
bootstrapResults$bootInt = count_int/r
bootstrapResults$bootSea = count_sea/r


write.table(bootstrapResults, file = "/SAVE_RESULTS_LOCATION/filename_Results.csv",  
            append=TRUE, sep = ",", row.names=FALSE, col.names=FALSE)

