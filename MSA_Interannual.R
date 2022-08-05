# Version in github

suppressPackageStartupMessages({
  library("mgcv")
  library("dplyr")
})

# Retrieving script input
args = commandArgs(TRUE)
s_input = args[1]
years_test = 2019
combos <- read.csv("/TEST_STRATA_LOCATION/strata_filename.csv",header=FALSE) # Read in strata to test
# Rows 1-20 are routine; Rows 21-31 are strata of interest; 32-36 are special interest; 44-47 are CAST, Rows 48-49 are for the paper
strata <- combos[s_input,]
strata = strata[, colSums(is.na(strata))==0]
numStrat = ncol(strata)

# Initialization - Importing data
data <- read.csv("/DATA_LOCATION/data_filename.csv") # Read in data
calfin <- data %>% select(,c("year","julian", "region","new_strata","old_strata","lat","lon","calfin_100m3")) # Select relevant covariate variables
input$calfin_100m3 <- log((input$calfin_100m3)/100 + 1)
input <- na.omit(input)

# Alternative models and plotting

knots <- list(julian = c(0.5, 366.5))
allYearDat = data.frame()
plotDat = data.frame(1977:years_test)
names(plotDat)[1] <- "Year"
newData <- data.frame(year = 1977:years_test, julian = 0, lon = 0, lat = 0)
RSS1 = data.frame()
maxData = 0
maxStrat = 0

title = "Strata"
leg = c()

for (i in 1:numStrat) {
  leg[i] = paste("Strata ", strata[1,i], " s(Year)")
  if (i == numStrat) {
    title = paste(title, strata[1,i])
  } else {
    title = paste(title, strata[1,i], '&')
  }
}

leg[i+1] = "Common s(Year)"
filename = paste("SAVE_PLOT_LOCATION/", title,".jpeg",sep="")
bootstrapResults <- data.frame(title)
bootstrapResults$pVal <- rep(NA,length(1))

for (i in 1:numStrat) {
  if (strata[1,i] == 48) {
    Xvar <- input %>% filter(old_strata == 38 | old_strata == 39) 
  } else {
    Xvar <- input %>% filter(old_strata == strata[1,i]) 
  }
  Xvar <- na.omit(Xvar)
  
  curGam <- bam(calfin_100m3 ~
                  s(julian, k = 12, bs = 'cc') + 
                  s(year, k = 10) +
                  s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                  ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                  ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                     m = list(c(1, 0.5), NA), k = c(25, 12)), 
                data = Xvar, method = 'fREML', knots = knots, 
                nthreads = c(4,1), discrete = TRUE)
  ndata = summary(curGam)$n
  
  if (ndata > maxData) {
    maxData = ndata
    maxStrat = i
  }
  
  RSS1 = rbind(RSS1,deviance(curGam))
  
  
  Xvar_year = Xvar
  cov = predict(curGam, type = "terms")
  sYear = Xvar_year$calfin_100m3 - coef(curGam)[1] - apply(cov,1,sum) + cov[,2] + cov[,4] # adding ti1 as well
  Xvar_year$calfin_100m3 = sYear
  
  allYearDat = rbind(allYearDat,Xvar_year)
  
  sYears = predict(curGam, newData,type="terms",exclude = c('s(julian)','s(lon,lat)','ti(julian,year)','ti(julian,lon,lat)','ti(year,lon,lat)'), newdata.guaranteed=TRUE)
  
  plotDat = cbind(plotDat,sYears)
  names(plotDat)[ncol(plotDat)] <- paste("Strat_",strata[1,i])
  
  if (i == 1) {
    jpeg(file = filename, 
         width = 2000,
         height = 2000,
         res = 300,
         quality = 100)
    plot(1977:years_test,sYears, type = "l", col=i, ylab = 's(Year)',xlab = 'Year',main=title,ylim=c(-1,1))
  } else {
    lines(1977:years_test,sYears,col=i)
  }
  
  
}

# Combined sYear term
allYearGAM <- bam(calfin_100m3 ~
                    s(year, k = 10)+
                    ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)), # adding ti1 as a regressor, 
                  data = allYearDat, method = 'fREML', knots = knots, 
                  nthreads = c(4,1), discrete = TRUE)

newData <- data.frame(year = 1977:years_test, julian = 0)
sYear_comb <- predict(allYearGAM, newData,type="terms", exculde = 'ti(julian,year)') #only extract the common s(Year) term for plotting, no intercept (+ coef(allDoYGAM)[1])
plotDat = cbind(plotDat,sYear_comb)

lines(1977:years_test,sYear_comb[,1],col=i+1)
legend( x = "topright",
        legend = leg,
        col = 0:i+1, lwd=1, lty=c(1,1,1),
        pch=c(NA,NA), cex = 0.5 )
dev.off()

# For Plotting Purposes
names(plotDat)[ncol(plotDat)] <- paste("Combined")
write.csv(plotDat, paste("SAVE_DATA_LOCATION/",title,".csv"))


# Fitting the null models
RSS0 = data.frame()
resAll = (input %>% filter(old_strata == strata[1,maxStrat]))$calfin_100m3
resAll[1:maxData] = NA
for (i in 1:(numStrat-1)) {
  resAll = cbind(resAll,NA)
}
lnNs_gam = resAll
bootCalfin = lnNs_gam

newData <- data.frame(year = rep(1977:years_test,each=365), julian = rep(1:365))
cov_allYear_comb = predict(allYearGAM, newData, type = "terms") # Does not include intercept
sDoY_ti1_comb <- cbind(newData,cov_allYear_comb) 
names(sDoY_ti1_comb)[names(sDoY_ti1_comb) == "s(year)"] <- "s_year"
names(sDoY_ti1_comb)[names(sDoY_ti1_comb) == "ti(julian,year)"] <- "ti_julianYear"
#HERE

for (i in 1:numStrat) {
  if (strata[1,i] == 48) {
    Xvar <- input %>% filter(old_strata == 38 | old_strata == 39) 
  } else {
    Xvar <- input %>% filter(old_strata == strata[1,i]) 
  }
  Xvar <- na.omit(Xvar)
  XvarNoYear <- Xvar
  
  for (y in 1977:years_test) {
    for (d in 1:365) {
      XvarNoYear$calfin_100m3[XvarNoYear$year == y & XvarNoYear$julian == d] <- 
        XvarNoYear$calfin_100m3[XvarNoYear$year == y & XvarNoYear$julian == d] - 
        sDoY_ti1_comb$s_year[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] - 
        sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] - 
        coef(allYearGAM)[1] 
    }
  }
  
  # Fitting new data back to GAMs without s(Year) or ti1 covariate term
  calfinGam_null0 <- bam(calfin_100m3 ~
                           s(julian, k = 12, bs = 'cc') + 
                           s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                           ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                              m = list(c(1, 0.5), NA), k = c(25, 12)), 
                         data = XvarNoYear, method = 'fREML', knots = knots, 
                         nthreads = c(4,1), discrete = TRUE)
  cov_null <- predict(calfinGam_null0, type = "terms")
  
  lnN = Xvar
  lnN$calfin_100m3 = NA
  
  for (y in 1977:years_test) {
    for (d in 1:365) {
      if (NROW(lnN$calfin_100m3[lnN$julian == d & lnN$year == y]) > 1) {
        lnN$calfin_100m3[lnN$year == y & lnN$julian == d] <- 
          apply(cov_null[Xvar$julian == d & Xvar$year == y,],1,sum) + # fitted estimates w/out sDoY or ti1
          coef(calfinGam_null0)[1] + 
          sDoY_ti1_comb$s_year[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] + 
          sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] +
          coef(allYearGAM)[1]  
      } else if (NROW(lnN$calfin_100m3[lnN$julian == d & lnN$year == y]) == 1) {
        lnN$calfin_100m3[lnN$year == y & lnN$julian == d] <- 
          sum(cov_null[Xvar$julian == d & Xvar$year == y,]) + 
          coef(calfinGam_null0)[1] + 
          sDoY_ti1_comb$s_year[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] + 
          sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] +
          coef(allYearGAM)[1]  
      }
    }
  }
  
  RSS0 <- rbind(RSS0,sum((Xvar$calfin_100m3 - lnN$calfin_100m3)^2))
  
  # For bootstrap later
  res = Xvar$calfin_100m3 - lnN$calfin_100m3
  resAll[1:length(res),i] = res
  lnNs_gam[1:length(res),i] = lnN$calfin_100m3
  
}

LRS = log(sum(RSS0)/sum(RSS1))


# Bootstrap Method
rep <- 1:100
count = 0
RSS1Boot = data.frame()
allYearDatBoot = data.frame()
RSS0Boot = data.frame()

for (r in rep) {
  
  if (count <= 10) {
    for (i in 1:numStrat) {
      bootSamp = sample(na.omit(resAll[,i]), length(na.omit(resAll[,i])), replace = TRUE)
      if (strata[1,i] == 48) {
        bootData <- input %>% filter(old_strata == 38 | old_strata == 39) 
      } else {
        bootData <- input %>% filter(old_strata == strata[1,i]) 
      }
      
      bootData$calfin_100m3 = bootSamp + na.omit(lnNs_gam[,i])
      bootCalfin[1:length(bootSamp),i] = bootData$calfin_100m3 #Saving this version of the bootstrapped data for use outside this for loop
      
      calfinGamBoot <- bam(calfin_100m3 ~
                             s(julian, k = 12, bs = 'cc') + 
                             s(year, k = 10) +
                             s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                             ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                             ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                                m = list(c(1, 0.5), NA), k = c(25, 12)), 
                           data = bootData, method = 'fREML', knots = knots, 
                           nthreads = c(4,1), discrete = TRUE)
      RSS1Boot = rbind(RSS1Boot,deviance(calfinGamBoot))
      
      XvarBoot_year = bootData
      covBoot = predict(calfinGamBoot, type = "terms")
      sYearBoot = XvarBoot_year$calfin_100m3 - coef(calfinGamBoot)[1] - apply(covBoot,1,sum) + covBoot[,2] + covBoot[,4] #adding ti1 as well
      XvarBoot_year$calfin_100m3 = sYearBoot
      
      allYearDatBoot = rbind(allYearDatBoot,XvarBoot_year)
    }
    
    allYearGAM_Boot <- bam(calfin_100m3 ~
                             s(year, k = 10)+
                             ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)), #adding ti1 as a regressor, 
                           data = allYearDatBoot, method = 'fREML', knots = knots, 
                           nthreads = c(4,1), discrete = TRUE)
    
    newData <- data.frame(year = rep(1977:years_test,each=365), julian = rep(1:365))
    covBoot_allYear_comb = predict(allYearGAM_Boot, newData, type = "terms")
    sYearBoot_ti1_comb <- cbind(newData,covBoot_allYear_comb) 
    names(sYearBoot_ti1_comb)[names(sYearBoot_ti1_comb) == "s(year)"] <- "s_year"
    names(sYearBoot_ti1_comb)[names(sYearBoot_ti1_comb) == "ti(julian,year)"] <- "ti_julianYear"
    
    for (i in 1:numStrat) {
      if (strata[1,i] == 48) {
        nullData_Boot <- input %>% filter(old_strata == 38 | old_strata == 39) 
      } else {
        nullData_Boot <- input %>% filter(old_strata == strata[1,i]) 
      }
      nullData_Boot <- na.omit(nullData_Boot)
      nullData_Boot$calfin_100m3 = na.omit(bootCalfin[,i])
      
      for (y in 1977:years_test) {
        for (d in 1:365) {
          nullData_Boot$calfin_100m3[nullData_Boot$year == y & nullData_Boot$julian == d] <- 
            nullData_Boot$calfin_100m3[nullData_Boot$year == y & nullData_Boot$julian == d] - 
            sYearBoot_ti1_comb$s_year[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] - 
            sYearBoot_ti1_comb$ti_julianYear[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] -
            coef(allYearGAM_Boot)[1]
        }
      }
      
      # Fitting new data back to GAMs without s(Year) and ti1 covariate terms
      calfinGam_nullBoot <- bam(calfin_100m3 ~
                                  s(julian, k = 12, bs = 'cc') + 
                                  s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                  ti(lon, lat, julian, d = c(2,1), bs = c('ds','cc'), 
                                     m = list(c(1, 0.5), NA), k = c(25, 12)), 
                                data = nullData_Boot, method = 'fREML', knots = knots, 
                                nthreads = c(4,1), discrete = TRUE)
      cov_nullBoot <- predict(calfinGam_nullBoot, type = "terms")
      
      lnNBoot = nullData_Boot
      lnNBoot$calfin_100m3 = NA
      
      for (y in 1977:years_test) {
        for (d in 1:365) {
          if (NROW(lnNBoot$calfin_100m3[lnNBoot$julian == d & lnNBoot$year == y]) > 1) {
            lnNBoot$calfin_100m3[lnNBoot$year == y & lnNBoot$julian == d] <- 
              apply(cov_nullBoot[lnNBoot$julian == d & lnNBoot$year == y,],1,sum) + 
              coef(calfinGam_nullBoot)[1] + 
              sYearBoot_ti1_comb$s_year[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] + 
              sYearBoot_ti1_comb$ti_julianYear[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] +
              coef(allYearGAM_Boot)[1] 
          } else if (NROW(lnNBoot$calfin_100m3[lnNBoot$julian == d & lnNBoot$year == y]) == 1) {
            lnNBoot$calfin_100m3[lnNBoot$year == y & lnNBoot$julian == d] <- 
              sum(cov_nullBoot[lnNBoot$julian == d & lnNBoot$year == y,]) + 
              coef(calfinGam_nullBoot)[1] + 
              sYearBoot_ti1_comb$s_year[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] + 
              sYearBoot_ti1_comb$ti_julianYear[sYearBoot_ti1_comb$year == y & sYearBoot_ti1_comb$julian == d] +
              coef(allYearGAM_Boot)[1] 
          }
        }
      }
      
      
      RSS0Boot <- rbind(RSS0Boot,sum((na.omit(bootCalfin[,i]) - lnNBoot$calfin_100m3)^2))
      #Change nullData_Boot$calfin_100m3 to na.omit(bootCalfin[,j])
      
    }
    
    LRS_Boot = log(sum(RSS0Boot)/sum(RSS1Boot))
    
    if (LRS_Boot > LRS) {
      count = count+1
    }
  }
}

bootstrapResults$pVal <- count/r
write.table(bootstrapResults, file = "/SAVE_RESULTS_LOCATION/filename_Results.csv",  
            append=TRUE, sep = ",", row.names=FALSE, col.names=FALSE)






