# Version in github

suppressPackageStartupMessages({
  library("mgcv")
  library("dplyr")
})

# Retrieving script input
args = commandArgs(TRUE)
s_input = args[1]
combos <- read.csv("/TEST_STRATA_LOCATION/strata_filename.csv",header=FALSE) # Read in strata to test
strata <- combos[s_input,]
strata = strata[, colSums(is.na(strata))==0]
numStrat = ncol(strata)

# Initialization - Importing data
data <- read.csv("/DATA_LOCATION/data_filename.csv") # Read in data
calfin <- data %>% select(,c("year","julian", "region","new_strata","old_strata","lat","lon","calfin_100m3")) # Select relevant covariate variables
input$calfin_100m3 <- log((input$calfin_100m3)/100 + 1)
input <- na.omit(input)


knots <- list(julian = c(0.5, 366.5))
allDoYDat = data.frame()
plotDat = data.frame(1:365)
names(plotDat)[1] <- "DoY"
newData <- data.frame(year = 0, julian = 1:365, lon = 0, lat = 0)
RSS1 = data.frame()
maxData = 0
maxStrat = 0

# For Plotting Climatology Curves

title = "Strata"
leg = c()

for (i in 1:numStrat) {
  leg[i] = paste("Strata ", strata[1,i], " s(DoY)")
  if (i == numStrat) {
    title = paste(title, strata[1,i])
  } else {
    title = paste(title, strata[1,i], '&')
  }
}
leg[i+1] = "Common s(DoY)"
filename = paste("SAVE_PLOT_LOCATION/", title,".jpeg",sep="")
bootstrapResults <- data.frame(title)
bootstrapResults$pVal <- rep(NA,length(1))

for (i in 1:numStrat) {
  
  if (strata[1,i] == 48) {
    Xvar <- input %>% filter(old_strata == 38 | old_strata == 39) # Combining data in Georges Basin
  } else {
    Xvar <- input %>% filter(old_strata == strata[1,i]) 
  }
  
  Xvar <- na.omit(Xvar)
  
  curGam <- bam(calfin_100m3 ~
                  s(julian, k = 12, bs = 'cc') + 
                  s(year, k = 10) +
                  s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                  ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                  ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                     m = list(c(1, 0.5), NA), k = c(25, 15)), 
                data = Xvar, method = 'fREML', knots = knots, 
                nthreads = c(4,1), discrete = TRUE)
  ndata = summary(curGam)$n
  
  if (ndata > maxData) {
    maxData = ndata
    maxStrat = i
  }
  
  RSS1 = rbind(RSS1,deviance(curGam))
  
  
  Xvar_DoY = Xvar
  cov = predict(curGam, type = "terms")
  sDoY = Xvar_DoY$calfin_100m3 - coef(curGam)[1] - apply(cov,1,sum) + cov[,1] + cov[,4] # adding ti1 as well
  Xvar_DoY$calfin_100m3 = sDoY # contains ti1
  
  allDoYDat = rbind(allDoYDat,Xvar_DoY)
  
  sDoYs = predict(curGam, newData,type="terms",exclude = c('s(year)','s(lon,lat)','ti(julian,year)','ti(julian,lon,lat)','ti(year,lon,lat)'), newdata.guaranteed=TRUE)
  
  plotDat = cbind(plotDat,sDoYs)
  names(plotDat)[ncol(plotDat)] <- paste("Strat_",strata[1,i])
  
  if (i == 1) {
    jpeg(file = filename, 
         width = 2000,
         height = 2000,
         res = 300,
         quality = 100)
    plot(1:365,sDoYs, type = "l", col=i, ylab = 's(DoY)',xlab = 'Day of Year',main=title,ylim=c(-3,3))
  } else {
    lines(1:365,sDoYs,col=i)
  }
  
}


# Combined sDoY and ti1 terms: ECF
allDoYGAM <- bam(calfin_100m3 ~
                   s(julian, k = 12) +
                   ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)),
                 data = allDoYDat, method = 'fREML', knots = knots, 
                 nthreads = c(4,1), discrete = TRUE)

newData <- data.frame(year = 0, julian = 1:365)
sDoY_comb <- predict(allDoYGAM, newData,type="terms", exculde = 'ti(julian,year)') 

plotDat = cbind(plotDat,sDoY_comb[,1])

lines(1:365,sDoY_comb[,1],col=i+1)
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
resAll = na.omit(resAll)
resAll[1:maxData] = NA
for (i in 1:(numStrat-1)) {
  resAll = cbind(resAll,NA)
}
lnNs_gam = resAll
bootCalfin = lnNs_gam

newData <- data.frame(year = rep(1977:2019,each=365), julian = rep(1:365))
cov_allDoY_comb = predict(allDoYGAM, newData, type = "terms") # Does not include intercept
sDoY_ti1_comb <- cbind(newData,cov_allDoY_comb) 
names(sDoY_ti1_comb)[names(sDoY_ti1_comb) == "s(julian)"] <- "s_julian"
names(sDoY_ti1_comb)[names(sDoY_ti1_comb) == "ti(julian,year)"] <- "ti_julianYear"


for (i in 1:numStrat) { 
  if (strata[1,i] == 48) {
    Xvar <- input %>% filter(old_strata == 38 | old_strata == 39) 
  } else {
    Xvar <- input %>% filter(old_strata == strata[1,i]) 
  }
  Xvar <- na.omit(Xvar)
  XvarNoDoY <- Xvar
  
  for (y in 1977:2019) {
    for (d in 1:365) {
      XvarNoDoY$calfin_100m3[XvarNoDoY$year == y & XvarNoDoY$julian == d] <- 
        XvarNoDoY$calfin_100m3[XvarNoDoY$year == y & XvarNoDoY$julian == d] - 
        sDoY_ti1_comb$s_julian[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] - 
        sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] - 
        coef(allDoYGAM)[1] 
    }
  }
  
  # Fitting new data back to GAMs without s(DoY) or ti1 covariate terms
  calfinGam_null0 <- bam(calfin_100m3 ~
                           s(year, k = 10) + 
                           s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                           ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                              m = list(c(1, 0.5), NA), k = c(25, 15)), 
                         data = XvarNoDoY, method = 'fREML', knots = knots, 
                         nthreads = c(4,1), discrete = TRUE)
  cov_null <- predict(calfinGam_null0, type = "terms")
  
  lnN = Xvar
  lnN$calfin_100m3 = NA
  
  
  for (y in 1977:2019) {
    for (d in 1:365) {
      if (NROW(lnN$calfin_100m3[lnN$julian == d & lnN$year == y]) > 1) {
        lnN$calfin_100m3[lnN$year == y & lnN$julian == d] <- 
          apply(cov_null[Xvar$julian == d & Xvar$year == y,],1,sum) + 
          coef(calfinGam_null0)[1] + 
          sDoY_ti1_comb$s_julian[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] + 
          sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] +
          coef(allDoYGAM)[1]  
      } else if (NROW(lnN$calfin_100m3[lnN$julian == d & lnN$year == y]) == 1) {
        lnN$calfin_100m3[lnN$year == y & lnN$julian == d] <- 
          sum(cov_null[Xvar$julian == d & Xvar$year == y,]) + 
          coef(calfinGam_null0)[1] + 
          sDoY_ti1_comb$s_julian[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] + 
          sDoY_ti1_comb$ti_julianYear[sDoY_ti1_comb$year == y & sDoY_ti1_comb$julian == d] +
          coef(allDoYGAM)[1]  
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
allDoYDatBoot = data.frame() 
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
      bootData <- na.omit(bootData)
      bootData$calfin_100m3 = bootSamp + na.omit(lnNs_gam[,i])
      bootCalfin[1:length(bootSamp),i] = bootData$calfin_100m3 #Saving this version of the bootstrapped 
      # data for use outside this for loop
      
      calfinGamBoot <- bam(calfin_100m3 ~
                             s(julian, k = 12, bs = 'cc') + 
                             s(year, k = 10) +
                             s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                             ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)) +
                             ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                                m = list(c(1, 0.5), NA), k = c(25, 15)), 
                           data = bootData, method = 'fREML', knots = knots, 
                           nthreads = c(4,1), discrete = TRUE)
      RSS1Boot = rbind(RSS1Boot,deviance(calfinGamBoot))
      
      XvarBoot_DoY = bootData
      covBoot = predict(calfinGamBoot, type = "terms")
      sDoYBoot = XvarBoot_DoY$calfin_100m3 - coef(calfinGamBoot)[1] - apply(covBoot,1,sum) + covBoot[,1] + covBoot[,4] #adding ti1 as well
      XvarBoot_DoY$calfin_100m3 = sDoYBoot
      
      allDoYDatBoot = rbind(allDoYDatBoot,XvarBoot_DoY)
    }
    
    allDoYGAM_Boot <- bam(calfin_100m3 ~
                            s(julian, k = 12) +
                            ti(julian, year, bs = c('cc', 'tp'), k = c(12, 15)), #adding ti1 as a regressor
                          data = allDoYDatBoot, method = 'fREML', knots = knots, 
                          nthreads = c(4,1), discrete = TRUE)
    
    newData <- data.frame(year = rep(1977:2019,each=365), julian = rep(1:365))
    covBoot_allDoY_comb = predict(allDoYGAM_Boot, newData, type = "terms")
    sDoYBoot_ti1_comb <- cbind(newData,covBoot_allDoY_comb) 
    names(sDoYBoot_ti1_comb)[names(sDoYBoot_ti1_comb) == "s(julian)"] <- "s_julian"
    names(sDoYBoot_ti1_comb)[names(sDoYBoot_ti1_comb) == "ti(julian,year)"] <- "ti_julianYear"
    
    for (i in 1:numStrat) {
      if (strata[1,i] == 48) {
        nullData_Boot <- input %>% filter(old_strata == 38 | old_strata == 39) 
      } else {
        nullData_Boot <- input %>% filter(old_strata == strata[1,i]) 
      }
      nullData_Boot <- na.omit(nullData_Boot)
      nullData_Boot$calfin_100m3 = na.omit(bootCalfin[,i])
      
      for (y in 1977:2019) {
        for (d in 1:365) {
          nullData_Boot$calfin_100m3[nullData_Boot$year == y & nullData_Boot$julian == d] <- 
            nullData_Boot$calfin_100m3[nullData_Boot$year == y & nullData_Boot$julian == d] - 
            sDoYBoot_ti1_comb$s_julian[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] - 
            sDoYBoot_ti1_comb$ti_julianYear[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] -
            coef(allDoYGAM_Boot)[1]
        }
      }
      
      
      # Fitting new data back to GAMs without s(DoY) or ti1 covariate terms
      calfinGam_nullBoot <- bam(calfin_100m3 ~
                                  s(year, k = 10) +
                                  s(lon, lat, k = 30, bs = 'ds', m = c(1, 0.5)) +
                                  ti(lon, lat, year, d = c(2,1), bs = c('ds','tp'),
                                     m = list(c(1, 0.5), NA), k = c(25, 15)),
                                data = nullData_Boot, method = 'fREML', knots = knots, 
                                nthreads = c(4,1), discrete = TRUE)
      cov_nullBoot <- predict(calfinGam_nullBoot, type = "terms")
      
      lnNBoot = nullData_Boot
      lnNBoot$calfin_100m3 = NA
      
      for (y in 1977:2019) {
        for (d in 1:365) {
          if (NROW(lnNBoot$calfin_100m3[lnNBoot$julian == d & lnNBoot$year == y]) > 1) {
            lnNBoot$calfin_100m3[lnNBoot$year == y & lnNBoot$julian == d] <- 
              apply(cov_nullBoot[lnNBoot$julian == d & lnNBoot$year == y,],1,sum) + 
              coef(calfinGam_nullBoot)[1] + 
              sDoYBoot_ti1_comb$s_julian[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] + 
              sDoYBoot_ti1_comb$ti_julianYear[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] +
              coef(allDoYGAM_Boot)[1] 
          } else if (NROW(lnNBoot$calfin_100m3[lnNBoot$julian == d & lnNBoot$year == y]) == 1) {
            lnNBoot$calfin_100m3[lnNBoot$year == y & lnNBoot$julian == d] <- 
              sum(cov_nullBoot[lnNBoot$julian == d & lnNBoot$year == y,]) + 
              coef(calfinGam_nullBoot)[1] + 
              sDoYBoot_ti1_comb$s_julian[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] + 
              sDoYBoot_ti1_comb$ti_julianYear[sDoYBoot_ti1_comb$year == y & sDoYBoot_ti1_comb$julian == d] +
              coef(allDoYGAM_Boot)[1] 
          }
        }
      }
      
      
      RSS0Boot <- rbind(RSS0Boot,sum((na.omit(bootCalfin[,i]) - lnNBoot$calfin_100m3)^2)) 
      
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

