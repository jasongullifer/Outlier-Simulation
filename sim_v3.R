# Assumptions:

# Need to set mean, sd, and N as to have appropriate power

# Outliers come from sampling the tail ends of the distribution also
# can specify a percentage of externally generated outliers, default is
# 5%

# Variance the same in each condition


popGen <- function(meanCog=700,meanNcog=700,sd=10,N=1000000){
  #Population generator
  cogs<-rnorm(N,meanCog,sd=100)
  ncogs<-rnorm(N,meanNcog,sd=100)
  return(as.data.frame(cbind(cogs,ncogs)))
}

sampler <- function(data,nsamples,nsims=10000,numExternalOutliers=0.05){
  #Sample nsamples from an input population and add an amount of
  #externally generated outliers. Conduct t-tests on the sample for 3
  #outlier removal procedures: treating cogs and ncogs as one
  #distribution, treating as two distributions, and no outlier
  #corrections.  Repeat this nsmims times.

  simdata<-data.frame()

  cogs.sample<-numeric()
  ncogs.sample<-numeric()

  for(simn in 1:nsims){
    cogs.sample<-sample(data$cogs,nsamples)
    ncogs.sample<-sample(data$ncogs,nsamples)

    RT <- c(cogs.sample,ncogs.sample)
    CogStat <- c(rep("cog",nsamples), rep("ncog", nsamples))
    sample <- data.frame(CogStat, RT)

    #If you want to include a % of random outliers
    if(numExternalOutliers > 0){
      cogrows <- sample(1:length(cogs.sample), numExternalOutliers * length(cogs.sample))
      ncogrows <- sample(1:length(ncogs.sample), numExternalOutliers * length(ncogs.sample))

      #Replace some cognates
      for(x in cogrows){
        if(x %% 2 == 0){ #If even row, make it a long outlier
          cogs.sample[x] <- rnorm(1, mean(cogs.sample)) + 4 * sd(cogs.sample)
        }
        else{ #If odd row, make it a short outlier
          cogs.sample[x] <- rnorm(1, mean(cogs.sample)) - 4 * sd(cogs.sample)
        }
      }

      #Replace some noncognates
      for(x in ncogrows){
        if(x %% 2 == 0){ #If even row, make it a long outlier
          ncogs.sample[x] <- rnorm(1,mean(ncogs.sample)) + 4 * sd(ncogs.sample)
        }
        else{ #If odd row, make it a short outlier
          ncogs.sample[x] <- rnorm(1,mean(ncogs.sample)) - 4 * sd(ncogs.sample)
        }
      }    
    }

    #No outlier treatment (t-test and effect size)
    ttest.nocorr <- t.test(sample$RT[sample$CogStat=="cog"],sample$RT[sample$CogStat=="ncog"])
    es.nocorr <- abs( (mean(sample$RT[sample$CogStat=="cog"]) - mean(sample$RT[sample$CogStat=="ncog"])) / sd(sample$RT))
 
    #Calculate 2.5SD above and below the mean of the sample
    sample<-ddply(sample,.(CogStat),transform,cuthigh=mean(RT) + 2.5*sd(RT), cutlow=mean(RT)-2.5*sd(RT))

    #Treatment as different distributions (t-test and effect size)
    sample.no.dif <- sample[sample$RT > sample$cutlow & sample$RT < sample$cuthigh,]
    ttest.no.dif<-t.test(sample.no.dif$RT[sample.no.dif$CogStat=="cog"],
                         sample.no.dif$RT[sample.no.dif$CogStat=="ncog"])
    es.no.dif <- abs((mean(sample.no.dif$RT[sample.no.dif$CogStat == "cog"])
                      - mean(sample.no.dif$RT[sample.no.dif$CogStat == "ncog"])) / sd(sample.no.dif$RT))
    
    #Treatment as same distribution (t-test and effect size)
    sample.no.same<-sample[sample$RT > mean(sample$RT) - 2.5 * sd(sample$RT)
                           & sample$RT < mean(sample$RT) + 2.5*sd(sample$RT), ]
    ttest.no.same<-t.test(sample.no.same$RT[sample.no.same$CogStat=="cog"],
                          sample.no.same$RT[sample.no.same$CogStat=="ncog"])   
    es.no.same <- abs((mean(sample.no.same$RT[sample.no.same$CogStat == "cog"])
                       - mean(sample.no.same$RT[sample.no.same$CogStat == "ncog"])) / sd(sample.no.same$RT))

    #Assemble the output data frame of t-test results and effect sizes
    curdata<-data.frame(simN=simn,
                        t.dif=abs(as.numeric(ttest.no.dif[1])),
                        p.dif=as.numeric(ttest.no.dif[3]),
                        es.dif=as.numeric(es.no.dif),
                        
                        t.same=abs(as.numeric(ttest.no.same[1])),
                        p.same=as.numeric(ttest.no.same[3]),
                        es.same=as.numeric(es.no.same),
                        
                        t.nocorr=abs(as.numeric(ttest.nocorr[1])),
                        p.nocorr=as.numeric(ttest.nocorr[3]),
                        es.nocorr=as.numeric(es.nocorr)
                        )
    simdata<-rbind(simdata,curdata) #Add to the stack of simulation datasets
  }
  return(simdata)
}

  pop.same<-popGen(700,700) #generate a known "same" population
  pop.dif<-popGen(600,710) #generate a known "difference" population
  
  simdata.same<-sampler(pop.same,50,1000)  #run simulations on same data
  simdata.dif<-sampler(pop.dif,50,1000)    #run simulations on different data
  
  
#################### Rejections of Null#######################
#Same pop #Incorrect rejection of the null
  print("Rejections of Null")
  print("Same pop, treated differently")
  nrow(simdata.same[simdata.same$p.dif < 0.05,])
  mean(simdata.same$t.dif)
  mean(simdata.same$es.dif)
  mean(simdata.same$t.dif[simdata.same$p.dif < 0.05])
  mean(simdata.same$es.dif[simdata.same$p.dif < 0.05])
  
  print("Same pop, treated same")
  nrow(simdata.same[simdata.same$p.same < 0.05,])
  mean(simdata.same$t.same)
  mean(simdata.same$es.same)
  mean(simdata.same$t.same[simdata.same$p.same < 0.05])
  mean(simdata.same$es.same[simdata.same$p.same < 0.05])
  
  print("Same pop, no outlier corrections")
  nrow(simdata.same[simdata.same$p.nocorr < 0.05,])
  mean(simdata.same$t.nocorr)
  mean(simdata.same$es.nocorr)
  mean(simdata.same$t.nocorr[simdata.same$p.nocorr < 0.05])
  mean(simdata.same$es.nocorr[simdata.same$p.nocorr < 0.05])
  
#Dif pop #Correct  rejection of the null
  print("Different pop, treated differently")
  nrow(simdata.dif[simdata.dif$p.dif < 0.05 ,])
  mean(simdata.dif$t.dif)
  mean(simdata.dif$es.dif)
  mean(simdata.dif$t.dif[simdata.dif$p.dif < 0.05])
  mean(simdata.dif$es.dif[simdata.dif$p.dif < 0.05])
  
  print("Different pop, treated same")
  nrow(simdata.dif[simdata.dif$p.same < 0.05,])
  mean(simdata.dif$t.same)
  mean(simdata.dif$es.same)
  mean(simdata.dif$t.same[simdata.dif$p.same < 0.05])
  mean(simdata.dif$es.same[simdata.dif$p.same < 0.05])
  
  print("Different pop, no outlier corrections")
  nrow(simdata.dif[simdata.dif$p.nocorr < 0.05,])
  mean(simdata.dif$t.nocorr)
  mean(simdata.dif$es.nocorr)
  mean(simdata.dif$t.nocorr[simdata.dif$p.nocorr < 0.05])
  mean(simdata.dif$es.nocorr[simdata.dif$p.nocorr < 0.05])
  

#################### Acceptance of Null#######################
#Same pop #Correct acceptance of the null
  print("Acceptance of Null")
  print("Same pop, treated differently")
  nrow(simdata.same[simdata.same$p.dif >= 0.05,])
  mean(simdata.same$t.dif)
  mean(simdata.same$es.dif)
  mean(simdata.same$t.dif[simdata.same$p.dif >= 0.05])
  mean(simdata.same$es.dif[simdata.same$p.dif >= 0.05])
  
  print("Same pop, treated same")
  nrow(simdata.same[simdata.same$p.same >= 0.05,])
  mean(simdata.same$t.same)
  mean(simdata.same$es.same)
  mean(simdata.same$t.same[simdata.same$p.same >= 0.05])
  mean(simdata.same$es.same[simdata.same$p.same >= 0.05])
  
  print("Same pop, no outlier corrections")
  nrow(simdata.same[simdata.same$p.nocorr >= 0.05,])
  mean(simdata.same$t.nocorr)
  mean(simdata.same$es.nocorr)
  mean(simdata.same$t.nocorr[simdata.same$p.nocorr >= 0.05])
  mean(simdata.same$es.nocorr[simdata.same$p.nocorr >= 0.05])

#Dif pop #Incorrect acceptance of the null
  print("Different pop, treated differently")
  nrow(simdata.dif[simdata.dif$p.dif >= 0.05 ,])
  mean(simdata.dif$t.dif)
  mean(simdata.dif$es.dif)
  mean(simdata.dif$t.dif[simdata.dif$p.dif >= 0.05])
  mean(simdata.dif$es.dif[simdata.dif$p.dif >= 0.05])
  
  print("Different pop, treated same")
  nrow(simdata.dif[simdata.dif$p.same >= 0.05,])
  mean(simdata.dif$t.same)
  mean(simdata.dif$es.same)
  mean(simdata.dif$t.same[simdata.dif$p.same >= 0.05])
  mean(simdata.dif$es.same[simdata.dif$p.same >= 0.05])
  
  print("Different pop, no outlier corrections")
  nrow(simdata.dif[simdata.dif$p.nocorr >= 0.05,])
  mean(simdata.dif$t.nocorr)
  mean(simdata.dif$es.nocorr)
  mean(simdata.dif$t.nocorr[simdata.dif$p.nocorr >= 0.05])
  mean(simdata.dif$es.nocorr[simdata.dif$p.nocorr >= 0.05])
