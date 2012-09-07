###Read in command line arguments###
args <- commandArgs(TRUE)
csv_output <- paste("/var/www/upload/files", args[4], sep = "/")
xml_output <- paste("/var/www/upload/files", args[3], sep = "/")
location_input <- paste("/var/www/upload/files", args[2], sep = "/")
data_input <- paste("/var/www/upload/files", args[1], sep = "/")

###IBI Calculator
SoCal_IBI <- function(locationinfo, data, DistinctCode=F, Grid=F, SampleDate=F, FieldReplicate=F){
  ###pull in other scripts###
  source("IBIlocation.r")
  #source("IBImap.r")
  source("IBIname_match.r")
  options(warn = -1)
  load("../Data/ibiv3.RData")
  ibi <- ibiv3
  starttime <- proc.time()
  data <- IBIname_match(data=data, DistinctCode=DistinctCode)
  colnames(data)[which(colnames(data) == "FinalID")] <- "Taxa"
  colnames(data)[which(colnames(data) == "BAResult")] <- "Result"
  ###Calculate total count###
  total_count <- tapply(data$Result, data$SampleID, sum)
  ###Create sample count flag###
  sample_count_flag <- rep(NA, length(total_count))
  names(sample_count_flag) <- names(total_count)
  sample_count_flag[(which(total_count>=500))] <- "Adequate"
  sample_count_flag[(which(total_count < 500 & total_count >= 450))] <- "Within specifications"
  sample_count_flag[(which(total_count < 450))] <- "Inadequate"
  ###Subsample down to 500###
  datalength <- length(data)
  rarifydown <- function(data){unlist(sapply(unique(data$SampleID), function(sample){
    v <- data[data$SampleID==sample, "Result"]
    
    if(sum(v)>=500){rrarefy(v, 500)} else
    {v}
  }
  )
  )
  }
  
  require(doParallel)
  require(vegan)
  registerDoParallel()
  rarificationresult <- foreach(i=1:20, .combine=cbind, .packages="vegan") %dopar% {
    rarifydown(data)
  }
  data <- cbind(data, rarificationresult)
  colnames(data)[(datalength + 1):(datalength + 20)]<- paste("Replicate", 1:20)
  
  ###Metrics set up###
  metrics <- as.data.frame(matrix(NA, nrow = length(unique(data$SampleID)), ncol = 140))
  samplenames <- names(tapply(data$SAFIT, data$SampleID, length))
  data$SAFIT <- as.character(data$SAFIT)
  ###Merge revelant ibi data into data table###
  data$MaxTol <- ibi$MaxTol[match(data$Taxa, ibi$FinalID)]
  data$MaxTol <- as.numeric(data$MaxTol)
  data$Class <- ibi$Class[match(data$Taxa, ibi$FinalID)]
  data$Order <- ibi$Order[match(data$Taxa, ibi$FinalID)]
  data$FunctionalFeedingGroup <- ibi$FunctionalFeedingGroup[match(data$Taxa, ibi$FinalID)]
  write.csv(data, file="midstream_test.csv")
  ###Number of Coleoptera taxa###
  for(i in 1:20){
    metrics[[i]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                          function(d)length(unique(d$SAFIT[d$Order=="Coleoptera"])))[, 2]			   
    
    ###Numer of EPT taxa
    metrics[[i+20]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$SAFIT[d$Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")])))[, 2]
    
    ###Number of predator taxa###
    metrics[[i+40]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$SAFIT[which(d$FunctionalFeedingGroup == "P")])))[, 2]
    
    ###Percent Non-Insect taxa###
    metrics[[i+60]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d)100*length(unique(d$SAFIT[which(d$Class != "Insecta")]))/length(unique(d$SAFIT)))[, 2]
    
    ###Percent tolerant taxa###
    metrics[[i+80]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                             function(d){
                               100*length(unique(d$SAFIT[d$MaxTol >= 8]))/length(unique(d$SAFIT))
                             })[, 2]
    
    ###Percent Intolerant###
    metrics[[i+100]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                              function(d){
                                100*sum(d$Result[d$MaxTol <= 2], na.rm=T)/sum(d$Result, na.rm=T)
                              })[, 2]
    
    ###Percent CF + CG individuals###
    metrics[[i+120]] <- ddply(data[data$distinct == "Distinct" & data[[datalength + i]]>0, ], "SampleID",
                              function(d){
                                100*sum(d$Result[which(d$FunctionalFeedingGroup %in% c("CF", "CG"))])/sum(d$Result)
                              })[, 2]	
  }
  
  for(i in 101:120){
    metrics[which(is.na(metrics[[i]])), i] <- 0
  }
  
  ###Identify Ecoregion###
  Ecoregion <- IBIlocation(locationinfo)
  data$Ecoregion <- Ecoregion[match(data$StationCode, names(Ecoregion))]
  stationlocation <- unique(data[, c("SampleID", "Ecoregion")])
  chap <- which(stationlocation[[2]]=="Chaparral")
  mount <- which(stationlocation[[2]]=="Mountains")
  missing <- which(is.na(stationlocation[[2]]))
  ###Convert metrics to scores###
  scores <- metrics
  ###Coleoptera scores###
  
  for(i in 1:20){
    scores[scores[[i]]>=6, i] <- 10
    scores[scores[[i]]==5, i] <- 8
    scores[scores[[i]]==4, i] <- 7
    scores[scores[[i]]==3, i] <- 5
    scores[scores[[i]]==2, i] <- 4  
    scores[missing, i] <- NA
    
    ###EPT scores###
    
    
    scores[intersect(chap, which(scores[[i+20]]<=1)), i+20] <- 0
    scores[intersect(chap, which(scores[[i+20]]>1 & scores[[i+20]]<=3)), i+20] <- 1
    scores[intersect(chap, which(scores[[i+20]]==4)), i+20] <- 2
    scores[intersect(chap, which(scores[[i+20]]>=5 & scores[[i+20]]<=6)), i+20] <- 3
    scores[intersect(chap, which(scores[[i+20]]>=7 & scores[[i+20]]<=8)), i+20] <- 4
    scores[intersect(chap, which(scores[[i+20]]>=9 & scores[[i+20]]<=10)), i+20] <- 5
    scores[intersect(chap, which(scores[[i+20]]>=11 & scores[[i+20]]<=12)), i+20] <- 6
    scores[intersect(chap, which(scores[[i+20]]>=13 & scores[[i+20]]<=14)), i+20] <- 7
    scores[intersect(chap, which(scores[[i+20]]==15)), i+20] <- 8
    scores[intersect(chap, which(scores[[i+20]]>=16 & scores[[i+20]]<=17)), i+20] <- 9
    scores[intersect(chap, which(scores[[i+20]]>=18)), i+20] <- 10
    
    scores[intersect(mount, which(scores[[i+20]]<=4)), 1+20] <- 0
    scores[intersect(mount, which(scores[[i+20]]>=5 & scores[[i+20]]<=6)), i+20] <- 1
    scores[intersect(mount, which(scores[[i+20]]==7)), i+20] <- 2
    scores[intersect(mount, which(scores[[i+20]]>=8 & scores[[i+20]]<=9)), i+20] <- 3
    scores[intersect(mount, which(scores[[i+20]]==10)), i+20] <- 4
    scores[intersect(mount, which(scores[[i+20]]>=11 & scores[[i+20]]<=12)), i+20] <- 5
    scores[intersect(mount, which(scores[[i+20]]>=13)), i+20] <- 6
    scores[intersect(mount, which(scores[[i+20]]>=14 & scores[[i+20]]<=15)), i+20] <- 7
    scores[intersect(mount, which(scores[[i+20]]==16)), i+20] <- 8
    scores[intersect(mount, which(scores[[i+20]]>=17 & scores[[i+20]]<=18)), i+20] <- 9
    scores[intersect(mount, which(scores[[i+20]]>=19)), i+20] <- 10
    
    scores[missing, i+20] <- NA
    
    ###Predators###
    
    scores[scores[[i+40]]<=3, i+40] <- 0
    scores[scores[[i+40]]==4, i+40] <- 1
    scores[scores[[i+40]]==5, i+40] <- 2
    scores[scores[[i+40]]==6, i+40] <- 3
    scores[scores[[i+40]]==7, i+40] <- 4
    scores[scores[[i+40]]==8, i+40] <- 5
    scores[scores[[i+40]]==9, i+40] <- 6
    scores[scores[[i+40]]==10, i+40] <- 7
    scores[scores[[i+40]]==11, i+40] <- 8
    scores[scores[[i+40]]==12, i+40] <- 9
    scores[scores[[i+40]]>=13, i+40] <- 10 
    
    scores[missing, i+40] <- NA
    
    ###Non-insect###
    
    scores[scores[[i+60]]>=9 & scores[[i+60]]<13, i+60] <- 9
    scores[scores[[i+60]]<=8, i+60] <- 10
    scores[scores[[i+60]]>= 47, i+60] <- 0
    scores[scores[[i+60]]>= 43 & scores[[i+60]]<47, i+60] <- 1
    scores[scores[[i+60]]>= 39 & scores[[i+60]]<43, i+60] <- 2
    scores[scores[[i+60]]>=35 & scores[[i+60]]<40, i+60] <- 3
    scores[scores[[i+60]]>=30 & scores[[i+60]]<35, i+60] <- 4
    scores[scores[[i+60]]>=26 & scores[[i+60]]<30, i+60] <- 5
    scores[scores[[i+60]]>=22 & scores[[i+60]]<26, i+60] <- 6
    scores[scores[[i+60]]>=18 & scores[[i+60]]<22, i+60] <- 7
    scores[scores[[i+60]]>=13 & scores[[i+60]]<18, i+60] <- 8
    
    scores[missing, i+60] <- NA
    
    ###Percent tolerant taxa###
    
    scores[scores[[i+80]]>=5 & scores[[i+80]]<9, i+80] <- 200
    scores[scores[[i+80]]>=9 & scores[[i+80]]<13, i+80] <- 8
    scores[scores[[i+80]]==200, i+80] <- 9
    scores[scores[[i+80]]<=4, i+80] <- 10
    scores[scores[[i+80]]>=13 & scores[[i+80]]<17, i+80] <- 7
    scores[scores[[i+80]]>= 38, i+80] <- 0
    scores[scores[[i+80]]>= 34 & scores[[i+80]]<38, i+80] <- 1
    scores[scores[[i+80]]>= 30 & scores[[i+80]]<34, i+80] <- 2
    scores[scores[[i+80]]>=26 & scores[[i+80]]<30, i+80] <- 3
    scores[scores[[i+80]]>=23 & scores[[i+80]]<26, i+80] <- 4
    scores[scores[[i+80]]>=20 & scores[[i+80]]<23, i+80] <- 5
    scores[scores[[i+80]]>=17 & scores[[i+80]]<19, i+80] <- 6
    
    
    scores[missing, i+80] <- NA
    
    
    ###Percent Intolerant###
    
    
    scores[intersect(chap, which(scores[[i+100]]==0)), i+100] <- 0
    scores[intersect(chap, which(scores[[i+100]]>=1 & scores[[i+100]]<4)), i+100] <- 1
    scores[intersect(chap, which(scores[[i+100]]>=4 & scores[[i+100]]<7)), i+100] <- 2
    scores[intersect(chap, which(scores[[i+100]]>=7 & scores[[i+100]]<10)), i+100] <- 3
    scores[intersect(chap, which(scores[[i+100]]>=10 & scores[[i+100]]<12)), i+100] <- 4
    scores[intersect(chap, which(scores[[i+100]]>=13 & scores[[i+100]]<15)), i+100] <- 5
    scores[intersect(chap, which(scores[[i+100]]>=16 & scores[[i+100]]<18)), i+100] <- 6
    scores[intersect(chap, which(scores[[i+100]]>=19 & scores[[i+100]]<20)), i+100] <- 7
    scores[intersect(chap, which(scores[[i+100]]>=21 & scores[[i+100]]<22)), i+100] <- 8
    scores[intersect(chap, which(scores[[i+100]]>=23 & scores[[i+100]]<24)), i+100] <- 9
    scores[intersect(chap, which(scores[[i+100]]>=25)), i+100] <- 10
    
    scores[intersect(mount, which(scores[[i+100]]<=1)), 1+100] <- 0
    scores[intersect(mount, which(scores[[i+100]]>=2 & scores[[i+100]]<6)), i+100] <- 1
    scores[intersect(mount, which(scores[[i+100]]>=6 & scores[[i+100]]<10)), i+100] <- 2
    scores[intersect(mount, which(scores[[i+100]]>=10 & scores[[i+100]]<14)), i+100] <- 3
    scores[intersect(mount, which(scores[[i+100]]>=14 & scores[[i+100]]<19)), i+100] <- 4
    scores[intersect(mount, which(scores[[i+100]]>=19 & scores[[i+100]]<23)), i+100] <- 5
    scores[intersect(mount, which(scores[[i+100]]>=23 & scores[[i+100]]<26)), i+100] <- 6
    scores[intersect(mount, which(scores[[i+100]]>=27 & scores[[i+100]]<32)), i+100] <- 7
    scores[intersect(mount, which(scores[[i+100]]>=32 & scores[[i+100]]<37)), i+100] <- 8
    scores[intersect(mount, which(scores[[i+100]]>=37 & scores[[i+100]]<42)), i+100] <- 9
    scores[intersect(mount, which(scores[[i+100]]>=42)), i+100] <- 10
    
    scores[missing, i+100] <- NA
    
    
    ###CF+CG###
    
    scores[intersect(chap, which(scores[[i+120]]<60)), i+120] <- 10
    scores[intersect(chap, which(scores[[i+120]]>=97)), i+120] <- 0
    scores[intersect(chap, which(scores[[i+120]]>=93 & scores[[i+120]]<97)), i+120] <- 1
    scores[intersect(chap, which(scores[[i+120]]>=89 & scores[[i+120]]<93)), i+120] <- 2
    scores[intersect(chap, which(scores[[i+120]]>=85 & scores[[i+120]]<89)), i+120] <- 3
    scores[intersect(chap, which(scores[[i+120]]>=81 & scores[[i+120]]<85)), i+120] <- 4
    scores[intersect(chap, which(scores[[i+120]]>=76 & scores[[i+120]]<81)), i+120] <- 5
    scores[intersect(chap, which(scores[[i+120]]>=72 & scores[[i+120]]<76)), i+120] <- 6
    scores[intersect(chap, which(scores[[i+120]]>=68 & scores[[i+120]]<72)), i+120] <- 7
    scores[intersect(chap, which(scores[[i+120]]>=64 & scores[[i+120]]<68)), i+120] <- 8
    scores[intersect(chap, which(scores[[i+120]]>=60 & scores[[i+120]]<64)), i+120] <- 9
    
    scores[intersect(mount, which(scores[[i+120]]<=39)), i+120] <- 10
    scores[intersect(mount, which(scores[[i+120]]>=95)), i+120] <- 0
    scores[intersect(mount, which(scores[[i+120]]>=89 & scores[[i+120]]<95)), i+120] <- 1
    scores[intersect(mount, which(scores[[i+120]]>=83 & scores[[i+120]]<89)), i+120] <- 2
    scores[intersect(mount, which(scores[[i+120]]>=77 & scores[[i+120]]<82)), i+120] <- 3
    scores[intersect(mount, which(scores[[i+120]]>=71 & scores[[i+120]]<77)), i+120] <- 4
    scores[intersect(mount, which(scores[[i+120]]>=65 & scores[[i+120]]<71)), i+120] <- 5
    scores[intersect(mount, which(scores[[i+120]]>=59 & scores[[i+120]]<65)), i+120] <- 6
    scores[intersect(mount, which(scores[[i+120]]>=53 & scores[[i+120]]<59)), i+120] <- 7
    scores[intersect(mount, which(scores[[i+120]]>=47 & scores[[i+120]]<53)), i+120] <- 8
    scores[intersect(mount, which(scores[[i+120]]>=40 & scores[[i+120]]<47)), i+120] <- 9
    
    scores[missing, i+120] <- NA
    
  }
  ###SCIBI###
  for(i in 0:19){
    scores[[i+141]] <- sapply(1:length(unique(data$SampleID)), function(j)(10/7)*(sum(c(scores[j, 1+i], scores[j, 20+i], 
                                                                                        scores[j, 40+i], scores[j, 60+i], scores[j, 80+i], scores[j, 100+i], scores[j, 120+i]))))
  }
  ###Calculate means for metrics and scores###
  means <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 15))
  for(i in 1:7){
    means[[i]] <- apply(metrics[, (((i-1)*20)+1):(20*i)], 1, function(d)sum(d)/20)
  }  
  for(i in 1:8){
    means[[i+7]] <- apply(scores[, (((i-1)*20)+1):(20*i)], 1, mean)
  }
  ###Construct output frame###
  results <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 21))
  results[[1]] <- unique(data[, c("StationCode", "SampleID")])$StationCode
  results[[2]] <- unique(data$SampleID)
  results <- results[match(samplenames, as.character(results[[2]])),]
  results[[3]] <- unique(data[, c("SampleID", "Ecoregion")])$Ecoregion
  results[[4]] <- total_count[!is.na(total_count)]
  results[[5]] <- sample_count_flag[!is.na(sample_count_flag)]
  results[[6]] <- rep(20, times=length(unique(data$SampleID)))
  results[[6]][which(results[[4]] < 500)] <- 1
  results[, 7:21] <- means
  results[[21]] <- round(results[[21]], digits=2)
  colnames(results) <- c("StationCode", "SampleID", "Ecoregion", "Total Count", "Count Flag", "Number of Iteration", 
                         "Number of Coleoptera Taxa", "Number of EPT Taxa", "Number of Predator Taxa", 
                         "Percent Non-Insect Taxa", "Percent Tolerant Taxa", "Percent Intolerant", 
                         "Percent CF + CG Individuals", "Coleoptera Score", "EPT Score", "Predator Taxa Score",
                         "Non-Insect Taxa Score", "Tolerant Taxa Score", "Intolerant Score",
                         "CF + CG Score", "SCIBI")
  ###Representativeness flag###
  if(Grid==T){
    data$Representativeness_Flag <- rep(NA, length(data$StationCode))
    data[which(is.na(data$TotalGrids)), "TotalGrids"] <- "Missing"
    flag <- sapply(1:length(data$StationCode), function(i)
      if(data$TotalGrids[i] == "Missing"){NA}else 
        if(as.numeric(data$TotalGrids[i]) >= 3 | data$GridsAnalyzed[i] > 2 | (data$GridsVolumeAnalyzed[i]/
          as.numeric(data$TotalGrids[i])) >= .25){"Representative"}else{"Potentially nonrepresentative"})
    data$Representativeness_Flag <- flag[!is.na(flag)]
  }
  ###Optional fields###
  options(warn = -1)
  if(Grid==T){
    results$"Representativeness Flag" <- data$"Representativeness_Flag"[match(results$SampleID, data$SampleID)]
  }
  if(FieldReplicate==T){
    results$FieldReplicate <- data$FieldReplicate[match(results$SampleID, data$SampleID)]
  }
  if(SampleDate==T){
    results$SampleDate <- data$SampleDate[match(results$SampleID, data$SampleID)]
  }
  extrastuff <- sum(c(Grid, FieldReplicate, SampleDate))
  if(extrastuff>0){
    results <- results[, c(1:5, 22:(21+extrastuff), 6:21)]  
  }
  ###Merge Lat/Long###
  results$Longitude <- locationinfo$Longitude[match(results$StationCode, locationinfo$StationCode)]
  results$Latitude <- locationinfo$Latitude[match(results$StationCode, locationinfo$StationCode)]
  ###Return results###
  return(results)
}

###Compute IBI score###
loc <- read.csv(location_input)
dat <- read.csv(data_input)
results <- SoCal_IBI(loc, dat, DistinctCode=T)

###Write results to csv###
write.csv(results, file=csv_output)

###XML writer###
writeXML <- function(results, output_file){
  sink(file=output_file)
  cat('<?xml version="1.0"?>', "\n")
  cat('<doc>', "\n")
  cat(' <markers>', "\n")
  for(i in 1:nrow(results)){
    cat(paste("  <marker stationid=", '"', results$StationCode[i], '"', " sampleid=",
              '"', results$SampleID[i], '"', " latitude=", '"', 
              results$Latitude[i], '"', " longitude=", '"', results$Longitude[i], '"', " SCIBI=", '"', 
              results$SCIBI[i], '"'," coleoptera_taxa=", '"', results$"Number of Coleoptera Taxa"[i], '"',
              " ept_taxa=", '"', results$"Number of EPT Taxa"[i], '"', " predator_taxa=",
              '"', results$"Number of Predator Taxa"[i], '"', " percent_non_insect=", '"',
              results$"Percent Non-Insect Taxa"[i], '"', " percent_tolerant_taxa=", '"',
              results$"Percent Tolerant Taxa"[i], '"', " percent_intolerant_ind=",  '"',
              results$"Percent Intolerant"[i], '"', " cf_cg=", '"', 
              results$"Percent CF + CG Individuals"[i], '"',  '/>', sep=""), "\n")
  }
  cat(' </markers>', "\n")
  cat('</doc>')
  sink()
}

###Write results to XML###
writeXML(results, xml_output)
