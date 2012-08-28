###Read in command line arguments###
args <- commandArgs(TRUE)
csv_output <- paste("/var/www/upload/files", args[4], sep = "/")
xml_output <- paste("/var/www/upload/files", args[3], sep = "/")
location_input <- paste("/var/www/upload/files", args[2], sep = "/")
data_input <- paste("/var/www/upload/files", args[1], sep = "/")

loc <- read.csv(location_input)
dat <- read.csv(data_input)

###Location determination###
N_or_S <- function(points){ ###"points must be a data frame with columns StationCode, Latitude, and Longitude
  load("Data/map.RData")
  load("Data/map2.RData")
  require(rgdal)
  require(plyr)
  coordinates(points) <- ~Longitude + Latitude
  points@proj4string <- map@proj4string
  socal_zone <- points[map,]
  socal_zone$X <- rep("socal", length(socal_zone$X))
  norcal_zone <- points[map2,]
  norcal_zone$X <- rep("norcal", length(norcal_zone$X))
  zone <- rbind.fill(socal_zone, norcal_zone)
  zone
}

zone <- N_or_S(loc)
socal <- zone$StationCode[which(zone$X == "socal")]
norcal <- zone$StationCode[which(zone$X == "norcal")]

###Compute SoCalIBI score###
if(length(socal) != 0){
source("SoCal_IBI.r")
sloc <- loc[which(loc$StationCode %in% socal), ]
sdat <- dat[which(dat$StationCode %in% socal),]
sresults <- SoCal_IBI(sloc, sdat)
}

###Compute NorCalIBI score###
if(length(norcal) != 0){
source(NorCal_IBI.r"")
nloc <- loc[which(loc$StationCode %in% norcal), ]
ndat <- dat[which(dat$StationCode %in% norcal),]
nresults <- NorCal_IBI(nloc, ndat)
}

###Combine results###
if(is.null(sresults)){
  results <- nresults
}
else(is.null(nresults)){
  results <- sresults
}
else{
  library(plyr)
  results <- rbind.fill(sresults, nresults)
}

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