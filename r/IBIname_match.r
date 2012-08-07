IBIname_match <- function(data, DistinctCode=F){
  colnames(data)[which(colnames(data) == "FinalID")] <- "Taxa"
  colnames(data)[which(colnames(data) == "BAResult")] <- "Result"
  data <- data[which(!is.na(data$Result)), ]
  load(system.file("data", "ibi.RData", package ="ibiscore"))
  require(plyr)
  ibi <- idata.frame(ibi)
  load(system.file("data", "taxonomy.RData", package ="ibiscore"))
  taxonomy <- idata.frame(taxonomy)
  ###Convert FinalID to SAFIT1###
  data$SAFIT <- rep(NA, length(data$Taxa))
  data$SAFIT <- ibi$SAFIT1[match(data$Taxa, ibi$FinalID)]
  data$SAFIT[is.na(data$SAFIT)] <- "Missing"
  ###Fix extra spaces
  data$Taxa <- as.character(data$Taxa)
  complex <- grep("Group", data$Taxa[data$SAFIT == "Missing"])
  extraspace <- data$Taxa[intersect(which(data$SAFIT == "Missing"), which(!(data$Taxa %in% (data$Taxa[data$SAFIT == "Missing"][complex]))))]
  data$Taxa[intersect(which(data$SAFIT == "Missing"), which(!(data$Taxa %in% (data$Taxa[data$SAFIT == "Missing"][complex]))))] <- 
    gsub("(\\w+)\\s+$", "\\1", extraspace)
  data$SAFIT <- ibi$SAFIT1[match(data$Taxa, ibi$FinalID)]
  data$SAFIT[is.na(data$SAFIT)] <- "Missing"
  ###Fix extra caps###
  cap1 <- gsub("(^\\w)[[:alnum:][:space:]]+", "\\1", data$Taxa[data$SAFIT == "Missing"])
  cap2 <- gsub("(^\\w)(\\w+)", "\\2", data$Taxa[data$SAFIT == "Missing"])
  data$Taxa[data$SAFIT == "Missing"] <- paste0(cap1, tolower(cap2))
  data$SAFIT <- ibi$SAFIT1[match(data$Taxa, ibi$FinalID)]
  data$SAFIT[is.na(data$SAFIT)] <- "Missing"
  ###Label FinalIDs that match SAFIT1 as distinct###
  data$distinct <- rep(NA, length(data$Taxa))
  #data$distinct[which(data$Taxa == data$SAFIT)] <- "Distinct"
  data$distinct[data$SAFIT == "Missing"] <- "Missing"
  ###Determine whether the rest of the FinalIDs are distinct###
  todetermine <- as.character(data$SAFIT[which(is.na(data$distinct))])
  todetermine <- as.data.frame(cbind(as.character(data$SampleID[which(is.na(data$distinct))]), todetermine))
  todeterminebystation <- tapply(as.character(todetermine$todetermine), todetermine$V1, list)
  for(j in 1:length(todeterminebystation)){
    data$distinct[intersect(which(data$SampleID %in% names(todeterminebystation[j])), which(data$SAFIT %in% todeterminebystation[[j]]))] <- sapply(1:length(todeterminebystation[[j]]), function(i, determine, tax){
      index <- which(tax$FinalID == todeterminebystation[[j]][i])
      level <- as.numeric(tax[index, "TaxonomicLevelCode"])
      levelname <- tax[index, "TaxonomicLevelName"]
      if(level >= 60){"Distinct"} else{
        criteron1 <- which(taxonomy$FinalID %in% todeterminebystation[[j]][which(!(unlist(todeterminebystation[j]) %in% unlist(todeterminebystation[[j]][i])))])
        criteron2 <- which(tax[criteron1, levelname] == tax[index, levelname])
        criteron3 <- which(tax[criteron2, "TaxonomicLevelCode"] > tax[index, "TaxonomicLevelCode"])
        if(length(criteron3) > 0){"Non-distinct "} else
        {"Distinct"}}}, determine=todeterminebystation[[j]], tax=taxonomy)}
  if(DistinctCode == T){
    data[data$distinct == "Not Distinct" & data$DistinctCode == "Yes", "distinct"] <- "Distinct"
  }
  return(data)
}