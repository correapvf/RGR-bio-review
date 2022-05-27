# standardize tables and join them, and generate tables for analysis
# use table_get.R first to generate gbif.csv and obis.csv

library(data.table)

# Retrive tables from HD
GBIF = data.table(read.csv("gbif.csv.gz", stringsAsFactors = F))
OBIS = data.table(read.csv("obis.csv.gz", stringsAsFactors = F))

OBIS$eventDate = as.Date(OBIS$eventDate)
GBIF$eventDate = as.Date(GBIF$eventDate)

OBIS2 = OBIS #backup
GBIF2 = GBIF

# return false where statment is false and true for NA (so we do not exclude where is NA)
comp <- function(x) {
  x[is.na(x)] = TRUE
  return(x)
}

# Filter records of OBIS table
OBIS2 = OBIS2[comp(OBIS$datasetName != "Vulnerable marine ecosystems in the South Pacific Ocean region")] # these are not from the Atlantic. See also waterBody collum

# Change basisOfRecord in OBIS to fossils of those that I sure
OBIS2$basisOfRecord[OBIS2$datasetName == "PANGAEA - Data from the Deep Sea Drilling Project (DSDP)"] = "FossilSpecimen"


# Filter records of GBIF table
GBIF2 = GBIF2[comp(GBIF2$phylum != "Tracheophyta")] # these are terrestrial  plants
GBIF2 = GBIF2[comp(GBIF2$locality != "Bushman's River, Estcourt")] # these are insects and fleshwater animals from south africa


# Correct Verbatim Depth in locality of GBIF
verbatimDepth <- function(x) {
  a = b = c = as.numeric(NA)
  
  temp = substr(x,17,nchar(x)-1)
  index = temp=="surface" | temp=="Surface"
  a[index] = 0
  
  temp2 = strsplit(temp[!index], "- ")
  temp2 = lapply(temp2, as.numeric)
  
  a[!index] = sapply(temp2, mean)
  b[!index] = sapply(temp2, function(x) x[1])
  c[!index] = sapply(temp2, function(x) x[2])
    
  list(a,b,c)
}

GBIF2[substr(locality,1,16) == "[VERBATIM DEPTH:",
      c("depth","minimumDepthInMeters","maximumDepthInMeters") := verbatimDepth(locality)]

# These are also fossils
GBIF2$basisOfRecord[GBIF2$identifiedBy == "Benson, Richard H" & GBIF2$class == "Ostracoda"] = "FossilSpecimen"
GBIF2$basisOfRecord[GBIF2$institutionCode == "NHMUK" & GBIF2$collectionCode == "PAL"] = "FossilSpecimen"

# add georeferenceRemarks to occurrenceRemarks
GBIF2[!is.na(georeferenceRemarks), "occurrenceRemarks" := paste(occurrenceRemarks,georeferenceRemarks,sep="; ")]

# add elevation to depth
GBIF2[is.na(depth) & !is.na(elevation), depth := elevation]

GBIF2$scientificName[GBIF2$scientificName=="Tremoctopus delle"] = "Tremoctopus" #delle is the author
GBIF2$genus[GBIF2$genus=="Tremoctopus delle"] = "Tremoctopus"

GBIF2$eventDate[GBIF2$scientificName=="Mirounga leonina" & is.na(GBIF2$eventDate)] = "1997-01-01" #correct this date based on OBIS

# Aggregate records for each table
# OBIS
#colnames(OBIS2)
cols = c("kingdom","phylum","class","order","family","genus",
  "infraspecificEpithet","scientificNameWithAuthor","taxonRank","basisOfRecord",
  "dynamicProperties","institutionCode","collectionCode","bibliographicCitation",
  "occurrenceRemarks","minimumDepthInMeters","maximumDepthInMeters",
  "fieldNumber","locality","recordedBy","identifiedBy")

OBIS3 = OBIS2[, lapply(.SD, function(x) paste(unique(x[!is.na(x)]),collapse="/")),
          by=.(scientificName,decimalLongitude,decimalLatitude,depth,eventDate),
          .SDcols = cols]

OBIS4 = OBIS2[, .(Noccurrences=.N,individualCount=sum(individualCount)),
              by=.(scientificName,decimalLongitude,decimalLatitude,depth,eventDate)]

OBIS3$samplingProtocol = NA
OBIS3$db = "OBIS"
OBIS3$Noccurrences = OBIS4$Noccurrences
OBIS3$individualCount = OBIS4$individualCount

# GBIF
#colnames(GBIF2)
GBIF3 = GBIF2[, lapply(.SD, function(x) paste(unique(x[!is.na(x)]),collapse="/")),
              by=.(scientificName,decimalLongitude,decimalLatitude,depth,eventDate),
              .SDcols = c(cols,"samplingProtocol")]

GBIF4 = GBIF2[, .(Noccurrences=.N, individualCount=sum(individualCount)),
              by=.(scientificName,decimalLongitude,decimalLatitude,depth,eventDate)]

GBIF3$db = "GBIF"
GBIF3$Noccurrences = GBIF4$Noccurrences
GBIF3$individualCount = GBIF4$individualCount

# Clear some space
OBIS2 = OBIS3
GBIF2 = GBIF3
rm(OBIS3,OBIS4,GBIF3,GBIF4)

# Bind tables
final_table = rbind(OBIS2,GBIF2)

# Search for minimum decimal places to aggregate table - keep for future reference
decimalplaces <- function(x) {
  y = x[!is.na(x)]
  if (length(y) == 0) {
    return(NA)
  }
  temp = vector("numeric", length(y))
  for (i in 1:length(y)) {
    if ((y[i] %% 1) != 0) {
      temp[i] = nchar(strsplit(sub('0+$', '', as.character(y[i])), ".", fixed=TRUE)[[1]][[2]])

    } else {
      temp[i] = 0
    }
  }
  return(temp)
}

temp = final_table[, .(longD=min(decimalplaces(decimalLongitude)),latD=min(decimalplaces(decimalLatitude))), by=.(scientificName)] #by=institutionCode

temp$longD[temp$longD < 2] = 2
temp$latD[temp$latD < 2] = 2
# Do NOT Round longitude and latitude less than 3 decimal places - ~ 100m precision

final_table = merge(final_table, temp, by = c("scientificName"), all.x = T) #by=institutionCode
final_table$longR = round(final_table$decimalLongitude, final_table$longD)
final_table$latR = round(final_table$decimalLatitude, final_table$latD)

# seach for depth, so considerer where depth is NA if same coordinates
depthBind <- function(x) {
  if (any(is.na(x)) & sum(!is.na(x))<=1) {
    return("any")
  } else {
    return("all")
  }
}

depthBind2 <- function(x,y) {
  ifelse(y=="any",NA,x)
}

temp = final_table[,.(depthbinded1=depthBind(depth)),by=.(scientificName,longR,latR)]
final_table = merge(final_table, temp, by = c("scientificName","longR","latR"), all.x = T)
final_table$depthbinded2 = depthBind2(final_table$depth,final_table$depthbinded1)



# search for dates
dateBind <- function(x) {
  if (any(is.na(x)) & sum(!is.na(x))<=1) {
    return("any")
  } else if (any(is.na(x))) {
    return("all")
  } else if (any(month(x)==1) & any(mday(x)==1)) {
    return("year")
  } else if (sum(mday(x)==1)>=2) {
    return("month")
  } else {
    return("all")
  }
}

dateBind2 <- function(x,y) {
  ifelse(y=="any",NA_character_,ifelse(y=="year",as.character(year(x)),ifelse(y=="month",paste0(year(x),'-',month(x)),as.character(x))))
}

temp = final_table[,.(datebinded1=dateBind(eventDate)),by=.(scientificName,longR,latR,depthbinded2)]
final_table = merge(final_table, temp, by = c("scientificName","longR","latR","depthbinded2"), all.x = T)
final_table$datebinded2 = dateBind2(final_table$eventDate,final_table$datebinded1)

# Aggregate records to remove duplicates
max2 <- function(x) {
  ifelse(all(is.na(x)),NA_real_,max(x, na.rm=T))
}

final_table2 = final_table[, lapply(.SD, function(x) paste(unique(x[!is.na(x)]),collapse="/")),
                           by = list(scientificName,longR,latR,depthbinded2,datebinded2),
                           .SDcols = c(cols,"samplingProtocol","db")]

final_table3 = final_table[, .(longitude=first(longR[which(longD==max(longD))]), latitude=first(latR[which(latD==max(latD))]),Noccurrences=max(Noccurrences),
                               individualCount=max(individualCount),depth=max(depth,na.rm=T), eventDate=max(eventDate,na.rm=T)),
                           by = list(scientificName,longR,latR,depthbinded2,datebinded2)]


# Extract GEBCO depth
library("raster")
rasterFile = "G:\\Meu Drive\\ArcGis\\gebco_08a_Subset0.img"
rasterFile = raster(rasterFile)
GEBCO_depth = extract(rasterFile, final_table2[,c("longR","latR")])

# Bind tables and reorder columns
final_table2$longR = final_table3$longitude
final_table2$latR = final_table3$latitude
final_table2$Noccurrences = final_table3$Noccurrences
final_table2$individualCount = final_table3$individualCount
final_table2$depthbinded2 = ifelse(is.infinite(final_table3$depth),NA,final_table3$depth)
final_table2$datebinded2 = as.Date(as.character(final_table3$eventDate))
final_table2$GEBCO_depth = GEBCO_depth
rm(final_table3, temp)

# clean up "/"
for (i in names(final_table2)) {
  if(class(final_table2[[i]]) == "character"){
    final_table2[[i]] = gsub("^/|/$", "", final_table2[[i]])
    final_table2[[i]] = tryCatch(as.numeric(final_table2[[i]]), warning = function(w) final_table2[[i]])
  }
}




# Create outputs and export
# # RGR_fossils
index = final_table2$basisOfRecord == "FossilSpecimen/HumanObservation" | final_table2$basisOfRecord == "FossilSpecimen" | final_table2$kingdom == "incertae sedis"
# RGR_fossils = final_table2[index]
# write.csv(RGR_fossils, "RGR_fossils.csv", row.names = FALSE)
# 
# # RGR_notAnimalia
# RGR_notAnimalia = final_table2[!index & kingdom != "Animalia"]
# write.csv(RGR_notAnimalia, "RGR_notAnimalia.csv", row.names = FALSE)


# RGR_Animalia and get qualityControl file
RGR_Animalia = final_table2[!index & kingdom == "Animalia"]
quality = read.csv("Animalia_qualityControl.csv", stringsAsFactors=FALSE)

RGR_Animalia2 = merge(RGR_Animalia, quality, by = "scientificName")
index_animalia = RGR_Animalia2$remove | RGR_Animalia2$dubious | RGR_Animalia2$marine.birds | RGR_Animalia2$lavae

RGR_Animalia2 = RGR_Animalia2[!index_animalia]
RGR_Animalia2$habitat = ifelse(RGR_Animalia2$benthic, "benthic", "pelagic")

write.csv(RGR_Animalia2[,c(1:31,39)], "RGR_Animalia.csv", row.names = FALSE)


# # RGR_Animalia_points
# x = RGR_Animalia2[, .(Nspecies=length(unique(scientificName)), Noccurrences=sum(Noccurrences),
#                     minDate=min(datebinded2), maxDate=max(datebinded2),
#                     minDepth=min(depthbinded2), maxDepth=max(depthbinded2)),
#                  by= .(longR,latR,habitat)]
# 
# rangeDate = ifelse(is.na(as.character(x$minDate)),paste(x$maxDate),ifelse(x$minDate==x$maxDate, paste(x$maxDate), paste(x$minDate,x$maxDate,sep=" to ")))
# rangeDepth = ifelse(x$minDepth==x$maxDepth, paste(x$maxDepth), paste(x$minDepth,x$maxDepth,sep=" to "))
# 
# RGR_Animalia_points =  cbind(x[,1:5],rangeDate,rangeDepth)
# write.csv(RGR_Animalia_points, "RGR_Animalia_points.csv", row.names = FALSE)
# 
# 
# # NOT_Animalia Points
# RGR_notAnimalia_points = RGR_notAnimalia[, .(Nspecies=length(unique(scientificName)), Noccurrences=sum(Noccurrences)), by= .(longR,latR)]
# write.csv(RGR_notAnimalia_points, "RGR_notAnimalia_points.csv", row.names = FALSE)


# # Total  number os species and records per phylum and class
RGR_Animalia2$individualCount2 = ifelse(is.na(RGR_Animalia2$individualCount) | RGR_Animalia2$individualCount==0,RGR_Animalia2$Noccurrences,RGR_Animalia2$individualCount)

x = RGR_Animalia2[, .(Nspecies=length(unique(scientificName)), Nrecords=.N, Noccurrences=sum(Noccurrences), Ninds=sum(individualCount2)), by= .(phylum, class, habitat)]
setorder(x, habitat, phylum, class)
# write.csv(x, "Nrecords_per_phylum.csv", row.names = FALSE)

