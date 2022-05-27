# Get tables from GBIF and OBIS and make a pre-clean up

#library(devtools)
#install_github("ropensci/rgbif")
#install_github("iobis/robis")
library(rgbif)
library(robis)
library(data.table)
library(taxize)

newdata = FALSE # If TRUE, retrieve new data from OBIS and GBIF. Setting this to TRUE will likely break this script, as records will differs!
# In this case, you must check line by line of the code and clean the new data accordingly.

newtaxaMatch = !file.exists("taxaMatch.RData") # file to store taxa data from WORMS. Can be reused, as it takes a long time to process.
wkt = "POLYGON((-39.5 -27.5,-39.5 -35.5,-27.5 -35.5,-27.5 -27.5,-39.5 -27.5))"
check_wkt(wkt)
options(gbif_user="",gbif_pwd="",gbif_email="") # insert user data here
ranks = c("Kingdom","Phylum","Class","Order","Family","Genus","Species") # keep basic taxon ranks


### GBIF
if (newdata) {
  d = occ_download(pred("geometry", wkt))
  # OR occ_download_list()
  
  # Get information about the download
  occ_download_meta(d)
  occ_download_get(d) # only run after finished
}

# Load taxaMatch
if (!newtaxaMatch) {
  load("taxaMatch.RData")
}

x = unz("0101012-160910150852091.zip", "occurrence.txt") # Change name of the file
gbif = read.delim(x, encoding="UTF-8", stringsAsFactors=FALSE, na.strings=c(""," ","NULL"))

# check occurreces "incertae sedis"
gbif = data.table(gbif)
index = gbif$scientificName == "incertae sedis" & !is.na(gbif$higherClassification)
index[is.na(index)] = F
x = gbif[index]
unique(x$higherClassification) # Check possible values
unique(x$taxonRank) # Check possible values

# Now correct some of them
gbif$kingdom[index] = "Foraminifera"
gbif$scientificName[index] = "Foraminifera"

# and without any taxonomic information
x = gbif[is.na(gbif$scientificName) & !is.na(gbif$higherClassification)]
x # these have none information on higherClassification, so ignore them

# now clear the data
gbif2 = gbif[!is.na(gbif$scientificName)]
gbif2 = gbif2[!scientificName == "incertae sedis"]
# total 378 (4744) occurrences deleted without taxonomic information

# check which collums have no information
x = colnames(gbif2)
index = vector("logical", length(x))
for (i in seq_along(x)) {
  if (all(is.na(gbif2[,get(x[i])]))) {
    index[i] = T
  }
}
x[index] # collums with no information

# and remove them
gbif3 = gbif2[,.SD,.SDcols = x[!index]]

# Correct some utf-8 to ISO conversion problems
gbif3$locality[substr(gbif3$locality, 12, 24) == "DO RIO GRANDE"] = "Elevação do Rio Grande" # to correct ELEVA????O DO RIO GRANDE
gbif3$locality[gbif3$locality == "ElevaÃ§Ã£o de Rio Grande"] = "Elevação do Rio Grande"

# Standardize fields with OBIS
gbif3$taxonRank = tolower(gbif3$taxonRank)

gbif3$basisOfRecord[gbif3$basisOfRecord == "HUMAN_OBSERVATION"] = "HumanObservation"
gbif3$basisOfRecord[gbif3$basisOfRecord == "OBSERVATION"] = "Observation"
gbif3$basisOfRecord[gbif3$basisOfRecord == "PRESERVED_SPECIMEN"] = "PreservedSpecimen"
gbif3$basisOfRecord[gbif3$basisOfRecord == "UNKNOWN"] = "Unknown"
gbif3$basisOfRecord[gbif3$basisOfRecord == "FOSSIL_SPECIMEN"] = "FossilSpecimen"

# i found this error in the GBIF table - some occurrences with kindom animalia without phylum are incorrect
if (newtaxaMatch) {
  x = unique(gbif3$species[gbif3$kingdom == "Animalia" & is.na(gbif3$phylum)])
  x = x[!is.na(x)]
  x
  
  taxaMatch = classification(x, db="worms", return_id=F)
  taxaMatch2 = rbind(taxaMatch)
  taxaMatch2 = reshape(taxaMatch2[1:3], idvar = "query", timevar = "rank", direction = "wide")
  
  taxaMatch2 = taxaMatch2[c("query",paste0("name.",ranks))]
  
  taxaMatch # Add species Ericsonia fenestrata manually - deleted from http://www.gbif.org/species/6878998
  taxaMatch2 = rbind(taxaMatch2, c("Ericsonia fenestrata","Chromista",rep(NA, 6)))
  colnames(taxaMatch2)[2:8] = ranks
}

#update table
gbif4 = gbif3
gbif4[taxaMatch2, tolower(ranks)[-7] := 
              .(Kingdom,Phylum,Class,Order,Family,Genus),
              on = c(species="query")]

if (newtaxaMatch) {
  # Also, some with no classification
  x2 = unique(gbif3$scientificName[gbif3$taxonRank == "species" & is.na(gbif3$species)])
  x2 # since there is only one, i'll do it manually
  x2 = classification("Pleuroskelidion unda", db="worms", return_id=F)
  x2 = cbind(x2)[1,]
}

#update table
gbif4[gbif4$scientificName=="Pleuroskelidion unda Patterson",
      tolower(ranks) := x2[,tolower(ranks)]]


# there is a Chauliodus fish and a Chauliodus aves
gbif4$genus[gbif4$class == "Aves" & gbif4$genus == "Chauliodus"] = "Chauliodus (aves)"

# Calculate a column for GBIF table with the lowest taxon name, without author
gbif4$scientificNameWithAuthor = gbif4$scientificName

# Calculate a column for GBIF table with the lowest taxon name
temp = gbif4$taxonRank
temp[temp == "subspecies"] = "species"
x = vector("character", length(temp))
for (i in 1:length(temp)) {
    x[i] = gbif4[[i,temp[i]]]
}
sum(is.na(x)) # this value must be 0
gbif4$scientificName = x

# save table
write.csv(gbif4, gzfile("gbif.csv.gz"), row.names = FALSE)



### OBIS
if (newdata) {
  obis = occurrence(geometry = wkt)
  write.csv(obis, gzfile("obis_bak.csv.gz"), row.names = FALSE) # save original table
} else {
  obis = read.csv("obis_bak.csv.gz", stringsAsFactors = F)
}

if (newtaxaMatch) {
  # Get classification - the OBIS table lacks kingdom, taxonRank and some species don't have classification
  x = unique(obis$aphiaID) # Get from the ones with aphiaID
  x = x[!is.na(x)]
  taxaMatch3 = classification(x, db = "worms", return_id = FALSE)
  
  # Get the rank
  trank = vector("character", length(taxaMatch3))
  for (i in 1:length(taxaMatch3)) {
    trank[i] = tail(taxaMatch3[[i]][[2]], 1)
  }
  
  # I'll just keep the basics ranks
  unique(trank) # run to see which ones to remove
  trank[trank == "Infraorder"] = "Order"
  trank[trank == "Infraphylum"] = "Phylum"
  trank[trank == "Variety"] = "Species"
  
  # Transform data in a nice table
  taxaMatch4 = rbind(taxaMatch3)
  taxaMatch4 = reshape(taxaMatch4[1:3], idvar = "query", timevar = "rank", direction = "wide")
  taxaMatch4 = taxaMatch4[c("query", paste0("name.", ranks))]
  taxaMatch4$trank = trank # Add taxonRank
  taxaMatch4$query = as.integer(taxaMatch4$query) # id must be integral
}

# Update the OBIS table
obis2 = data.table(obis)

# Update only kingdom and taxonRank, as others ranks are already with correct values
obis2[taxaMatch4, c("kingdom","taxonRank") := .(name.Kingdom,trank), on = c(aphiaID="query")]

if (newtaxaMatch) {
  ## Now to the species without aphiaID
  x=unique(obis2$scientificName[is.na(obis2$aphiaID)])
  
  taxaMatch5 = classification(x, db="gbif", return_id=F) #, ask=F) # gbif is better to find these records
  
  # Same process as above
  trank = vector("character", length(taxaMatch5))
  for (i in 1:length(taxaMatch5)){
    if (all(is.na(taxaMatch5[[i]]))) {next}
    trank[i] = tail(taxaMatch5[[i]][[2]], 1)
  }
  unique(trank) # No need to modified anything
  
  trank = trank[trank!=""] # remove blank vectors
  
  taxaMatch6 = rbind(taxaMatch5)
  taxaMatch6 = reshape(taxaMatch6[,1:3], idvar = "query", timevar = "rank", direction = "wide")
  taxaMatch6$trank = trank # Add taxonRank
  colnames(taxaMatch6)[2:8] = ranks # change colnames for easier syntax
  
  # Also check for "incertae sedis"
  taxaMatch6 = data.table(taxaMatch6)
  taxaMatch6[taxaMatch6$Kingdom == "incertae sedis"]
  taxaMatch6 = taxaMatch6[!taxaMatch6$Kingdom == "incertae sedis"] # and remove
}

# Update
obis3 = obis2
obis3[taxaMatch6, c(tolower(ranks),"taxonRank") := 
              .(Kingdom,Phylum,Class,Order,Family,Genus,Species,trank),
              on = c(scientificName="query")]


### Now check occurrences still without rank and kingdom
index = is.na(obis3$kingdom)
if (newtaxaMatch) {
  x = obis3[index]
  View(x)
  
  # Their classfication is more complete, so lets check which kingdom they are
  x = unique(obis3$scientificName[index])
  x # All taxonRank are species
  
  x = unique(gsub(" .*$", "", x)) # check genus only
  
  x3 = classification(x, db="gbif", return_id=F)
  x3 # all are from kingdom Chromista
}

# Update table based on above
obis3$kingdom[index] = "Chromista"
obis3$taxonRank[index] = "species"

#save data
if (newtaxaMatch) save(taxaMatch2, taxaMatch4, taxaMatch6, x2, x3, file = "taxaMatch.RData")

## Standardize data with gbif
# unique(obis3$basisOfRecord)
obis3$basisOfRecord[is.na(obis3$basisOfRecord)] = "Unknown"
obis3$basisOfRecord[obis3$basisOfRecord == "S"] = "PreservedSpecimen"
# Also: D = Unknown; O = Observation

obis3$taxonRank = tolower(obis3$taxonRank)

# some ocurrences of foraminifera are negative. when I investigated why, I found out that some records
# from the pangea (at least) are with the depth value wrong. CORRECT this if we are going to use this data
obis3$depth = abs(obis3$depth)

# Create a field "Name with authors" in OBIS
obis3$scientificNameWithAuthor = ifelse(is.na(obis3$scientificNameAuthorship), obis3$scientificName,
                              paste(obis3$scientificName, obis3$scientificNameAuthorship))

# remove subgenus and subspecies from OBIS
temp = obis3$taxonRank
obis3$infraspecificEpithet = NA
x = vector("character", length(temp))
for (i in 1:length(temp)) {
  if (temp[i] == "subspecies") {
    x[i] = obis3[[i,"species"]]
    obis3$infraspecificEpithet[i] = gsub("^.* ", "", obis3[[i,"scientificName"]])
  } else {
    x[i] = obis3[[i,temp[i]]]
    
    # get rid of subgenus
    x2 = strsplit(x[i], " ")[[1]]
    if (length(x2) == 3) {x[i] = paste(x2[c(1,3)], collapse = " ")}
    
    # when species is NA, go to scientificName
    if (is.na(x[i])) {x[i] = obis3[[i,"scientificName"]]} 
  }
}
sum(is.na(x)) # this value
obis3$scientificName = x

## Finally, save table
write.csv(obis3, gzfile("obis.csv.gz"), row.names = FALSE)
