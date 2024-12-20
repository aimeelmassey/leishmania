setwd("/Users/aimeemassey/DOCUMENTS/RESEARCH/BRAZIL_phd/RESULTS/Bug1 - Carrion Flies/")
library(dplyr)

### Final script as of 9/11/2021 ###

### 1. import master data files ###
###Master file
cf12s <- read.csv("Run11_CF.csv", stringsAsFactors = FALSE)
### make sample column
cf12s$Counts <- as.numeric(as.character(cf12s$Counts))
cf12s$Sample <- paste(cf12s$site, cf12s$pool)

## pertenient data
cf12s <- cf12s %>% select(site, Sample, rep, Counts, Scientific_Name, Taxa_ID, Percent_Identical_Matches, Query_Sequence_Matching_Part)


### 2. Remove blanks and NTC
blanks <- cf12s %>% filter(str_detect(Sample, "BLANK|Blank|NTC"))
cf12s <- cf12s %>% filter(!str_detect(Sample,"BLANK|Blank|NTC"))



### 3. Remove reads <1% of sample total, and remove obs. with %identity <96%
cf12s <- cf12s %>% group_by(Sample, Query_Sequence_Matching_Part) %>%
  mutate(OTUsum = sum(Counts)) %>%
  ungroup()

cf12s <- cf12s %>% group_by(Sample) %>%
  mutate(totalsum = sum(Counts)) %>%
  ungroup()

## %totalsum
cf12s <- cf12s %>% filter(OTUsum>0.01*totalsum)
## %identity
cf12s <- cf12s %>% group_by(Sample) %>%
  filter(Percent_Identical_Matches>95.99)



### 4. Create new data frame with only pertinent information
## sppsum = sum of reads for a particular NAME in a sample
cf12s <-cf12s %>% select(site, Sample, rep, Scientific_Name, Percent_Identical_Matches, Query_Sequence_Matching_Part, OTUsum)



### 5. Take out species that don't appear in BOTH reps per sample ###
cf12s1 <- cf12s %>% select(site, Sample, rep, Scientific_Name)
cf12s1 <- unique(cf12s1)


cf12s1$rep<-gsub("a", "2", cf12s1$rep, fixed=TRUE)
cf12s1$rep<-gsub("b", "3", cf12s1$rep, fixed=TRUE)
cf12s1$rep<-as.numeric(cf12s1$rep)

cf12s1 <- cf12s1 %>% 
  group_by(Sample, Scientific_Name) %>% 
  mutate(repsum = sum(rep)) %>%
  filter(repsum >= 5)

cf12s1$rep<-gsub("2", "a", cf12s1$rep, fixed=TRUE)
cf12s1$rep<-gsub("3", "b", cf12s1$rep, fixed=TRUE)

cf12s <- merge(cf12s, cf12s1, by=c("site","Sample", "rep", "Scientific_Name"))





### 6. Contaminant Removals ###
##humans
cf12s<-cf12s[!grepl("Catopuma", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Chordata", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Macaca", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Pan", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Pongo", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Trachypithecus", cf12s$Scientific_Name),]

##other
cf12s<-cf12s[!grepl("vertebrate", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Ursus", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Acomys", cf12s$Scientific_Name),] 
cf12s<-cf12s[!grepl("Cervus elaphus", cf12s$Scientific_Name),] 
cf12s<-cf12s[!grepl("Zaedyus pichiy", cf12s$Scientific_Name),] 

### REMOVALS after manuscript revision 9/11/2021
## 98-100% matches
cf12s<-cf12s[!grepl("Dicamptodon tenebrosus", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Melospiza melodia", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Podilymbus podiceps", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Spilogale putorius", cf12s$Scientific_Name),]
## 90-97.99% matches
cf12s<-cf12s[!grepl("Sus scrofa", cf12s$Scientific_Name),]
cf12s<-cf12s[!grepl("Lophuromys sikapusi", cf12s$Scientific_Name),]



### REPLACEMENTS 5/25/24 (based on leishmania sandfly  data)
cf12s$Scientific_Name[cf12s$Scientific_Name == "Ateles belzebuth"] = "Ateles marginatus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Canis latrans"] = "Canis lupus familiaris"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Coccyzus erythropthalmus"] = "Piaya cayana"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Dasyprocta punctata"] = "Dasyprocta spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Dasypus yepesi"] = "Dasypus novemcinctus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Didelphis marsupialis"] = "Didelphis spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Felis margarita"] = "Felis catus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Felis silvestris"] = "Felis catus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Gymnogyps californianus"] = "Coragyps atratus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Phoenicopterus roseus"] = "Coragyps atratus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Fratercula cirrhata"] = "Coragyps atratus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Lycalopex sechurae"] = "Cerdocyon thous"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Plecturocebus cupreus"] = "Plecturocebus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Sapajus apella"] = "Sapajus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Tapirus pinchaque"] = "Tapirus terrestris"

cf12s$Scientific_Name[cf12s$Scientific_Name == "Alouatta seniculus macconnelli"] = "Mico spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Amazona aestiva"] = "Amazona spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Anthus hodgsoni"] = "Anthus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Anthus novaeseelandiae"] = "Anthus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Brycon sp. UFRGS11377"] = "Leporinus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Buteo lagopus"] = "Buteo spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Buteo polyosoma"] = "Buteo spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Fringilla polatzeki"] = "Thraupidae (family)"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Gracilinanus emiliae"] = "Gracilinanus spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Lepidocolaptes wagleri"] = "Lepidocolaptes spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Nothoprocta ornata"] = "Nothoprocta spp."
cf12s$Scientific_Name[cf12s$Scientific_Name == "Nyctiphrynus mcleodii"] = "Nyctiphrynus ocellatus"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Ovis ammon darwini"] = "Ovis aries"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Tiaris olivaceus"] = "Thraupidae (family)"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Vanellus vanellus"] = "Tinamidae (family)"
cf12s$Scientific_Name[cf12s$Scientific_Name == "Vultur gryphus"] = "Cathartes aura"



### select only pertinent information and sum counts by Scientific Name
cf12s <- cf12s %>% select(site, Sample, OTUsum, Scientific_Name)

cf12s <- cf12s %>% group_by(Sample, Scientific_Name) %>% mutate(sppsum = sum(OTUsum))
cf12s <- cf12s %>% select(site, Sample, Scientific_Name, sppsum)
cf12s <- unique(cf12s)


## how many samples in final data
length(unique(cf12s[["site"]])) 
length(unique(cf12s[["Scientific_Name"]])) 
list(unique(cf12s[["Scientific_Name"]]))
length(unique(cf12s[["Sample"]]))





