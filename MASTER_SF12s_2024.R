####################################
### Final SF 12s quality control ###
####################################
setwd("/Users/aimeemassey/DOCUMENTS/RESEARCH/BRAZIL_phd/RESULTS/Bug2 - Sandflies/SF_12s/")
library(dplyr)
library(stringr)

### 1. import csv files
run10sf12s <- read.csv("Run10_Miseq_SF_Results_Master.csv", header = T, stringsAsFactors = FALSE)
run11sf12s <- read.csv("Run11_SF_Results_Master.csv", header = T, stringsAsFactors = FALSE)
run12sf12s <- read.csv("Run12_SF_Results_Master.csv", header = T, stringsAsFactors = FALSE)


### 2. collate master files
sf12s <- rbind(run10sf12s, run11sf12s, run12sf12s)
sf12s$Counts <- as.numeric(as.character(sf12s$Counts))
sf12s$Sample <- paste(sf12s$site, sf12s$pool)
sf12s <- as.data.frame(sf12s, stringsAsFactors = FALSE)

## pertenient data
sf12s <- sf12s %>% select(site, Sample, rep, Counts, Scientific_Name, Taxa_ID, Percent_Identical_Matches, Query_Sequence_Matching_Part)


### 3. Remove blanks and NTC
blanks <- sf12s %>% filter(str_detect(Sample, "BLANK|Blank|NTC"))
sf12s <- sf12s %>% filter(!str_detect(Sample,"BLANK|Blank|NTC"))




### 4. Remove reads <1% of sample total, and remove obs. with %identity <96%
sf12s <- sf12s %>% group_by(Sample, Query_Sequence_Matching_Part) %>%
  mutate(OTUsum = sum(Counts)) %>%
  ungroup()

sf12s <- sf12s %>% group_by(Sample) %>%
  mutate(totalsum = sum(Counts)) %>%
  ungroup()

## %totalsum
sf12s <- sf12s %>% filter(OTUsum>0.01*totalsum)
## %identity
sf12s <- sf12s %>% group_by(Sample) %>%
  filter(Percent_Identical_Matches>95.99)




### 6. Create new data frame with only pertinent information
## sppsum = sum of reads for a particular NAME in a sample
sf12s <-sf12s %>% select(site, Sample, rep, Scientific_Name, Percent_Identical_Matches, Query_Sequence_Matching_Part, OTUsum)

### 7. Take out species that don't appear in BOTH reps per sample ###
sf12s1 <- sf12s %>% select(site, Sample, rep, Scientific_Name)
sf12s1 <- unique(sf12s1)


sf12s1$rep<-gsub("a", "2", sf12s1$rep, fixed=TRUE)
sf12s1$rep<-gsub("b", "3", sf12s1$rep, fixed=TRUE)
sf12s1$rep<-as.numeric(sf12s1$rep)

sf12s1 <- sf12s1 %>% 
  group_by(Sample, Scientific_Name) %>% 
  mutate(repsum = sum(rep)) %>%
  filter(repsum >= 5)

sf12s1$rep<-gsub("2", "a", sf12s1$rep, fixed=TRUE)
sf12s1$rep<-gsub("3", "b", sf12s1$rep, fixed=TRUE)

sf12s <- merge(sf12s, sf12s1, by=c("site","Sample", "rep", "Scientific_Name"))




### Contaminants REMOVALS###
### humans
sf12s<-sf12s[!grepl("Macaca mulatta", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Pan paniscus", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Homo heidelbergensis", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Chordata environmental sample", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Catopuma temminckii", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Pongo abelii", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Pan troglodytes troglodytes", sf12s$Scientific_Name),]
sf12s<-sf12s[!grepl("Pan troglodytes verus", sf12s$Scientific_Name),]


## other = PNW lab contaminants; or no regional match in Brazil
sf12s<-sf12s[!grepl("Acomys cahirinus", sf12s$Scientific_Name),] 
sf12s<-sf12s[!grepl("Callospermophilus lateralis", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Cervus elaphus", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Erethizon dorsatum", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Glaucomys sabrinus", sf12s$Scientific_Name),] 
sf12s<-sf12s[!grepl("Junco hyemalis", sf12s$Scientific_Name),] 
sf12s<-sf12s[!grepl("Marmota flaviventris", sf12s$Scientific_Name),] ##super low count after cleaning (<100)
sf12s<-sf12s[!grepl("Martes", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Microtus", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Myodes gapperi", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Regulus regulus", sf12s$Scientific_Name),] ##unknown
sf12s<-sf12s[!grepl("Tamias", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Urocitellus columbianus", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Ursus americanus", sf12s$Scientific_Name),] ##lab
sf12s<-sf12s[!grepl("Ursus spelaeus", sf12s$Scientific_Name),] ##lab - likely brown bear
sf12s<-sf12s[!grepl("Zaedyus pichiy", sf12s$Scientific_Name),] ##very low read number
sf12s<-sf12s[!grepl("vertebrate environmental sample", sf12s$Scientific_Name),]


## also remove canis reads below blank threshold from blanks SF15 and SF13
sf12s <- sf12s %>% filter(!(Sample %in% c("B5 36" , "B5 38" , "B5 41" , "B5 42" , "B5 47" , "B5 48" , "B5 49" , "B5 51") &
                              Scientific_Name == "Canis lupus familiaris"))



### REPLACEMENTS AFTER MANUAL BLASTING
sf12s$Scientific_Name[sf12s$Scientific_Name == "Ateles belzebuth"] = "Ateles marginatus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Canis latrans"] = "Canis lupus familiaris"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Coccyzus erythropthalmus"] = "Piaya cayana"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Coendou melanurus"] = "Coendou prehensilis"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Cuniculus taczanowskii"] = "Cuniculus paca"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Dasypus pilosus"] = "Dasypus novemcinctus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Dasypus yepesi"] = "Dasypus novemcinctus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Dasyprocta punctata"] = "Dasyprocta spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Didelphis marsupialis"] = "Didelphis spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Didelphis imperfecta"] = "Didelphis spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Felis margarita"] = "Felis catus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Felis silvestris"] = "Felis catus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Fratercula cirrhata"] = "Coragyps atratus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Gymnogyps californianus"] = "Coragyps atratus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Lepus americanus"] = "Leporidae (family)"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Lepus californicus"] = "Leporidae (family)"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Lycalopex sechurae"] = "Cerdocyon thous" ##better match
sf12s$Scientific_Name[sf12s$Scientific_Name == "Lynx rufus"] = "Puma concolor"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Meleagris gallopavo"] = "Gallus gallus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Odocoileus spp."] = "Mazama americana"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Oecomys concolor"] = "Oecomys spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Oecomys cf. rex AK-2017"] = "Oecomys spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Phoenicopterus roseus"] = "Coragyps atratus"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Plecturocebus cupreus"] = "Plecturocebus spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Proechimys cuvieri"] = "Proechimys spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Sapajus apella"] = "Sapajus spp."
sf12s$Scientific_Name[sf12s$Scientific_Name == "Sylvilagus audubonii"] = "Leporidae (family)"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Sylvilagus floridanus"] = "Leporidae (family)"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Tamandua mexicana"] = "Tamandua tetradactyla"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Tapirus pinchaque"] = "Tapirus terrestris"
sf12s$Scientific_Name[sf12s$Scientific_Name == "Ramphastos cuvieri"] = "Ramphastos tucanus" 
sf12s$Scientific_Name[sf12s$Scientific_Name == "Ramphastos ambiguus ambiguus"] = "Ramphastos tucanus" 
sf12s$Scientific_Name[sf12s$Scientific_Name == "Ramphastos swainsonii"] = "Ramphastos tucanus" 
sf12s$Scientific_Name[sf12s$Scientific_Name == "Thamnophilus nigrocinereus"] = "Thamnophilidae (family)"










### select only pertinent information and sum counts by Scientific Name
sf12s <- sf12s %>% select(site, Sample, OTUsum, Scientific_Name)

sf12s <- sf12s %>% group_by(Sample, Scientific_Name) %>% mutate(sppsum = sum(OTUsum))
sf12s <- sf12s %>% select(site, Sample, Scientific_Name, sppsum)
sf12s <- unique(sf12s)













## how many samples in final data
length(unique(sf12s[["site"]]))
length(unique(sf12s[["Scientific_Name"]]))
list(unique(sf12s[["Scientific_Name"]])) 










