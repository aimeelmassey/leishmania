setwd("/Users/aimeemassey/DOCUMENTS/RESEARCH/BRAZIL_phd/RESULTS/Bug2 - Sandflies/SFspp_COI/")
library(dplyr)
library(tidyr)


### 1. import master data files ###
#Master files
Msfspp12 <- read.csv("Run12_SFspp_COI_Midori.csv", stringsAsFactors = FALSE)
MsfsppA <- read.csv("Run15_LibA_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppB <- read.csv("Run15_LibB_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppC <- read.csv("Run15_LibC_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppD <- read.csv("Run15_LibD_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppE <- read.csv("Run15_LibE_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppF <- read.csv("Run15_LibF_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppG <- read.csv("Run15_LibG_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)
MsfsppH <- read.csv("Run15_LibH_sfCOI_MIDORI.csv", stringsAsFactors = FALSE)

Msfspp <- rbind(Msfspp12,MsfsppA, MsfsppB, MsfsppC, MsfsppD, MsfsppE, MsfsppF, MsfsppG, MsfsppH)

Msfspp$Counts <- as.numeric(as.character(Msfspp$Counts))
Msfspp$Sample <- paste(Msfspp$site, Msfspp$pool)
Msfspp <- as.data.frame(Msfspp, stringsAsFactors = FALSE)
Msfspp$Scientific_Name <- Msfspp$X.6
Msfspp$Family <- Msfspp$X.4

Msfspp <- Msfspp %>% select(site, Sample, rep, Counts, Family,Scientific_Name, Percent_Identical_Matches, Query_Sequence_Matching_Part)

### 1a. Remove blanks and NTC ###
Msfspp<-Msfspp[!grepl("BLANK", Msfspp$Sample),]
Msfspp<-Msfspp[!grepl("Blank", Msfspp$Sample),]
Msfspp<-Msfspp[!grepl("NTC", Msfspp$Sample),]


### 2. Remove sample replicates below 500 read threshold per sample ###
Msfspp <- Msfspp %>% group_by(Sample) %>% 
  mutate(totalsum =sum(Counts)) %>% 
  ungroup()

Msfspp <- Msfspp %>% group_by(Sample) %>% 
  filter(totalsum>500) %>%
  ungroup()


### 3. select only SF species - from Family Psychodidae
SF <- Msfspp %>% filter(Family=="Psychodidae")
SF <- SF %>% separate(Scientific_Name, c("Genus","Species"), sep = "([_])")
## remove non-SF genus (likely moth species since they share Psychodidae family)
SF <- SF[!grepl("Parabazarella", SF$Genus),] ## moth genus
SF <- SF[!grepl("Parajungiella", SF$Genus),] ## moth genus
SF <- SF[!grepl("Psychoda", SF$Genus),] ## moth genus


### 4. Remove reads <1% of sample total, and remove obs. with %identity <96%
SF <- SF %>% group_by(Sample, Query_Sequence_Matching_Part) %>%
  mutate(OTUsum = sum(Counts)) %>%
  ungroup()

SF <- SF %>% group_by(Sample) %>%
  mutate(totalsum = sum(Counts)) %>%
  ungroup()

## %totalsum
SF <- SF %>% filter(OTUsum>0.01*totalsum)
## %identity
SF <- SF %>% group_by(Sample) %>%
  filter(Percent_Identical_Matches>95.99)




### 5. Create new data frame with only pertinent information
SF <- SF %>% select(site, Sample, rep, Genus, Species, Counts, Query_Sequence_Matching_Part, Percent_Identical_Matches)

SF <- SF %>% group_by(Sample, Query_Sequence_Matching_Part) %>%
  mutate(sppsum = sum(Counts)) %>%
  ungroup()

SF <- SF %>% select(site, Sample, rep, Genus, Species, sppsum, Query_Sequence_Matching_Part, Percent_Identical_Matches)
SF <- unique(SF)


##Take out species that don't appear in BOTH reps per sample ###
SF1 <- SF %>% select(site, Sample, rep, Genus, Species)
SF1 <- unique(SF1)

SF1$rep<-gsub("a", "2", SF1$rep, fixed=TRUE)
SF1$rep<-gsub("b", "3", SF1$rep, fixed=TRUE)
SF1$rep<-as.numeric(SF1$rep)
SF1<-unique(SF1)

SF1 <- SF1 %>% 
  group_by(Sample, Genus, Species) %>% 
  mutate(repsum = sum(rep)) %>%
  filter(repsum >= 5)

SF1$rep<-gsub("2", "a", SF1$rep, fixed=TRUE)
SF1$rep<-gsub("3", "b", SF1$rep, fixed=TRUE)

SF <- merge(SF, SF1, by=c("site","Sample", "rep", "Genus", "Species"))


### REPLACEMENTS
SF$Scientific_Name = paste(SF$Genus, SF$Species, sep="_")
SF$Scientific_Name[SF$Scientific_Name == "Lutzomyia_walkeri"] = "Evandromyia_walkeri"
SF$Scientific_Name[SF$Scientific_Name == "Trichophoromyia_ininii"] = "Trichophoromyia_auraensis"

SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_umbratilis"] = "Nyssomyia_spp."
SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_whitmani"] = "Nyssomyia_spp."
SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_anduzei"] = "Nyssomyia_spp."
SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_flaviscutellata"] = "Nyssomyia_spp."
SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_intermedia"] = "Nyssomyia_spp."
SF$Scientific_Name[SF$Scientific_Name == "Nyssomyia_yuilli"] = "Nyssomyia_spp."






### LAST STEP. Remove rep information and make new totalsum variable ###
SF <- subset(SF, select = -c(rep, repsum))
SF <- unique(SF)
SF <- SF %>% group_by(Sample, Scientific_Name) %>% mutate(sppsum = sum(sppsum))
SF <- subset(SF, select = -c(Query_Sequence_Matching_Part, Percent_Identical_Matches))
SF <- unique(SF)


SF <- SF %>% group_by(Sample) %>%
  mutate(totalsumQC = sum(sppsum))


## how many samples in final data
length(unique(SF[["site"]]))
length(unique(SF[["Scientific_Name"]]))
list(unique(SF[["Scientific_Name"]])) 



