## Global wildlife diversity in team sport organisations
## Corey Bradshaw & Ugo Arbieu
## April 2025

rm(list=ls()) # remove all objects

## load required libraries
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(lme4)
library(Hmisc)
library(boot)
library(ggpubr)
library(performance)
library(sjPlot)
library(ggplot2)

# source files for statistical analyses
source("new_lmer_AIC_tables3.R") 
source("r.squared.R") 

# load data
dat <- read.csv("dat.csv", header=T, sep=";")
meta <- read.csv("metadat.csv", header=T, sep=";")

colnames(meta) <- c("Gender","Continent","Country","Sport","TotTeams","AnimTeams")
meta$Percent <- meta$AnimTeams/meta$TotTeams*100
head(meta)

##############################################
##############################################
## 1. PRELIMINARY DESCRIPTION OF THE DATA 
##############################################
##############################################

    # number of professional leagues
    metaduce <- meta[complete.cases(meta),] # 230 professional leagues
    
    # total number of teams represented in those leagues
    sum(metaduce$TotTeams) # 2970
    length(which(metaduce$AnimTeams==0)) # 25 leagues without any animal >> 205 prof leagues with an animal
    
    # Exploring the number of Men and Women teams
    nteams <- length(unique(dat$Club)) #834 clubs in total in the dataset
    nteams
    length(unique(dat$Club[which(dat$Gender=="Men")])) #682
    length(unique(dat$Club[which(dat$Gender=="Women")])) #174
    
    # Number of unique taxons
    sum(table(dat$Classification_General)) #879
    
    # Number of teams that have a wild animal as a representation
    ndom <- length(which(dat$Trophic != "Domestic")) # 747 occurrences of wild animal (note that "occurrences" are different from "teams": some teams have several animals depicted in the logo)
    ndom
    nmendom <- length(unique(dat$Club[which(dat$Gender=="Men" & dat$Trophic != "Domestic")])) # 579 Men's teams have a wild animal
    nmendom
    nwomdom <- length(unique(dat$Club[which(dat$Gender=="Women" & dat$Trophic != "Domestic")])) # 149 Women's teams have a wild animal
    nwomdom
    
    nteams-(nmendom+nwomdom) # 106 teams use domestic forms
    
    # the number of occurrences of wild taxa are identifiable at the species level
    length(which(dat$Classification_General=="Species" & dat$Trophic!="Domestic")) # 445 occurrences
    length(which(dat$Classification_General=="Species" & dat$Gender=="Men" & dat$Trophic!="Domestic")) #353 men's occurrences
    length(which(dat$Classification_General=="Species" & dat$Gender=="Women" & dat$Trophic!="Domestic")) #92 women's occurrences
    
    # Dometic species concerned (n=14)
    unique(dat$Classification_Specific[which(dat$Trophic=="Domestic")])


#############################################################
# Percentage of wildlife teams in each league, across sports
#############################################################

# Note that we have not removed domestic species from the dataset yet

for (i in 1:dim(meta)[1]){
  focus_all <- subset(dat, Gender==meta$Gender[i] & Continent==meta$Continent[i] & Country==meta$Country[i] & Sport==meta$Sport[i])
  meta$naf_all[i] <- length(unique(focus_all$Classification_Specific)) # number of all animal forms (wild + domestic) in the league
  focus_wild <- subset(dat, Gender==meta$Gender[i] & Continent==meta$Continent[i] & Country==meta$Country[i] & Sport==meta$Sport[i] & Trophic!="Domestic")
  meta$naf_wild[i] <- length(unique(focus_wild$Classification_Specific)) # number of all wild animal forms in the league
  meta$AnimTeamsWild[i] <- length(unique(focus_wild$Club))
}

meta$Percent2 <- meta$AnimTeamsWild/meta$TotTeams*100

meta$Sport=factor(meta$Sport, levels=c("Rugby XIII", "Cricket", "Baseball", "American Football", "Rugby XV", "Ice Hockey", "Football", "Basketball", "Volleyball", "Handball"))

# Plot of the percentage of teams using only WILD "animal forms" in each sport
panel_a <- ggplot(meta, aes(y=Percent2, x=Sport, colour=Sport)) + 
  geom_boxplot(outlier.shape=NA)+
  theme_classic()+
  geom_jitter(width = 0.2, aes(colour=Sport))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(size=15),
        legend.position = "none")+
  scale_colour_manual(values=c("#94F2D1", "#D59404", "#8A8C8A", "#A0082B", "#C4657A","#99CCF6","#027017", "#A684BE","#B8B208","#060D89"))

panel_a

# Exploration of the number of leagues concerned
metared <- meta[is.na(meta$TotTeams)==F,]
dim(metared)[1] # 230 professional leagues investigated
length(which(metared$Gender=="Men")) #163 men's leagues
length(which(metared$Gender=="Women")) #67 women's leagues
dim(metared)[1]-length(which(metared$Percent2==0)) # 197 leagues with wild animals
dim(metared)[1]-length(which(metared$Percent==0)) # 205 leagues with wild+domestic animals


# Storing the number of wildlife teams across sports, genders, continents and countries
sport_ta <- aggregate(metared$AnimTeamsWild, by=list(Sport=metared$Sport), FUN=sum)
gender_ta <- aggregate(metared$AnimTeamsWild, by=list(Gender=metared$Gender), FUN=sum)
continent_ta <- aggregate(metared$AnimTeamsWild, by=list(Continent=metared$Continent), FUN=sum)
country_ta <- aggregate(metared$AnimTeamsWild, by=list(Country=metared$Country), FUN=sum)

########################################################################################################
# Ranking and Distributions among classes / orders/ families / species
########################################################################################################

# We now remove the domestic animal forms from the dataset
dat <- subset(dat, dat$Trophic!="Domestic")

# ranking of classes in the dataset
tab_class <- table(dat$Classification_Class)
datclass <- data.frame(sort(tab_class))
datclass$id  <- 1:nrow(datclass)
colnames(datclass) <- c("Classification_Class","Freq","id")

# ranking of orders in the dataset
tab_order <- table(dat$Classification_Order)
datorder <- data.frame(sort(tab_order))
datorder$id  <- 1:nrow(datorder)
colnames(datorder) <- c("Classification_Order","Freq","id")

# ranking of families in the dataset
tab_family <- table(dat$Classification_Family)
datfamily <- data.frame(sort(tab_family))
datfamily$id  <- 1:nrow(datfamily)
colnames(datfamily) <- c("Classification_Family","Freq","id")

#ranking of species in the dataset
tab_species <- table(dat$Classification_Species)
datspecies <- data.frame(sort(tab_species))
datspecies$id  <- 1:nrow(datspecies) #48 teams are the only ones using an identifiable species.
colnames(datspecies) <- c("Classification_Species","Freq","id") # 85 species in total

####################################################
# Ranking of the different groups
####################################################

# ranking of animal representations
tab_forms <- table(dat$Classification_Specific)
datforms <- data.frame(sort(tab_forms))
datforms$id  <- 1:nrow(datforms)
dim(datforms)[1] # number of distinct taxa
colnames(datforms) <- c("Classification_Specific","Freq","id")
length(which(datforms$Freq==1)) # 79 teams are the only ones to use a unique animal taxon as an emblem

# creating files to store the distributions among classes / orders / families / species
dat_class <- data.frame(dat$Classification_Class)
dat_u_class <- table(dat_class)
#write.csv(dat_u_class, "unique_class_list.csv")

dat_oc_comb <- data.frame(dat$Classification_Order,dat$Classification_Class)
dat_oc_comb <- unique(dat_oc_comb)
colnames(dat_oc_comb) <- c("Classification_Order","Classification_Class")
datm_oc <- merge(datorder,dat_oc_comb, by="Classification_Order",all.x=TRUE)
#write.csv(datm_oc, "unique_order_list.csv")
datm_oc <- datm_oc[order(datm_oc$id),]

dat_fc_comb <- data.frame(dat$Classification_Family,dat$Classification_Class)
dat_fc_comb <- unique(dat_fc_comb)
colnames(dat_fc_comb) <- c("Classification_Family","Classification_Class")
datm_fc <- merge(datfamily,dat_fc_comb, by="Classification_Family",all.x=TRUE)
#write.csv(datm_fc, "unique_family_list.csv")
datm_fc <- datm_fc[order(datm_fc$id),]

dat_sc_comb <- data.frame(dat$Classification_Species,dat$Classification_Class)
dat_sc_comb <- unique(dat_sc_comb)
colnames(dat_sc_comb) <- c("Classification_Species","Classification_Class")

datm_sc <- merge(datspecies,dat_sc_comb, by="Classification_Species",all.x=TRUE)
#write.csv(datm_sc, "unique_species_list.csv")
datm_sc <- datm_sc[order(datm_sc$id),]


#######################################################
#######################################################
## 3. PREPARATION OF FIGURE 2 - histogram of wild taxa
#######################################################
#######################################################

# Classes
colors_class <- c("#4abec2",   "#ce465a",   "#7cb843",   "#7e65cc",   "#a5903e",   "#c554b2",   "#529c61",  "#bd698f",    "#7089cc",   "#ce703d")
plot_class <- ggplot(datclass, aes(x = Classification_Class, y = Freq, fill=as.factor(Classification_Class))) + 
  geom_bar(stat = "identity", width=.8) +
  labs(x = "\nclass", y = "frequency\n") +
  theme_classic()+
  theme(axis.text = element_text(size=15))+
  scale_fill_manual(values=colors_class)+
  theme(
    legend.position = "none"
  )+
  coord_flip()
plot_class

# Orders
colors_order <- c("#529c61","#4abec2", "#7e65cc","#7089cc","#c554b2","#bd698f","#7cb843","#ce703d","#a5903e")
plot_order <- ggplot(datm_oc, aes(x = Classification_Order, y = Freq, fill=as.factor(Classification_Class))) + 
  geom_bar(stat = "identity", width=.8) +
  labs(x = "\norder", y = "frequency\n") +
  theme_classic()+
  theme(axis.text = element_text(size=15))+
  scale_fill_manual(values=colors_order)+
  theme(
    legend.position = "none"
  )+
  coord_flip()
plot_order

colors_family <- c("#529c61", "#7e65cc","#7089cc","#c554b2","#bd698f","#ce703d","#a5903e")
plot_family <- ggplot(datm_fc, aes(x = Classification_Family, y = Freq, fill=as.factor(Classification_Class))) + 
  geom_bar(stat = "identity", width=.8) +
  labs(x = "\nfamily", y = "frequency\n") +
  theme_classic()+
  scale_fill_manual(values=colors_family)+
  theme(axis.text = element_text(size=15))+
  theme(
    legend.position = "none"
  )+
  coord_flip()
plot_family

colors_species <- c("#529c61", "#7e65cc","#7089cc","#c554b2","#bd698f","#ce703d","#a5903e")
plot_species <- ggplot(datm_sc, aes(x = Classification_Species, y = Freq, fill=as.factor(Classification_Class))) + 
  geom_bar(stat = "identity", width=.8) +
  labs(x = "\nspecies", y = "frequency\n") +
  theme_classic()+
  theme(axis.text = element_text(face='italic', size=15))+
  scale_fill_manual(values=colors_species)+
  theme(
    #legend.position = "none"
  )+
  coord_flip()

colors_species_15 <- c("#7089cc","#ce703d")
datm_sc_15 <- tail(datm_sc,15)
plot_species_15 <- ggplot(datm_sc_15, aes(x = Classification_Species, y = Freq, fill=as.factor(Classification_Class))) + 
  geom_bar(stat = "identity", width=.8) +
  labs(x = "\nspecies", y = "frequency\n") +
  theme_classic()+
  theme(axis.text = element_text(face='italic', size=15))+
  scale_fill_manual(values=colors_species_15)+
  theme(
    legend.position = "none"
  )+
  coord_flip()
plot_species_15

grid.arrange(plot_class,plot_species_15, nrow = 1)


##############################################
##############################################
## 3. PREPARATION OF SUPPLEMENTARY TABLE
##############################################
##############################################

record <- unique(dat$Classification_General)
# species list
dats <- subset(dat,Classification_General==record[1])
uniques <- table(dats$Classification_Specific)
#write.csv(uniques, "uniquespecies.csv")

# Family list
datf <- subset(dat,Classification_General==record[2])
uniquef <- table(datf$Classification_Specific)
#write.csv(uniquef, "uniquefamilies.csv")

# Class list
datc <- subset(dat,Classification_General==record[3])
uniquec <- table(datc$Classification_Specific)
#write.csv(uniquec, "uniqueclasses.csv")

# Order list
dato <- subset(dat,Classification_General==record[4])
uniqueo <- table(dato$Classification_Specific)
#write.csv(uniqueo, "uniqueorders.csv")          

# Infraorder list
dati <- subset(dat,Classification_General==record[5])
uniquei <- table(dati$Classification_Specific)
#write.csv(uniquei, "uniqueinfraorders.csv")

# Suborder list
datsub <- subset(dat,Classification_General==record[6])
uniquesub <- table(datsub$Classification_Specific)
#write.csv(uniquesub, "uniquesuborders.csv")

# Genus list
datg <- subset(dat,Classification_General==record[7])
uniqueg <- table(datg$Classification_Specific)
#write.csv(uniqueg, "uniquegenus.csv")

# Subfamily list
datsf <- subset(dat,Classification_General==record[8])
uniquesf <- table(datsf$Classification_Specific)
#write.csv(uniquesf, "uniquesubfamily.csv")

# Superorder list
datso <- subset(dat,Classification_General==record[9])
uniqueso <- table(datso$Classification_Specific)
#write.csv(uniqueso, "uniquesuperorder.csv")

##############################################################
##############################################################
## 4. PREPARATION OF FIGURE S1 (=differences across genders)
##############################################################
##############################################################

# Calculating the raw number of wild taxa per gender
nsp3=c()
for (i in 1:length(unique(dat$Gender))){
  cntsp <- subset(dat,Gender==unique(dat$Gender)[i])
  nsp3[i] <- length(unique(cntsp$Classification_Specific))
}

gender <- unique(dat$Gender)
datgender <- data.frame(gender,nsp3)
datgender <- datgender[order(datgender$gender),]
datgender$relative <- datgender$nsp3/gender_ta$x

panel_d_relative <- ggplot(datgender, aes(y=relative, x=gender, fill=gender)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text=element_text(size=15),
        legend.position = "none")+
  scale_fill_manual(values=c("#A4118E","#11A427"))

# plot of % per gender and continent, for only wild taxa
panel_c <- ggplot(meta, aes(y=Percent2, x=Continent,colour=Gender), group = interaction(Continent, Gender)) + 
  #geom_bar(position="dodge", stat="identity")+
  #geom_violin(trim=FALSE)+
  geom_boxplot(outlier.shape=NA)+
  theme_classic()+
  geom_jitter(aes(colour = Gender), position = position_jitterdodge(.1)) +
  scale_fill_manual(values=c("#A4118E","#11A427"))+
  scale_colour_manual(values=c("#A4118E","#11A427"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none"
  )

# average percentage of teams having a wild animal as an emblem, per gender
aggregate(metared$Percent2, list(metared$Gender), FUN=mean)

grid.arrange(panel_d_relative,panel_c, nrow = 1)


#####################################################
#####################################################
## 4. PREPARATION OF FIGURE 3 (=biplot of countries)
#####################################################
#####################################################

# Calculating the raw number of wild taxa per continent
nsp2=c()
for (i in 1:length(unique(dat$Continent))){
  cntsp <- subset(dat,Continent==unique(dat$Continent)[i])
  nsp2[i] <- length(unique(cntsp$Classification_Specific))
}

contin <- unique(dat$Continent)
datcontinent <- data.frame(contin,nsp2)
datcontinent <- datcontinent[order(datcontinent$contin),]
datcontinent$relative <- datcontinent$nsp2/continent_ta$x


# Calculating the raw number of wild taxa per country and continent
nsp4=c()
continentcountry <- c()
for (i in 1:length(unique(dat$Country))){
  cntsp <- subset(dat,Country==unique(dat$Country)[i])
  nsp4[i] <- length(unique(cntsp$Classification_Specific))
  continentcountry[i] <- cntsp$Continent[1]
}
country <- unique(dat$Country)
datcountry <- data.frame(country,continentcountry,nsp4)
datcountry <- datcountry[order(datcountry$country),]
country_ta <- country_ta[-which(country_ta$Country=="French Polynesia"),] #removing French Polynesia that has no team represented
datcountry$NWildTeams <- country_ta$x
datcountry$relative <- datcountry$nsp4/country_ta$x
datcountry <- datcountry[order(datcountry$relative),]
datcountry$country <- factor(datcountry$country, levels=datcountry$country)

datcountry <- datcountry[order(datcountry$nsp4),]
datcountry$country <- factor(datcountry$country, levels=datcountry$country)

panel_e_combined <- ggplot(datcountry, aes(x=nsp4, y=relative, size=NWildTeams, label=country)) + 
  geom_point(aes(color=continentcountry))+
  geom_text_repel(size=3)+
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  scale_color_manual(values=c("#1DA961","#E16B15","#158BE1","#A91D65","#A9A71D"))
panel_e_combined



##############################################
##############################################
## 5. PREPARATION OF FIGURE S3 (=sports biplot)
##############################################
##############################################
nsp=c()
for (i in 1:length(unique(dat$Sport))){
  ctsp <- subset(dat,Sport==unique(dat$Sport)[i])
  nsp[i] <- length(unique(ctsp$Classification_Specific))
}
nsp
sports <- unique(dat$Sport)
datsport <- data.frame(sports,nsp)
colnames(datsport)=c("Sport", "nsp")
datsports <- merge(datsport,sport_ta, by="Sport",all.x=TRUE)
#datsport <- datsport[order(datsport$sports),]
datsports$relative <- datsports$nsp/datsports$x # relative to the number of teams that have a wild taxon

#ordering the labels
datsports$Sport=factor(datsports$Sport, levels=c("Rugby XIII", "Cricket", "Baseball", "American Football", "Rugby XV", "Ice Hockey", "Football", "Basketball", "Volleyball", "Handball"))

panel_sportrawrel <- ggplot(datsports, aes(x=nsp, y=relative,label=Sport)) + 
  geom_point(aes(color=Sport))+
  geom_text_repel(size=3)+
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  scale_color_manual(values=c("#94F2D1", "#D59404", "#8A8C8A", "#A0082B", "#C4657A","#99CCF6","#027017", "#A684BE","#B8B208","#060D89"))
panel_sportrawrel


##################################################################################################################
##################################################################################################################
## 7. PREPARATION OF FIGURE 5 (=presence/absence of species for each country)
## Focusing on the 83 species (85 species - 2 species absent from the IUCN database)
##################################################################################################################
##################################################################################################################

datpresence <- subset(dat,Classification_General=="Species")
length(which(datpresence$N_Species_In_Country==0)) # 229 instances of species not present in the country
length(which(datpresence$N_Species_In_Country==1)) # 214 instances of species present in the country
length(which(is.na(datpresence$N_Species_In_Country)==TRUE)) # 2 NAs (absent from the IUCN database)

# Calculating the  number of species per country
pres=c()
abs=c()
Nspecies=c()
continentcountry=c()
country=c()
for (i in 1:length(unique(datpresence$Country))){
  cntsp <- subset(datpresence,Country==unique(datpresence$Country)[i])
  pres[i] <- length(which(cntsp$N_Species_In_Country==1))
  abs[i] <- length(which(cntsp$N_Species_In_Country==0))
  Nspecies[i] <- length(unique(cntsp$Classification_Specific))
  country[i] <- cntsp$Country[1]
  continentcountry[i] <- cntsp$Continent[1]
}

datpresabs <- data.frame(country,continentcountry,pres,abs,Nspecies)
datpresabs$continentcountry <-  factor(datpresabs$continentcountry, levels=c("Europe", "Americas", "Oceania", "Asia", "Africa"))
datpresabs <- datpresabs[order(datpresabs$continentcountry, datpresabs$pres),]

# transforming to long format
# Reshape data from wide to long format
long_df <- melt(datpresabs, id.vars = c("country", "continentcountry","Nspecies"), variable.name = "Presabs", value.name = "Value")
long_df$Presabs <- factor(long_df$Presabs, levels=c("abs","pres"))
long_df$country <- as.factor(long_df$country)
ord <- long_df[1:43,]
ord <- ord[order(ord$continentcountry, ord$Value),]
long_df$country <- factor(long_df$country, levels=ord$country)

presabsplot <- ggplot(long_df, aes(x = country, y = Value, fill = Presabs)) + 
  geom_hline(yintercept=c(10,20,30), color="darkgrey", lty=2)+
  geom_bar(position="stack", stat = "identity")+
  coord_flip()+
  theme_classic()+
  theme(axis.text=element_text(size=12))+
  scale_fill_manual(values=c("#d3dc85","#0598a8"))
presabsplot

##################################################################################################################################
##################################################################################################################################
## 8. PREPARATION OF FIGURE 3 (=conservation status)

## Focusing on the 83 species (85 species - 2 species absent from the IUCN database)
##################################################################################################################################
##################################################################################################################################
datconserv <- subset(dat,Classification_General=="Species")
unique(datconserv$Club[which(datconserv$Gender=="Men")]) #349
unique(datconserv$Club[which(datconserv$Gender=="Women")]) #91
91+349 #440 teams are considered in this analysis

datconserv$Classification_Species[is.na(datconserv$N_Species_In_Country)] ### species for which we could not retrieve IUCN data

# plotting the representation of classes among the 83 species
classifsp <- data.frame(unique(datconserv$Classification_Species))
colnames(classifsp) <- "Sp"
classifgen <- data.frame(datconserv$Classification_Species, datconserv$Classification_Class)
colnames(classifgen)=c("Sp","Class")
classifgen <- classifgen[!duplicated(classifgen),]

tab_class <- table(classifgen$Class)
datclass <- data.frame(sort(tab_class))
datclass$id  <- 1:nrow(datclass)
colnames(datclass) <- c("Class","Freq","id")

colors_class <- c("#7e65cc", "#c554b2", "#bd698f", "#529c61", "#a5903e", "#7089cc", "#ce703d")
# plotting the absolute number of species found for each class in sport emblems
class_abs <- ggplot(datclass, aes(x="", y=Freq, fill=Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_class)
class_abs

# plotting the number of times these species are found across all sport organizations
tab_class <- table(datconserv$Classification_Class)
datclass2 <- data.frame(sort(tab_class))
datclass2$id  <- 1:nrow(datclass2)
colnames(datclass2) <- c("Classification_Class","Freq","id")

class_rel <- ggplot(datclass2, aes(x="", y=Freq, fill=Classification_Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_class)
class_rel

# plotting the representation of conservation statuses
consstatus <- data.frame(datconserv[,24:36],datconserv$Classification_Species)
consstatus <- consstatus[!duplicated(consstatus),]
apply(consstatus[1:13], 2,sum, na.rm=T)

datcv <- data.frame(sort(apply(consstatus[1:9], 2,sum, na.rm=T)))
datcv$id  <- 1:nrow(datcv)
datcv$Category <- rownames(datcv)
colnames(datcv) <- c("Freq","id","Category")

colors_conserv <- c("#D81E05", "#D1D1C6", "#FC7F3F", "#542344", "#000000", "#60C659", "#FFFFFF","#CCE226", "#F9E814")
# plotting the number of species in each threat category of the IUCN Red List across the 83 species
conserv_abs <- ggplot(datcv, aes(x="", y=Freq, fill=Category)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_conserv)
conserv_abs

# raw number of conservation status accounting for club use of species
consstatus <- data.frame(datconserv[,24:36],datconserv$Classification_Species)
datcv2 <- data.frame(sort(apply(consstatus[1:9], 2,sum, na.rm=T)))
datcv2$id  <- 1:nrow(datcv2)
datcv2$Category <- rownames(datcv2)
colnames(datcv2) <- c("Freq","id","Category")
conserv_rel <- ggplot(datcv2, aes(x="", y=Freq, fill=Category)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_conserv)
conserv_rel

# plotting the representation of population trends across the 83 species
poptrend <- data.frame(datconserv[,33:36],datconserv$Classification_Species)
poptrend <- poptrend[!duplicated(poptrend),]

datpt <- data.frame(sort(apply(poptrend[1:4], 2,sum, na.rm=T)))
datpt$id  <- 1:nrow(datpt)
datpt$Trend <- rownames(datpt)
colnames(datpt) <- c("Freq","id","Trend")
colors_trend <- c("#d0270e", "#03af16", "#058dc7", "#dfdb02")

poptrend_abs <- ggplot(datpt, aes(x="", y=Freq, fill=Trend)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_trend)
poptrend_abs

# population trends across all sport organizations
poptrend <- data.frame(datconserv[,33:36],datconserv$Classification_Species)
datpt2 <- data.frame(sort(apply(poptrend[1:4], 2,sum, na.rm=T)))
datpt2$id  <- 1:nrow(datpt)
datpt2$Trend <- rownames(datpt)
colnames(datpt2) <- c("Freq","id","Trend")

poptrend_rel <- ggplot(datpt2, aes(x="", y=Freq, fill=Trend)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(
    legend.position = "none"
  )+
  scale_fill_manual(values=colors_trend)
poptrend_rel

grid.arrange(class_abs, class_rel, conserv_abs, conserv_rel, poptrend_abs,poptrend_rel, nrow = 3)


##############################################
##############################################
## 9. ISOLATING DATA FOR STATISTICAL ANALYSIS
##############################################
##############################################
data_Analysis_1 <- metared
data_Analysis_2 <- datconserv[c("ID", "Gender", "Sport", "Continent", "Country", "Club", "Classification_Specific", "Classification_General", "Classification_Class", "Classification_Order", "Classification_Family", "Classification_Genus", "Classification_Species", "N_Species_In_Country", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD", "NE", "Decreasing", "Increasing", "Stable", "Unknown")]

##############################################
##############################################
## 10. STATISTICAL ANALYSIS
##############################################
##############################################

# functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}


## Professional leagues in each country. To compare women’s and men’s differences in 'Percent2'
# TotTeams = total number of teams in the league
# AnimTeams = number of teams with animal (domesticated + wild) emblems in the league
# Percent = TotTeams / AnimTeams
# naf_all = total number of animal organisms found in the league
# naf_wild = total number of wildlife organisms found in the league
# AnimTeamsWild = number of teams with wildlife emblems in the league
# Percent2 = AnimTeamsWild / TotTeams

# load the data
league.dat <- data_Analysis_1
head(league.dat)

# scale & transform
hist(league.dat$Percent2, breaks=20, col="grey", xlab="% of teams with wildlife emblems", main="")
twmwldembllgtcor1 <- ifelse(league.dat$Percent2==0, 0.05, league.dat$Percent2/100)
twmwldembllgtcor <- ifelse(twmwldembllgtcor1==1, 0.95, twmwldembllgtcor1)
league.dat$tmwldembl.sc <- scale(logit(twmwldembllgtcor), scale=T, center=T)
hist(league.dat$tmwldembl.sc, col="grey", xlab="scaled logit prop of teams with wildlife emblems", main="")

# number of leagues with wildlife emblems x gender
table(league.dat$Gender)
# number of leagues with wildlife emblems x continent
table(league.dat$Continent)
xtabs(league.dat$tmwldembl.sc ~ league.dat$Continent) / table(league.dat$Continent)

ggdensity(league.dat$tmwldembl.sc, xlab = "scaled logit prop of teams with wildlife emblems")
ggqqplot(league.dat$tmwldembl.sc)
shapiro.test(league.dat$tmwldembl.sc)

###########################################################################################
hist(league.dat$AnimTeamsWild)
hist(league.dat$TotTeams)

## generalised linear mixed-effects models
# model set
m1 <- "AnimTeamsWild ~ TotTeams + Gender + TotTeams*Gender + Sport + Gender*Sport + (1|Continent)"
m2 <- "AnimTeamsWild ~ TotTeams + Gender + Sport + Gender*Sport + (1|Continent)"
m3 <- "AnimTeamsWild ~ TotTeams + Gender + Sport + (1|Continent)"
m4 <- "AnimTeamsWild ~ TotTeams + Gender + TotTeams*Gender + (1|Continent)"
m5 <- "AnimTeamsWild ~ TotTeams + Gender + (1|Continent)"
m6 <- "AnimTeamsWild ~ TotTeams + Sport + (1|Continent)"
m7 <- "AnimTeamsWild ~ TotTeams + (1|Continent)"
m8 <- "AnimTeamsWild ~ 1 + (1|Continent)"

# relevelling the sport factor to get football as the reference (most number of leagues)
league.dat$Sport=factor(league.dat$Sport, levels=c("Football", "Rugby XIII", "Cricket", "Baseball", "American Football", "Rugby XV", "Ice Hockey", "Basketball", "Volleyball", "Handball"))

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:8,1),sep="")))

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table

# here the best model contains TotTeams & Sport (Continent is a random factor) = model 6
summary(mod.list[[1]])
check_model(mod.list[[1]])
# this plot only shows the random effects, because the model only contains the random factor of continents
plot_model(mod.list[[1]], show.values=T, type="re", vline.color = "purple")
plot_model(mod.list[[1]], show.values=T, type="est", vline.color = "purple")
plot_model(mod.list[[6]], show.values=T, type="re", vline.color = "purple")
plot_model(mod.list[[6]], show.values=T, type="est", vline.color = "purple")
plot_model(mod.list[[5]], show.values=T, type="re", vline.color = "purple")
plot_model(mod.list[[5]], show.values=T, type="est", vline.color = "purple")


mod6.out <- plot_model(mod.list[[6]], show.values=T, type="est", vline.color = "purple")
str(mod6.out)

mod6.outtab <- data.frame(var=mod6.out$data$term, est=mod6.out$data$estimate, up=mod6.out$data$conf.high,
           lo=mod6.out$data$conf.low, se=mod6.out$data$std.error)
mod6.outtab

# Preparing Figure S2
ggplot(mod6.outtab, aes(x = var, y = est)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.2) +
  coord_flip() +
  labs(x = "Variable", y = "Estimate", title = "Estimates with Confidence Intervals") +
  theme_classic()+
  theme(axis.text=element_text(size=12))


########################################
## Database focusing on the taxa that could be identified at the species level
## “Classification” columns = finest possible level of identification
# “N_Species_In_Country” =  number of corresponding species that occur in the respective country,
#   according to the IUCN database
# “EX” … “NE” = IUCN threat levels &  corresponding number of species in that country
## “Decreasing” … “Unknown” = population trends according to the IUCN database

## “N_Species_In_Country”  used for analysis of presence/absence (0=absent, 1=present) in the country

# load data
spp.dat <- data_Analysis_2
head(spp.dat)


## differences in class representation
table(spp.dat$Continent)
table(spp.dat$Country)
table(spp.dat$Classification_Class)
table(spp.dat$N_Species_In_Country)
table(spp.dat$N_Species_In_Country, spp.dat$Continent, spp.dat$Classification_Class)
# regrouping to avoid divergence issues in the model
spp.dat$ClassComb <- ifelse(spp.dat$Classification_Class == "Actinopterygii" | 
                              spp.dat$Classification_Class == "Arachnida" |
                              spp.dat$Classification_Class == "Chondrichthyes" |
                              spp.dat$Classification_Class == "Reptilia" |
                              spp.dat$Classification_Class == "Insecta", "other", 
                            spp.dat$Classification_Class)
table(spp.dat$ClassComb)
table(spp.dat$Continent, spp.dat$ClassComb)

## simulation test of independence
class.tab <- table(spp.dat$Continent, spp.dat$ClassComb)
class.tab
contSums <- rowSums(class.tab)
classSums <- colSums(class.tab)

# exact binomial test
binom.test(classSums[2], sum(classSums), 1/3, alternative="greater")


exp.tab <- class.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- classSums[1] * (contSums[i]/sum(class.tab))
  exp.tab[i,2] <- classSums[2] * (contSums[i]/sum(class.tab))
  exp.tab[i,3] <- classSums[3] * (contSums[i]/sum(class.tab))
}
exp.tab
chi.tab <- (class.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chisq.test(class.tab, simulate.p.value=T, B=10000)
# significant differences of animal group selection across continents

iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(0,iter)
comp.vec <- rep(0,iter)

for (c in 1:iter) {
  rnd.tab <- class.tab
  chi.tab.adj <- chi.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],classSums[1]/sum(classSums))
  rnd.tab[1,2] <- rbinom(1,contSums[1],classSums[2]/sum(classSums))
  rnd.tab[1,3] <- contSums[1] - sum(rnd.tab[1,1], rnd.tab[1,2])
  mam.less <- ifelse(rnd.tab[1,3] < 0, rbinom(1, abs(rnd.tab[1,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[1,3]), abs(rnd.tab[1,3]) - mam.less, 0)
  rnd.tab[1,2] <- rnd.tab[1,2] - mam.less
  rnd.tab[1,1] <- ifelse(rnd.tab[1,3] > 0, rnd.tab[1,1], rnd.tab[1,1] - ave.less)
  rnd.tab[1,3] <- ifelse(rnd.tab[1,3] < 0, 0, rnd.tab[1,3])
  
  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],classSums[1]/sum(classSums))
  rnd.tab[2,2] <- rbinom(1,contSums[2],classSums[2]/sum(classSums))
  rnd.tab[2,3] <- contSums[2] - sum(rnd.tab[2,1], rnd.tab[2,2])
  mam.less <- ifelse(rnd.tab[2,3] < 0, rbinom(1, abs(rnd.tab[2,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[2,3]), abs(rnd.tab[2,3]) - mam.less, 0)
  rnd.tab[2,2] <- rnd.tab[2,2] - mam.less
  rnd.tab[2,1] <- ifelse(rnd.tab[2,3] > 0, rnd.tab[2,1], rnd.tab[2,1] - ave.less)
  rnd.tab[2,3] <- ifelse(rnd.tab[2,3] < 0, 0, rnd.tab[2,3])
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],classSums[1]/sum(classSums))
  rnd.tab[3,2] <- rbinom(1,contSums[3],classSums[2]/sum(classSums))
  rnd.tab[3,3] <- contSums[3] - sum(rnd.tab[3,1], rnd.tab[3,2])
  mam.less <- ifelse(rnd.tab[3,3] < 0, rbinom(1, abs(rnd.tab[3,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[3,3]), abs(rnd.tab[3,3]) - mam.less, 0)
  rnd.tab[3,2] <- rnd.tab[3,2] - mam.less
  rnd.tab[3,1] <- ifelse(rnd.tab[3,3] > 0, rnd.tab[3,1], rnd.tab[3,1] - ave.less)
  rnd.tab[3,3] <- ifelse(rnd.tab[3,3] < 0, 0, rnd.tab[3,3])
  
  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],classSums[1]/sum(classSums))
  rnd.tab[4,2] <- rbinom(1,contSums[4],classSums[2]/sum(classSums))
  rnd.tab[4,3] <- contSums[4] - sum(rnd.tab[4,1], rnd.tab[4,2])
  mam.less <- ifelse(rnd.tab[4,3] < 0, rbinom(1, abs(rnd.tab[4,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[4,3]) - mam.less, 0)
  rnd.tab[4,2] <- rnd.tab[4,2] - mam.less
  rnd.tab[4,1] <- ifelse(rnd.tab[4,3] > 0, rnd.tab[4,1], rnd.tab[4,1] - ave.less)
  rnd.tab[4,3] <- ifelse(rnd.tab[4,3] < 0, 0, rnd.tab[4,3])
  
  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],classSums[1]/sum(classSums))
  rnd.tab[5,2] <- rbinom(1,contSums[5],classSums[2]/sum(classSums))
  rnd.tab[5,3] <- contSums[5] - sum(rnd.tab[5,1], rnd.tab[5,2])
  mam.less <- ifelse(rnd.tab[5,3] < 0, rbinom(1, abs(rnd.tab[5,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[5,3]) - mam.less, 0)
  rnd.tab[5,2] <- rnd.tab[5,2] - mam.less
  rnd.tab[5,1] <- ifelse(rnd.tab[5,3] > 0, rnd.tab[5,1], rnd.tab[5,1] - ave.less)
  rnd.tab[5,3] <- ifelse(rnd.tab[5,3] < 0, 0, rnd.tab[5,3])
  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  inf.sub <- which(is.infinite(chi.tab.rnd)==T, arr.ind=T)
  if (length(inf.sub) > 0) {
    chi.tab.rnd[inf.sub] <- 0
    chi.tab.adj[inf.sub] <- 0
  } # end if  
  
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)
  comp.vec[c] <- ifelse(chi.rnd[c] >= sum(chi.tab.adj, na.rm=T), 1, 0) 
  
  if (c %% itdiv==0) print(c)
  
} # end c
pr <- sum(comp.vec)/iter
pr # high type I error for hypothesis that distribution of species by major class differs by continent


## conservation status
## simulation test of independence
head(spp.dat)
threat.sub <- spp.dat[,which(colnames(spp.dat) %in% c("EX","EW","CR","EN","VU","NT","LC","DD","NE"))]
dim(threat.sub)
threat.stat <- which(threat.sub == 1, arr.ind=T)
nothreat.sub <- which(rowSums(threat.sub, na.rm=T) == 0)
threat.stat.ord <- threat.stat[order(threat.stat[,1],decreasing=F),]
threat <- ifelse(threat.stat.ord[,2] == 1, "EX", 
                 ifelse(threat.stat.ord[,2] == 2, "EW", 
                        ifelse(threat.stat.ord[,2] == 3, "CR", 
                               ifelse(threat.stat.ord[,2] == 4, "EN", 
                                      ifelse(threat.stat.ord[,2] == 5, "VU", 
                                             ifelse(threat.stat.ord[,2] == 6, "NT", 
                                                    ifelse(threat.stat.ord[,2] == 7, "LC", 
                                                           ifelse(threat.stat.ord[,2] == 8, "DD",
                                                                  ifelse(threat.stat.ord[,2] == 9, "NE", NA)))))))))
dim(spp.dat)
length(threat)
spp.red <- spp.dat[-nothreat.sub,]
spp.red$IUCN <- threat
spp.red$threat <- ifelse(spp.red$IUCN == "CR" | spp.red$IUCN == "EN" | spp.red$IUCN == "VU", 1, 0)
threat.tab <- table(spp.red$Continent, spp.red$threat)
threat.tab
threat.tab[,1]/rowSums(threat.tab)
contSums <- rowSums(threat.tab)
threatSums <- colSums(threat.tab)

# exact binomial test
binom.test(threatSums[2], sum(threatSums), 1/2, alternative="two.sided")

exp.tab <- threat.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- threatSums[1] * (contSums[i]/sum(threat.tab))
  exp.tab[i,2] <- threatSums[2] * (contSums[i]/sum(threat.tab))
}
exp.tab
rowSums(exp.tab)
chi.tab <- (threat.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chisq.test(threat.tab, simulate.p.value=T, B=10000)

iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(NA,iter)

for (c in 1:iter) {
  rnd.tab <- threat.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],threatSums[1]/sum(threatSums))
  rnd.tab[1,2] <- contSums[1] - rnd.tab[1,1]
  
  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],threatSums[1]/sum(threatSums))
  rnd.tab[2,2] <- contSums[2] - rnd.tab[2,1]
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],threatSums[1]/sum(threatSums))
  rnd.tab[3,2] <- contSums[3] - rnd.tab[3,1]
  
  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],threatSums[1]/sum(threatSums))
  rnd.tab[4,2] <- contSums[4] - rnd.tab[4,1]
  
  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],threatSums[1]/sum(threatSums))
  rnd.tab[5,2] <- contSums[5] - rnd.tab[5,1]
  
  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)
  
  if (c %% itdiv==0) print(c)
  
} # end c
pr <- length(which(chi.rnd >= chi.obs))/iter
pr # low type I error for hypothesis that distribution of species by major threat category (threatened or not) differs by continent


## population trend
head(spp.dat)
trend.sub <- spp.dat[,which(colnames(spp.dat) %in% c("Decreasing","Increasing","Stable","Unknown"))]
head(trend.sub)
dim(trend.sub)
notrend.sub <- which(rowSums(trend.sub, na.rm=T) == 0)
trend.stat <- which(trend.sub == 1, arr.ind=T)
trend.stat.ord <- trend.stat[order(trend.stat[,1],decreasing=F),]
trend <- ifelse(trend.stat.ord[,2] == 1, "decr", 
                ifelse(trend.stat.ord[,2] == 2, "incr", 
                       ifelse(trend.stat.ord[,2] == 3, "stab", 
                              ifelse(trend.stat.ord[,2] == 4, "unkn", NA))))
trend
dim(spp.dat)
length(trend)
spp.red <- spp.dat[-notrend.sub,]
spp.red$trend <- trend
spp.red2 <- spp.red[spp.red$trend != "unkn", ]
dim(spp.red2)
trend.tab <- table(spp.red2$Continent, spp.red2$trend)
trend.tab
contSums <- rowSums(trend.tab)
trendSums <- colSums(trend.tab)

# exact binomial test
binom.test(trendSums[1], sum(threatSums), 1/3, alternative="greater")

exp.tab <- trend.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- trendSums[1] * (contSums[i]/sum(trend.tab))
  exp.tab[i,2] <- trendSums[2] * (contSums[i]/sum(trend.tab))
  exp.tab[i,3] <- trendSums[3] * (contSums[i]/sum(trend.tab))
}
exp.tab
trend.tab
chi.tab <- (trend.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chi.obs
chisq.test(trend.tab, simulate.p.value=T, B=10000)


iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(0,iter)
comp.vec <- rep(0,iter)

for (c in 1:iter) {
  rnd.tab <- trend.tab
  chi.tab.adj <- chi.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],trendSums[1]/sum(trendSums))
  rnd.tab[1,2] <- rbinom(1,contSums[1],trendSums[2]/sum(trendSums))
  rnd.tab[1,3] <- contSums[1] - sum(rnd.tab[1,1], rnd.tab[1,2])
  mam.less <- ifelse(rnd.tab[1,3] < 0, rbinom(1, abs(rnd.tab[1,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[1,3]), abs(rnd.tab[1,3]) - mam.less, 0)
  rnd.tab[1,2] <- rnd.tab[1,2] - mam.less
  rnd.tab[1,1] <- ifelse(rnd.tab[1,3] > 0, rnd.tab[1,1], rnd.tab[1,1] - ave.less)
  rnd.tab[1,3] <- ifelse(rnd.tab[1,3] < 0, 0, rnd.tab[1,3])
  
  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],trendSums[1]/sum(trendSums))
  rnd.tab[2,2] <- rbinom(1,contSums[2],trendSums[2]/sum(trendSums))
  rnd.tab[2,3] <- contSums[2] - sum(rnd.tab[2,1], rnd.tab[2,2])
  mam.less <- ifelse(rnd.tab[2,3] < 0, rbinom(1, abs(rnd.tab[2,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[2,3]), abs(rnd.tab[2,3]) - mam.less, 0)
  rnd.tab[2,2] <- rnd.tab[2,2] - mam.less
  rnd.tab[2,1] <- ifelse(rnd.tab[2,3] > 0, rnd.tab[2,1], rnd.tab[2,1] - ave.less)
  rnd.tab[2,3] <- ifelse(rnd.tab[2,3] < 0, 0, rnd.tab[2,3])
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],trendSums[1]/sum(trendSums))
  rnd.tab[3,2] <- rbinom(1,contSums[3],trendSums[2]/sum(trendSums))
  rnd.tab[3,3] <- contSums[3] - sum(rnd.tab[3,1], rnd.tab[3,2])
  mam.less <- ifelse(rnd.tab[3,3] < 0, rbinom(1, abs(rnd.tab[3,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[3,3]), abs(rnd.tab[3,3]) - mam.less, 0)
  rnd.tab[3,2] <- rnd.tab[3,2] - mam.less
  rnd.tab[3,1] <- ifelse(rnd.tab[3,3] > 0, rnd.tab[3,1], rnd.tab[3,1] - ave.less)
  rnd.tab[3,3] <- ifelse(rnd.tab[3,3] < 0, 0, rnd.tab[3,3])
  
  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],trendSums[1]/sum(trendSums))
  rnd.tab[4,2] <- rbinom(1,contSums[4],trendSums[2]/sum(trendSums))
  rnd.tab[4,3] <- contSums[4] - sum(rnd.tab[4,1], rnd.tab[4,2])
  mam.less <- ifelse(rnd.tab[4,3] < 0, rbinom(1, abs(rnd.tab[4,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[4,3]) - mam.less, 0)
  rnd.tab[4,2] <- rnd.tab[4,2] - mam.less
  rnd.tab[4,1] <- ifelse(rnd.tab[4,3] > 0, rnd.tab[4,1], rnd.tab[4,1] - ave.less)
  rnd.tab[4,3] <- ifelse(rnd.tab[4,3] < 0, 0, rnd.tab[4,3])
  
  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],trendSums[1]/sum(trendSums))
  rnd.tab[5,2] <- rbinom(1,contSums[5],trendSums[2]/sum(trendSums))
  rnd.tab[5,3] <- contSums[5] - sum(rnd.tab[5,1], rnd.tab[5,2])
  mam.less <- ifelse(rnd.tab[5,3] < 0, rbinom(1, abs(rnd.tab[5,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[5,3]) - mam.less, 0)
  rnd.tab[5,2] <- rnd.tab[5,2] - mam.less
  rnd.tab[5,1] <- ifelse(rnd.tab[5,3] > 0, rnd.tab[5,1], rnd.tab[5,1] - ave.less)
  rnd.tab[5,3] <- ifelse(rnd.tab[5,3] < 0, 0, rnd.tab[5,3])
  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  inf.sub <- which(is.infinite(chi.tab.rnd)==T, arr.ind=T)
  if (length(inf.sub) > 0) {
    chi.tab.rnd[inf.sub] <- 0
    chi.tab.adj[inf.sub] <- 0
  } # end if  
  
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)
  comp.vec[c] <- ifelse(chi.rnd[c] >= sum(chi.tab.adj, na.rm=T), 1, 0)
  
  if (c %% itdiv==0) print(c)
  
} # end c
pr <- sum(comp.vec)/iter
pr # low type I error for hypothesis that distribution of species by major trend category differs by continent


## spatial congruence
## generalised linear mixed-effects models
# model set
m1 <- "N_Species_In_Country ~ ClassComb + (1|Continent)"
m2 <- "N_Species_In_Country ~ 1 + (1|Continent)"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:2,1),sep="")))

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),data=spp.dat, family=binomial(link = "logit"), na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

# Table S2
sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table

## evidence for variation in proportion of species in country by class and continent

plot_model(mod.list[[1]], show.values=T, type="re", vline.color = "purple")
plot_model(mod.list[[1]], show.values=T, type="est", vline.color = "purple")
summary(mod.list[[1]])
probs <- plogis(predict(mod.list[[1]]))
missdatsub <- which(is.na(spp.dat$N_Species_In_Country))

# predicted probabilities by category
probs.dat <- data.frame("cont"=spp.dat$Continent[-missdatsub], "Class"=spp.dat$Classification_Class[-missdatsub],
                        "ClassComb"=spp.dat$ClassComb[-missdatsub], "probs"=probs)

# by continent
EUR.mn <- mean(probs.dat$probs[probs.dat$cont=="Europe"], na.rm=T)
EUR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Europe"], probs=0.025, na.rm = T))
EUR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Europe"], probs=0.975, na.rm = T))
EUR.sd <- sd(probs.dat$probs[probs.dat$cont=="Europe"], na.rm=T)
print(c(EUR.lo, EUR.mn, EUR.up))

AMR.mn <- mean(probs.dat$probs[probs.dat$cont=="Americas"], na.rm=T)
AMR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Americas"], probs=0.025, na.rm = T))
AMR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Americas"], probs=0.975, na.rm = T))
AMR.sd <- sd(probs.dat$probs[probs.dat$cont=="Americas"], na.rm=T)
print(c(AMR.lo, AMR.mn, AMR.up))

OCN.mn <- mean(probs.dat$probs[probs.dat$cont=="Oceania"], na.rm=T)
OCN.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Oceania"], probs=0.025, na.rm = T))
OCN.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Oceania"], probs=0.975, na.rm = T))
OCN.sd <- sd(probs.dat$probs[probs.dat$cont=="Oceania"], na.rm=T)
print(c(OCN.lo, OCN.mn, OCN.up))

ASIA.mn <- mean(probs.dat$probs[probs.dat$cont=="Asia"], na.rm=T)
ASIA.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Asia"], probs=0.025, na.rm = T))
ASIA.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Asia"], probs=0.975, na.rm = T))
ASIA.sd <- sd(probs.dat$probs[probs.dat$cont=="Asia"], na.rm=T)
print(c(ASIA.lo, ASIA.mn, ASIA.up))

AFR.mn <- mean(probs.dat$probs[probs.dat$cont=="Africa"], na.rm=T)
AFR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Africa"], probs=0.025, na.rm = T))
AFR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Africa"], probs=0.975, na.rm = T))
AFR.sd <- sd(probs.dat$probs[probs.dat$cont=="Africa"], na.rm=T) 
print(c(AFR.lo, AFR.mn, AFR.up))

prob.mns <- c(EUR.mn, AMR.mn, OCN.mn, ASIA.mn, AFR.mn)
prob.los <- c(EUR.lo, AMR.lo, OCN.lo, ASIA.lo, AFR.lo)
prob.ups <- c(EUR.up, AMR.up, OCN.up, ASIA.up, AFR.up)
prob.sds <- c(EUR.sd, AMR.sd, OCN.sd, ASIA.sd, AFR.sd)

prob.cont.summ <- data.frame("mean"=prob.mns, "sd"=prob.sds, "lo"=prob.los, "up"=prob.ups)
rownames(prob.cont.summ) <- c("EUR", "AMR", "OCN", "ASIA", "AFR")
prob.cont.summ

prob.cont.plot <- ggplot(prob.cont.summ) +
  geom_bar(aes(x=reorder(row.names(prob.cont.summ), mean), y=mean), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(prob.cont.summ), ymin=mean-sd, ymax=mean+sd),
                 linewidth=0.4, colour="black", alpha=0.9, size=0.3)
prob.cont.plot + coord_flip() +
  xlab("continent") + ylab("probability")

## probability of species chosen being in country increasing EUR < AMR < OCN < ASIA < AFR

## To complement Figure 4
prob.cont.plot <- ggplot(prob.cont.summ) +
  geom_point(aes(x=reorder(row.names(prob.cont.summ), mean), y=mean), stat="identity", fill="blue") +
  geom_errorbar( aes(x=row.names(prob.cont.summ), ymin=mean-sd, ymax=mean+sd),
                 linewidth=0.4, colour="black", alpha=0.9, size=0.3, width=0)+
  coord_flip()+
  theme_classic()+
  theme(axis.text=element_text(size=12))+
  geom_hline(yintercept=c(0.3,0.4,0.5,0.6), color="darkgrey", lty=2)+
  xlab("continent") + ylab("probability")

grid.arrange(presabsplot, prob.cont.plot, nrow=1, widths = c(3,1))

# by Class
MAM.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], na.rm=T)
MAM.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], probs=0.025, na.rm = T))
MAM.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], probs=0.975, na.rm = T))
MAM.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], na.rm=T) 
print(c(MAM.lo, MAM.mn, MAM.up))

AVE.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="Aves"], na.rm=T)
AVE.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Aves"], probs=0.025, na.rm = T))
AVE.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Aves"], probs=0.975, na.rm = T))
AVE.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="Aves"], na.rm=T) 
print(c(AVE.lo, AVE.mn, AVE.up))

OTH.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="other"], na.rm=T)
OTH.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="other"], probs=0.025, na.rm = T))
OTH.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="other"], probs=0.975, na.rm = T))
OTH.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="other"], na.rm=T) 
print(c(OTH.lo, OTH.mn, OTH.up))

probcl.mns <- c(MAM.mn, AVE.mn, OTH.mn)
probcl.los <- c(MAM.lo, AVE.lo, OTH.lo)
probcl.ups <- c(MAM.up, AVE.up, OTH.up)
probcl.sds <- c(MAM.sd, AVE.sd, OTH.sd)

prob.class.summ <- data.frame("mean"=probcl.mns, "sd"=probcl.sds, "lo"=probcl.los, "up"=probcl.ups)
rownames(prob.class.summ) <- c("MAM", "AVE", "OTH")
prob.class.summ

prob.cont.plot <- ggplot(prob.class.summ) +
  geom_bar(aes(x=reorder(row.names(prob.class.summ), mean), y=mean), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(prob.class.summ), ymin=mean-sd, ymax=mean+sd),
                 linewidth=0.4, colour="black", alpha=0.9, size=0.3)
prob.cont.plot + coord_flip() +
  xlab("class") + ylab("probability")

## lower representation of mammals within country than Aves or 'other'

