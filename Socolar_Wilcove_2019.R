# This code performs the analysis and produces the figures reported in Socolar & Wilcove (2019) INSERT CITATION
# The script can be run in RStudio 'as is' with no modification. If not using RStudio, omit lines 19-20 and set 
# the working directory to the folder on your computer where the .zip file was unarchived (i.e. the parent
# directory of the "code" folder where this script resides).

# HMC fitting in rstanarm is stochastic, and re-running the script will produce very minor variation in model
# results. To analyze and plot the exact data in the paper, omit the `Model Fitting` and `Cross Validation` 
# sections (lines 332-400) and include (un-comment) lines 404-405 and lines 433-434. This will load stable 
# versions of the fitted model and cross-validation objects from the fData directory. Even if run in full, this 
# script will not overwrite the data objects residing in the fData directory.

# WARNING: as written, the model fitting and cross-validation will execute parallel processes on 4 cores. The
# total runtime may exceed an hour, depending on the computer. 

# Note: this code uses the term "forest-present" as a synonym of the main text's term "forest-user".

library(rstanarm)
library(doBy)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

##### User-defined functions #####
'%ni%' <- Negate('%in%')

##### Data import, cleaning and manipulation #####
# Load data from Socolar et al 2019, available at https://datadryad.org/bitstream/handle/10255/dryad.220673/birds.csv?sequence=1
download.file("https://datadryad.org/bitstream/handle/10255/dryad.220673/birds.csv?sequence=1",
              destfile = 'dData/birds.csv')
iquitos_data <- read.csv('dData/birds.csv', stringsAsFactors = F)

# Load data from Parker et al 1996.
download.file("https://ndownloader.figshare.com/files/16378250", destfile = 'dData/Parker_Stotz_Fitzpatrick_1996.zip')
unzip('dData/Parker_Stotz_Fitzpatrick_1996.zip', exdir = 'dData/Parker_Stotz_Fitzpatrick_1996')

parker <- read.csv('dData/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv')
parker.nb <- read.csv('dData/Parker_Stotz_Fitzpatrick_1996/databases/cdata.csv')

# Load data file for taxa that do not appear in Parker et al (splits and newly-described species)
extras <- read.csv("Data/parker_additions.csv")

# Assemble single dataframe
parker <- gtools::smartbind(parker, parker.nb)
parker <- gtools::smartbind(parker, extras)
parker$GENUS <- as.character(parker$GENUS)

# Update taxonomy
parker$GENUS[parker$GENUS == "Columba"] <- "Patagioenas"
parker$sp <- paste(parker$GENUS, parker$SPECIES, sep = "_")

parker$sp[grep("Capito_niger", parker$sp)] <- "Capito_auratus"
parker$sp[grep("Pitylus_grossus", parker$sp)] <- "Saltator_grossus"
parker$sp[grep("Xiphorhynchus_\\(guttatus\\) guttatus", parker$sp)] <- "Xiphorhynchus_guttatus"
parker$sp[grep("Momotus_\\(momota\\) momota", parker$sp)] <- "Momotus_momota"
parker$sp[grep("Myrmotherula_\\(brachyura\\) ignota", parker$sp)] <- "Myrmotherula_ignota"
parker$sp[grep("Pionites_melanocephala", parker$sp)] <- "Pionites_melanocephalus"
parker$sp[grep("Pipra_coronota", parker$sp)] <- "Lepidothrix_coronata"
parker$sp[grep("Hypocnemis_cantator", parker$sp)] <- "Hypocnemis_peruviana"
parker$sp[grep("Ortalis_\\(motmot\\) guttata", parker$sp)] <- "Ortalis_guttata"
parker$sp[grep("Hylophylax_poecilinota", parker$sp)] <- "Willisornis_poecilinotus"
parker$sp[grep("Tangara_\\(velia\\) velia", parker$sp)] <- "Tangara_velia"
parker$sp[grep("Psarocolius_latirostris", parker$sp)] <- "Cacicus_latirostris"
parker$sp[grep("Myrmotherula_hauxwelli", parker$sp)] <- "Isleria_hauxwelli"
parker$sp[grep("Schiffornis_turdinus", parker$sp)] <- "Schiffornis_turdina"
parker$sp[grep("Tolmomyias_\\(assimilis\\) assimilis", parker$sp)] <- "Tolmomyias_assimilis"
parker$sp[grep("Pipra_pipra", parker$sp)] <- "Dixiphia_pipra"
parker$sp[grep("Oryzoborus_\\(angolensis\\) angolensis", parker$sp)] <- "Oryzoborus_angolensis"
parker$sp[grep("Thamnophilus_\\(doliatus\\) doliatus", parker$sp)] <- "Thamnophilus_doliatus"
parker$sp[grep("Brotogeris_\\(versicolurus\\) versicolurus", parker$sp)] <- "Brotogeris_versicolurus"
parker$sp[grep("Todirostrum_latirostre", parker$sp)] <- "Poecilotriccus_latirostris"
parker$sp[grep("Buteo_magnirostris", parker$sp)] <- "Rupornis_magnirostris"
parker$sp[grep("Xiphorhynchus_picus", parker$sp)] <- "Dendroplex_picus"
parker$sp[grep("Glaucis_hirsuta", parker$sp)] <- "Glaucis_hirsutus"
parker$sp[grep("Trogon_violaceus", parker$sp)] <- "Trogon_ramonianus"
parker$sp[grep("Jacamerops_aurea", parker$sp)] <- "Jacamerops_aureus"
parker$sp[grep("Psarocolius_oseryi", parker$sp)] <- "Clypicterus_oseryi"
parker$sp[grep("Leptotila_\\(rufaxilla\\) rufaxilla", parker$sp)] <- "Leptotila_rufaxilla"
parker$sp[grep("Hylophylax_punctulata", parker$sp)] <- "Hylophylax_punctulatus"
parker$sp[grep("Tolmomyias_\\(flaviventris\\) flaviventris", parker$sp)] <- "Tolmomyias_flaviventris"
parker$sp[grep("Myrmotherula_\\(brachyura\\) brachyura", parker$sp)] <- "Myrmotherula_brachyura"
parker$sp[grep("Donacobius_atricapillus", parker$sp)] <- "Donacobius_atricapilla"
parker$sp[grep("Nonnula_\\(ruficapilla\\) ruficapilla", parker$sp)] <- "Nonnula_ruficapilla"
parker$sp[grep("Amazona_\\(ochrocephala\\) ochrocephala", parker$sp)] <- "Amazona_ochrocephala"
parker$sp[grep("Thryothorus_\\(genibarbis\\) genibarbis", parker$sp)] <- "Thryothorus_genibarbis"
parker$sp[grep("Polioptila_\\(plumbea\\) plumbea", parker$sp)] <- "Polioptila_plumbea"
parker$sp[grep("Columbina_\\(talpacoti\\) talpacoti", parker$sp)] <- "Columbina_talpacoti"
parker$sp[grep("Sporophila_\\(lineola\\) bouvronides", parker$sp)] <- "Sporophila_buvronides"
parker$sp[grep("Leucopternis_schistacea", parker$sp)] <- "Buteogallus_schistaceus"
parker$sp[grep("Vireo_\\(olivaceus\\) chivi", parker$sp)] <- "Vireo_chivi"
parker$sp[grep("Phylloscartes_flaveolus", parker$sp)] <- "Capsiempis_flaveola"
parker$sp[grep("Xiphorhynchus_\\(spixii\\) elegans", parker$sp)] <- "Xiphorhynchus_elegans"
parker$sp[grep("Frederickena_unduligera", parker$sp)] <- "Frederickena_unduliger"
parker$sp[grep("Hylophylax_naevia", parker$sp)] <- "Hylophylax_naevius"
parker$sp[grep("Conopias_\\(albovittata\\) parva", parker$sp)] <- "Conopias_parvus"
parker$sp[grep("Daptrius_americanus", parker$sp)] <- "Ibycter_americanus"
parker$sp[grep("Myrmotherula_\\(haematonota\\) haematonota", parker$sp)] <- "Epinecrophylla_haematonota"
parker$sp[grep("Percnostola_leucostigma", parker$sp)] <- "Schistocichla_leucostigma"
parker$sp[grep("Troglodytes_\\(aedon\\) aedon", parker$sp)] <- "Troglodytes_aedon"
parker$sp[grep("Dendrocincla_\\(fuliginosa\\) fuliginosa", parker$sp)] <- "Dendrocincla_fuliginosa"
parker$sp[grep("Percnostola_schistacea", parker$sp)] <- "Schistocichla_schistacea"
parker$sp[grep("Herpsilochmus_\\(sticturus\\) dugandi", parker$sp)] <- "Herpsilochmus_dugandi"
parker$sp[grep("Notharchus_macrorhynchos", parker$sp)] <- "Notharchus_hyperrhynchus"
parker$sp[grep("Ceryle_torquata", parker$sp)] <- "Megaceryle_torquata"
parker$sp[grep("Philydor_erythropterus", parker$sp)] <- "Philydor_erythropterum"
parker[nrow(parker)+1, ] <- parker[which(parker$sp == "Frederickena_unduliger"), ]
parker$sp[nrow(parker)] <- "Frederickena_fulva"
parker$sp[grep("Philydor_erythropterus", parker$sp)] <- "Philydor_erythropterum"
parker$sp[grep("barrabandi", parker$sp)] <- "Pyrilia_barrabandi"
parker$sp[grep("Pyrrhura_\\(melanura\\) melanura", parker$sp)] <- "Pyrrhura_melanura"
parker$sp[grep("Patagioenas_livia", parker$sp)] <- "Columba_livia"
parker$sp[grep("Chondrohierax_\\(uncinatus\\) uncinatus", parker$sp)] <- "Chondrohierax_uncinatus"
parker$sp[grep("Butorides_\\(striatus\\) striatus", parker$sp)] <- "Butorides_striatus"
parker$sp[grep("barrabandi", parker$sp)] <- "Pyrilia_barrabandi"
parker$sp[grep("Rostrhamus_hamatus", parker$sp)] <- "Helicolestes_hamatus"
parker$sp[grep("Furnarius_\\(leucopus\\) leucopus", parker$sp)] <- "Furnarius_leucopus"
parker$sp[grep("Forpus_crassirostris", parker$sp)] <- "Forpus_xanthopterygius"
parker$sp[grep("Xiphorhynchus_necopinus", parker$sp)] <- "Dendroplex_kienerii"
parker$sp[grep("Aratinga_leucophthalmus", parker$sp)] <- "Psittacara_leucophthalma"
parker$sp[grep("Lepidocolaptes_albolineatus", parker$sp)] <- "Lepidocolaptes_fatimalimae"
parker[nrow(parker)+1, ] <- parker[which(parker$sp == "Lepidocolaptes_fatimalimae"), ]
parker$sp[nrow(parker)] <- "Lepidocolaptes_duidae"
parker$sp[grep("Dendrocolaptes_\\(certhia\\) certhia", parker$sp)] <- "Dendrocolaptes_certhia"
parker$sp[grep("Tangara_\\(mexicana\\) mexicana", parker$sp)] <- "Tangara_mexicana"
parker$sp[grep("Pyrrhura_\\(picta\\) picta", parker$sp)] <- "Pyrrhura_roseifrons"
parker[nrow(parker)+1, ] <- parker[which(parker$sp == "Pyrrhura_roseifrons"), ]
parker$sp[nrow(parker)] <- "Pyrrhura_lucianni"
parker$sp[grep("Onychorhynchus_\\(coronatus\\) coronatus", parker$sp)] <- "Onychorhynchus_coronatus"
parker$sp[grep("Gymnopithys_lunulata", parker$sp)] <- "Gymnopithys_lunulatus"
parker$sp[grep("Passerina_cyanoides", parker$sp)] <- "Cyanocompsa_cyanoides"
parker$sp[grep("Agelaius_icterocephalus", parker$sp)] <- "Chrysomus_icterocephalus"
parker$sp[grep("Porphyrula_martinica", parker$sp)] <- "Porphyrio_martinica"
parker$sp[grep("Chiroxiphia_\\(pareola\\) pareola", parker$sp)] <- "Chiroxiphia_pareola"
parker$sp[grep("Myrmotherula_erythrura", parker$sp)] <- "Epinecrophylla_erythrura"
parker$sp[grep("Sclerurus_mexicanus", parker$sp)] <- "Sclerurus_obscurior"
parker$sp[grep("Ardeola_ibis", parker$sp)] <- "Bubulcus_ibis"
parker$sp[grep("Scaphidura_oryzivora", parker$sp)] <- "Molothrus_oryzivorus"
parker$sp[grep("Chlorestes_notatus", parker$sp)] <- "Chlorestes_notata"
parker$sp[grep("Xenops_milleri", parker$sp)] <- "Microxenops_milleri"
parker$sp[grep("Synallaxis_\\(gujanensis\\) gujanensis", parker$sp)] <- "Synallaxis_gujanensis"
parker$sp[grep("Certhiaxis_mustelina", parker$sp)] <- "Certhiaxis_mustelinus"
parker$sp[grep("Metopothrix_aurantiacus", parker$sp)] <- "Metopothrix_aurantiaca"
parker$sp[grep("Phaethornis_longuemareus", parker$sp)] <- "Phaethornis_atrimentalis"
parker$sp[grep("Xiphocolaptes_\\(promeropirhynchus\\) orenocensis", parker$sp)] <- "Xiphocolaptes_promeropirhynchus"
parker$sp[grep("Sporophila_\\(lineola\\) lineola", parker$sp)] <- "Sporophila_lineola"
parker$sp[grep("Myiodynastes_\\(maculatus\\) solitarius", parker$sp)] <- "Myiodynastes_maculatus"

# correct iquitos_data taxonomy where problematic:
# remove Trogon melanurus and Hylophylax naevius
iquitos_data <- iquitos_data[-grep("Trogon_melanurus", iquitos_data$Species), ]
iquitos_data <- iquitos_data[-grep("Hylophylax_naevius", iquitos_data$Species), ]
# correct spellings
iquitos_data$Species[which(iquitos_data$Species == "Odonotophorus_gujanensis")] <- "Odontophorus_gujanensis"
iquitos_data$Species[which(iquitos_data$Species == "Brotogeris_versicolorus")] <- "Brotogeris_versicolurus"
iquitos_data$Species[which(iquitos_data$Species == "Megarhynchus_pitangua")] <- "Megarynchus_pitangua"
iquitos_data$Species[which(iquitos_data$Species == "Vireo_olivaceus")] <- "Vireo_chivi"
iquitos_data$Species[which(iquitos_data$Species == "Galbula_flavirostris")] <- "Galbula_albirostris"
iquitos_data$Species[which(iquitos_data$Species == "Ancistrops_strigulatus")] <- "Ancistrops_strigilatus"
iquitos_data$Species[which(iquitos_data$Species == "Gymnoderis_foetidus")] <- "Gymnoderus_foetidus"
iquitos_data$Species[which(iquitos_data$Species == "Todirostrum_chrysocephalum")] <- "Todirostrum_chrysocrotaphum"
iquitos_data$Species[which(iquitos_data$Species == "Dromococcyx_phaisanellus")] <- "Dromococcyx_phasianellus"
iquitos_data$Species[which(iquitos_data$Species == "Loreto_antwren")] <- "Herpsilochmus_loreto"
iquitos_data$Species[which(iquitos_data$Species == "Attila_spadecius")] <- "Attila_spadiceus"

write.csv(parker, file = "Data/parker_iquitos_data_taxonomy.csv")
write.csv(iquitos_data, file = "Data/iquitos_data_updated.csv")

# Aggregate habitats in Parker et al (1996) database
parker$flood.for <- parker$F2 %in% c("Y", "Q") | parker$F3 %in% c("Y", "Q") | parker$F13 %in% c("Y", "Q")
parker$flood.nfor <- parker$N11 %in% c("Y", "Q") | parker$N12 %in% c("Y", "Q") | parker$A1 %in% c("Y", "Q") | parker$A5 %in% c("Y", "Q") | 
  parker$A6 %in% c("Y", "Q") | parker$A8 %in% c("Y", "Q") | parker$A9 %in% c("Y", "Q") | parker$F2 == "E" | parker$F3 == "E"
parker$flood.all <- parker$flood.for | parker$flood.nfor
parker$tf <- parker$F1 %in% c("Y", "Q") | parker$F12 %in% c("Y", "Q")
parker$forest.present <- (parker$flood.for == 1) | (parker$tf == 1)
parker$upland <- parker$F1 %in% c("Y", "Q")

# Get lists of specialist species
# forest-user species
FBsp <- as.character(parker$sp[parker$forest.present])
cFBsp <- FBsp[FBsp %in% iquitos_data$Species]

# forest-specialist species
FSsp <- as.character(parker$sp[parker$forest.present & !parker$flood.nfor])
cFSsp <- FSsp[FSsp %in% iquitos_data$Species]

# forest-based floodplain specialists
FFsp <- as.character(parker$sp[parker$flood.for & !parker$tf])
cFFsp <- FFsp[FFsp %in% iquitos_data$Species]

# floodplain specialists that occur in non-forest (may overlap with forest-based species)
FNsp <- as.character(parker$sp[parker$flood.nfor & !parker$tf])
cFNsp <- FNsp[FNsp %in% iquitos_data$Species]

# all floodplain specialists
allfloodsp <- unique(c(cFFsp, cFNsp))

# terra firme specialists
TFsp <- as.character(parker$sp[parker$tf & !(parker$flood.for | parker$flood.nfor)])
cTFsp <- TFsp[TFsp %in% iquitos_data$Species]

# Species classified as terra firme specialists based on Parker et al (1996) but
# reclassified here as non-specialists
tfexclude <- read.csv("Data/tf_exclude.csv", header = F)
cTFsp <- cTFsp[cTFsp %ni% tfexclude[,1]]

# Species not classified as floodplain specialists based on Parker et al (1996) but
# reclassified here as floodplain specialists.
floodextras <- read.csv("Data/flood_specialist_extras.csv")

# rich-soil, poor-soil, river-limitation, and migratory designations, compiled from
# Pomara et al (2012), Alvarez Alonso et al (2013), and Schulenberg et al (2010)
spec <- read.csv("Data/specialists.csv")
spec[is.na(spec)] <- 0

##### build dataframe for modeling
specialists <- data.frame(species=spec$Species, poor=((spec$WhiteSand.Alvarez+spec$Poorsoil.Pomara+spec$Poorsoil.extra)>0),
                          rich=((spec$Richsoil.Alvarez + spec$Richsoil.Pomara)>0), river = spec$RiverLimit, migratory = spec$Mig)
specialists$flood <- 0
specialists$flood[specialists$species %in% allfloodsp] <- 1
specialists$flood[specialists$species %in% floodextras$Extras] <- 1
specialists$tf <- 0
specialists$tf[specialists$species %in% cTFsp] <- 1
specialists$forest.present <- 0
specialists$forest.present[specialists$species %in% cFBsp] <- 1
specialists$forest.specialist <- 0
specialists$forest.specialist[specialists$species %in% cFSsp] <- 1

##### Import Eltontraits data
download.file("https://ndownloader.figshare.com/files/5631081",
              destfile = 'dData/elton.txt')
traits <- read.delim(paste0('dData/elton.txt'), header=T, stringsAsFactors = F)
traits <- traits[1:9993,]
traits$sciName <- gsub(" ", "_", traits$Scientific)

# Handle taxonomic discrepancies 
traits$sciName[traits$sciName == "Ocyalus_latirostris"] <- "Cacicus_latirostris"
traits$sciName[traits$sciName == "Myrmotherula_hauxwelli"] <- "Isleria_hauxwelli"
traits$sciName[traits$sciName == "Pipra_pipra"] <- "Dixiphia_pipra"
traits$sciName[traits$sciName == "Buteo_magnirostris"] <- "Rupornis_magnirostris"
traits$sciName[traits$sciName == "Cissopis_leverianus"] <- "Cissopis_leveriana"
traits$sciName[traits$sciName == "Trogon_violaceus"] <- "Trogon_ramonianus"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Xiphorhynchus_ocellatus", ]
traits$sciName[nrow(traits)] <- "Xiphorhynchus_chunchotambo"
traits$sciName[traits$sciName == "Corythopis_torquatus"] <- "Corythopis_torquata"
traits$sciName[traits$sciName == "Sporophila_bouvronides"] <- "Sporophila_buvronides"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Turdus_hauxwelli", ]
traits$sciName[nrow(traits)] <- "Turdus_sanchezorum"
traits$sciName[traits$sciName == "Leucopternis_schistaceus"] <- "Buteogallus_schistaceus"
traits$sciName[traits$sciName == "Vireo_olivaceus"] <- "Vireo_chivi"
traits$sciName[traits$sciName == "Frederickena_unduligera"] <- "Frederickena_unduliger"
traits$sciName[traits$sciName == "Anurolimnas_fasciatus"] <- "Laterallus_fasciatus"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Frederickena_unduliger", ]
traits$sciName[nrow(traits)] <- "Frederickena_fulva"
traits$sciName[traits$sciName == "Butorides_striata"] <- "Butorides_striatus"
traits$sciName[traits$sciName == "Aratinga_leucophthalma"] <- "Psittacara_leucophthalma"
traits$sciName[traits$sciName == "Lepidocolaptes_albolineatus"] <- "Lepidocolaptes_fatimalimae"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Lepidocolaptes_fatimalimae", ]
traits$sciName[nrow(traits)] <- "Lepidocolaptes_duidae"
traits$sciName[traits$sciName == "Pyrrhura_picta"] <- "Pyrrhura_roseifrons"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Pyrrhura_roseifrons", ]
traits$sciName[nrow(traits)] <- "Pyrrhura_lucianni"
traits$sciName[traits$sciName == "Sclerurus_mexicanus"] <- "Sclerurus_obscurior"
traits$sciName[traits$sciName == "Xenops_milleri"] <- "Microxenops_milleri"
traits$sciName[traits$sciName == "Xenops_milleri"] <- "Microxenops_milleri"
traits[nrow(traits) + 1, ] <- traits[traits$sciName == "Herpsilochmus_gentryi", ]
traits$sciName[nrow(traits)] <- "Herpsilochmus_loreto"

# Update dataframe for modeling with trait data
specialists$bodymass <- specialists$stratum <- specialists$diet <- specialists$ForStrat.wataroundsurf <-
  specialists$ForStrat.watbelowsurf <- specialists$ForStrat.ground <- specialists$ForStrat.understory <-
  specialists$ForStrat.midhigh <- specialists$ForStrat.canopy <- specialists$ForStrat.aerial <- NA

for(i in 1:nrow(specialists)){
  sprow <- traits[traits$sciName == specialists$species[i], ]
  specialists$bodymass[i] <- sprow$BodyMass.Value
  specialists$diet[i] <- sprow$Diet.5Cat
  specialists$ForStrat.wataroundsurf[i] <- sprow$ForStrat.wataroundsurf
  specialists$ForStrat.watbelowsurf[i] <- sprow$ForStrat.watbelowsurf
  specialists$ForStrat.ground[i] <- sprow$ForStrat.ground
  specialists$ForStrat.understory[i] <- sprow$ForStrat.understory
  specialists$ForStrat.midhigh[i] <- sprow$ForStrat.midhigh
  specialists$ForStrat.canopy[i] <- sprow$ForStrat.canopy
  specialists$ForStrat.aerial[i] <- sprow$ForStrat.aerial
  if((sprow$ForStrat.ground + sprow$ForStrat.understory) > 50){
    if(sprow$ForStrat.ground > sprow$ForStrat.understory){specialists$stratum[i] <- "G"}else{
      specialists$stratum[i] <- "U"
    }
  }else if((sprow$ForStrat.midhigh + sprow$ForStrat.canopy) > 50){
    if(sprow$ForStrat.midhigh > sprow$ForStrat.canopy){specialists$stratum[i] <- "M"}else{
      specialists$stratum[i] <- "C"
    }
  }else if((sprow$ForStrat.understory + sprow$ForStrat.midhigh) > 50){
    if(sprow$ForStrat.understory > sprow$ForStrat.midhigh){specialists$stratum[i] <- "U"}else{
      specialists$stratum[i] <- "M"
    }
  }
  else{specialists$stratum[i] <- "O"}
}

i1 <- summaryBy(Count ~ Transect + Point + Visit + Species, 
                data = iquitos_data, FUN = sum) # count of each species on each visit
i22 <- summaryBy(Count.sum ~ Transect + Point + Species, 
                data = i1, FUN = max)  # obtain maximum count
i22$PC <- with(i22, interaction(Transect, Point))
i22$Dis <- NA
iquitos_data$PC <- paste(iquitos_data$Transect, iquitos_data$Point, sep = ".")
for(i in 1:nrow(i22)){
  i22$Dis[i] <- unique(iquitos_data$Dis[iquitos_data$PC == i22$PC[i]])
}


specialists$abun.p <- 0
specialists$abun.d <- 0
for(i in 1:nrow(specialists)){
  cp <- i22$Count.sum.max[i22$Species == specialists$species[i] & i22$Dis == "P"]
  cd <- i22$Count.sum.max[i22$Species == specialists$species[i] & i22$Dis == "D"]
  if(length(cp) > 0){specialists$abun.p[i] <- sum(cp)}
  if(length(cd) > 0){specialists$abun.d[i] <- sum(cd)}
}

specialists$abun.tot <- specialists$abun.d + specialists$abun.p
specialists$abun.double <- cbind(specialists$abun.d, specialists$abun.p)
specialists$scale.mass <- scale(specialists$bodymass)

sum(c(specialists$ForStrat.wataroundsurf, specialists$ForStrat.watbelowsurf, specialists$ForStrat.aerial))/
  sum(c(specialists$ForStrat.wataroundsurf, specialists$ForStrat.watbelowsurf, specialists$ForStrat.aerial,
        specialists$ForStrat.canopy, specialists$ForStrat.ground, specialists$ForStrat.midhigh, specialists$ForStrat.understory))

specialists$StanStrat.ground <- scale(specialists$ForStrat.ground)
specialists$StanStrat.understory <- scale(specialists$ForStrat.understory)
specialists$StanStrat.midstory <- scale(specialists$ForStrat.midhigh)


##### Model fitting #####
model.fits <- list()
set.seed(13)
model.fits$global_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         StanStrat.ground + StanStrat.understory + StanStrat.midstory + diet + scale.mass + migratory + #traits
                         river + flood + tf + poor + rich, #specialization
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

model.fits$null_stanMER <-
  rstanarm::stan_glmer(formula = abun.double ~ (1|species),  #null
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

model.fits$coarse_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist, #forest
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

model.fits$specialization_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         river + flood + tf + poor + rich, #specialization
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

model.fits$traits_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         StanStrat.ground + StanStrat.understory + StanStrat.midstory + diet + scale.mass + migratory, #traits
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

model.fits$traits2_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         diet + scale.mass + migratory, #traits
                       family = 'binomial', data = specialists, seed = 8, cores = 4,
                       prior_covariance = rstanarm::decov(regularization = 1, concentration = 1, shape = 1, scale = 100),
                       prior_intercept = rstanarm::normal(location = 0, scale = 10, autoscale = F),
                       prior = rstanarm::normal(location = 0, scale = 5, autoscale = F))

save(model.fits, file = "Data/model_fits.Rdata")


##### Cross-validation #####
load("Data/model_fits.Rdata")
Xvalid <- list()
Xvalid$global <- rstanarm::kfold(model.fits$global_stanMER)
Xvalid$null <- rstanarm::kfold(model.fits$null_stanMER)
Xvalid$coarse <- rstanarm::kfold(model.fits$coarse_stanMER)
Xvalid$specialization <- rstanarm::kfold(model.fits$specialization_stanMER)
Xvalid$traits <- rstanarm::kfold(model.fits$traits_stanMER)
Xvalid$traits2 <- rstanarm::kfold(model.fits$traits2_stanMER)
save(Xvalid, file = "Data/Xvalid.Rdata")

##### Model-comparison #####
load("Data/model_fits.Rdata")
load("Data/Xvalid.Rdata")
#load("fData/model_fits.Rdata")
#load("fData/Xvalid.Rdata")

rstanarm::compare_models(Xvalid$global, Xvalid$null, Xvalid$coarse, Xvalid$specialization, Xvalid$traits, Xvalid$traits2)
rstanarm::compare_models(Xvalid$global, Xvalid$specialization)
rstanarm::compare_models(Xvalid$specialization, Xvalid$traits2)
rstanarm::compare_models(Xvalid$specialization, Xvalid$traits)
rstanarm::compare_models(Xvalid$traits2, Xvalid$traits)
rstanarm::compare_models(Xvalid$traits2, Xvalid$coarse)
rstanarm::compare_models(Xvalid$traits, Xvalid$coarse)

# Percent (latent) variance explained
re_variances_df <- data.frame(global = as.data.frame(model.fits$global_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              null = as.data.frame(model.fits$null_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              coarse = as.data.frame(model.fits$coarse_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              specialization = as.data.frame(model.fits$specialization_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              traits = as.data.frame(model.fits$traits_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              traits2 = as.data.frame(model.fits$traits2_stanMER)$`Sigma[species:(Intercept),(Intercept)]`)
var_explained <- function(name, df = re_variances_df){return((df["null"] - df[name])/df["null"])}
summary(var_explained("global")[,1])["Mean"]
summary(var_explained("specialization")[,1])["Mean"]
summary(var_explained("traits")[,1])["Mean"]
summary(var_explained("traits2")[,1])["Mean"]
summary(var_explained("coarse")[,1])["Mean"]


###### Plotting #####
load("Data/model_fits.Rdata")
load("Data/Xvalid.Rdata")
#load("fData/model_fits.Rdata")
#load("fData/Xvalid.Rdata")

re_variances_df <- data.frame(global = as.data.frame(model.fits$global_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              null = as.data.frame(model.fits$null_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              coarse = as.data.frame(model.fits$coarse_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              specialization = as.data.frame(model.fits$specialization_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              traits = as.data.frame(model.fits$traits_stanMER)$`Sigma[species:(Intercept),(Intercept)]`,
                              traits2 = as.data.frame(model.fits$traits2_stanMER)$`Sigma[species:(Intercept),(Intercept)]`)
var_explained <- function(name, df = re_variances_df){return((df["null"] - df[name])/df["null"])}

# Figure 2: Model comparison
mcdf <- data.frame(model = c("global", "specialization", "traits", "traits2", "coarse"),
                   elpd = c(Xvalid$global$elpd_kfold,
                            Xvalid$specialization$elpd_kfold,
                            Xvalid$traits$elpd_kfold,
                            Xvalid$traits2$elpd_kfold,
                            Xvalid$coarse$elpd_kfold),
                   r2 = c(summary(var_explained("global")[,1])["Mean"],
                          summary(var_explained("specialization")[,1])["Mean"],
                          summary(var_explained("traits")[,1])["Mean"],
                          summary(var_explained("traits2")[,1])["Mean"],
                          summary(var_explained("coarse")[,1])["Mean"])
)
mcdf$elpd_diff <- mcdf$elpd - Xvalid$null$elpd_kfold
mcdf$r2scaled <- mcdf$r2 * 160

# Not run (figure not in final paper):
# pdf(file = "Figures/model_comparison.pdf", width = 8, height = 5)
#     barplot(height = t(as.matrix(mcdf[ ,c("elpd_diff", "r2scaled")])), beside = T, ylim = c(0,160), 
#             names.arg = c("global", "specialization", "traits", "traits2", "coarse"),
#             col = c("gray60", "gray90"), main = "model comparison", 
#             xlab = "model", legend = c("ELPD improvement", "R-squared"),
#             axes = F)
#     axis(side = 2, at = seq(0,160,40))
#     axis(side = 4, at = seq(0,160,16), labels=as.character(seq(0,1,.1)))
# dev.off()

pdf(file = "Figures/f2a.pdf", width = 6, height = 5)
barplot(height = mcdf$elpd_diff, ylim = c(0,160), 
        names.arg = c("global", "specialization", "traits", "traits2", "coarse"),
        col = c("gray60"), main = "ELPD", 
        xlab = "",
        axes = F)
axis(side = 2, at = seq(0,160,40))
dev.off()

pdf(file = "Figures/f2b.pdf", width = 6, height = 5)
barplot(height = mcdf$r2scaled, ylim = c(0,140), 
        names.arg = c("global", "specialization", "traits", "traits2", "coarse"),
        col = c("gray90"), main = "R-squared", 
        xlab = "",
        axes = F)
axis(side = 2, at = seq(0,160,32), labels=as.character(seq(0,1,.2)))
dev.off()

# Figure 3: Parameter estimates from GLMs
pn <- vector()
pnames <- vector()
pmu <- vector()
plci <- vector()
puci <- vector()
for(i in 1:length(model.fits)){
  pn[i] <- which(rownames(summary(model.fits[[i]])) == "b[(Intercept) species:Amazilia_fimbriata]") - 1
  pnames <- c(pnames, rownames(summary(model.fits[[i]]))[1:pn[i]])
  pmu <- c(pmu, summary(model.fits[[i]])[1:pn[i], 1])
  plci <- c(plci, summary(model.fits[[i]])[1:pn[i], 4])
  puci <- c(puci, summary(model.fits[[i]])[1:pn[i], 8])
}
param_est <- data.frame(model = rep(names(model.fits), pn), parameter = pnames, pMean = pmu,
                        L95 = plci, U95 = puci)
coef.names <- unique(param_est$parameter)
n.coef <- length(coef.names)
model.names <- unique(param_est$model)
n.model <- length(model.names)
colors <- viridis::viridis(n.model, begin = .1)
param_est$color <- NA
param_est$color[param_est$model == "global_stanMER"] <- colors[1]
param_est$color[param_est$model == "null_stanMER"] <- colors[2]
param_est$color[param_est$model == "coarse_stanMER"] <- colors[3]
param_est$color[param_est$model == "specialization_stanMER"] <- colors[4]
param_est$color[param_est$model == "traits_stanMER"] <- colors[5]
param_est$color[param_est$model == "traits2_stanMER"] <- colors[6]
pchs <- c(24,22,21,23,20,25)
plotting_height <- 0
ph <- vector()
counter1 <- 0
for(i in 1:n.coef){
  plotting_height <- plotting_height + .5
  counter <- 0
  counter1 <- counter1 + 1
  ph[counter1] <- plotting_height
  for(j in 1:n.model){
    k <- which(param_est$model == model.names[n.model + 1 - j] & param_est$parameter == coef.names[n.coef + 1 - i])
    if(length(k) == 1){
      counter <- counter + 1
      plotting_height <- plotting_height + .25
    }
  }
}

pdf(file = "Figures/f3.pdf", width = 10.5/.84, height = 10)
    plot(c(1,1), t='n', axes=F, xlab = 'effect size', ylab = "", xlim = c(-12,13), ylim = c(0,22))
    rect(xleft = -12, xright = 13, ybottom = ph[5] - .25, ytop = ph[6] - .1, col = 'gray94', border = NA)
    rect(xleft = -12, xright = 13, ybottom = ph[15] - .25, ytop = ph[17] - .1, col = 'gray94', border = NA)
    lines(c(0,0), c(0,22), col = 'gray')
    axis(side = 1, at = 2*c(-4:4))
    plotting_height <- 0
    for(i in 1:n.coef){
      plotting_height <- plotting_height + .5
      counter <- 0
      for(j in 1:n.model){
        k <- which(param_est$model == model.names[n.model + 1 - j] & param_est$parameter == coef.names[n.coef + 1 - i])
        if(length(k) == 1){
          counter <- counter + 1
          plotting_height <- plotting_height + .25
          lines(c(param_est$L95[k], param_est$U95[k]), rep(plotting_height, 2), col = param_est$color[k])
          points(x = param_est$pMean[k], y = plotting_height, col = param_est$color[k], pch = pchs[j], bg = param_est$color[k])
        }
      }
    }
dev.off()

# Figure 4: posterior predictions
specialists2 <- specialists[order(specialists$species), ]
post <- as.data.frame(model.fits$global_stanMER)
lp_post <- as.data.frame(matrix(data=NA, nrow=4000, ncol=451))
names(lp_post) <- specialists2$species
for(j in 1:451){
  lp_post[,j] <- post[, 1] + post[, 17 + j] + post[, 2]*specialists2$forest.present[j] +
    post[, 3]*specialists2$forest.specialist[j] + post[, 4]*specialists2$StanStrat.ground[j] +
    post[, 5]*specialists2$StanStrat.understory[j] + post[, 6]*specialists2$StanStrat.midstory[j] +
    post[, 7]*(specialists2$diet[j] == "Invertebrate") + post[, 8]*(specialists2$diet[j] == "Omnivore") +
    post[, 9]*(specialists2$diet[j] == "PlantSeed") + post[, 10]*(specialists2$diet[j] == "VertFishScav") +
    post[, 11]*specialists2$scale.mass[j] + post[, 12]*specialists2$migratory[j] +
    post[, 13]*specialists2$river[j] + post[, 14]*specialists2$flood[j] + 
    post[, 15]*specialists2$tf[j] + post[, 16]*specialists2$poor[j] +
    post[, 17]*specialists2$rich[j]
}
post_intervals <- data.frame(species = specialists2$species, mean=NA, lci=NA, uci=NA)
for(i in 1:451){
  post_intervals$mean[i] <- exp(mean(lp_post[,i]))
  post_intervals$lci[i] <- exp(quantile(lp_post[,i], .025))
  post_intervals$uci[i] <- exp(quantile(lp_post[,i], .975))
}
post_intervals2 <- cbind(post_intervals, specialists2)
post_intervals2 <- post_intervals2[order(post_intervals2$mean), ]

pdf(file = "Figures/f4.pdf", width = 9, height = 6)
par(mfrow = c(3,2), mar = c(1,3,1,0.1))

ifexpr <- c("T", "post_intervals2$forest.specialist[i]", "post_intervals2$forest.present[i] & post_intervals2$flood[i]",
            "post_intervals2$poor[i]", "post_intervals2$river[i]", "post_intervals2$diet[i] == \"Invertebrate\"")
            
for(k in 1:6){
  plot(post_intervals2$mean, #col = pdcolors[(fold_changes_sig$col)],
       xaxt = 'n', xlab="", ylab = "", ylim=c(10^-5, 10^5), 
       yaxt="n", log="y", pch=16, cex=.5, cex.axis=0.7,
       col = "white", bty = 'n')
  box(col = "gray80")
  for(i in 1:nrow(post_intervals2)){
    if(eval(parse(text = ifexpr[k]))){
      segments(i, post_intervals2$lci[i], y1=post_intervals2$uci[i], col = "gray50")
    }
  }
  for(i in 1:nrow(post_intervals2)){
    if(eval(parse(text = ifexpr[k]))){
      points(i, post_intervals2$mean[i], col = "gray25", pch=16)
    }
  }
  abline(h=1)
  at.y <- 10^(-5:5)
  lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
  if(k %in% c(1,3,5)){axis(2, at=at.y, labels=lab.y, las=1, cex.axis=.8)}
}
dev.off()

# Figure 5: floodplain specialist proliferation
iquitos_data$tpt <- paste(iquitos_data$Transect, iquitos_data$Point, sep = "_")
allpoints <- iquitos_data[!duplicated(iquitos_data$tpt), c('Transect', 'Hab', 'Dis', 'Point', 'tpt')]
allpoints$nfloodsp <- 0
allpoints$nfloodindiv <- 0
fdata <- iquitos_data[iquitos_data$Species %in% allfloodsp, ]
for(i in 1:nrow(allpoints)){
  fdata_pt <- fdata[fdata$tpt == allpoints$tpt[i], ]
  allpoints$nfloodsp[i] <- length(unique(fdata_pt$Species))
  counter <- 0
  flood_abuns <- vector()
  for(j in allfloodsp){
    fsppt <- fdata_pt[fdata_pt$Species == j, ]
    if(nrow(fsppt) > 0){
      counter <- counter + 1
      flood_abuns[counter] <- max(fsppt$Count)
    }
  }
  
  allpoints$nfloodindiv[i] <- sum(flood_abuns)
}

se <- function(x){sqrt(var(x)/length(x))}

barplotdata <- data.frame(group = c("f.p", "f.d", "t.p", "t.d"),
                          mean.rich = c(mean(allpoints$nfloodsp[allpoints$Hab == "V" & allpoints$Dis == "P"]),
                                        mean(allpoints$nfloodsp[allpoints$Hab == "V" & allpoints$Dis == "D"]),
                                        mean(allpoints$nfloodsp[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "P"]),
                                        mean(allpoints$nfloodsp[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "D"])),
                          se.rich = c(se(allpoints$nfloodsp[allpoints$Hab == "V" & allpoints$Dis == "P"]),
                                      se(allpoints$nfloodsp[allpoints$Hab == "V" & allpoints$Dis == "D"]),
                                      se(allpoints$nfloodsp[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "P"]),
                                      se(allpoints$nfloodsp[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "D"])),
                          mean.abun = c(mean(allpoints$nfloodindiv[allpoints$Hab == "V" & allpoints$Dis == "P"]),
                                        mean(allpoints$nfloodindiv[allpoints$Hab == "V" & allpoints$Dis == "D"]),
                                        mean(allpoints$nfloodindiv[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "P"]),
                                        mean(allpoints$nfloodindiv[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "D"])),
                          se.abun = c(se(allpoints$nfloodindiv[allpoints$Hab == "V" & allpoints$Dis == "P"]),
                                      se(allpoints$nfloodindiv[allpoints$Hab == "V" & allpoints$Dis == "D"]),
                                      se(allpoints$nfloodindiv[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "P"]),
                                      se(allpoints$nfloodindiv[allpoints$Hab %in% c("W", "U") & allpoints$Dis == "D"]))
)
barplotdata$lse.rich <- barplotdata$mean.rich - barplotdata$se.rich
barplotdata$use.rich <- barplotdata$mean.rich + barplotdata$se.rich
barplotdata$lse.abun <- barplotdata$mean.abun - barplotdata$se.abun
barplotdata$use.abun <- barplotdata$mean.abun + barplotdata$se.abun

pdf(file = "Figures/f5a.pdf", width = 5.2, height = 5)
    barplot(height = matrix(barplotdata$mean.rich, nrow = 2), beside = T, ylim = c(0,22), names.arg = c("floodplain", "terra firme"),
            col = c("gray60", "gray90"), main = "floodplain specialist richness", xlab = "forest-type", ylab = "# species", legend = c("primary forest", "agriculture"))
    for(i in 1:4){
      lines(c(.5 + i + (i > 2), .5 + i + (i > 2)), c(barplotdata$lse.rich[i], barplotdata$use.rich[i]), col = "black", lwd = 2)
    }
dev.off()

pdf(file = "Figures/f5b.pdf", width = 5.2, height = 5)
    barplot(height = matrix(barplotdata$mean.abun, nrow = 2), beside = T, ylim = c(0,37), names.arg = c("floodplain", "terra firme"),
            col = c("gray60", "gray90"), main = "floodplain specialist abundance", xlab = "forest-type", ylab = "# individuals", legend = c("primary forest", "agriculture"))
    for(i in 1:4){
      lines(c(.5 + i + (i > 2), .5 + i + (i > 2)), c(barplotdata$lse.abun[i], barplotdata$use.abun[i]), col = "black", lwd = 2)
    }
dev.off()

# Figure S1: rich vs poor species
i22$Hab <- NA
for(i in 1:nrow(i22)){
  i22$Hab[i] <- unique(iquitos_data$Hab[iquitos_data$PC == i22$PC[i]])
}

iq2 <- i22[i22$Hab %in% c("U","W") & i22$Dis == "P",]


iq2$tpt <- paste(iq2$Transect, iq2$Point, sep = "_")
iqp <- iq2[iq2$Species %in% specialists$species[specialists$poor == 1], ]
iqr <- iq2[iq2$Species %in% specialists$species[specialists$rich == 1], ]

mdf <- data.frame(pt = unique(iq2$tpt), nr = 0, np = 0)
for(i in 1:nrow(mdf)){
  mdf$nr[i] <- length(unique(iqr$Species[iqr$tpt == mdf$pt[i]]))
  mdf$np[i] <- length(unique(iqp$Species[iqp$tpt == mdf$pt[i]]))
}

pdf(file = "Figures/S1.pdf", width = 5, height = 5)
  plot(jitter(mdf$np) ~ jitter(mdf$nr), xlab = "rich-soil specialist richness", ylab = "poor-soil specialist richness")
dev.off()

# Figure S2: model checking: PPC
pdf(file = "Figures/S2.pdf", width = 5, height = 3)
rstanarm::pp_check(model.fits$global_stanMER)
dev.off()

# Figure S3: model checking--mixed
link.resids <- as.data.frame(model.fits$global_stanMER)[, 18:468]
specialists2 <- specialists[order(specialists$species), ]
post <- as.data.frame(model.fits$global_stanMER)
lp_post_noRintercepts <- as.data.frame(matrix(data=NA, nrow=4000, ncol=451))
names(lp_post_noRintercepts) <- specialists2$species
for(j in 1:451){
  lp_post_noRintercepts[,j] <- post[, 1] + post[, 2]*specialists2$forest.present[j] +
    post[, 3]*specialists2$forest.specialist[j] + post[, 4]*specialists2$StanStrat.ground[j] +
    post[, 5]*specialists2$StanStrat.understory[j] + post[, 6]*specialists2$StanStrat.midstory[j] +
    post[, 7]*(specialists2$diet[j] == "Invertebrate") + post[, 8]*(specialists2$diet[j] == "Omnivore") +
    post[, 9]*(specialists2$diet[j] == "PlantSeed") + post[, 10]*(specialists2$diet[j] == "VertFishScav") +
    post[, 11]*specialists2$scale.mass[j] + post[, 12]*specialists2$migratory[j] +
    post[, 13]*specialists2$river[j] + post[, 14]*specialists2$flood[j] + 
    post[, 15]*specialists2$tf[j] + post[, 16]*specialists2$poor[j] +
    post[, 17]*specialists2$rich[j]
}
dev.off()
pdf(file = "Figures/S3.pdf", width = 8, height = 8)
par(mar = c(2,2,0.1,0.1), mfrow = c(5,4))
for(i2 in 1:20){
  i <- 200*i2
  xaxt <- NULL
  yaxt <- NULL
  if(!(i2 %in% c(1,5,9,13,17))){yaxt <- 'n'}
  if(!(i2 %in% 17:20)){xaxt <- 'n'}
  plot(as.numeric(link.resids[i, ]) ~ as.numeric(lp_post_noRintercepts[i, ]), xlab = "",
       ylab = "", ylim = c(-12,12), xlim = c(-10,12), xaxt=xaxt, yaxt=yaxt)
}
dev.off()

# Figure S4: forest cover near points
download.file("https://datadryad.org/bitstream/handle/10255/dryad.220655/Census_Points.csv?sequence=1",
              destfile = 'dData/veg.csv')
veg_data <- read.csv('dData/veg.csv', fileEncoding = 'latin1')
veg_data[is.na(veg_data)] <- 0
veg_data <- veg_data[-which(veg_data$Discard == 1), ]
veg_d <- veg_data[veg_data$Disturbance == "D", ]
veg_p <- veg_data[veg_data$Disturbance == "U", ]
veg_dd <- veg_d
veg_dd$Habitat <- as.character(veg_dd$Habitat)
veg_dd$Habitat[veg_dd$Habitat == "Arena Blanca"] <- "White-sand points"
veg_dd$Habitat[veg_dd$Habitat == "Tierra Firme"] <- "Upland points"
veg_dd$Habitat[veg_dd$Habitat == "Varzea"] <- "Floodplain points"
veg_dd$Habitat <- as.factor(veg_dd$Habitat)
veg_dd$Habitat <- factor(veg_dd$Habitat, levels=c("Upland points","Floodplain points","White-sand points"))
veg_dd$riverNS <- NA
veg_dd$riverNS[veg_dd$Location == "Miraflores"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Km77"] <- "N"
veg_dd$riverNS[veg_dd$Location == "ExplorNapo"] <- "N"
veg_dd$riverNS[veg_dd$Location == "El Varillal"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Moronacocha"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Nueva Esperanza"] <- "N"
veg_dd$riverNS[veg_dd$Location == "El Dorado"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Indiana"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Km35"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Sabalillo"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Santa Clotilde"] <- "N"
veg_dd$riverNS[veg_dd$Location == "Jenaro Herrera"] <- "S"
veg_dd$riverNS[veg_dd$Location == "Requena"] <- "S"
veg_dd$riverNS[veg_dd$Location == "Chino"] <- "S"
veg_dd$riverNS[veg_dd$Location == "Tamshiyacu"] <- "S"
veg_dd$riverNS[veg_dd$Location == "Madre Selva"] <- "S"
veg_dd$CCF <- 100*(veg_dd$TallCCF + veg_dd$ShortCCF.including.old.gap.)/(100-veg_dd$Lake.stream.15m)
vegNS <- data.frame(group = c("N", "S"),
                    CCFmean = c(mean(veg_dd$CCF[veg_dd$riverNS == "N"]),
                                mean(veg_dd$CCF[veg_dd$riverNS == "S"])),
                    CCFse = c(se(veg_dd$CCF[veg_dd$riverNS == "N"]),
                              se(veg_dd$CCF[veg_dd$riverNS == "S"])))

vegNS$lse <- vegNS$CCFmean - vegNS$CCFse
vegNS$use <- vegNS$CCFmean + vegNS$CCFse
vegHab <- data.frame(group = c("W", "U", "F"),
                     CCFmean = c(mean(veg_dd$CCF[veg_dd$Habitat == "White-sand points"]),
                                 mean(veg_dd$CCF[veg_dd$Habitat == "Upland points"]),
                                 mean(veg_dd$CCF[veg_dd$Habitat == "Floodplain points"])),
                     CCFse = c(se(veg_dd$CCF[veg_dd$Habitat == "White-sand points"]),
                               se(veg_dd$CCF[veg_dd$Habitat == "Upland points"]),
                               se(veg_dd$CCF[veg_dd$Habitat == "Floodplain points"])))

vegHab$lse <- vegHab$CCFmean - vegHab$CCFse
vegHab$use <- vegHab$CCFmean + vegHab$CCFse
vegHab <- vegHab[3:1,]

pdf(file = "Figures/S4a.pdf", width = 5.2, height = 5)
barplot(height = vegNS$CCFmean, ylim = c(0,45), names.arg = c("north-bank", "south-bank"),
        col = c("gray70", "gray70"), main = "across the Amazon", ylab = "closed-canopy forest (% cover)")
for(i in 1:2){
  lines(c(i - .3 + 0.2*(i > 1), i - .3 + 0.2*(i > 1)), c(vegNS$lse[i], vegNS$use[i]), col = "black", lwd = 2)
}
dev.off()

pdf(file = "Figures/S4b.pdf", width = 5.2, height = 5)
barplot(height = vegHab$CCFmean, ylim = c(0,50), names.arg = c("floodplain", "terra firme", "white-sands"),
        col = rep("gray70", 3), main = "across forest-types", ylab = "closed-canopy forest (% cover)")
for(i in 1:3){
  lines(c(i - .3 + 0.2*(i - 1), i - .3 + 0.2*(i - 1)), c(vegHab$lse[i], vegHab$use[i]), col = "black", lwd = 2)
}
dev.off()

# Table S4: Random intercepts (species-specific intercept terms from global model)
summary(model.fits$global_stanMER)
