user <- "JacobSocolar"
#user <- "Jacob"

setwd(paste0("/Users/",user,"/Dropbox/Work/Iquitos/Data"))
'%ni%' <- Negate('%in%')

cdata <- read.csv("cdata.csv", stringsAsFactors = F)
stotz <- read.csv(paste0("/Users/",user,"/Dropbox/Work/Useful_data/Stotz_et_al/stotz_adata.csv"))
extras <- read.csv("stotz_cdata_additions.csv")
stotz.nb <- read.csv(paste0("/Users/",user,"/Dropbox/Work/Useful_data/Stotz_et_al/Stotz_nearctic_nonbreeding.csv"))
stotz <- gtools::smartbind(stotz, stotz.nb)
stotz <- gtools::smartbind(stotz, extras)
stotz$GENUS <- as.character(stotz$GENUS)
stotz$GENUS[stotz$GENUS == "Columba"] <- "Patagioenas"
stotz$sp <- paste(stotz$GENUS, stotz$SPECIES, sep = "_")

# Update stotz taxonomy to reflect cdata taxonomy:
stotz$sp[grep("Capito_niger", stotz$sp)] <- "Capito_auratus"
stotz$sp[grep("Pitylus_grossus", stotz$sp)] <- "Saltator_grossus"
stotz$sp[grep("Xiphorhynchus_\\(guttatus\\) guttatus", stotz$sp)] <- "Xiphorhynchus_guttatus"
stotz$sp[grep("Momotus_\\(momota\\) momota", stotz$sp)] <- "Momotus_momota"
stotz$sp[grep("Myrmotherula_\\(brachyura\\) ignota", stotz$sp)] <- "Myrmotherula_ignota"
stotz$sp[grep("Pionites_melanocephala", stotz$sp)] <- "Pionites_melanocephalus"
stotz$sp[grep("Pipra_coronota", stotz$sp)] <- "Lepidothrix_coronata"
stotz$sp[grep("Hypocnemis_cantator", stotz$sp)] <- "Hypocnemis_peruviana"
stotz$sp[grep("Ortalis_\\(motmot\\) guttata", stotz$sp)] <- "Ortalis_guttata"
stotz$sp[grep("Hylophylax_poecilinota", stotz$sp)] <- "Willisornis_poecilinotus"
stotz$sp[grep("Tangara_\\(velia\\) velia", stotz$sp)] <- "Tangara_velia"
stotz$sp[grep("Psarocolius_latirostris", stotz$sp)] <- "Cacicus_latirostris"
stotz$sp[grep("Myrmotherula_hauxwelli", stotz$sp)] <- "Isleria_hauxwelli"
stotz$sp[grep("Schiffornis_turdinus", stotz$sp)] <- "Schiffornis_turdina"
stotz$sp[grep("Tolmomyias_\\(assimilis\\) assimilis", stotz$sp)] <- "Tolmomyias_assimilis"
stotz$sp[grep("Pipra_pipra", stotz$sp)] <- "Dixiphia_pipra"
stotz$sp[grep("Oryzoborus_\\(angolensis\\) angolensis", stotz$sp)] <- "Oryzoborus_angolensis"
stotz$sp[grep("Thamnophilus_\\(doliatus\\) doliatus", stotz$sp)] <- "Thamnophilus_doliatus"
stotz$sp[grep("Brotogeris_\\(versicolurus\\) versicolurus", stotz$sp)] <- "Brotogeris_versicolurus"
stotz$sp[grep("Todirostrum_latirostre", stotz$sp)] <- "Poecilotriccus_latirostris"
stotz$sp[grep("Buteo_magnirostris", stotz$sp)] <- "Rupornis_magnirostris"
stotz$sp[grep("Xiphorhynchus_picus", stotz$sp)] <- "Dendroplex_picus"
stotz$sp[grep("Glaucis_hirsuta", stotz$sp)] <- "Glaucis_hirsutus"
stotz$sp[grep("Trogon_violaceus", stotz$sp)] <- "Trogon_ramonianus"
stotz$sp[grep("Jacamerops_aurea", stotz$sp)] <- "Jacamerops_aureus"
stotz$sp[grep("Psarocolius_oseryi", stotz$sp)] <- "Clypicterus_oseryi"
stotz$sp[grep("Leptotila_\\(rufaxilla\\) rufaxilla", stotz$sp)] <- "Leptotila_rufaxilla"
stotz$sp[grep("Hylophylax_punctulata", stotz$sp)] <- "Hylophylax_punctulatus"
stotz$sp[grep("Tolmomyias_\\(flaviventris\\) flaviventris", stotz$sp)] <- "Tolmomyias_flaviventris"
stotz$sp[grep("Myrmotherula_\\(brachyura\\) brachyura", stotz$sp)] <- "Myrmotherula_brachyura"
stotz$sp[grep("Donacobius_atricapillus", stotz$sp)] <- "Donacobius_atricapilla"
stotz$sp[grep("Nonnula_\\(ruficapilla\\) ruficapilla", stotz$sp)] <- "Nonnula_ruficapilla"
stotz$sp[grep("Amazona_\\(ochrocephala\\) ochrocephala", stotz$sp)] <- "Amazona_ochrocephala"
stotz$sp[grep("Thryothorus_\\(genibarbis\\) genibarbis", stotz$sp)] <- "Thryothorus_genibarbis"
stotz$sp[grep("Polioptila_\\(plumbea\\) plumbea", stotz$sp)] <- "Polioptila_plumbea"
stotz$sp[grep("Columbina_\\(talpacoti\\) talpacoti", stotz$sp)] <- "Columbina_talpacoti"
stotz$sp[grep("Sporophila_\\(lineola\\) bouvronides", stotz$sp)] <- "Sporophila_buvronides"
stotz$sp[grep("Leucopternis_schistacea", stotz$sp)] <- "Buteogallus_schistaceus"
stotz$sp[grep("Vireo_\\(olivaceus\\) chivi", stotz$sp)] <- "Vireo_chivi"
stotz$sp[grep("Phylloscartes_flaveolus", stotz$sp)] <- "Capsiempis_flaveola"
stotz$sp[grep("Xiphorhynchus_\\(spixii\\) elegans", stotz$sp)] <- "Xiphorhynchus_elegans"
stotz$sp[grep("Frederickena_unduligera", stotz$sp)] <- "Frederickena_unduliger"
stotz$sp[grep("Hylophylax_naevia", stotz$sp)] <- "Hylophylax_naevius"
stotz$sp[grep("Conopias_\\(albovittata\\) parva", stotz$sp)] <- "Conopias_parvus"
stotz$sp[grep("Daptrius_americanus", stotz$sp)] <- "Ibycter_americanus"
stotz$sp[grep("Myrmotherula_\\(haematonota\\) haematonota", stotz$sp)] <- "Epinecrophylla_haematonota"
stotz$sp[grep("Percnostola_leucostigma", stotz$sp)] <- "Schistocichla_leucostigma"
stotz$sp[grep("Troglodytes_\\(aedon\\) aedon", stotz$sp)] <- "Troglodytes_aedon"
stotz$sp[grep("Dendrocincla_\\(fuliginosa\\) fuliginosa", stotz$sp)] <- "Dendrocincla_fuliginosa"
stotz$sp[grep("Percnostola_schistacea", stotz$sp)] <- "Schistocichla_schistacea"
stotz$sp[grep("Herpsilochmus_\\(sticturus\\) dugandi", stotz$sp)] <- "Herpsilochmus_dugandi"
stotz$sp[grep("Notharchus_macrorhynchos", stotz$sp)] <- "Notharchus_hyperrhynchus"
stotz$sp[grep("Ceryle_torquata", stotz$sp)] <- "Megaceryle_torquata"
stotz$sp[grep("Philydor_erythropterus", stotz$sp)] <- "Philydor_erythropterum"
stotz[nrow(stotz)+1, ] <- stotz[which(stotz$sp == "Frederickena_unduliger"), ]
stotz$sp[nrow(stotz)] <- "Frederickena_fulva"
stotz$sp[grep("Philydor_erythropterus", stotz$sp)] <- "Philydor_erythropterum"
stotz$sp[grep("barrabandi", stotz$sp)] <- "Pyrilia_barrabandi"
stotz$sp[grep("Pyrrhura_\\(melanura\\) melanura", stotz$sp)] <- "Pyrrhura_melanura"
stotz$sp[grep("Patagioenas_livia", stotz$sp)] <- "Columba_livia"
stotz$sp[grep("Chondrohierax_\\(uncinatus\\) uncinatus", stotz$sp)] <- "Chondrohierax_uncinatus"
stotz$sp[grep("Butorides_\\(striatus\\) striatus", stotz$sp)] <- "Butorides_striatus"
stotz$sp[grep("barrabandi", stotz$sp)] <- "Pyrilia_barrabandi"
stotz$sp[grep("Rostrhamus_hamatus", stotz$sp)] <- "Helicolestes_hamatus"
stotz$sp[grep("Furnarius_\\(leucopus\\) leucopus", stotz$sp)] <- "Furnarius_leucopus"
stotz$sp[grep("Forpus_crassirostris", stotz$sp)] <- "Forpus_xanthopterygius"
stotz$sp[grep("Xiphorhynchus_necopinus", stotz$sp)] <- "Dendroplex_kienerii"
stotz$sp[grep("Aratinga_leucophthalmus", stotz$sp)] <- "Psittacara_leucophthalma"
stotz$sp[grep("Lepidocolaptes_albolineatus", stotz$sp)] <- "Lepidocolaptes_fatimalimae"
stotz[nrow(stotz)+1, ] <- stotz[which(stotz$sp == "Lepidocolaptes_fatimalimae"), ]
stotz$sp[nrow(stotz)] <- "Lepidocolaptes_duidae"
stotz$sp[grep("Dendrocolaptes_\\(certhia\\) certhia", stotz$sp)] <- "Dendrocolaptes_certhia"
stotz$sp[grep("Tangara_\\(mexicana\\) mexicana", stotz$sp)] <- "Tangara_mexicana"
stotz$sp[grep("Pyrrhura_\\(picta\\) picta", stotz$sp)] <- "Pyrrhura_roseifrons"
stotz[nrow(stotz)+1, ] <- stotz[which(stotz$sp == "Pyrrhura_roseifrons"), ]
stotz$sp[nrow(stotz)] <- "Pyrrhura_lucianni"
stotz$sp[grep("Onychorhynchus_\\(coronatus\\) coronatus", stotz$sp)] <- "Onychorhynchus_coronatus"
stotz$sp[grep("Gymnopithys_lunulata", stotz$sp)] <- "Gymnopithys_lunulatus"
stotz$sp[grep("Passerina_cyanoides", stotz$sp)] <- "Cyanocompsa_cyanoides"
stotz$sp[grep("Agelaius_icterocephalus", stotz$sp)] <- "Chrysomus_icterocephalus"
stotz$sp[grep("Porphyrula_martinica", stotz$sp)] <- "Porphyrio_martinica"
stotz$sp[grep("Chiroxiphia_\\(pareola\\) pareola", stotz$sp)] <- "Chiroxiphia_pareola"
stotz$sp[grep("Myrmotherula_erythrura", stotz$sp)] <- "Epinecrophylla_erythrura"
stotz$sp[grep("Sclerurus_mexicanus", stotz$sp)] <- "Sclerurus_obscurior"
stotz$sp[grep("Ardeola_ibis", stotz$sp)] <- "Bubulcus_ibis"
stotz$sp[grep("Scaphidura_oryzivora", stotz$sp)] <- "Molothrus_oryzivorus"
stotz$sp[grep("Chlorestes_notatus", stotz$sp)] <- "Chlorestes_notata"
stotz$sp[grep("Xenops_milleri", stotz$sp)] <- "Microxenops_milleri"
stotz$sp[grep("Synallaxis_\\(gujanensis\\) gujanensis", stotz$sp)] <- "Synallaxis_gujanensis"
stotz$sp[grep("Certhiaxis_mustelina", stotz$sp)] <- "Certhiaxis_mustelinus"
stotz$sp[grep("Metopothrix_aurantiacus", stotz$sp)] <- "Metopothrix_aurantiaca"
stotz$sp[grep("Phaethornis_longuemareus", stotz$sp)] <- "Phaethornis_atrimentalis"
stotz$sp[grep("Xiphocolaptes_\\(promeropirhynchus\\) orenocensis", stotz$sp)] <- "Xiphocolaptes_promeropirhynchus"
stotz$sp[grep("Sporophila_\\(lineola\\) lineola", stotz$sp)] <- "Sporophila_lineola"
stotz$sp[grep("Myiodynastes_\\(maculatus\\) solitarius", stotz$sp)] <- "Myiodynastes_maculatus"

# correct cdata taxonomy where problematic:
# remove Trogon melanurus and Hylophylax naevius
cdata <- cdata[-grep("Trogon_melanurus", cdata$Species), ]
cdata <- cdata[-grep("Hylophylax_naevius", cdata$Species), ]

# correct spellings
cdata$Species[which(cdata$Species == "Odonotophorus_gujanensis")] <- "Odontophorus_gujanensis"
cdata$Species[which(cdata$Species == "Brotogeris_versicolorus")] <- "Brotogeris_versicolurus"
cdata$Species[which(cdata$Species == "Megarhynchus_pitangua")] <- "Megarynchus_pitangua"
cdata$Species[which(cdata$Species == "Vireo_olivaceus")] <- "Vireo_chivi"
cdata$Species[which(cdata$Species == "Galbula_flavirostris")] <- "Galbula_albirostris"
cdata$Species[which(cdata$Species == "Ancistrops_strigulatus")] <- "Ancistrops_strigilatus"
cdata$Species[which(cdata$Species == "Gymnoderis_foetidus")] <- "Gymnoderus_foetidus"
cdata$Species[which(cdata$Species == "Todirostrum_chrysocephalum")] <- "Todirostrum_chrysocrotaphum"
cdata$Species[which(cdata$Species == "Dromococcyx_phaisanellus")] <- "Dromococcyx_phasianellus"
cdata$Species[which(cdata$Species == "Loreto_antwren")] <- "Herpsilochmus_loreto"

cdata$Species[cdata$Species %ni% stotz$sp]

write.csv(stotz, file = "stotz_cdata_taxonomy.csv")
write.csv(cdata, file = "cdata_updated.csv")

################################
stotz <- read.csv("stotz_cdata_taxonomy.csv")
cdata <- read.csv("cdata_updated.csv")

# habitat lumping
stotz$flood.for <- stotz$F2 %in% c("Y", "Q") | stotz$F3 %in% c("Y", "Q") | stotz$F13 %in% c("Y", "Q")
stotz$flood.nfor <- stotz$N11 %in% c("Y", "Q") | stotz$N12 %in% c("Y", "Q") | stotz$A1 %in% c("Y", "Q") | stotz$A5 %in% c("Y", "Q") | 
  stotz$A6 %in% c("Y", "Q") | stotz$A8 %in% c("Y", "Q") | stotz$A9 %in% c("Y", "Q")
stotz$flood.all <- stotz$flood.for | stotz$flood.nfor
stotz$tf <- stotz$F1 %in% c("Y", "Q") | stotz$F12 %in% c("Y", "Q")
stotz$forest.present <- (stotz$flood.for == 1) | (stotz$tf == 1)
stotz$upland <- stotz$F1 %in% c("Y", "Q")

# forest-present species
FBsp <- as.character(stotz$sp[stotz$forest.present])
cFBsp <- FBsp[FBsp %in% cdata$Species]

# forest-specialist species
FSsp <- as.character(stotz$sp[stotz$forest.present & !stotz$flood.nfor])
cFSsp <- FSsp[FSsp %in% cdata$Species]

# forest-based floodplain specialists
FFsp <- as.character(stotz$sp[stotz$flood.for & !stotz$tf])
cFFsp <- FFsp[FFsp %in% cdata$Species]

# non-forest-based floodplain specialists (may overlap with forest-based species)
FNsp <- as.character(stotz$sp[stotz$flood.nfor & !stotz$tf])
cFNsp <- FNsp[FNsp %in% cdata$Species]

# all floodplain specialists
allfloodsp <- unique(c(cFFsp, cFNsp))

# terra firme specialists
TFsp <- as.character(stotz$sp[stotz$tf & !(stotz$flood.for | stotz$flood.nfor)])
cTFsp <- TFsp[TFsp %in% cdata$Species]

tfexclude <- read.csv("tf_exclude.csv", header = F)
cTFsp <- cTFsp[cTFsp %ni% tfexclude[,1]]

spec <- read.csv("specialists.csv")
floodextras <- read.csv("flood_specialist_extras.csv")

spec[is.na(spec)] <- 0

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

traits <- read.delim(paste0("/Users/",user,"/Dropbox/Work/Useful_data/EltonTraits/BirdFuncDat.txt"), header=T, stringsAsFactors = F)   # traits data (includes diet)
traits$sciName <- gsub(" ", "_", traits$Scientific)
traits <- traits[1:9993,]

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

specialists$abun.p <- 0
specialists$abun.d <- 0
for(i in 1:nrow(specialists)){
  cp <- cdata$Count[cdata$Species == specialists$species[i] & cdata$Dis == "P"]
  cd <- cdata$Count[cdata$Species == specialists$species[i] & cdata$Dis == "D"]
  if(length(cp) > 0){specialists$abun.p[i] <- sum(cp)}
  if(length(cd) > 0){specialists$abun.d[i] <- sum(cd)}
}

specialists$abun.tot <- specialists$abun.d + specialists$abun.p
specialists$abun.double <- cbind(specialists$abun.p, specialists$abun.d)
specialists$scale.mass <- scale(specialists$bodymass)

sum(c(specialists$ForStrat.wataroundsurf, specialists$ForStrat.watbelowsurf, specialists$ForStrat.aerial))/
  sum(c(specialists$ForStrat.wataroundsurf, specialists$ForStrat.watbelowsurf, specialists$ForStrat.aerial,
        specialists$ForStrat.canopy, specialists$ForStrat.ground, specialists$ForStrat.midhigh, specialists$ForStrat.understory))

specialists$StanStrat.ground <- scale(specialists$ForStrat.ground)
specialists$StanStrat.understory <- scale(specialists$ForStrat.understory)
specialists$StanStrat.midstory <- scale(specialists$ForStrat.midhigh)

model.fits <- list()

model.fits$global_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         StanStrat.ground + StanStrat.understory + StanStrat.midstory + diet + scale.mass + migratory + #traits
                         river + flood + upland + poor + rich, #specialization
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$global_stanMER)[1:17,]

model.fits$null_stanMER <-
  rstanarm::stan_glmer(formula = abun.double ~ (1|species),  #null
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$null_stanMER)[1:17,]

model.fits$null2_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist, #forest
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$null2_stanMER)[1:17,]

model.fits$specialization_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         river + flood + upland + poor + rich, #specialization
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$specialization_stanMER)[1:17,]

model.fits$traits_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         StanStrat.ground + StanStrat.understory + StanStrat.midstory + diet + scale.mass + migratory, #traits
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$traits_stanMER)[1:17,]

model.fits$traits2_stanMER <- 
  rstanarm::stan_glmer(formula = abun.double ~ (1|species) +  #null
                         forest.present + forest.specialist + #forest
                         diet + scale.mass + migratory, #traits
                       family = 'binomial', data = specialists, seed = 8, cores = 4)
summary(model.fits$traits2_stanMER)[1:17,]

save(model.fits, file = "model_fits.Rdata")
load("model_fits.Rdata")

Xvalid <- list()
Xvalid$global <- rstanarm::kfold(model.fits$global_stanMER)
Xvalid$null <- rstanarm::kfold(model.fits$null_stanMER)
Xvalid$null2 <- rstanarm::kfold(model.fits$null2_stanMER)
Xvalid$specialization <- rstanarm::kfold(model.fits$specialization_stanMER)
Xvalid$traits <- rstanarm::kfold(model.fits$traits_stanMER)
Xvalid$traits2 <- rstanarm::kfold(model.fits$traits2_stanMER)

save(Xvalid, file = "Xvalid.Rdata")

rstanarm::compare_models(Xvalid$global, Xvalid$null, Xvalid$null2, Xvalid$specialization, Xvalid$traits, Xvalid$traits2)
rstanarm::compare_models(Xvalid$global, Xvalid$specialization)
rstanarm::compare_models(Xvalid$specialization, Xvalid$traits)

###### Percent (latent) variance explained #########
load("model_fits.Rdata")
re_variances_df <- data.frame(global = as.data.frame(model.fits$global_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2,
                              null = as.data.frame(model.fits$null_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2,
                              null2 = as.data.frame(model.fits$null2_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2,
                              specialization = as.data.frame(model.fits$specialization_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2,
                              traits = as.data.frame(model.fits$traits_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2,
                              traits2 = as.data.frame(model.fits$traits2_stanMER)$`Sigma[species:(Intercept),(Intercept)]`^2)

var_explained <- function(name, df = re_variances_df){return((df["null"] - df[name])/df["null"])}

summary2 <- function(x){return(list(mean = mean(x), lci = quantile(x, .025), uci = quantile(x, .975)))}

summary(var_explained("global")[,1])

###### Extract parameter estimates for plotting ###########
# How many total parameters?
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
param_est$color[param_est$model == "null2_stanMER"] <- colors[3]
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


pdf(file = paste0("/Users/", user, "/Dropbox/Work/Iquitos/Specialist_responses/GLM_coefs.pdf"), width = 10.5, height = 10)

plot(c(1,1), t='n', axes=F, xlab = 'effect size', ylab = "", xlim = c(-8,13), ylim = c(0,22))
rect(xleft = -8, xright = 13, ybottom = ph[1] - .25, ytop = ph[6] - .1, col = 'gray94', border = NA)
rect(xleft = -8, xright = 13, ybottom = ph[15] - .25, ytop = ph[17] - .1, col = 'gray94', border = NA)
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

###### Figure 2: floodplain birds in disturbed and intact floodplains and uplands #####
cdata$tpt <- paste(cdata$Transect, cdata$Point, sep = "_")
allpoints <- cdata[!duplicated(cdata$tpt), c('Transect', 'Hab', 'Dis', 'Point', 'tpt')]
allpoints$nfloodsp <- 0
allpoints$nfloodindiv <- 0
fdata <- cdata[cdata$Species %in% allfloodsp, ]
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


pdf(file = paste0("/Users/", user, "/Dropbox/Work/Iquitos/Specialist_responses/flood_specialist_rich.pdf"), width = 5.2, height = 5)
barplot(height = matrix(barplotdata$mean.rich, nrow = 2), beside = T, ylim = c(0,22), names.arg = c("floodplain", "terra firme"),
        col = c("gray60", "gray90"), main = "floodplain specialist richness", xlab = "forest-type", legend = c("primary forest", "agriculture"))
for(i in 1:4){
  lines(c(.5 + i + (i > 2), .5 + i + (i > 2)), c(barplotdata$lse.rich[i], barplotdata$use.rich[i]), col = "black", lwd = 2)
}
dev.off()

pdf(file = paste0("/Users/", user, "/Dropbox/Work/Iquitos/Specialist_responses/flood_specialist_abun.pdf"), width = 5.2, height = 5)
barplot(height = matrix(barplotdata$mean.abun, nrow = 2), beside = T, ylim = c(0,37), names.arg = c("floodplain", "terra firme"),
        col = c("gray60", "gray90"), main = "floodplain specialist abundance", xlab = "forest-type", legend = c("primary forest", "agriculture"))
for(i in 1:4){
  lines(c(.5 + i + (i > 2), .5 + i + (i > 2)), c(barplotdata$lse.abun[i], barplotdata$use.abun[i]), col = "black", lwd = 2)
}
dev.off()


###### Plotting veg data ######
veg_data <- read.csv("/Users/JacobSocolar/Dropbox/Work/Iquitos/Data/Transect_locations.csv", fileEncoding = "latin1")
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
                    CCFmean = c(mean(veg_dd$CCF[veg_dd$riverNS == "N" & veg_dd$Habitat != "White-sand points"]),
                                mean(veg_dd$CCF[veg_dd$riverNS == "S" & veg_dd$Habitat != "White-sand points"])),
                    CCFse = c(se(veg_dd$CCF[veg_dd$riverNS == "N" & veg_dd$Habitat != "White-sand points"]),
                              se(veg_dd$CCF[veg_dd$riverNS == "S" & veg_dd$Habitat != "White-sand points"])))

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


pdf(file = paste0("/Users/", user, "/Dropbox/Work/Iquitos/Specialist_responses/NS_ccf_cover.pdf"), width = 5.2, height = 5)
barplot(height = vegNS$CCFmean, ylim = c(0,45), names.arg = c("north-bank", "south-bank"),
        col = c("gray70", "gray70"), main = "across the Amazon", ylab = "closed-canopy forest (% cover)")
for(i in 1:2){
  lines(c(i - .3 + 0.2*(i > 1), i - .3 + 0.2*(i > 1)), c(vegNS$lse[i], vegNS$use[i]), col = "black", lwd = 2)
}
dev.off()

pdf(file = paste0("/Users/", user, "/Dropbox/Work/Iquitos/Specialist_responses/Hab_ccf_cover.pdf"), width = 5.2, height = 5)
barplot(height = vegHab$CCFmean, ylim = c(0,50), names.arg = c("floodplain", "terra firme", "white-sands"),
        col = rep("gray70", 3), main = "across forest-types", ylab = "closed-canopy forest (% cover)")
for(i in 1:3){
  lines(c(i - .3 + 0.2*(i - 1), i - .3 + 0.2*(i - 1)), c(vegHab$lse[i], vegHab$use[i]), col = "black", lwd = 2)
}
dev.off()
