setwd("/Users/Jacob/Dropbox/Work/Iquitos/Data")
'%ni%' <- Negate('%in%')

cdata <- read.csv("cdata.csv", stringsAsFactors = F)
stotz <- read.csv("/Users/Jacob/Dropbox/Work/Useful_data/Stotz_et_al/stotz_adata.csv")
extras <- read.csv("stotz_cdata_additions.csv")
stotz.nb <- read.csv("/Users/Jacob/Dropbox/Work/Useful_data/Stotz_et_al/Stotz_nearctic_nonbreeding.csv")
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
stotz <- read.csv("stotz_cdata_taxonomy.csv")

# habitat lumping
stotz$flood.for <- stotz$F2 %in% c("Y") | stotz$F3 %in% c("Y") | stotz$F13 %in% c("Y")
stotz$flood.nfor <- stotz$N11 %in% c("Y") | stotz$N12 %in% c("Y") | stotz$A1 %in% c("Y") | stotz$A5 %in% c("Y") | 
  stotz$A6 %in% c("Y") | stotz$A8 %in% c("Y") | stotz$A9 %in% c("Y")
stotz$upland.for <- stotz$F1 %in% c("Y") | stotz$F12 %in% c("Y")
stotz$forest <- (stotz$flood.for == 1) | (stotz$upland.for == 1)
stotz$tf <- stotz$F1 %in% c("Y")

FFsp <- as.character(stotz$sp[stotz$flood.for & !stotz$upland.for])
cFFsp <- FFsp[FFsp %in% cdata$Species]

FNsp <- as.character(stotz$sp[stotz$flood.nfor & !stotz$upland.for])
cFNsp <- FSsp[FSsp %in% cdata$Species]

allfloodsp <- unique(c(cFFsp, cFNsp))

TFsp <- as.character(stotz$sp[stotz$tf & !(stotz$flood.for | stotz$flood.nfor)])
cTFsp <- TFsp[TFsp %in% cdata$Species]

tfexclude <- read.csv("tf_exclude.csv", header = F)
cTFsp <- cTFsp[cTFsp %ni% tfexclude[,1]]

Forestsp <- stotz$sp[stotz$forest]
cForestsp <- Forestsp[Forestsp %in% cdata$Species]

spec <- read.csv("specialists.csv")
extras <- read.csv("flood_specialist_extras.csv")

spec[is.na(spec)] <- 0

specialists <- data.frame(species=spec$Species, poor=((spec$WhiteSand.Alvarez+spec$Poorsoil.Pomara+spec$Poorsoil.extra)>0),
                          rich=((spec$Richsoil.Alvarez + spec$Richsoil.Pomara)>0), river = spec$RiverLimit)
specialists$flood <- 0
specialists$flood[specialists$species %in% allfloodsp] <- 1
specialists$flood[specialists$species %in% extras$Extras] <- 1
specialists$flood.forest <- 0
specialists$flood.forest[specialists$species %in% cFFsp] <- 1
specialists$flood.forest[specialists$species %in% extras$Extras[which(!is.na(extras$Forest))]] <- 1
specialists$tf <- 0
specialists$tf[specialists$species %in% cTFsp] <- 1
specialists$forest <- 0
specialists$forest[specialists$species %in% cForestsp] <- 1
specialists$flood.nonforest <- specialists$flood - specialists$flood.forest


traits <- read.delim("/Users/Jacob/Dropbox/Work/Useful_data/EltonTraits/BirdFuncDat.txt", header=T, stringsAsFactors = F)   # traits data (includes diet)
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

all_vars_stanMER <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                      stratum + diet + scale.mass + 
                                      forest + river + flood.forest + 
                                      flood.nonforest + tf + poor + 
                                      rich, family = 'binomial', data = specialists, seed = 349)
summary(all_vars_stanMER)[1:17,]

all_vars_stanMER_2 <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                           ForStrat.ground + ForStrat.understory + ForStrat.midhigh + 
                                           diet + scale.mass + 
                                           forest + river + flood.forest + 
                                           flood.nonforest + tf + poor + 
                                           rich, family = 'binomial', data = specialists, seed = 349)
summary(all_vars_stanMER_2)[1:17,]

functional_vars_stanMER <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                                  stratum + diet + scale.mass + 
                                                  forest, family = 'binomial', data = specialists, seed = 349)
summary(functional_vars_stanMER)[1:11,]

functional_vars_stanMER_2 <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                                  ForStrat.ground + ForStrat.understory + ForStrat.midhigh +
                                                  diet + scale.mass + forest, 
                                                family = 'binomial', data = specialists, seed = 349)
summary(functional_vars_stanMER_2)[1:11,]

hab_vars_stanMER <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                           forest + river + flood.forest + 
                                           flood.nonforest + tf + poor + 
                                           rich, family = 'binomial', data = specialists, seed = 349)
summary(hab_vars_stanMER)[1:8,]

functional_vars_stanMER_3 <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                                  stratum + diet + scale.mass, 
                                                  family = 'binomial', data = specialists, seed = 349)
summary(functional_vars_stanMER_3)[1:11,]

functional_vars_stanMER_4 <- rstanarm::stan_glmer(formula = abun.double ~ (1|species) + 
                                                    ForStrat.ground + ForStrat.understory + ForStrat.midhigh +
                                                    diet + scale.mass, 
                                                  family = 'binomial', data = specialists, seed = 349)
summary(functional_vars_stanMER_4)[1:11,]


kfold_all_vars <- rstanarm::kfold(all_vars_stanMER)
kfold_all_vars_2 <- rstanarm::kfold(all_vars_stanMER_2)
kfold_functional_vars <- rstanarm::kfold(functional_vars_stanMER)
kfold_functional_vars_2 <- rstanarm::kfold(functional_vars_stanMER_2)
kfold_functional_vars_3 <- rstanarm::kfold(functional_vars_stanMER_3)
kfold_functional_vars_4 <- rstanarm::kfold(functional_vars_stanMER_4)
kfold_hab_vars <- rstanarm::kfold(hab_vars_stanMER)

rstanarm::compare_models(kfold_all_vars, kfold_all_vars_2, kfold_functional_vars,
                         kfold_functional_vars_2, kfold_functional_vars_3, kfold_functional_vars_4,
                         kfold_hab_vars)
rstanarm::compare_models(kfold_all_vars, kfold_hab_vars)

summary(all_vars_stanMER)[c(185,193,194,263,265,374,446),]
