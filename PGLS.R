##1. Install Packages/Set Working Directory##
install.packages("ape")
install.packages("geiger")
install.packages("nlme")
install.packages("phytools")
install.packages("readr")
library(readr)
library(ape)
library(geiger)
library(nlme)
library(phytools)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue", force = TRUE)
install.packages("qvalue")
library(qvalue)
##2.Setup###
getwd()
orthogroups <- read.csv("../Results_Mar03/Orthogroups/Orthogroups.GeneCount.no.Y.csv") #Removed OGs mapped to the Y-chr.
orthogroups <- as.data.frame(orthogroups)
orthogroups$Variance <- apply(orthogroups[, 2:97], 1, var)
orthogroups_fixed <- subset(orthogroups, Variance > 0)
orthogroups <- orthogroups_fixed
traits <- read.csv("../Traits.csv") #Traits.csv is a reduced form of the sample information table. Only trait info.
traits$Habitat_cont <- rep(NA, 96)
#Change habitat from numerical to categorical.
for(i in 1:nrow(traits)){
  if(traits$Habitat[i] == "Terrestrial"){
    traits$Habitat_cont[i] <- 0
  } else if (traits$Habitat[i] == "Amphibious" ){
    traits$Habitat_cont[i] <- .5
  } else {traits$Habitat_cont[i] <- 1}
}

#tree_1 <- read.tree(text = "(Platypus:0.157723,(((African_elephant:0.0303445,(Yellow_spotted_hyrax:0.00416229,Hyrax:0.00571213)1:0.0557999)1:0.00356936,Dugong:0.0309277)1:0.0296162,((((Capybara_2:0.0860688,Hazel_dormouse:0.0672463)1:0.00625397,((American_beaver_2:0.0448253,Ords_kangaroo_rat:0.0793509)1:0.009123,(Lesser_egyptian_jerboa:0.0701534,(((Mice:0.0122298,Southern_african_pygmy_mouse_2:0.0143725)1:0.0145737,(European_harvest_mouse:0.0264934,Norway_rat:0.0251374)1:0.00331957)1:0.0218472,((Creeping_vole:0.00573036,Water_vole_arvicoloides:0.00676778)1:0.0283558,Golden_hamster:0.033128)1:0.0102822)1:0.0429343)1:0.0177497)1:0.00746546)1:0.0144452,(Philippine_flying_lemur:0.0454054,((Ring_tailed_lemur:0.0302977,(Small_eared_galago:0.00580689,Moholi_bushbaby:0.00844716)1:0.0408986)1:0.0135429,(Guinea_baboon:0.0147543,((Western_lowland_gorilla:0.00380359,(Pygmy_chimpanzee:0.00295066,Human:0.0031733)1:0.000858898)1:0.00326683,Bornean_orangutan:0.00708194)1:0.00497629)1:0.0315969)1:0.00606256)1:0.0023505)1:0.0080641,((European_hedgehog:0.0954306,(Etruscan_shrew:0.0763853,Common_shrew:0.0564348)1:0.0541364)1:0.0210407,((Large_flying_fox:0.0443224,(Brazillian_free_tailed_bat:0.0310468,Little_brown_bat:0.0484963)1:0.0176202)1:0.0126661,((Chinese_pangolin:0.0614831,(((Maned_wolf:0.00504278,(Arctic_fox:0.00204756,Corsac_fox:0.00151803)1:0.00330194)1:0.0281012,((Giant_panda:0.0100186,(Polar_bear:0.00233426,American_black_bear:0.0019886)1:0.00512402)1:0.0124601,(((Least_weasel:0.00282877,Domestic_ferret:0.0031188)1:0.00668264,(Northern_Sea_Otter:0.00441857,Eurasian_river_otter:0.00538625)1:0.00399938)1:0.0195228,(((Guadaulpe_fur_seal:0.00189954,California_sea_lion:0.00215946)1:0.00407239,Pacific_walrus:0.00480466)1:0.00551529,((Northern_elephant_seal:0.00447166,Weddell_seal:0.00875675)1:0.00286709,((Baikal_seal:0.00464198,Gray_seal:0.00294039)1:0.000137378,(Harbor_seal:0.00169868,Ringed_seal:0.00193619)0.997:0.000192505)1:0.00352331)1:0.00479884)1:0.00933463)1:0.00147568)1:0.00559052)1:0.00638111,(((Pallas_cat:0.00432203,(Black_footed_cat:0.0018024,Domestic_cat:0.00182584)1:0.000900561)1:0.00115431,(Tiger:0.00177703,(Jaguar:0.0023221,Lion:0.00148124)1:0.000381671)1:0.00234517)1:0.0131795,Meerkat:0.0316994)1:0.0134595)1:0.0164107)1:0.00373,(Horse:0.0398269,((Pig:0.0434115,((((NA_right_whale:0.00509841,(Minke_whale:0.00412006,((Grey_whale:0.0038828,(Fin_whale:0.00439345,Humpback_whale:0.00302962)1:0.000845165)0.623:0.000175204,(Blue_whale:0.00389468,Rice_whale:0.00286431)1:0.000692764)1:0.000417649)1:0.00265973)1:0.00449639,((Pygmy_sperm_whale:0.00946746,Sperm_whale:0.00843791)1:0.00576859,((Northern_bottle_nose_whale:0.00337571,(Blainville_beaked_whale:0.0033694,(True_beaked_whale:0.00294494,Sowerby_beaked_whale:0.00281306)0.994:0.000294155)1:0.000919245)1:0.0074585,(((Harbor_porpoise:0.00148882,Vaquita:0.00252016)1:0.0035693,(Beluga_whale:0.00194119,Narwhal:0.00210991)1:0.00236895)1:0.0016634,(Atlantic_white_sided_dolphin:0.00310169,(((Long_finned_pilot_whale:0.00134586,False_killer_whale:0.00134843)1:0.00108871,Common_bottlenose_dolphin:0.00222982)1:0.000484018,Orca:0.00274291)1:0.00053767)1:0.00292883)1:0.00586348)1:0.00179482)1:0.00152093)1:0.0158386,(Pygmy_hippopotamus:0.00291176,Hippo:0.00209541)1:0.0262509)1:0.00306419,((Pronghorn:0.0186311,(Okapi:0.00949471,Kordofan_giraffe:0.00645096)1:0.00712447)1:0.00113447,(((Sheep:0.00480104,Muskox:0.00471037)1:0.00807131,(African_buffalo:0.00887689,(Bison:0.00540638,Cow:0.00222812)1:0.00406645)1:0.00519074)1:0.00407506,((Chinese_water_deer:0.0111018,Eurasian_elk:0.00739978)1:0.00142005,(Reeves_muntjac:0.00747273,((Red_deer:0.00116039,Wapiti:0.00104633)1:0.000675627,Sika_deer:0.00208155)1:0.00417629)1:0.00149236)1:0.00824109)1:0.00138005)1:0.0319161)1:0.00615668)1:0.00363209,(Arabian_camel:0.00140975,(Bactrian_camel:0.00521,Wild_bactrian_camel:0.00174366)1:0.000727676)1:0.0359946)1:0.0139)1:0.00183385)1:0.00118264)1:0.00293668)1:0.00812629)1:0.01003)1:0.157723);")
tree_1 <- read.tree("../Results_Mar03/Species_Tree/SpeciesTree_rooted.txt")
##3. Tree scaling and restructuring###
plot(tree_1)
#Calibrating the tree at certain nodes
nodes<-c(findMRCA(tree_1,c("Platypus","Cow")),
         findMRCA(tree_1,c("African_elephant","Mice")),
         findMRCA(tree_1,c("Human","Mice")),
         findMRCA(tree_1,c("Cow","Common_shrew")),
         findMRCA(tree_1,c("Horse","Orca")),
         findMRCA(tree_1,c("Human","Pygmy_chimpanzee")),
         findMRCA(tree_1,c("Blue_whale","Orca")),
         findMRCA(tree_1,c("Lion","Polar_bear")),
         findMRCA(tree_1,c("Etruscan_shrew","Common_shrew")))
age.min=c(163.7,94.7,81.3,83.3,72.8, 6,32.3,52.9,23.8)
age.max=c(185.9,101.9,91,88,78.3,6.5,35.2,57.3, 46.6)
calibration<-makeChronosCalib(tree_1,node=nodes,
                              age.min=age.min,age.max=age.max)#Making the Tree ultra metric based on the calibration. 
#pl.tree<-chronos(tree_1,calibration=calibration)
#plot(pl.tree)
#tree_ultra_.01 <- chronos(tree_1, model = "relaxed", lambda = .01, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
#tree_ultra_.1 <- chronos(tree_1, model = "relaxed", lambda = .1, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
#tree_ultra_1 <- chronos(tree_1, model = "relaxed", lambda = 1, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
#tree_ultra_100 <- chronos(tree_1, model = "relaxed", lambda = 100, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
#tree_ultra_1000 <- chronos(tree_1, model = "relaxed", lambda = 1000, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
tree_ultra <- chronos(tree_1, model = "relaxed", lambda = 1, calibration = calibration, control = chronos.control(iter.max = 10000, eval.max = 10000))
write.tree(tree_ultra, file = "my_tree.tre")
plot(tree_ultra)

##4. PGLS
#Setting up the matrix that will hold the results.
mass_OGs <- matrix(nrow = nrow(orthogroups), ncol = 7) #nrow(remove_OGS_test)
OGs <- nrow(orthogroups)
colnames(mass_OGs) <- list("Orthogroup", "Intercept", "Slope_log_mass", "Slope_Habitat_cont", "P_value_log_mass", "P_value_Habitat_cont", "Alpha")
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(-5, 10), 
     ylim = c(0, 30))
#The loop that runs through all OGs.
for(i in 1:OGs){
  #Setting up the data for PGLS
  test1 <- orthogroups[i,]
  test1 <- test1[,c(2:97)]
  test2 <- t(test1)
  test3 <- as.data.frame(test2)
  test3[,2] <- row.names(test3)
  colnames(test3) <- c("counts", "Sample")    
  test4<- merge(test3, traits, by = "Sample")
  test4$Sample <- as.factor(test4$Sample)
  rownames(test4) <- test4$Sample
  match(tree_1$tip.label, test4$Sample)->VariablesSort
  test4[VariablesSort,]-> Variables
  Species <- Variables$Sample
  ###Strength of stabilizing section. Try multiple alpha values and compare the models AIC.
  alpha_starts <- c(0.01, 0.1, 1, 10, 100, 1000)
  gls_results <- list()
  for (a in alpha_starts) {
    
    # tryCatch allows the loop to keep going even if one model fails
    model_fit <- tryCatch({
      gls(counts ~ Log_Mass + Habitat_cont,
          data = Variables, 
          correlation = corMartins(a, tree_ultra, form = ~Species, fixed = TRUE),
          method = "ML") # Use ML to allow for AIC comparison
    }, error = function(e) {
      message(paste("  Error at alpha", a, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(model_fit)) {
      gls_results[[as.character(a)]] <- model_fit
    }
  }
  if (length(gls_results) > 0) {
    aic_table <- data.frame(
      Start_Alpha = names(gls_results),
      AIC = sapply(gls_results, AIC)
    )# Sort by best fit
  } else {
    stop("All alpha values failed to converge.")
  }
  
  a <- as.numeric(aic_table[order(aic_table$AIC), ][1,1])
###Using the alpha value, fit the PGLS
  
  fit1 <- gls(counts ~ Log_Mass + Habitat_cont,
              data = Variables, correlation=corMartins(a,tree_ultra, form = ~Species, fixed = TRUE))
  sum_model <- summary(fit1)
  p_value_log_mass <- coef(sum_model)[2,4]
  p_value_Habitat_cont <- coef(sum_model)[3,4]
  slope_log_mass <- coef(sum_model)[2,1]
  slope_Habitat_cont <- coef(sum_model)[3,1]
  intercept <- coef(sum_model)[1,1]
  mass_OGs[i,1] <- orthogroups[i,1]
  mass_OGs[i,2] <- intercept
  mass_OGs[i,3] <- slope_log_mass
  mass_OGs[i,4] <- slope_Habitat_cont
  mass_OGs[i,5] <- p_value_log_mass 
  mass_OGs[i,6] <- p_value_Habitat_cont
  mass_OGs[i,7] <- a
  abline(sum_model)
  print(i)}
##5. Q-value##
#Do q-value correction.
trait_OGs <- as.data.frame(mass_OGs)
mass_P_value <- trait_OGs$P_value_log_mass
Habitat_P_value <- trait_OGs$P_value_Habitat_cont
mass_P_value <- as.numeric(mass_P_value)
Habitat_P_value <- as.numeric(Habitat_P_value)
hist(mass_P_value)
hist(Habitat_P_value)
q_value_mass <- qvalue(p = mass_P_value)
q_value_Habitat <- qvalue(p = Habitat_P_value)
summary(q_value_mass)
summary(q_value_Habitat)
trait_OGs$q_value_mass  <-(q_value_mass$qvalues)
trait_OGs$q_value_Habitat  <-(q_value_Habitat$qvalues)


