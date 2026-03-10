##### IV - Merge aquamaps data and run SACs #####

library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(vegan)

load(file = "Data/fish_env_cleandf_230cells.rData")
fish_env <- fish_env2; rm(fish_env2)

load("Data/sti_values.RData")
load("Data/all_fish.RData")

sti <- sti %>%
  left_join(all_fish %>% select(Species, Genus), by = "Species")

spp_missing <- setdiff(unique(fish_env$Species), unique(sti$Species))

length(unique(all_fish$Species))
length(unique(fish_env$Species))

# Extract species names from all_fish df (from aquamaps) 
sp_names = sti %>% filter(Species %in% unique(fish_env$Species)) %>% select(SpeciesID, Species)
sp_names = sp_names[!duplicated(sp_names[,2]),]

sp_envelope <- sti %>% 
  filter(SpeciesID %in% sp_names$SpeciesID)

sp_envelope <- sp_envelope[,c("TempPrefMin", "TempPrefMax", "TempMin", "TempMax",
                              "SalinityPrefMin", "SalinityPrefMax", "SalinityMin", "SalinityMax",
                              "OxyPrefMin", "OxyPrefMax", "OxyMin", "OxyMax",
                              "Species", "SpeciesID")]

# Merge env preferences with fish data
fish_env2 <- left_join(fish_env, sp_envelope, by = "Species")

# check how many species have data
fish_env2 %>%
  dplyr::filter(!is.na(TempPrefMin) & !is.na(TempPrefMax)) %>%
  dplyr::summarise(count_species = n_distinct(Species))
length(unique(fish_env2$Species))

# calculate  width per specie
fish_env2$PrefSalwidth <- fish_env2$SalinityPrefMax - fish_env2$SalinityPrefMin
fish_env2$PrefTempwidth <- fish_env2$TempPrefMax - fish_env2$TempPrefMin
fish_env2$Prefo2width <- fish_env2$OxyPrefMax - fish_env2$OxyPrefMin

fish_env2$Salwidth <- fish_env2$SalinityMax - fish_env2$SalinityMin
fish_env2$Tempwidth <- fish_env2$TempMax - fish_env2$TempMin
fish_env2$o2width <- fish_env2$OxyMax - fish_env2$OxyMin

save(fish_env2, file = "Data/fish_aquamaps_env_230cells.rData")

# transform negative values of depth
fish_env2$depth <- abs(fish_env2$depth)
fish_env2$cell_decade <- paste(fish_env2$cell, fish_env2$decade, sep = "_")

# exclude cells with less than 2 hauls in a decade
filtered_data <- fish_env2 %>%
  dplyr::group_by(decade, cell) %>%
  dplyr::filter(n_distinct(haul_id) >= 2) %>%
  dplyr::ungroup()

#### 1. Fit SACs #####
source("R/SAC_functions.r")

samples <- filtered_data[,c("haul_id","Species","wgt_cpua","cell","year","cell_decade")]
summary(samples)
names(samples)
names(samples)[c(2,3)] <- c("species","biomass")

#### 2. create objects and fit SAC models ####
cells <- sort(unique(samples$cell_decade))
cells_sub <- cells[]

model=c('gleason','gitay','arrhenius','michaelis-menten')
yprim=vector()
identity=vector()
nb.hauls=vector()
nb.hauls.min=vector()
keep=vector()
sr=vector()
chao=vector()
limit0.8 <- limit0.75 <- limit0.65 <- vector()
asympt <- vector()

sac2 <- list()
model92 <- list()
x2 <- list()
y3 <- list()

for(i in 1:length(cells)){
  
  print(i)
  
  sub.surv <- subset(samples, samples$cell_decade %in% cells[i])
  if(length(unique(sub.surv$haul_id))>=3){
    tryCatch({
      
      # Transform into a matrix for SAC function
      sample.mat <- dcast(sub.surv, haul_id~species, sum, value.var='biomass')
      sample.mat <- sample.mat[!is.na(sample.mat$haul_id), ]
      rownames(sample.mat) <- sample.mat[,1]
      sample.mat[,1] <- NULL
      sample.mat[is.na(sample.mat)] <- 0
      
      # Get the unobserved total number of species
      unobs <- specpool(sample.mat)
      sample.mat <- round(sample.mat)
      sac <- specaccum(sample.mat, method = "random", permutations=10)
      sac2[i] <- list(sac)
      
      #rare <- specaccum(sample.mat, method='rarefaction', xvar='individuals')
      model9 <- tryCatch({
        fitspecaccum(sac, model = 'michaelis-menten', method = "random")
      }, error = function(e) {
        message("Michaelis-Menten failed, trying logis")
        tryCatch({
          fitspecaccum(sac, model = 'logis', method = "random")
        }, error = function(e) {
          message("Logis failed, trying Gleason")
          tryCatch({
            fitspecaccum(sac, model = 'gleason', method = "random")
          }, error = function(e) {
            message("Gleason failed, trying Arrhenius")
            tryCatch({
              fitspecaccum(sac, model = 'arrhenius', method = "random")
            }, error = function(e) {
              message("All models failed.")
              NULL  # Retorna NULL se nenhum modelo for ajustado
            })
          })
        })
      })
      
      
      model92[i] <- list(model9)
      
      fitted9 <- fitted(model9)
      aic9 <- sapply(model9$models, AIC)
      num9 <- which(aic9==min(aic9))
      result9 <- fitted9[,num9]
      means9 <- apply(fitted9, 1, FUN=mean)
      w <- coef(model9)
      v <- apply(w, 1, FUN=mean)
      x <- seq(from=1, to=nrow(fitted9), by=1)
      x2[i] <- list(x)
      
      
      ### Get fit function
      y <- curves(v[1], v[2], x, model[u])
      y3[i] <- list(y)
      
      ### Get first derivative at max nhauls
      y2 <- first.der(v[1],v[2],x, model[u])
      
      plot(sac, xlab='', ylab='')
      plot(model9, add=TRUE, col='grey')
      plot(sac, xlab='', ylab='', add=TRUE)
      points(y~x, type='l', lwd=5, col='orange')
      
      ### Get all values
      if(length(y2)>0){
        yprim[length(yprim)+1] <- y2[length(y2)]
        identity[length(identity)+1] <- as.character(cells[i])
        nb.hauls[length(nb.hauls)+1] <- nrow(sample.mat)
        chao[length(chao)+1] <- round(unobs$chao)
        sr[length(sr)+1] <- unobs$Species
        asympt[length(asympt)+1] <- v[1]
        
        ### Get nb of hauls necessary to reach % of chao 
        limit <- 0.8*round(v[1])
        obj <- which(y>limit)
        
        if(length(obj)>0){limit0.8[length(limit0.8)+1] <- obj[1]}
        else{limit0.8[length(limit0.8)+1] <- NA}
        
        limit <- 0.75*round(v[1])
        obj <- which(y>limit)
        
        if(length(obj)>0){limit0.75[length(limit0.75)+1] <- obj[1]}
        else{limit0.75[length(limit0.75)+1] <- NA}
        
        limit <- 0.65*round(v[1])
        obj <- which(y>limit)
        
        if(length(obj)>0){limit0.65[length(limit0.65)+1] <- obj[1]}
        else{limit0.65[length(limit0.65)+1] <- NA}
        
        #Keep or not if y'<0.65 meaning decreasing trend
        x2 <- seq(from=1, to=100, by=1)
        y0 <- curves(v[1], v[2], x2, model[u])
        yprim2 <- first.der(v[1],v[2],x2, model[u])
        lets <- which(yprim2<0.65)
        if(length(lets)>0){
          if(lets[1] <= length(x)){keep[length(keep)+1] <- 'YES'}
          if(lets[1] > length(x)){keep[length(keep)+1] <- 'NO'}
          
          #Min number of hauls necessary
          nb.hauls.min[length(nb.hauls.min)+1]  <-lets[1]                          
          
        }
        
      }
      
      if(length(lets)==0){keep[length(keep)+1] <- 'NO'
      nb.hauls.min[length(nb.hauls.min)+1] <- NA}
      rm(y2,sample.mat)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})                    
  }
}

mat.deriv <- as.data.frame(cbind(identity, yprim, nb.hauls,nb.hauls.min, keep,sr,chao,limit0.8, limit0.75, limit0.65, asympt))
mat.deriv$yprim <- as.numeric(as.vector(mat.deriv$yprim))
mat.deriv$nb.hauls <- as.numeric(as.vector(mat.deriv$nb.hauls))
mat.deriv$nb.hauls.min <- as.numeric(as.vector(mat.deriv$nb.hauls.min))
mat.deriv$limit0.75 <- as.numeric(as.vector(mat.deriv$limit0.75))

save(mat.deriv, file = "Data/mat.deriv.230cell_decade.RData")

mat.deriv75 <- mat.deriv[!is.na(mat.deriv$limit0.75),]

# Filter mat.deriv75 to keep only rows where "keep" is "YES"
mat.deriv75 <- subset(mat.deriv75, keep == "YES")

# Join fish_env2 with the filtered mat.deriv75 by matching "id" in fish_env2 to "identity" in mat.deriv75
filtered_data2 <- filtered_data[filtered_data$cell_decade %in% mat.deriv75$identity, ]

# make a vector with the number of hauls needed for each grid cell to reach 75% of asymptotic species richness
nb.hauls.min <- mat.deriv75$limit0.75

# Empty vector to store first loop results
testlist <- vector(mode="list")

cells2 <- sort(unique(filtered_data2$cell_decade))

# For loop randomly selecting X number of hauls 99 times for each grid cell
for(i in 1:length(cells2)){
  print(i)
  sub.surv <- subset(filtered_data2, filtered_data2$cell_decade %in% cells2[i])
  length(unique(sub.surv$haul_id))
  unique_hauls <- unique(sub.surv$haul_id)
  
  y <- replicate(999,sample(unique_hauls, nb.hauls.min[i], replace=FALSE))
  duplicated(y)
  as.array(y)
  testlist[[i]] <- y
  
}

# converting each list into a dataframe structure
N <- seq(1:length(unique(filtered_data2$cell_decade)))
cell_names <- data.frame(matrix(ncol = 999, nrow = 0))

for (i in 1:length(N)){
  print(i)
  
  y <- testlist[i]
  gh <- data.frame(y)
  
  if(ncol(gh) == 1){
    gh <- t(gh)
    colnames(gh) <- colnames(cell_names)[1:ncol(gh)]
  }

  cell_names <- rbind(cell_names, gh)
}

test <- sac2
test2 <- model92

##### 3. plot all SAC curves ####
allcurves <- list()

for(m in seq_along(test)){
  
  df <- data.frame(sites = test[[m]]$sites)
  df$richness <- test[[m]]$richness
  df$sd <- test[[m]]$sd
  df$id <- rep(m, nrow(df))
  allcurves[m] <- list(df)
  
  print(m)
}

allmodels <- list()

for(n in seq_along(test2)){
  
  df2 <- data.frame(sites = test2[[n]]$sites)
  df2$richness <- test2[[n]]$richness
  df2$sd <- test2[[n]]$sd
  df2$id <- rep(n, nrow(df2))
  allmodels[n] <- list(df2)
  
  print(n)
} 

b <- bind_rows(allcurves)
b$sac <- "sac"
c <- bind_rows(allmodels)
c$sac <- "model"
d <- rbind(b,c)

cells2 <- as.data.frame(cells2)
cells2$id <- seq(1:590)

d <- merge(d, cells2, by = "id")

ggplot(d, aes(x = sites, y = richness, group = interaction(cells2, sac))) +
  geom_line() +
  geom_point(aes(col = as.factor(cells2))) +
  geom_ribbon(aes(
    ymin = richness - 2 * sd,
    ymax = richness + 2 * sd,
    group = interaction(cells2, sac)
  ), alpha = 0.05) +
  theme_bw() +
  labs(x = "Number of hauls", y = "Richness", color = "Sample unit") +
  theme(legend.position = "none")

save(cell_names, file = "Data/cell_names_hex_230cells.rData")
save(filtered_data2, file = "Data/filtered_data2_sac_hex_230cells.rData")
