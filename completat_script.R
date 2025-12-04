
# Definim matriu M (ls SNPs)
M <- matrix(c(
  2,0,1,1,0,0,0,2,1,2,
  1,0,0,0,0,2,0,2,1,0,
  1,1,2,1,1,0,0,2,1,2,
  0,0,2,1,0,1,0,2,2,1,
  0,1,1,2,0,0,0,2,1,2,
  1,1,0,1,0,2,0,2,2,1,
  0,0,1,1,0,2,0,2,2,0,
  0,1,1,0,0,1,0,2,2,0,
  2,0,0,0,0,1,2,2,1,2,
  0,0,0,1,1,2,0,2,0,0,
  0,1,1,0,0,1,0,2,2,1,
  1,0,0,0,1,1,0,2,0,0,
  0,0,0,1,1,2,0,2,1,0,
  1,0,1,1,0,2,0,1,0,0), nrow=14, byrow=T)


# 2. Definim funció run_BLUP_GBLUP
run_BLUP_GBLUP <- function(M){
  
  # Pedigrí i fenotip
  PED <- matrix(c(
    1,0,0, 2,0,0, 3,0,0, 4,0,0,
    5,0,0, 6,0,0, 7,0,0, 8,0,0, 9,0,0, 10,0,0,
    11,0,0, 12,0,0, 13,0,0, 14,0,0, 15,13,4,
    16,15,2, 17,15,5, 18,14,6, 19,14,9, 20,14,9,
    21,1,3, 22,14,8, 23,14,11, 24,14,10,
    25,14,7, 26,14,12), nrow=26, byrow=TRUE)
  
  FAT <- c(9,13.4,12.7,15.4,5.9,7.7,10.2,
           4.8,7.6,8.8,9.8,9.2,11.5,13.3)
  
  vara <- 35.241
  vare <- 245
  
  library(AGHmatrix)
  A <- Amatrix(PED)
  
  # 2a. SNP-BLUP (SR) ####################
  Z_sr <- M[-c(9:14), ]
  n_snp <- ncol(M)
  
  freq_sr <- colSums(Z_sr)/(2*nrow(Z_sr))
  for(j in 1:n_snp) Z_sr[,j] <- Z_sr[,j] - 2*freq_sr[j]
  
  het_sr <- sum(2*freq_sr*(1-freq_sr))
  
  X_sr <- matrix(1, 8, 1)
  I_sr <- diag(n_snp)
  y_sr <- FAT[1:8]
  alpha_sr <- het_sr*(vare/vara)
  
  LHS_sr <- crossprod(cbind(X_sr,Z_sr))
  neq <- ncol(LHS_sr)
  LHS_sr[2:neq, 2:neq] <- LHS_sr[2:neq, 2:neq] + I_sr*alpha_sr
  RHS_sr <- t(cbind(X_sr,Z_sr)) %*% y_sr
  
  sol_sr <- solve(LHS_sr, RHS_sr)
  
  # Valors millorants SNP-BLUP
  Z_ss <- M[-c(1:8), ]
  freq_ss <- colSums(Z_ss)/(2*nrow(Z_ss))
  for(j in 1:n_snp) Z_ss[,j] <- Z_ss[,j] - 2*freq_ss[j]
  a_ss <- Z_ss %*% sol_sr[2:(n_snp+1)]
  
  # 2b. GBLUP #########################
  Z_g <- M
  freq_g <- colSums(Z_g)/(2*nrow(Z_g))
  for(j in 1:n_snp) Z_g[,j] <- Z_g[,j] - 2*freq_g[j]
  
  het_g <- sum(2*freq_g*(1-freq_g))
  G_g <- (Z_g %*% t(Z_g)) / het_g
  G_gm <- G_g + diag(nrow(M))*0.01
  G_gmi <- solve(G_gm)
  
  X_g <- matrix(1,14,1)
  W_g <- diag(14)
  y_g <- FAT
  alpha_g <- vare/vara
  
  LHS_g <- crossprod(cbind(X_g,W_g))
  neq <- ncol(LHS_g)
  LHS_g[2:neq,2:neq] <- LHS_g[2:neq,2:neq] + G_gmi*alpha_g
  RHS_g <- t(cbind(X_g,W_g)) %*% y_g
  sol_g <- solve(LHS_g, RHS_g)
  
  sol_g_last6 <- tail(sol_g,6)
  
  # 2c. SSD ###########
  FAT_last6 <- tail(FAT,6)
  
  a_BLUP  <- a_ss + sol_sr[1]
  a_GBLUP <- sol_g_last6 + sol_g[1]
  
  SSD_BLUP  <- sum((FAT_last6 - a_BLUP)^2)
  SSD_GBLUP <- sum((FAT_last6 - a_GBLUP)^2)
  
  # OUTPUT #########
  list(SSD_BLUP=SSD_BLUP, SSD_GBLUP=SSD_GBLUP)
}

# Generem totes les combinacions d'eliminació de columnes (un total de n!)
ncol_mat <- ncol(M)
llista_resultats <- list()
k <- 1

for(mida in 1:ncol_mat){  # pots canviar 1:3 per 1:total SNPs si vols totes
  combs <- combn(1:ncol_mat, mida, simplify = FALSE)
  for(cols in combs){
    llista_resultats[[k]] <- M[,-cols,drop=FALSE]
    names(llista_resultats)[k] <- paste0("elim_", paste(cols, collapse="_")) # donar noms als possibles combinacions
    k <- k + 1
  }
}

# Calcul BLUP_GBLUP per cada combinacio de SNPs ############
# Calcul BLUP_GBLUP per cada combinacio de SNPs ############
# Create vectors to store results for CSV
combination_names <- c()
ssd_blup_values <- c()
ssd_gblup_values <- c()

for(i in seq_along(llista_resultats)){
  
  result <- tryCatch({
    run_BLUP_GBLUP(llista_resultats[[i]])
  }, error = function(e){
    return(NULL)
  })
  
  if(!is.null(result)){
    combination_names <- c(combination_names, names(llista_resultats)[i])
    ssd_blup_values <- c(ssd_blup_values, result$SSD_BLUP)
    ssd_gblup_values <- c(ssd_gblup_values, result$SSD_GBLUP)
  }
}
head(resultats_BLUP_GBLUP, 10)  # comprobar q ha funcionat solo mirar el top 10 files (sino hi ha ~3.000.000)

# escribim tot calculat al CSV
results_df <- data.frame(
  Combination = combination_names,
  SSD_BLUP = ssd_blup_values,
  SSD_GBLUP = ssd_gblup_values
)
write.csv(results_df, "blup_gblup_results.csv", row.names = FALSE)

# mirem els resultats per encontrar el millor BLUP, i el millor GLUP (els valors min)
element_BLUP  <- Inf
element_GBLUP <- Inf
index_BLUP  <- NA
index_GBLUP <- NA

for(i in seq_along(combination_names)){
  
  if(ssd_blup_values[i] < element_BLUP){
    element_BLUP <- ssd_blup_values[i]
    index_BLUP <- i
  }
  
  if(ssd_gblup_values[i] < element_GBLUP){
    element_GBLUP <- ssd_gblup_values[i]
    index_GBLUP <- i
  }
}



# Resultat finals
cat("Millor SSD BLUP:", element_BLUP, "\n")
cat("Columnes eliminades BLUP:", names(llista_resultats)[index_BLUP], "\n\n")

cat("Millor SSD GBLUP:", element_GBLUP, "\n")
cat("Columnes eliminades GBLUP:", names(llista_resultats)[index_GBLUP], "\n")
