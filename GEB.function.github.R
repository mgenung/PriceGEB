# input is a matrix object with species as rows and sites as columns
# if you are using this version of the Price equation as it was used in the GEB paper,
# then cell values should be measurements of ecosystem function

#BCI2 = t(as.matrix(BCI))[,1:6]

geb.price = function(x) {
  nplots    <- c() # empty vector for number of spatial or temporal reps in each dataset

  # number of species
  nsp    <- nrow(x)
  
  # number of 'plots' (or years etc.)
  nplots <- ncol(x)
  
  # empty Price terms
  srel <- matrix(NA, nplots, nplots) # richness loss
  sreg <- matrix(NA, nplots, nplots) # richness gain
  scel <- matrix(NA, nplots, nplots) # composition loss
  sceg <- matrix(NA, nplots, nplots) # composition gain
  CDE <- matrix(NA, nplots, nplots) # context dependence effect

  # put columns in order of lowest to highest function bc that's how we define baseline
  x <- x[,rev(order(colSums(x)))]
  
  # i+j loop does all pairwise comparisons
  for(i in 1:(nplots-1)){
    
    for(j in ((i+1):nplots)){
      
      z        <- x[,i] # EF by species at site 1
      z_p      <- x[,j] # EF by species at site 2
      zbar     <- mean(x[which(x[,i] != 0),i]) # mean species EF at site 1
      zbar_p   <- mean(x[which(x[,j] != 0),j]) # mean species EF at site 2
      sc       <- length(which(x[,i] > 0 & x[,j] > 0)) # no. of shared species
      s        <- length(which(x[,i] > 0)) # richness at site 1
      s_p      <- length(which(x[,j] > 0)) # richness at site 2
      w        <- x[,i]; w[which(w > 0)] <- 1 # pres/abs for site 1
      w_p      <- x[,j]; w_p[which(w_p > 0)] <- 1 # pres/abs for site 2
      wbar     <- mean(w[which(x[,i] != 0 | x[,j] != 0)]) # mean pres/abs for site 1
      wbar_p   <- mean(w_p[which(x[,j] != 0 | x[,i] != 0)]) # mean pres/abs for site 2
      deltaz   <- z-z_p # changes in species' EF between sites
      both_pa  <- w*w_p # pres/abs vector, only pres if at both sites
      zbarc    <- mean(x[which(both_pa == 1),i]) # mean EF of site 1 sp if they are pres at site 2
      zbarc_p  <- mean(x[which(both_pa == 1),j]) # mean EF of site 2 sp if they are pres at site 1
      
      srel[i,j] <- zbar*(sc-s)/sum(z) # Richness-Loss
      sreg[i,j] <- zbar_p*(s_p-sc)/sum(z) # Richness-Gain
      scel[i,j] <- sc*(zbarc-zbar)/sum(z) # Composition-Loss
      sceg[i,j] <- -(sc*(zbarc_p-zbar_p))/sum(z) # Composition-Gain
      CDE[i,j] <- sc*(zbarc_p-zbarc)/sum(z) # Abundance
    }
  }
  
  # consolidated output
  main.out <- data.frame("1_RL" = as.vector(srel),
                         "2_RG" = as.vector(sreg),
                         "3_CL" = as.vector(scel),
                         "4_CG" = as.vector(sceg),
                         "5_R+C" = as.vector(srel+sreg+scel+sceg),
                         "6_AB" = as.vector(CDE)
  )
  
  return(main.out)
}

geb.price(BCI2)
