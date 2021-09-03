## function to generate model matrix manually 
require(stringr)

generate.model.matrix <- function(data, scale.X = TRUE) {
  # lets figure out who the potential competitors are
  focals <- as.character(unique(data$Focal))
  non.X.colnames <- c("Plot", "Year", "Tag.tree", "Focal", "DBH_mi", "G")
  competitors <- colnames(data)[!colnames(data) %in% non.X.colnames]

  # variables for alpha coefficients
  all.alphas <- sort(competitors)
  
  # variables for beta coefficients
  # betas between conspecific neighbors
  possible.intrabetas <- all.alphas
  possible.intrabetas <-
    unlist(lapply(possible.intrabetas, function(x) {
      # paste0("I((", x, "^2)/2)")
      paste0("I((", x, "^2))")
    }))
  # betas between heterospecific neighbors
  possible.interbetas <- combn(all.alphas, 2)
  possible.interbetas <-
    apply(possible.interbetas, 2, paste, collapse = ":")
  # combine all betas together into a single variable
  all.betas <- c(possible.intrabetas, possible.interbetas)

  # model formula (only the "fixed" part)
  # this is not for model fitting, but for data preparation only
  alpha.part <- paste(all.alphas, collapse = "+")
  beta.part  <- paste(all.betas, collapse = "+")
  formula.tmp <- as.formula(paste0("G ~ 0 +", alpha.part, "+", beta.part))
  
  # generate model matrix
  mm <- model.matrix(formula.tmp, data = data)
  # merge model matrix to original data.frame
  out <- cbind(data[, non.X.colnames], mm)
  
  # for each focal species, check for columns that sums to zero
  # these are neighbour/neighbour-combinations that never co-occured
  neighbour.check <- matrix(NA, nrow = ncol(mm), ncol = length(focals))
  colnames(neighbour.check) <- focals
  for (i in focals) {
    out.sub <- droplevels(subset(out, Focal == i))
    # neighbour.check[, i] <- colSums(out.sub[, colnames(mm)])
    # neighbour.check[, i] <- colSums(out.sub[, colnames(mm)] > 0)
    neighbour.check[, i] <-
      apply(out.sub[, colnames(mm)], 2, function(x) {
        diff(range(x[x>0]))
      })
  }
  if (all(neighbour.check > 0)) {
    cat("Neighbour co-occurrence check OK!")
  } else {
    cat("Some neighbours / neighbour combinations did not co-occur with the focal; do not proceed!")
  }
  
  # change X colnames to ensure regression function does not
  # self-generate model matrix again
  all.betas.newnames <- str_replace(all.betas, ":", "_")
  all.betas.newnames <- str_remove(all.betas.newnames, "I\\(\\(")
  # all.betas.newnames <- str_remove(all.betas.newnames, "\\)/2\\)")
  all.betas.newnames <- str_remove(all.betas.newnames, "\\)\\)")
  all.betas.newnames <- str_replace(all.betas.newnames, "\\^", "_")
  names(out) <- c(non.X.colnames, all.alphas, all.betas.newnames)
  
  # now mannually scale X
  if (scale.X) {
    # mm.scale <- scale(out[, !colnames(out) %in% non.X.colnames])
    alpha.center <- colMeans(mm[, all.alphas])
    intrabeta.center <- alpha.center ^ 2
    interbeta.center <- apply(combn(alpha.center, 2), 2, prod)
    all.betas.center <- c(intrabeta.center, interbeta.center)
    names(all.betas.center) <- all.betas.newnames
    all.center <- c(alpha.center, all.betas.center)
    mm.center <- t(apply(out[, names(all.center)], 1, function(x) x - all.center))
    all.scale <- apply(mm.center, 2, sd)
    mm.scale <- t(apply(mm.center, 1, function(x) x / all.scale))
    
    out[, names(all.center)] <- mm.scale
    return(list(in.dat = out, center = all.center, scale = all.scale))
  } else {
    return(list(in.dat = out, center = NULL, scale = NULL))
  }
}


# function to expand new data
expand.newdat <- function(data, scale.X = TRUE) {
  # lets figure out who the potential competitors are
  non.X.colnames <- c("Pred.type", "Focal", "DBH_mi")
  competitors <- colnames(data)[!colnames(data) %in% non.X.colnames]
  
  # variables for alpha coefficients
  all.alphas <- sort(competitors)
  
  # variables for beta coefficients
  # betas between conspecific neighbors
  possible.intrabetas <- all.alphas
  possible.intrabetas <-
    unlist(lapply(possible.intrabetas, function(x) {
      # paste0("I((", x, "^2)/2)")
      paste0("I((", x, "^2))")
    }))
  # betas between heterospecific neighbors
  possible.interbetas <- combn(all.alphas, 2)
  possible.interbetas <-
    apply(possible.interbetas, 2, paste, collapse = ":")
  # combine all betas together into a single variable
  all.betas <- c(possible.intrabetas, possible.interbetas)
  
  # model formula (only the "fixed" part)
  # this is not for model fitting, but for data preparation only
  alpha.part <- paste(all.alphas, collapse = "+")
  beta.part  <- paste(all.betas, collapse = "+")
  formula.tmp <- as.formula(paste0("~ 0 +", alpha.part, "+", beta.part))
  
  # generate model matrix
  mm <- model.matrix(formula.tmp, data = data)
  # change X colnames to ensure regression function does not
  # self-generate model matrix again
  all.betas.newnames <- str_replace(all.betas, ":", "_")
  all.betas.newnames <- str_remove(all.betas.newnames, "I\\(\\(")
  # all.betas.newnames <- str_remove(all.betas.newnames, "\\)/2\\)")
  all.betas.newnames <- str_remove(all.betas.newnames, "\\)\\)")
  all.betas.newnames <- str_replace(all.betas.newnames, "\\^", "_")
  colnames(mm) <- c(all.alphas, all.betas.newnames)
  
  # center and scale attributes
  x.center <- in.dat$center
  x.scale <- in.dat$scale
  dbh.scale <- sd(mandai$DBH_mi)
  
  # mannually scale X
  if (scale.X) {
    mm <- t(apply(mm, 1, function(x) (x - x.center[names(x)]) / x.scale[names(x)]))
    data$DBH_mi <- data$DBH_mi / dbh.scale
  } else {
    mm <- mm
    data$DBH_mi <- data$DBH_mi
  }
  
  # merge model matrix to original data.frame
  out <- cbind(data[, non.X.colnames], mm)
  return(out)
}
