# predict cummulative effect of alpha and beta on growth across sizes

growth.curve.cumm.pred <- function(model, dat, DBH.range = c("species", "same")) {
  
  focals <- as.character(sort(unique(mandai$Focal)))
  
  # neighbour attributes
  pred.by <- 0.1
  if (DBH.range == "species") {
    dbh_mi.range <- with(dat, tapply(DBH_mi, Focal, range))
  } else if (DBH.range == "same") {
    dbh_mi.range <- lapply(focals, function(x) c(1, 40))
    names(dbh_mi.range) <- focals
  }
  
  dbh_mi.seq   <- lapply(dbh_mi.range, function(x) seq(x[1], x[2], by = pred.by))
  Pred.type    <- c("Max",
                    "D only",
                    "D + HOI_intra", 
                    "D + HOI_intra + HOI_inter")
  
  # center and scale attributes
  x.center <- in.dat$center
  x.scale <- in.dat$scale
  dbh.scale <- sd(mandai$DBH_mi)
  
  # predict at mean neighbour BA
  xv <- x.center[neigh]
  xv["OTHERS"] <- 0
  xv_2 <- xv ^2
  names(xv_2) <- paste0(names(xv), "_2")
  
  xv_xv_names <- setdiff(names(x.center), c(names(xv), names(xv_2)))
  xv_xv <- xv[word(xv_xv_names, 1, sep = "_")] * xv[word(xv_xv_names, 2, sep = "_")]
  names(xv_xv) <- xv_xv_names
  
  xv_all <- c(xv, xv_2, xv_xv)
  xv_all <- xv_all[names(x.center)] # make sure order is identical
  
  # create new predictor design matrix
  Xv1 <- Xv2 <- Xv3 <- Xv4 <- xv_all
  Xv1 <- Xv1 * 0
  Xv2[setdiff(names(Xv2), neigh)] <- 0
  Xv3[setdiff(names(Xv3), c(neigh, paste0(neigh, "_2")))] <- 0

  # expand new data
  newdat.org <- list()
  for (i in seq_len(length(dbh_mi.seq))) {
    newdat.org[[i]] <- 
      data.frame(
        Focal = names(dbh_mi.seq)[i],
        Plot = NA,
        Year = NA,
        DBH_mi = rep(dbh_mi.seq[[i]], length(Pred.type))
      ) %>% 
      bind_cols(
        bind_rows(
          data.frame(t(replicate(length(dbh_mi.seq[[i]]), Xv1))),
          data.frame(t(replicate(length(dbh_mi.seq[[i]]), Xv2))),
          data.frame(t(replicate(length(dbh_mi.seq[[i]]), Xv3))),
          data.frame(t(replicate(length(dbh_mi.seq[[i]]), Xv4)))
        ) %>% 
          mutate(Pred.type = 
                   rep(Pred.type, each = length(dbh_mi.seq[[i]])))
      )
  }
  newdat.org <- do.call(rbind, newdat.org)
  
  # scale new data
  newdat.scale <- newdat.org
  newdat.scale[, names(x.center)] <- 
    t(apply(newdat.scale[, names(x.center)], 1, function(x) (x - x.center[names(x)]) / x.scale[names(x)]))
  newdat.scale$DBH_mi <- newdat.scale$DBH_mi / dbh.scale 
  re_form <- NULL
  
  # Predict
  out <- t(apply(
    inv_modulus_trans(
      # fitted(
      predict(
        model, 
        newdata = newdat.scale,
        re_formula = re_form,
        summary = FALSE, 
        allow_new_levels = TRUE, 
        sample_new_levels = "gaussian",
        cores = 4
        ),
      0.55), 
    2, quantile, probs = c(0.05, 0.5, 0.95)))
  colnames(out) <- c("Q5", "Estimate", "Q95")
  out <- cbind(newdat.org, out)  
  out$Pred.type <- factor(out$Pred.type, 
                          levels = c("Max",
                                     "D only",
                                     "D + HOI_intra", 
                                     "D + HOI_intra + HOI_inter"))
  
  return(out)
}
