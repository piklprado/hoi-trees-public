# Function to numerically simulate diameter over time

simulateD_focal_garden_sep_mod <- 
  function(model.n, model.a, model.b, t.step = 100, t.max = 2, focals) {
    # time steps
    t.step <- t.step
    t.max <- t.max
    dt <- t.max / t.step
    
    # focal species or individuals
    focals <- focals
    focals.no <- table(focals)
    
    # initial DBH
    dbh.scale <- sd(mandai$DBH_mi)
    D0 <- 1
    
    # neighbour attributes
    Pred.type <- c("Null", "Direct", "HOI")
    
    # center and scale attributes
    x.center <- in.dat$center
    x.scale <- in.dat$scale
    
    # setup neighbourhood X's
    focal.D0 <- rep(pi * (D0/2)^2, length(focals))
    names(focal.D0) <- focals
    
    alpha0 <- tapply(focal.D0, names(focal.D0), sum)
    absentee.names <- setdiff(names(x.center)[1:10], names(alpha0))
    absentee <- numeric(length(absentee.names))
    names(absentee) <- absentee.names
    alpha0 <- c(alpha0, absentee)
    alpha0 <- alpha0[order(names(alpha0))]
    
    beta_intra0 <- alpha0 ^2
    names(beta_intra0) <- paste0(names(alpha0), "_2")
    
    beta_inter_names <- setdiff(names(x.center), c(names(alpha0), names(beta_intra0)))
    beta_inter0 <- alpha0[word(beta_inter_names, 1, sep = "_")] * alpha0[word(beta_inter_names, 2, sep = "_")]
    names(beta_inter0) <- beta_inter_names
    
    neigh0 <- c(alpha0, beta_intra0, beta_inter0)
    neigh0 <- t(replicate(length(Pred.type), neigh0))
    neigh0 <- 
      data.frame(Pred.type = Pred.type,
                 neigh0)
    # scale and center
    neigh0[, names(x.center)] <- t((t(neigh0[, names(x.center)]) - x.center) / x.scale)
    # bind neighbours with focals
    neigh <- 
      bind_cols(
        Focal = rep(focals, each = length(Pred.type)),
        do.call(rbind, replicate(length(focals), neigh0, simplify = FALSE))
      )
    
    # set up diameter matrix for all focals
    D <- matrix(NA, nrow = length(focals)*length(Pred.type), ncol = t.step+1)
    D[, 1] <- D0
    
    neigh_t <- 
      neigh %>% 
      mutate(DBH_mi = D[, 1])
    
    for (t in 1:t.step) {
      pb <- txtProgressBar(min = 1, max = t.step, style = 3)
      
      # bind neighbour with DBH_mi and then scale it before prediction
      newdat.tmp <- 
        neigh_t %>% 
        mutate(DBH_mi = DBH_mi / dbh.scale) %>% 
        # create unique individual ID for left-joining predictions later
        bind_cols(UID = rep(1:length(focals), each = length(Pred.type)))
      newdat.tmp.n <- newdat.tmp %>% filter(Pred.type == "Null")
      newdat.tmp.a <- newdat.tmp %>% filter(Pred.type == "Direct")
      newdat.tmp.b <- newdat.tmp %>% filter(Pred.type == "HOI")
      re_form <- NULL
      
      # prediction
      pred.tmp.n <-
        bind_cols(newdat.tmp.n %>% select(UID, Focal, Pred.type),
                  Estimate = 
                    as.numeric(posterior_predict(
                      model.n,
                      newdata = newdat.tmp.n,
                      re_formula = re_form,
                      allow_new_levels = TRUE,
                      nsamples = 1,
                      cores = 4
                    ))) %>% 
        # backtransform growth
        mutate(Estimate = inv_modulus_trans(Estimate, 0.55))
      pred.tmp.a <-
        bind_cols(newdat.tmp.a %>% select(UID, Focal, Pred.type),
                  Estimate = 
                    as.numeric(posterior_predict(
                      model.a,
                      newdata = newdat.tmp.a,
                      re_formula = re_form,
                      allow_new_levels = TRUE,
                      nsamples = 1,
                      cores = 4
                    ))) %>% 
        # backtransform growth
        mutate(Estimate = inv_modulus_trans(Estimate, 0.55))
      pred.tmp.b <-
        bind_cols(newdat.tmp.b %>% select(UID, Focal, Pred.type),
                  Estimate = 
                    as.numeric(posterior_predict(
                      model.b,
                      newdata = newdat.tmp.b,
                      re_formula = re_form,
                      allow_new_levels = TRUE,
                      nsamples = 1,
                      cores = 4
                    ))) %>% 
        # backtransform growth
        mutate(Estimate = inv_modulus_trans(Estimate, 0.55)) 
      pred.tmp <- bind_rows(pred.tmp.n, pred.tmp.a, pred.tmp.b)
      newdat.tmp <- left_join(newdat.tmp, pred.tmp,
                              by = c("Focal", "Pred.type", "UID"))
      D[, t+1] <- D[, t] + (pred.tmp[, "Estimate"] * dt)

      # Update neighbourhood BA
      alpha_t <- matrix(D[, t+1], 
                        nrow = length(Pred.type), 
                        ncol = length(focals), 
                        byrow = FALSE,
                        dimnames = list(Pred.type, focals))
      alpha_t <- pi * (alpha_t / 2)^2  # convert diameters to basal areas
      alpha_t <- t(apply(alpha_t, 1, function(x) tapply(x, names(x), sum)))
      # add missing species
      alpha_t_absentee <- matrix(0, nrow = length(Pred.type), ncol = length(absentee))
      colnames(alpha_t_absentee) <- absentee.names
      alpha_t <- cbind(alpha_t, alpha_t_absentee)
      alpha_t <- alpha_t[, order(colnames(alpha_t))]
      
      # remove focal BA from neighbour sum BA
      alpha_t_df <- 
        bind_cols(
          newdat.tmp[, 1:2],
          DBH_mi = D[, t+1],
          as.data.frame(do.call(rbind, replicate(length(focals), alpha_t, simplify = FALSE)))
        ) %>% 
        pivot_longer(cols = c(-Focal, -Pred.type, -DBH_mi),
                     names_to = "Neigh",
                     values_to = "BA_j") %>% 
        # remove focal BA from neighbour sum BA
        mutate(BA_j = ifelse(Focal == Neigh, BA_j - pi*(DBH_mi/2)^2, BA_j)) %>% 
        pivot_wider(names_from = Neigh,
                    values_from = BA_j)
      
      beta_intra_t <- t(apply(alpha_t_df[, -c(1:3)], 1, function(x) x^2))
      colnames(beta_intra_t) <- paste0(colnames(alpha_t), "_2")
      
      beta_inter_t <- 
        t(apply(alpha_t_df[, -c(1:3)], 1, function(x) {
          x[word(beta_inter_names, 1, sep = "_")] * x[word(beta_inter_names, 2, sep = "_")]
        }))
      colnames(beta_inter_t) <- beta_inter_names
      
      # combine all predictor terms
      neigh_t <- cbind(alpha_t_df, beta_intra_t, beta_inter_t)
      # scale and center
      neigh_t[, names(x.center)] <- t((t(neigh_t[, names(x.center)]) - x.center) / x.scale)

      setTxtProgressBar(pb, t)
    }
    close(pb)
    
    # tidy output
    out <- 
      cbind(neigh[, c("Focal", "Pred.type")], D) %>% 
      pivot_longer(cols = -c("Focal", "Pred.type"),
                   names_to = "t",
                   values_to = "D") %>% 
      mutate(t = ((as.numeric(t) - 1) / t.step) * t.max,
             Pred.type = 
               factor(Pred.type, levels = c("Null", "Direct", "HOI"))) %>% 
      # add individual ID
      bind_cols(Ind = rep(seq_len(length(focals)), each = ncol(D)*length(Pred.type)))
    
    return(out)
    
  }
