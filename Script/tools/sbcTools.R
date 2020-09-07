util <- new.env()
source('tools/stanTools.r', local=util)
sbcFitFile <- function(save_progress, stanmodel, modelName, S) {
  file.path(save_progress, paste0(modelName, '-', S, '.rda'))
}

sbc_new <- function(stanmodel, modelName, data, N, L, n_eff_reltol=0.2, ..., save_progress, load_incomplete=FALSE) {
  doSave <- !missing(save_progress)
  # parameter names
  stan_code <- stanmodel$code()
  stan_code <- scan(what = character(), sep = "\n", quiet = TRUE, text = stan_code)
  pars_lines <- grep("[[:space:]]*(pars_)|(pars_\\[.*\\])[[:space:]]*=", stan_code, value = TRUE)
  pars_lines <- pars_lines[!grepl("^[[:space:]]*vector", pars_lines) & 
                             !grepl("^[[:space:]]*real", pars_lines)]   
  pars_names <- trimws(sapply(strsplit(pars_lines, split = "=", fixed = TRUE), tail, n = 1))
  pars_names <- unique(sub("^([a-z,A-Z,0-9,_]*)_.*;", "\\1", pars_names))
  noUnderscore <- grepl(";", pars_names, fixed=TRUE)
  
  if (!load_incomplete) {
    todo <- as.integer(seq(from = 1, to = .Machine$integer.max, length.out = N))
  } else {
    mn <- modelName
    runs <- dir(save_progress)
    runs <- runs[grepl(paste0("^", mn,'-(\\d+).rda$'), runs)]
    if (length(runs) == 0) {
      stop(paste("No completed runs found in", dir,
                 "matching regular expression", paste0("^", mn,'-(\\d+).rda$'),
                 "\nDid you use sbc(..., save_progress='/path/to/results')?"))
    }
    todo <- as.integer(sub(paste0(mn,'-(\\d+).rda'), "\\1", runs))
  }
  post = list() 
  times = list()
  for(n in 1:N) {     
    S <- seq(1:N)[n]
    #load if exists
    if (doSave) {
      file <- sbcFitFile(save_progress, stanmodel,modelName, S) #TODO data chage should be reflected in fitname
      if (file.exists(file)) {
        got <- try(load(file), silent=TRUE)
        time <- out$time()['total']
        times <- cbind(times, time)
        post[[n]] <- out
        next
      }
    }
    #iter_sampling = M = 3 * (20n - 1) # iter_warmup = 500, 597,
    out <- stanmodel$sample(data, chains = 1, iter_warmup = 500, iter_sampling = L, parallel_chains = 1, adapt_delta = 0.8,  save_warmup = FALSE, thin = 1) 
    time <- out$time()['total']
    times <- cbind(times, time)
    print(out$summary())  
    post[[n]] = out
    #save
    if (doSave) {
      save(out, file=file)
    }
  }
  avg_fit_time = mean(as.numeric(times))
  write.csv(avg_fit_time, file = file.path(delivDir, paste0(modelName, "_times.csv", sep = "")))
  
  # prior predictive distribution
  Y <- sapply(post, FUN = function(p) {
    summary <- p$summary()
    summary %>%
      filter(str_detect(variable, "y_$|y_\\[.*\\]")) %>%
      pull(mean) # not mean all same value but for extract purpose
  })
  
  # realizations of parameter
  pars <- sapply(post, FUN = function(p) {
    summary <- p$summary()
    summary %>%
      filter(str_detect(variable, "pars_\\[[[:digit:]]+\\]")) %>%
      pull(mean)
  })
  
  # ranks: unthinned binary values for draw > true
  ranks <- lapply(post, FUN = function(p) {
    r <- subset_draws(p$draws(), variable = "ranks_") %>%
      as_draws_matrix()
    if (is.null(dim(r))) {
      r <- as.matrix(r)
    }
    colnames(r) <- pars_names
    r[] <- r > 0
    return(r)
  })
  
  # divergences
  # sampler_params <- lapply(post, FUN = function(p){ p$sampler_diagnostics() %>%
  #     as_draws_matrix()
  # })
  out <- list(ranks = ranks, Y = Y, pars = pars, fits = post, time = avg_fit_time) #, sampler_params = sampler_params
  return(out)
}
ppc.sbc <- function(x, modelName, data, thin = 3, ...){
  # thinner <- seq(from = 1, to = dim(x$Y)[1], by = thin)
  # yhat <- x$Y[thinner]
  N = length(x$Y) / length(data$ye)
  yhat_df<- data.frame(matrix(x$Y, ncol = N))
  names(yhat_df) <- str_replace(names(data.frame(matrix(x$Y, ncol = N))), "X", "~y")
  yhat_df$idu <- as.numeric(row.names(yhat_df))
  yhat_df.melt <- reshape2::melt(yhat_df, id.vars="idu")
  y_df <- data.frame(idu = yhat_df$idu, value = data$ye)
  hp_mu <- toString(lapply(c(data$alpha_mu_prior,data$rho_mu_prior), FUN = function(p){round(p, 1)}))
  hp_sd <- toString(lapply(c(data$alpha_sd_prior,data$rho_sd_prior), FUN = function(p){round(p, 3)}))
  if(N< 5){
    plot<- ggplot(yhat_df.melt, aes(idu, value, color = variable),show.legend = FALSE) + 
    geom_point()+
    geom_line(data = y_df, aes(idu, value, colour="ye")) + 
    ylab(paste0("y_sim, ye ")) + xlab("covariate index") + 
    annotate("text",x = 10, y = 1000, label = paste0(round(x$time,1), "s for hp" , hp_mu," , ", hp_sd))
  }else{
    plot<- ggplot(yhat_df.melt, aes(idu, value, color = variable),show.legend = FALSE) + 
      geom_point()+
      theme(legend.position = "none")+
      geom_line(data = y_df, aes(idu, value, colour="ye")) + 
      ylab(paste0("y_sim, ye ")) + xlab("covariate index") + 
      annotate("text",x = 10, y = 800, label = paste0(round(x$time,1), "s for hp" , hp_mu," , ", hp_sd))
  }
  # plot <- ggplot(yhat_df.melt, aes(idu, value),show.legend = FALSE) + 
  #   geom_point() +  theme(legend.position = "none") +
  #   geom_line(data = y_df, aes(idu, value), colour="orange") + 
  #   ylab(paste0("y_sim, ye ")) + xlab("covariate index") + 
  #   annotate("text",x = 10, y = 1000, label = paste0(round(x$time,1), "s for hp" , hp_mu," , ", hp_sd))
  ggsave(file = file.path(delivDir, paste0(modelName, "_pp.png")), width = 5, height = 5)
  return(plot)
}

plot.sbc <- function(x, modelName, bins = 20, thin = 3, ...) {
  thinner <- seq(from = 1, to = nrow(x$ranks[[1]]), by = thin)
  u <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(u), each = nrow(u)))
  d <- data.frame(u = c(u), parameter)
  d <- d %>% filter(d$parameter != "y") ## remove y
  L = length(thinner)
  N = length(x$fits)
  CI = qbinom(c(0.005,0.5,0.995), size= N, prob = 1/(bins))
  print(CI)
  binwidth <- (length(thinner) + 1) / bins
  sbc_plot <- suppressWarnings(ggplot2::ggplot(d, aes(x = u)) +
                                 geom_histogram(binwidth = binwidth, color = "black", fill = "#ffffe8", boundary = 0) +
                                 facet_wrap(vars(parameter),  scale = "free") + theme(panel.spacing.x = unit(2, "lines"))) + 
                  geom_polygon(data=data.frame(x=c(-5,0,-5,L+5,L,L+5,-5),y= c(CI[1],CI[2],CI[3],CI[3],CI[2],CI[1],CI[1])),aes(x=x,y=y),fill="grey45",color="grey25",alpha=0.5)
  ggsave(file = file.path(delivDir, paste0(modelName, ".png")), width = 5, height = 5)
  return(sbc_plot)
}  

MW1 <- function(bin_count){
  bins <- length(bin_count)
  unif <- rep(1/bins, bins)
  M <- sum(bin_count)
  tempf <- Vectorize(function(i)  abs(bin_count[i]/M  - unif[i]))
  val <- integrate(tempf,1,bins, rel.tol=.Machine$double.eps^.05)$value
  return(val)
}
MKM <- function(bin_count){
  bins <- length(bin_count)
  diff <- abs(mean(bin_count) - bin_count)
  val <- diff[which.max(diff)] / mean(bin_count)
  return(val)
}
MChisq <- function(bin_count){
  return(chisq.test(bin_count)$p.value)
}

uniformity.sbc <- function(x, modelName, bins = 20, thin = 3){
  samples <- nrow(x$ranks[[1]])
  thinner <- seq(from = 1, to = samples, by = thin)
  ranks <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(ranks), each = nrow(ranks)))
  num_params <- ncol(ranks)
  M = samples / thin
  max_rank = M + 1
  bin_size <- max_rank / bins
  pval <- rep(NA, num_params)
  integ_msr <- rep(NA, num_params)
  max_diff <- rep(NA, num_params)
  bin_counts <- data.frame() #rep(NA, num_params)
  for (i in 1:num_params){
    bin_count <- rep(0, bins)
    integ <- rep(0, bins)
    for (m in 1: length(ranks[,1])) {
      bin <- ceiling(ranks[m,i] / bin_size)
      bin_count[bin] <- bin_count[bin] + 1
    }
    # sum of bin_count is N (= total fit number)
    print(bin_count)
    pval[i] <- MChisq(bin_count)
    integ_msr[i] <- MW1(bin_count)
    max_diff[i] <- MKM(bin_count)
    bin_counts <- rbind(bin_counts, bin_count)
  }
  colnames(bin_counts) <- seq(1:bins)
  write.csv(list("pval" = pval,"bin_counts"=  bin_counts, "integ" = integ_msr), file = file.path(delivDir, paste0(modelName, "_counts.csv", sep = "")))
  return(list(bin_counts = bin_counts, pval = pval, integ_msr = integ_msr, max_diff = max_diff))
}
