library(ggplot2);
library(knitr); 
library(tidyverse);
library(rstan);
library(tufte);
library(parallel);
library(cmdstanr);
library(posterior);
library("bayesplot");
library("ggplot2");
library("rstanarm");
set_cmdstan_path("/Users/hyunjimoon/Dropbox/20_paper/charles/code/cmdstan")
options(digits = 2);  options(htmltools.dir.version = FALSE)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

scriptDir <- getwd()
modelDir <- file.path(scriptDir, "models")
dataDir <- file.path(scriptDir, "data")
delivDir <- file.path(scriptDir, "deliv", modelName)
nChains <- 4
parallel_chains <- min(nChains, detectCores())

data <- read_rdump(file.path(dataDir, dataFile))
println <- function(msg) cat(msg); cat("\n")
printf <- function(pattern, ...) println(sprintf(pattern, ...))
util <- new.env()
source('stan_utility.R', local=util)
source('gp_utility.R', local=util)
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

c_light_teal="#6B8E8E"
c_mid_teal="#487575"
c_dark_teal="#1D4F4F"

c_green_trans <- c("#00FF0080")
c_superfine <- c("#8F272705")

sbcFitFile <- function(save_progress, stanmodel, modelName, S) {
  file.path(save_progress, paste0(modelName, '-', S, '.rda'))
}

sbc <- function(stanmodel, modelName, data, N, M, n_eff_reltol=0.2, ..., save_progress, load_incomplete=FALSE) {
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
  for(n in 1:N) {     
    S <- seq(1:N)[n]
    #load if exists
    if (doSave) {
      file <- sbcFitFile(save_progress, stanmodel,modelName, S) #TODO data chage should be reflected in fitname
      if (file.exists(file)) {
        got <- try(load(file), silent=TRUE)
        post[[n]] <- out
        next
      }
    }
    #iter_sampling = M = 3 * (20n - 1) # iter_warmup = 500, 597,
    out <- stanmodel$sample(data, chains = 1, iter_warmup = 500, iter_sampling = 597, parallel_chains = 1, adapt_delta = 0.8,  save_warmup = FALSE, thin = 1) 
    print(out$summary())  
    post[[n]] = out
    #save
    if (doSave) {
      save(out, file=file)
    }
  }
  
  # prior predictive distribution
  Y <- sapply(post, FUN = function(p) {
    summary <- p$summary()
    summary %>%
      filter(str_detect(variable, "y_$|y_\\[.*\\]")) %>%
      pull(mean)
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
  out <- list(ranks = ranks, Y = Y, pars = pars, fits = post) #, sampler_params = sampler_params
  return(out)
}
ppc.sbc <- function(x, thin = 3, ...){
  thinner <- seq(from = 1, to = dim(x$Y)[1], by = thin)
  yhat <- x$Y[thinner,]
  yhat_df<- data.frame(matrix(yhat, ncol = N))
  yhat_df$idu <- as.numeric(row.names(yhat_df))
  yhat_df.melt <- reshape2::melt(yhat_df, id.vars="idu")
  ggplot(yhat_df.melt, aes(idu, value, col=variable)) + 
    geom_line()
}

plot.sbc <- function(x, thin = 3, ...) {
  thinner <- seq(from = 1, to = nrow(x$ranks[[1]]), by = thin)
  u <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(u), each = nrow(u)))
  d <- data.frame(u = c(u), parameter)
  bins <- 20
  binwidth <- (length(thinner) + 1) / bins
  suppressWarnings(ggplot2::ggplot(d, aes(x = u)) +
                     geom_histogram(binwidth = binwidth, color = "black",
                                    fill = "#ffffe8", boundary = 0) +
                     facet_wrap(vars(parameter)) +
                     theme(panel.spacing.x = unit(2, "lines"))) +
    geom_hline(yintercept= nrow(u)/bins)
}

uniformity.sbc <- function(x, bins = 20, thin = 3){
  samples <- nrow(x$ranks[[1]])
  thinner <- seq(from = 1, to = samples, by = thin)
  ranks <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(ranks), each = nrow(ranks)))
  num_params <- length(colnames(ranks))
  M = samples / thin
  max_rank = M + 1
  bin_size <- max_rank / bins
  pval <- rep(NA, num_params)
  
  for (i in 1:num_params){
    bin_count <- rep(0, bins)
    for (m in 1: length(ranks[,1])) {
      bin <- ceiling(ranks[m,i] / bin_size)
      bin_count[bin] <- bin_count[bin] + 1
    }
    # sum of bin_count is N (= total fit number)
    pval[i] <- chisq.test(bin_count)$p.value
    print("bin_count")
    print(bin_count)
    print("pval")
    print(pval)
  }
  #list(bin_count, pval)
}

set_Data_Model_Dir <- function(dataFile, modelName){
  scriptDir <- getwd()
  modelDir <- file.path(scriptDir, "models")
  dataDir <- file.path(scriptDir, "data")
  outDir <- file.path(scriptDir, "deliv", modelName)
  delivDir <- file.path(scriptDir, "deliv", modelName)
  dir.create(delivDir)
  data <- read_rdump(file.path(dataDir, dataFile))
  print(data)
  file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
  mod <- cmdstan_model(file)
  return(list(data = data, mod = mod, delivDir = delivDir)) 
}
# list[data, mod, delivDir] <- set_Data_Model_Dir("ovarian_reduced.data.R", "bernoulli_logit_glm_ela_ACself")

data = set_Data_Model_Dir("ovarian_reduced.data.R", "bernoulli_logit_glm_ela_ACself")$data
mod = set_Data_Model_Dir("ovarian_reduced.data.R", "bernoulli_logit_glm_ela_ACself")$mod
delivDir = set_Data_Model_Dir("ovarian_reduced.data.R", "bernoulli_logit_glm_ela_ACself")$delivDir
submodelName <- "bernoulli_logit_glm_ela_ACself"
N = 3 #(bins = 20)*N # num of sim and fit
M = 597 # 897 (6n -3) # num of post. draw - related to iter_sampling - not used currently
sbc_res_bern_lik <- sbc(mod, submodelName, data = data, N = N, M = M, save_progress = delivDir)
ppc.sbc(sbc_res_bern_lik)
uniformity.sbc(sbc_res_bern_lik)
plot.sbc(sbc_res_bern_lik)