library(ggplot2);
library(knitr); 
library(tidyverse);
library(rstan);
library(tufte);
library(parallel);
library(cmdstanr);
library(posterior)
set_cmdstan_path("/Users/hyunjimoon/Dropbox/20_paper/charles/code/cmdstan")
source(file.path("tools", "cmdStanTools.r"))
source(file.path("tools", "stanTools.r"))
options(digits = 2);  options(htmltools.dir.version = FALSE)

modelName <- "disease_map_ela_sbc"
dataFile <- "disease_data_100.r"
scriptDir <- getwd()

modelDir <- file.path(scriptDir, "models")
dataDir <- file.path(scriptDir, "data")
outDir <- file.path(scriptDir, "deliv", modelName)
delivDir <- file.path("deliv", modelName)
data <- read_rdump(file.path(dataDir, dataFile))

reduced <- FALSE
n_sample <- 100 
nChains <- 4
num_cores <- min(nChains, detectCores())

if (reduced) {
  bin <- floor(100 / n_sample)
  select_index <- rep(NA, n_sample)
  for (i in 1:n_sample) select_index[i] <- sample(((i - 1) * bin + 1):(i * bin), 1)
  data <- data[c("x", "n_obs", "n_covariates", "ye")]
  data$x <- data$x[select_index, ]
  data$ye <- data$ye[select_index]
  data$n_obs <- n_sample
} else{
  data <- data[c("x", "n_obs", "n_covariates", "ye")]
}

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
    todo <- as.integer(seq(from = 0, to = .Machine$integer.max, length.out = N))
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
    S <- seq(from = 1, to = .Machine$integer.max, length.out = N)[n]
    #load if exists
    if (doSave) {
      file <- sbcFitFile(save_progress, stanmodel,modelName, S) #TODO data chage should be reflected in fitname
      if (file.exists(file)) {
        got <- try(load(file), silent=TRUE)
        post[[n]] <- out
        next
      }
    }
    #iter_sampling = M = 3 * (20n - 1)
    out <- stanmodel$sample(data, chains = 1, iter_warmup = 500, iter_sampling = 597, parallel_chains = 1, save_warmup = FALSE, thin = 1) #seed = floor(S)
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
  out <- list(ranks = ranks, Y = Y, pars = pars) #, sampler_params = sampler_params
  return(out)
}

plot.sbc <- function(x, thin = 3, ...) {
  thinner <- seq(from = 1, to = nrow(x$ranks[[1]]), by = thin)
  u <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(u), each = nrow(u)))
  d <- data.frame(u = c(u), parameter)
  suppressWarnings(ggplot2::ggplot(d) +
                     ggplot2::geom_freqpoly(ggplot2::aes(x = u)) +
                     ggplot2::facet_wrap("parameter"))
}

uniformity.sbc <- function(x, bins = 20, thin = 3){
  samples <- nrow(x$ranks[[1]])
  thinner <- seq(from = 1, to = samples, by = thin)
  ranks <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(ranks), each = nrow(ranks)))
  num_params <- length(unique(parameter))
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
  print("uniformity pvalue")
  print(pval)
  }
}

data$rho_alpha_prior <- 2 #2.42393
data$rho_beta_prior <- 20 #14.8171
modelName <- "disease_map_ela_sbc_2_20_3rd"
file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
mod <- cmdstan_model(file, quiet = FALSE)

# modelName <- "bern_param_1_100_2"
# mod = cmdstan_model(stan_file = "models/sbc_test.stan")
# data = list(N = 10, a = 1, b = 100)
N = 10 #999 # num of sim and fit
M = 117 # 897 (6n -3) # num of post. draw - related to iter_sampling - not used currently

sbc_res <- sbc(mod, modelName, data = data, N = N, M = M, save_progress = outDir)
plot.sbc(sbc_res)
uniformity.sbc(sbc_res, thin = 3)

print.sbc <- function(x, ...) {
  divergences <- apply(x$sampler_params, MARGIN = 3, FUN = function(y) sum(y[,"divergent__"]))
  bad <- sum(divergences > 0L)
  cat(paste(bad, "chains had divergent transitions after warmup\n"))
  if (bad > 0L) cat(paste("there were a total of", sum(divergences), 
                          "divergent transitions across all chains\n"))
  if (length(x$pareto_k)) {
    cut_pareto_k <- cut(c(x$pareto_k), breaks = c(-Inf, 0.5, 0.7, 1, Inf))
    cat("Aggregate Pareto k estimates:\n")
    print(prop.table(table(cut_pareto_k)))
  }
  return(invisible(NULL))
}