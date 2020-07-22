library(ggplot2); library(knitr);  library(tidyverse); library(rstan); library(tufte); library(parallel); library(cmdstanr); library(posterior)
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
n_sample <- 100 # 911
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
data$rho_alpha_prior <- 2.42393
data$rho_beta_prior <- 14.8171

file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
mod <- cmdstan_model(file, quiet = FALSE)
out <- mod$sample(data, chains = 1,iter_warmup = 100, iter_sampling = 100, parallel_chains = 1, save_warmup = FALSE, thin = 1)
fit <- out$summary()

sbc <- function(stanmodel, modelName, data, M, ...) {
  # parameter names
  stan_code <- stanmodel$code()
  stan_code <- scan(what = character(), sep = "\n", quiet = TRUE, text = stan_code)
  pars_lines <- grep("[[:space:]]*(pars_)|(pars_\\[.*\\])[[:space:]]*=", stan_code, value = TRUE)
  pars_lines <- pars_lines[!grepl("^[[:space:]]*vector", pars_lines) & 
                             !grepl("^[[:space:]]*real", pars_lines)]   
  pars_names <- trimws(sapply(strsplit(pars_lines, split = "=", fixed = TRUE), tail, n = 1))
  pars_names <- unique(sub("^([a-z,A-Z,0-9,_]*)_.*;", "\\1", pars_names))
  noUnderscore <- grepl(";", pars_names, fixed=TRUE)
  
  ranks <- rep(NA, M)
  sampler_params <- rep(NA, M)
  m <- 1
  while (m <= M){
    S <- seq(from = 1, to = .Machine$integer.max, length.out = M)[m]
    out <- stanmodel$sample(data, chains = 1,iter_warmup = 100, iter_sampling = 100, seed = floor(S), parallel_chains = 1, save_warmup = FALSE, thin = 1)
    if (class(out)[1] != "CmdStanMCMC") next
    m = m + 1
    summary <- out$summary()
    Y <- summary %>% #as_draws_matrix()
      filter(str_detect(variable, "y_$|y_\\[.*\\]")) %>%
      pull(mean)
    
    means <-  summary %>%
      filter(str_detect(variable, "y_$|y_\\[.*\\]")) %>%
      pull(mean)
    mark <- summary %>%
      filter(str_detect(variable, "pars_\\[[[:digit:]]+\\]")) %>%
      pull(mean)
    pars <- means[mark]
    
    sampler_params <- out$sampler_diagnostics() %>%
      as_draws_matrix()
    
    r <- subset_draws(out$draws(), "ranks_") %>%
      as_draws_matrix()
    if (is.null(dim(r))) {
      r <- as.matrix(r)
    }
    colnames(r) <- pars_names
    r[] <- r > 0
    ranks[m] <- r
    #post[m] <- list( Y = Y, pars = pars, sampler_params = sampler_params, ranks = r)
  }
  return(ranks)
}
plot.sbc <- function(ranks, thin = 3, ...) {
  thinner <- seq(from = 1, to = nrow(ranks[[1]]), by = thin)
  u <- t(sapply(ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(u), each = nrow(u)))
  d <- data.frame(u = c(u), parameter)
  suppressWarnings(ggplot2::ggplot(d) + 
                     ggplot2::geom_freqpoly(ggplot2::aes(x = u), ...) + 
                     ggplot2::facet_wrap("parameter"))
}

# plot.sbc <- function(x, thin = 3, ...) {
#   thinner <- seq(from = 1, to = nrow(x$ranks[[1]]), by = thin)
#   u <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
#   parameter <- as.factor(rep(colnames(u), each = nrow(u)))
#   d <- data.frame(u = c(u), parameter)
#   suppressWarnings(ggplot2::ggplot(d) + 
#                      ggplot2::geom_freqpoly(ggplot2::aes(x = u), ...) + 
#                      ggplot2::facet_wrap("parameter"))
# }

print.sbc <- function(sampler_params, ...) {
  divergences <- apply(sampler_params, MARGIN = 3, FUN = function(y) sum(y[,"divergent__"]))
  bad <- sum(divergences > 0L)
  cat(paste(bad, "chains had divergent transitions after warmup\n"))
  if (bad > 0L) cat(paste("there were a total of", sum(divergences),
                          "divergent transitions across all chains\n"))
  return(invisible(NULL))
}

# file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
# mod <- cmdstan_model(file, quiet = FALSE)
sbc_res <- sbc(mod, modelName, data, M = 5) #, save_progress = outDir, refresh = 0
# pdf(file = file.path(delivDir, paste(modelName,"SBCPlots%03d.pdf", sep = "")), width = 6, height = 6, onefile = F)
plot(sbc_res)
# dev.off()
#print(sbc_res)
