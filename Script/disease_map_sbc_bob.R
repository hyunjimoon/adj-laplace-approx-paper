library(ggplot2); library(knitr);  library(tidyverse); library(rstan); library(tufte); library(parallel); library(cmdstanr); library(posterior)
set_cmdstan_path("/Users/hyunjimoon/Dropbox/20_paper/charles/code/cmdstan")
source(file.path("tools", "cmdStanTools.r"))
source(file.path("tools", "stanTools.r"))
options(digits = 2);  options(htmltools.dir.version = FALSE)

modelName <- "disease_map_ela_sbc_bob"
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
data$rho_alpha_prior <- 2.42393
data$rho_beta_prior <- 14.8171

file <- file.path(modelDir, modelName, paste0(modelName, ".stan"))
mod <- cmdstan_model(file, quiet = FALSE)


println <- function(msg) { cat(msg); cat("\n") }
printf <- function(pattern, ...) println(sprintf(pattern, ...))
print_file <- function(file) cat(paste(readLines(file), "\n", sep=""), sep="")
knitr::opts_chunk$set(
  include = TRUE,  cache = FALSE,  collapse = TRUE,  echo = TRUE,
  message = FALSE, tidy = FALSE,  warning = FALSE,   comment = "  ",
  dev = "png", dev.args = list(bg = '#FFFFF8'), dpi = 300,
  fig.align = "center",  fig.width = 7,  fig.asp = 0.618,  fig.show = "hold",
  out.width = "90%")

ggtheme_tufte <- function() {
  theme(plot.background =
          element_rect(fill = "#fffff8",
                       colour = "#fffff8",
                       size = 0.5,
                       linetype = "solid"),
        plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),
        panel.background =
          element_rect(fill = "#fffff8",
                       colour = "#fffff8",
                       size = 0.5,
                       linetype = "solid"),
        panel.grid.major = element_line(colour = "white",
                                        size = 1, linetype="dashed"),
        panel.grid.minor = element_blank(),
        legend.box.background =
          element_rect(fill = "#fffff8",
                       colour = "#fffff8",
                       linetype = "solid"),
        axis.ticks = element_blank(),
        axis.text = element_text(family = "Palatino", size = 14),
        axis.title.x = element_text(family = "Palatino", size = 16,
                                    margin = margin(t = 15,
                                                    r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family = "Palatino", size = 16,
                                    margin = margin(t = 0,
                                                    r = 15, b = 0, l = 0)),
        strip.background = element_rect(fill = "#fffff8",
                                        colour = "#fffff8",
                                        linetype = "solid"),
        strip.text = element_text(family = "Palatino", size = 14),
        legend.text = element_text(family = "Palatino", size = 14),
        legend.title = element_text(family = "Palatino", size = 16,
                                    margin = margin(b = 5)),
        legend.background = element_rect(fill = "#fffff8",
                                         colour = "#fffff8",
                                         linetype = "solid"),
        legend.key = element_rect(fill = "#fffff8",
                                  colour = "#fffff8",
                                  linetype = "solid")
  )
}

## helper functions - sbc

# @param y:  sequence of ranks in 1:max_rank
# @param max_rank: maximum rank of data in y
# @param bins (default 20):  bins to use for chi-square test
# @error return NA if max rank not divisible by number of bins
# @return p-value for chi-square test that data is evenly
#         distributed among the bins
test_uniform_ranks <- function(y, max_rank, bins = 4) { #20
  if (max_rank / bins != floor(max_rank / bins)) {
    printf("ERROR in test_uniform_ranks")
    printf("  max rank must be divisible by bins.")
    printf("  found max rank = %d;  bins = %d", max_rank, bins)
    return(NA)
  }
  bin_size <- max_rank / bins
  bin_count <- rep(0, bins)
  N <- length(y)
  for (n in 1:N) {
    bin <- ceiling(y[n] / bin_size)
    bin_count[bin] <- bin_count[bin] + 1
  }
  chisq.test(bin_count)$p.value
}



# Run one iteration of the sampler, then parse the output column names
# to determine number of monitored parameters.
# @param model: cmdstan_model object
# @param data: data for model
# @return number of entries in array I_lt_sim
num_monitored_params <- function(model, data) {
  fit <- model$sample(data=data, chains=1, iter_warmup=20, iter_sampling=1, refresh=0)
  metadata <- fit$metadata()
  sum(str_detect(metadata$model_params,"I_lt_sim"))
}

# @param model: cmdstan_model object
# @param data: data for model
# @param sbc_sims:  number of total simulation runs for SBC
# @param stan_sims: number of posterior draws per Stan simulation
# @param init_thin: initial thinning (doubles thereafter up to max)
# @param max_thin: max thinning level
# @param seed: PRNG seed to use for Stan program to generate data
# @param target_n_eff: target effective sample size (should be 80%
#                      or 90% of stan_sims to stand a chance)
# @return list with keys (rank, p_value, thin) for 2D array of ranks
#         and 1D array of p-values, and 1D array of thinning rates
sbc <- function(model, data,
                sbc_sims = 1000, stan_sims = 999,
                init_thin = 4, max_thin = 64,
                target_n_eff = 0.8 * stan_sims) {
  num_params <- num_monitored_params(model, data)
  ranks <- matrix(nrow = sbc_sims, ncol = num_params)
  thins <- rep(NA, sbc_sims)
  for (n in 1:sbc_sims) {
    n_eff <- 0
    thin <- init_thin
    while (TRUE) {
      fit <- model$sample(
        data = data,
        chains = 1,
        iter_sampling = thin * stan_sims,
        thin = thin,
        adapt_delta = 0.99,
        refresh = 0)
      #stanfit <- read_stan_csv(fit$output_files()) 
      #n_eff <- parameterTable(stanfit, c("lp__"))$n_eff
      n_eff <- fit$summary("lp__") %>% pull("ess_bulk")
      print(fit$summary("lp__"))

      if (n_eff >= target_n_eff || (2 * thin) > max_thin) break; 
      thin <- 2 * thin
    }
    thins[n] <- thin
    for (i in 1:num_params) {
      col = paste(c("I_lt_sim[",i,"]"), collapse="")
      ranks[n, i] <- sum(fit$draws(col)) + 1
    }
  }
  pval <- rep(NA, num_params)
  for (i in 1:num_params)
    pval[i] <- test_uniform_ranks(ranks[ , i],
                                  max_rank = stan_sims + 1)
  list(rank = ranks, p_value = pval, thin = thins)
}

# working example:
# model <- cmdstan_model("models3/normal-sbc.stan")
# result <- sbc(model, data = NULL, sbc_sims = 5, stan_sims = 7,
#               max_thin = 32)

result <- sbc(mod, data = data, sbc_sims = 3, stan_sims = 7,
              max_thin = 64)
table(result$thin)
result$p_value
A <- dim(result$rank)[1]
rank_df <-
  rbind(data.frame(parameter = rep("mu", A),
                   y = result$rank[ , 1]),
        data.frame(parameter = rep("sigma", A),
                   y = result$rank[, 2]))
rank_plot <-
  ggplot(rank_df, aes(x = y)) +
  geom_histogram(binwidth = 1, color = "black", #binwidth = 50 when 1000 sbc_sims
                 fill = "#ffffe8", boundary = 0) +
  facet_wrap(vars(parameter)) +
  ggtheme_tufte() +
  theme(panel.spacing.x = unit(2, "lines"))
rank_plot


