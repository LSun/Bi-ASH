#######
## Simulations corresponding to the Statistica Sinica revisions requiring
## misspecification of negative controls.
#######

one_rep <- function(new_params, current_params) {
  source("~/GitHub/Bi-ASH/code/biashr.R")
  source("~/GitHub/Bi-ASH/code/nc_adjustment_methods.R")
  args_val <- append(current_params, new_params)
  set.seed(new_params$current_seed)

  d_out <- seqgendiff::poisthin(mat = args_val$mat,
                                nsamp = args_val$Nsamp,
                                ngene = args_val$Ngene + args_val$ncontrols,
                                skip_gene = args_val$skip_gene,
                                signal_params = list(mean = 0, sd = args_val$log2foldsd),
                                gselect = "random",
                                prop_null = (args_val$nullpi * args_val$Ngene + args_val$ncontrols) / (args_val$Ngene + args_val$ncontrols)
                                )

  which_null <- abs(d_out$beta) < 10 ^ -6
  nnull         <- sum(which_null)
  control_genes <- which_null
  control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrols)] <- FALSE

  stopifnot(sum(control_genes) == args_val$ncontrols)
  stopifnot(sum(control_genes & !which_null) == round((1 - args_val$prop_cont) * args_val$ncontrols))

  beta_true <- d_out$beta

  X <- d_out$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(d_out$Y + 1)

  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

  start.time <- proc.time()

  ## KEY:
  ## o = original
  ## m = MAD inflation
  ## c = control-gene calibrated
  ## l = limma-adjusted
  ## lb = limma-adjusted Before gls
  ## la = limma-adjusted After gls
  ## d = delta-adjustment from CATE package

  method_list        <- list()
  
  ## OLS ---------------------------------------------------------------------
  method_list$ols  <- ols(Y = Y, X = X)
  # method_list$ols_m  <- mad_adjust(method_list$ols_o)
  # method_list$ols_c  <- ctl_adjust(method_list$ols_o, control_genes = control_genes)
  # method_list$ols_l  <- limma_adjust(method_list$ols_o)
  # method_list$ols_lm <- mad_adjust(method_list$ols_l)
  # method_list$ols_lc <- ctl_adjust(method_list$ols_l, control_genes = control_genes)

  ## limma ---------------------------------------------------------------------  
  method_list$limma <- limma_simp(Y = d_out$Y, X = X)
  
  ## RUV2 --------------------------------------------------------------------
  method_list$ruv2 <- limma_adjust(ruv2_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes))
  # method_list$ruv2_o  <- ruv2_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv2_m  <- mad_adjust(method_list$ruv2_o)
  # method_list$ruv2_c  <- ctl_adjust(method_list$ruv2_o, control_genes = control_genes)
  # method_list$ruv2_l  <- limma_adjust(method_list$ruv2_o)
  # method_list$ruv2_lm <- mad_adjust(method_list$ruv2_l)
  # method_list$ruv2_lc <- ctl_adjust(method_list$ruv2_l, control_genes = control_genes)

  ## RUV3 --------------------------------------------------------------------
  method_list$ruv3 <- ruv3_limma_post(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_o  <- ruv3_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_m  <- mad_adjust(method_list$ruv3_o)
  # method_list$ruv3_c  <- ruv3_ctl_adjust(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_lb <- ruv3_limma_pre(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_lbm <- mad_adjust(method_list$ruv3_lb)
  # method_list$ruv3_lbc <- ruv3_limma_pre_adjust(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_la  <- ruv3_limma_post(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv3_lam <- mad_adjust(method_list$ruv3_la)
  # method_list$ruv3_lac <- ruv3_limma_post_adjust(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)

  ## RUV4 (not CATE) ---------------------------------------------------------
  # method_list$ruv4_o  <- ruv4_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$ruv4_m  <- mad_adjust(method_list$ruv4_o)
  # method_list$ruv4_c  <- ctl_adjust(method_list$ruv4_o, control_genes = control_genes)
  # method_list$ruv4_l  <- limma_adjust(method_list$ruv4_o)
  # method_list$ruv4_lm <- mad_adjust(method_list$ruv4_l)
  # method_list$ruv4_lc <- ctl_adjust(method_list$ruv4_l, control_genes = control_genes)

  ## CATE -------------------------------------------------------------------
  method_list$cate  <- limma_adjust(cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes))
  # method_list$cate_o   <- cate_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$cate_m   <- mad_adjust(method_list$cate_o)
  # method_list$cate_c   <- ctl_adjust(method_list$cate_o, control_genes = control_genes)
  # method_list$cate_lb  <- cate_limma(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$cate_lbm <- mad_adjust(method_list$cate_lb)
  # method_list$cate_lbc <- ctl_adjust(method_list$cate_lb, control_genes = control_genes)
  # method_list$cate_la  <- limma_adjust(method_list$cate_o)
  # method_list$cate_lam <- mad_adjust(method_list$cate_la)
  # method_list$cate_lac <- ctl_adjust(method_list$cate_la, control_genes = control_genes)
  # method_list$cate_d   <- cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  # method_list$cate_dm  <- mad_adjust(method_list$cate_d)
  # method_list$cate_dc  <- ctl_adjust(method_list$cate_d, control_genes = control_genes)
  # method_list$cate_dl  <- limma_adjust(method_list$cate_d)
  # method_list$cate_dlm <- mad_adjust(method_list$cate_dl)
  # method_list$cate_dlc <- ctl_adjust(method_list$cate_dl, control_genes = control_genes)

  ## t-approximation version of RUVB -------------------------------------
  # method_list$ruvbn    <- ruvb_bfa_gs_linked_se(Y = Y, X = X, num_sv = num_sv,
  #                                               control_genes = control_genes)
  # method_list$ruvbnl   <- limma_adjust(method_list$ruvbn)

  ## Normal approximation version of RUVB --------------------------------
  # method_list$ruvbnn   <- method_list$ruvbn
  # method_list$ruvbnn$df <- Inf
  
  ## ashr  -------------------------------------------------------------------
  method_list$ashr <- ashr_nc(Y = d_out$Y, X = X, control_genes = control_genes)
  # method_list$ashro <- ashr_o(Y = Y, X = X, control_genes = control_genes)

  ## biashr -------------------------------------------------------------------
  method_list$biashr <- biashr_nc(Y = d_out$Y, X = X, control_genes = control_genes)
  # method_list$biashro <- biashr_o(Y = Y, X = X, control_genes = control_genes)
  
  ## Get CI's and p-values ----------------------------------------------------
  pci_list <- lapply(X = method_list, FUN = calc_ci_p)
  pci_list$ashr <- list(betahat = method_list$ashr$betahat, pvalues = method_list$ashr$lfsr, lower = pci_list$ashr$lower, upper = pci_list$ashr$upper)
  pci_list$biashr <- list(betahat = method_list$biashr$betahat, pvalues = method_list$biashr$lfsr, lower = pci_list$biashr$lower, upper = pci_list$biashr$upper)
  # pci_list$ashro <- list(betahat = method_list$ashro$betahat, pvalues = method_list$ashro$lfsr, lower = pci_list$ashro$lower, upper = pci_list$ashro$upper)
  # pci_list$biashro <- list(betahat = method_list$biashro$betahat, pvalues = method_list$biashro$lfsr, lower = pci_list$biashro$lower, upper = pci_list$biashro$upper)
  
  ## Get summary quantities --------------------------------------------------
  get_mse <- function(args, beta_true, control_genes) {
    mean((args$betahat[!control_genes] - beta_true[!control_genes]) ^ 2)
  }

  get_auc <- function(args, which_null, control_genes) {
    if (sum(which_null) == length(which_null)) {
      return(NA)
    }
    pROC::roc(predictor = args$pvalues[!control_genes],
              response = which_null[!control_genes])$auc
  }

  get_coverage <- function(args, beta_true, control_genes) {
    mean(args$lower[!control_genes] < beta_true[!control_genes] &
           args$upper[!control_genes] > beta_true[!control_genes])
  }

  mse_vec <- sapply(pci_list, get_mse, beta_true = beta_true,
                    control_genes = control_genes)
  auc_vec <- sapply(pci_list, get_auc, which_null = which_null,
                    control_genes = control_genes)
# cov_vec <- sapply(pci_list, get_coverage, beta_true = beta_true,
#                   control_genes = control_genes)

  ## Save parameters
  par_settings <- c(log2foldsd = current_params$log2foldsd,
                    Ngene = current_params$Ngene,
                    log2foldmean = current_params$log2foldmean,
                    skip_gene = current_params$skip_gene,
                    current_seed = new_params$current_seed,
                    nullpi = new_params$nullpi,
                    Nsamp = new_params$Nsamp,
                    ncontrols = new_params$ncontrols,
                    prop_cont = args_val$prop_cont,
                    poisthin = new_params$poisthin)

  return_vec <- c(par_settings, mse_vec, auc_vec)
  
  if (abs(method_list$ashr$g.pi0 - 1) < 1e-10) {return_vec[length(return_vec) - 1] = NA}
  if (abs(method_list$biashr$g.pi0 - 1) < 1e-10) {return_vec[length(return_vec)] = NA}
  
  xtot.time <- proc.time() - start.time
  return(return_vec)
}

itermax <- 1000 ## itermax should be 500
seed_start <- 777

## these change
nullpi_seq   <- c(0.1, 0.5, 0.9)
Nsamp_seq    <- c(10)
ncontrol_seq <- c(100)
prop_control <- c(1)

par_vals <- expand.grid(current_seed = (1 + seed_start):(itermax + seed_start),
                        nullpi       = nullpi_seq,
                        Nsamp        = Nsamp_seq,
                        ncontrols    = ncontrol_seq,
                        prop_cont    = prop_control)
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

## Re order par_vals to get even computation time on all notes
par_vals <- par_vals[sample(seq_len(nrow(par_vals))), ]

suppressPackageStartupMessages(library(tidyverse))
par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    for (inner_list_index in 1:ncol(par_vals)) {
        par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
        names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
    }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 0.8
args_val$Ntotal       <- 1e4
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## Create muscle_mat with most expressed genes
mat <- t(as.matrix(readRDS("~/GitHub/Bi-ASH/data/muscle.rds")))
args_val$mat <- mat[, order(apply(mat, 2, median), decreasing = TRUE)[1:args_val$Ntotal]]
rm(mat)

## try-catch
safe_one_rep <- safely(one_rep)

# oout <- safe_one_rep(par_list[[1000]], args_val)
# oout

## ## If on your own computer, use this
suppressPackageStartupMessages(library(snow))
library(parallel)
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

saveRDS(object = sout, file = "~/GitHub/Bi-ASH/output/nc_sims.rds")
