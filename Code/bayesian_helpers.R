#### Bayesian Approach Helper Functions ####

#### GBPV function ####
GBPV = function(MCMC_draws){
  c(B, t) %<-% dim(MCMC_draws)
  sd_P <- apply(MCMC_draws, 2, sd)
  mean_P <- apply(MCMC_draws, 2, mean)
  z_P <- sapply(1:B, function(x) max(abs(MCMC_draws[x,]- mean_P)/sd_P))
  levels <- c(seq(0.5, 0.01, -0.01), seq(0.009, 0.0001, -0.001), seq(0.0009, 0.00001, -0.0001))
  cb1 <- quantile(z_P, prob = 1-levels, type = 7)
  MAPS = 0.5*rep(1L, t)
  for(k in 2:length(levels)){
    temp_ind = (mean_P - cb1[k]*sd_P) * (mean_P+cb1[k]*sd_P)
    MAPS[which(temp_ind>0)] = levels[k]
  }
  return(min(MAPS))
}

#### Bayes Coefficient extraction function ####
extract_coef_from_stanfit <- function(fit, formula, data, family = "gaussian",
                                      joint_FPCA = NULL, intercept = TRUE) {
  if (inherits(family, "family")) {
    family <- family$family
  }
  formula_use <- refundBayes:::refundBayesformula(formula)
  if (is.null(joint_FPCA)) {
    if (length(formula_use$func_var) == 0) {
      func_comp <- NULL
    } else {
      func_comp <- rep(FALSE, length(formula_use$func_var))
    }
  } else {
    func_comp <- joint_FPCA
  }
  datacode_use <- refundBayes:::refundBayesdata(formula_use, data, family, func_comp, intercept)
  if (!is.null(func_comp)) {
    func_effect_basis <- refundBayes:::extract_basis(formula_use, data, func_comp)
  } else {
    func_effect_basis <- NA
  }
  fit.samp <- rstan::extract(fit, permuted = TRUE)
  if (!is.null(func_comp)) {
    func_coef <- refundBayes:::recon_fun_coef(
      basis     = func_effect_basis,
      fit.samp  = fit.samp,
      func_comp = func_comp,
      trans.mat = datacode_use$trans.mat
    )
  } else {
    func_coef <- NA
  }
  if (length(formula_use$scalar_var) > 0) {
    scalar_coef <- fit.samp$gamma
    colnames(scalar_coef) <- colnames(datacode_use$X_scalar)
  } else {
    scalar_coef <- NULL
  }
  if (intercept == TRUE) {
    int <- fit.samp$eta_0
  } else {
    int <- NULL
  }
  res.fit <- list(
    stanfit      = fit,
    spline_basis = func_effect_basis,
    standata     = datacode_use$standata,
    int          = int,
    scalar_coef  = scalar_coef,
    func_coef    = func_coef,
    family       = family
  )
  class(res.fit) <- "refundBayes"
  return(res.fit)
}
