#### Dependent libraries ####
library(splines)
library(Matrix)
library(numbers)
library(numbers)
library(mvnfast)
library(zeallot)
library(rstan)
library(parallel)
library(refundBayes)
library(refund)
library(mgcv)

#### Helper functions ####
source("vsof.R")
source("bayesian_helpers.R")

#### settings ####
N       <- c(100, 200, 500, 1000)
dots    <- 50
up      <- 1000
maxIt   <- 2000
xi      <- seq(0, 1, by = 0.1) 
B       <- 5000
b1      <- 0
dfs     <- 8 
sig     <- 0.05 # 0.5
rho     <- 0.5
S       <- 100 # 50, 100
s       <- seq(0, 1, length = S)
CovStr  <- matrix(0, nrow = S, ncol = S)
for(r in 1:S){
  for(c in 1:S){
    CovStr[r,c] <- sig*rho^(abs(r - c))
  }
}

#### Bayesian model setup ####
rstan_options(auto_write = TRUE)
my_formula <- y ~ s(tmat, by = lmat*wmat, bs = "ps", k = 8)
# Generate the Stan code once, using the original Bayesian variable names.
# Each parallel worker compiles it once at startup.
set.seed(10312019)
i <- 1; n <- 1
mu1   <- xi[i] * sin((4*s - 1)*pi) * (s + 1222/10000)
fx    <- mvnfast::rmvn(n = N[n], mu = rep(0, length(s)), sigma = CovStr)
e1    <- rnorm(N[1])
y1    <- b1 + fx %*% mu1 + e1
nt    <- ncol(fx)
n_num <- N[1]
tind  <- seq(0, 1, length.out = nt)

data.SoFR <- data.frame(
  wmat = I(fx),
  lmat = I(matrix(1, ncol = nt, nrow = n_num)),
  tmat = I(matrix(tind, ncol = nt, nrow = n_num, byrow = TRUE))
)
data.SoFR$y <- as.numeric(y1)

cat("=== Generating Stan code ===\n")
temp <- refundBayes::sofr_bayes(
  my_formula,
  data = data.SoFR,
  family = gaussian(),
  runStan = FALSE
)
model_code <- temp$stancode

onep  <- ones <- oned <- oneb <- onem <- onet <- array(NA, dim = c(B, length(xi), length(N)))
twop  <- twos <- twod <- twob <- twom <- twot <- array(NA, dim = c(B, length(xi), length(N)))
threep <- threeb <- threem <- threet <- array(NA, dim = c(B, length(xi), length(N)))


#### parallel setup ####
ncores_parallel <- 50
cat(sprintf("=== Setting up %d parallel workers ===\n", ncores_parallel))
cl <- makeCluster(ncores_parallel)

clusterExport(cl, "model_code", envir = environment())
clusterEvalQ(cl, {
  library(splines)
  library(Matrix)
  library(numbers)
  library(mvnfast)
  library(zeallot)
  library(rstan)
  library(refundBayes)
  library(refund)
  library(mgcv)
  rstan_options(auto_write = TRUE)
  worker_model <- rstan::stan_model(model_code = model_code)
  NULL
})
cat("=== All workers compiled ===\n")

clusterExport(
  cl, c("vsof", "construct_fe", "construct_P", "eval_logp_gauss",
        "construct_se", "eval_logp_binpr", "zls_test", "zls_test.vsof",
        "extract_coef_from_stanfit", "GBPV", "my_formula",
        "N", "b1", "dfs", "maxIt", "s", "S", "nt", "tind"),
  envir = environment()
)

#### run simulation ####
iter  <- 0

for(i in 1:length(xi)){
  # i <- 11
  for(n in 1:length(N)){
    # n <- 4
    set.seed(10312019)
    mu1   <- xi[i]*sin((4*s - 1)*pi)*(s+1222/10000)
    fx    <- mvnfast::rmvn(n = N[n], mu = rep(0, length(s)), sigma = CovStr)
    # sx2   <- xi*sx1
    n_num <- N[n]
    clusterExport(cl, c("mu1", "fx", "n", "n_num"), envir = environment())
    
    cat(sprintf("=== Starting i=%d, n=%d: %d repetitions across %d workers ===\n",
                i, n, B, ncores_parallel))
    results <- parSapply(cl, seq_len(B), function(b) {
      # b <- 1
      set.seed(b)
      
      ##### one sample #####
      e1        <- rnorm(N[n])
      y1        <- b1 + fx%*%mu1 + e1
      datab     <- data.frame(y = y1)
      
      stv  <- proc.time()
      vsofb  <- vsof(fixed = y ~ 1, fe = ~ fx, data = datab, Ts = list(fcca = ncol(fx)),
                   fe.spline = list(df = dfs, degree = 3, xi = 0.01), # 12 or 14 or 16, 3, 0.001
                   maxIter = maxIt, verbose = FALSE)
      etv <- as.numeric(proc.time() - stv)
      
      zlsb  <- zls_test(vsofb)
      
      onep   <- zlsb$p.value
      ones   <- zlsb$statistic
      oned   <- zlsb$parameters$df
      oneb   <- sum(vsofb$funcEffect[[1]] - mu1)/S
      onem   <- sum((vsofb$funcEffect[[1]] - mu1)^2)/S
      onet   <- etv[3]
      
      stp         <- proc.time()
      pfrb        <- pfr(y ~ lf(fx, argvals = s*S, k = dfs), data = datab)
      etp         <- as.numeric(proc.time() - stp)
      
      sumb        <- summary(pfrb)
      
      twop   <- sumb$s.table[4]
      twos   <- sumb$s.table[3]
      twod   <- sumb$s.table[1]
      twob   <- sum(as.numeric(coef(pfrb)$value) - mu1)/S
      twom   <- sum((as.numeric(coef(pfrb)$value) - mu1)^2)/S
      twot   <- etp[3]
      
      ##### Bayesian approach #####
      # Use the same response y1 and predictor curves fx as VSOF and PFR.
      stb <- proc.time()
      data.SoFR.iter <- data.frame(
        wmat = I(fx),
        lmat = I(matrix(1, ncol = nt, nrow = n_num)),
        tmat = I(matrix(tind, ncol = nt, nrow = n_num, byrow = TRUE))
      )
      data.SoFR.iter$y <- as.numeric(y1)
      
      formula_use <- refundBayes:::refundBayesformula(my_formula)
      func_comp <- rep(FALSE, length(formula_use$func_var))
      datacode <- refundBayes:::refundBayesdata(
        formula_use, data.SoFR.iter, "gaussian",
        func_comp, intercept = TRUE
      )
      
      fit <- rstan::sampling(
        worker_model,
        data = datacode$standata,
        iter = 4000,
        warmup = 2000,
        chains = 2,
        cores = 1,
        refresh = 0,
        seed = b
      )
      etb <- as.numeric(proc.time() - stb)
      
      result <- extract_coef_from_stanfit(
        fit, my_formula, data.SoFR.iter, family = "gaussian"
      )
      bias <- sum(colMeans(result$func_coef[[1]]) - mu1)/S
      mise <- sum((colMeans(result$func_coef[[1]]) - mu1)^2)/S
      gbpv <- GBPV(result$func_coef[[1]])
      
      threep <- gbpv
      threeb <- bias
      threem <- mise
      threet <- etb[3]
      
      # Return one result vector per replication; workers do not write arrays.
      c(
        onep = unname(onep),
        ones = unname(ones),
        oned = unname(oned),
        oneb = unname(oneb),
        onem = unname(onem),
        onet = unname(onet),
        twop = unname(twop),
        twos = unname(twos),
        twod = unname(twod),
        twob = unname(twob),
        twom = unname(twom),
        twot = unname(twot),
        threep = unname(threep),
        threeb = unname(threeb),
        threem = unname(threem),
        threet = unname(threet)
      )
    })
    
    # Store all worker results in the original arrays in the main session.
    onep[,i,n] <- results["onep", ]
    ones[,i,n] <- results["ones", ]
    oned[,i,n] <- results["oned", ]
    oneb[,i,n] <- results["oneb", ]
    onem[,i,n] <- results["onem", ]
    onet[,i,n] <- results["onet", ]
    twop[,i,n] <- results["twop", ]
    twos[,i,n] <- results["twos", ]
    twod[,i,n] <- results["twod", ]
    twob[,i,n] <- results["twob", ]
    twom[,i,n] <- results["twom", ]
    twot[,i,n] <- results["twot", ]
    threep[,i,n] <- results["threep", ]
    threeb[,i,n] <- results["threeb", ]
    threem[,i,n] <- results["threem", ]
    threet[,i,n] <- results["threet", ]
    
    
    for(b in seq_len(B)){
      ##### simulation controls #####
      iter <- iter + 1
      
      if(mod(iter, dots) == 0){
        cat('.')
      }
      
      if(mod(iter, up) == 0){
        cat(paste("\n",iter,"datasets completed\n"))
      }           
    }
  }
}

stopCluster(cl)

# par(mfrow = c(2,4))
saveRDS(onep, file = 'fs_power_vsof.rds')
saveRDS(twop, file = 'fs_power_pfr.rds')
saveRDS(threep, file = 'fs_power_bayes.rds')

saveRDS(oneb, file = 'fs_bias_vsof.rds')
saveRDS(twob, file = 'fs_bias_pfr.rds')
saveRDS(threeb, file = 'fs_bias_bayes.rds')

saveRDS(onem, file = 'fs_mse_vsof.rds')
saveRDS(twom, file = 'fs_mse_pfr.rds')
saveRDS(threem, file = 'fs_mse_bayes.rds')

saveRDS(onet, file = 'fs_time_vsof.rds')
saveRDS(twot, file = 'fs_time_pfr.rds')
saveRDS(threet, file = 'fs_time_bayes.rds')

onep  <- readRDS('fs_power_vsof.rds')
twop  <- readRDS('fs_power_pfr.rds')
threep <- readRDS('fs_power_bayes.rds')

oneb  <- readRDS('fs_bias_vsof.rds')
twob  <- readRDS('fs_bias_pfr.rds')
threeb <- readRDS('fs_bias_bayes.rds')

onem  <- readRDS('fs_mse_vsof.rds')
twom  <- readRDS('fs_mse_pfr.rds')
threem <- readRDS('fs_mse_bayes.rds')

onet  <- readRDS('fs_time_vsof.rds')
twot  <- readRDS('fs_time_pfr.rds')
threet <- readRDS('fs_time_bayes.rds')


#### nominal alpha ####
pow_one   <- apply(onep < 0.05, c(2,3), mean)
pow_two   <- apply(twop < 0.05, c(2,3), mean)
pow_three <- apply(threep < 0.05, c(2,3), mean)
fscale    <- 1.5

##### size #####
pow_one[1,]
pow_two[1,]
pow_three[1,]

##### power #####
matplot(xi, pow_one, type = 'n', xlab = '', ylab = 'Power, Nominal', main = paste('One Sample, Func.\nS = ', S, ', K = ', dfs, sep = ''), 
        cex.axis = fscale, cex.lab = fscale, cex.main = fscale, ylim = c(0,1))
grid()
matplot(xi, pow_one, add = TRUE, type = 'b', lwd = 3, col = 4, lty = 4:1, pch = 16)
matplot(xi, pow_two, add = TRUE, type = 'b', lwd = 3, col = 2, lty = 4:1, pch = 17)
matplot(xi, pow_three, add = TRUE, type = 'b', lwd = 3, col = 3, lty = 4:1, pch = 15)
legend('bottomright', c('VSOF', 'PFR', 'Bayesian'),
       col = c(4, 2, 3), lty = 1, pch = c(16, 17, 15), bty = 'n')


#### smaller alpha ####
pow_one   <- apply(onep < 0.01, c(2,3), mean)
pow_two   <- apply(twop < 0.01, c(2,3), mean)
pow_three <- apply(threep < 0.01, c(2,3), mean)
fscale    <- 1.5

##### size #####
pow_one[1,]
pow_two[1,]
pow_three[1,]

##### power #####
matplot(xi, pow_one, type = 'n', xlab = '', ylab = 'Power, 0.01', main = paste('One Sample, Func.\nS = ', S, ', K = ', dfs, sep = ''), 
        cex.axis = fscale, cex.lab = fscale, cex.main = fscale, ylim = c(0,1))
grid()
matplot(xi, pow_one, add = TRUE, type = 'b', lwd = 3, col = 4, lty = 4:1, pch = 16)
matplot(xi, pow_two, add = TRUE, type = 'b', lwd = 3, col = 2, lty = 4:1, pch = 17)
matplot(xi, pow_three, add = TRUE, type = 'b', lwd = 3, col = 3, lty = 4:1, pch = 15)
legend('bottomright', c('VSOF', 'PFR', 'Bayesian'),
       col = c(4, 2, 3), lty = 1, pch = c(16, 17, 15), bty = 'n')


#### smallest alpha ####
pow_one   <- apply(onep < 0.001, c(2,3), mean)
pow_two   <- apply(twop < 0.001, c(2,3), mean)
pow_three <- apply(threep < 0.001, c(2,3), mean)
fscale    <- 1.5

##### size #####
pow_one[1,]
pow_two[1,]
pow_three[1,]

##### power #####
matplot(xi, pow_one, type = 'n', xlab = '', ylab = 'Power, 0.001', main = paste('One Sample, Func.\nS = ', S, ', K = ', dfs, sep = ''), 
        cex.axis = fscale, cex.lab = fscale, cex.main = fscale, ylim = c(0,1))
grid()
matplot(xi, pow_one, add = TRUE, type = 'b', lwd = 3, col = 4, lty = 4:1, pch = 16)
matplot(xi, pow_two, add = TRUE, type = 'b', lwd = 3, col = 2, lty = 4:1, pch = 17)
matplot(xi, pow_three, add = TRUE, type = 'b', lwd = 3, col = 3, lty = 4:1, pch = 15)
legend('bottomright', c('VSOF', 'PFR', 'Bayesian'),
       col = c(4, 2, 3), lty = 1, pch = c(16, 17, 15), bty = 'n')


#### bias ####
vsof_bias  <- apply(oneb, c(1,3), mean)
pfr_bias  <- apply(twob, c(1,3), mean)
bayes_bias <- apply(threeb, c(1,3), mean)

bias_mat <- cbind(vsof_bias[,1], pfr_bias[,1], bayes_bias[,1],
                  vsof_bias[,2], pfr_bias[,2], bayes_bias[,2],
                  vsof_bias[,3], pfr_bias[,3], bayes_bias[,3],
                  vsof_bias[,4], pfr_bias[,4], bayes_bias[,4])
colnames(bias_mat) <- paste(rep(N, each = 3),
                            rep(c('VSOF', 'PFR', 'Bayes'), length(N)), sep = '\n')
boxplot(bias_mat, col = rep(c(4,2,3), length(N)), xlab = 'N', ylab = 'Bias')
abline(h = 0, lty = 3)


#### mse ####
vsof_mse   <- apply(onem, c(1,3), mean)
pfr_mse   <- apply(twom, c(1,3), mean)
bayes_mse <- apply(threem, c(1,3), mean)

mse_mat <- cbind(vsof_mse[,1], pfr_mse[,1], bayes_mse[,1],
                 vsof_mse[,2], pfr_mse[,2], bayes_mse[,2],
                 vsof_mse[,3], pfr_mse[,3], bayes_mse[,3],
                 vsof_mse[,4], pfr_mse[,4], bayes_mse[,4])
colnames(mse_mat) <- paste(rep(N, each = 3),
                           rep(c('VSOF', 'PFR', 'Bayes'), length(N)), sep = '\n')
boxplot(mse_mat, col = rep(c(4,2,3), length(N)), xlab = 'N', ylab = 'MSE')


#### time ####
vsof_time  <- apply(onet, c(1,3), mean)
pfr_time  <- apply(twot, c(1,3), mean)
bayes_time <- apply(threet, c(1,3), mean)

time_mat <- cbind(vsof_time[,1], pfr_time[,1], bayes_time[,1],
                  vsof_time[,2], pfr_time[,2], bayes_time[,2],
                  vsof_time[,3], pfr_time[,3], bayes_time[,3],
                  vsof_time[,4], pfr_time[,4], bayes_time[,4])
colnames(time_mat) <- paste(rep(N, each = 3),
                            rep(c('VSOF', 'PFR', 'Bayes'), length(N)), sep = '\n')
boxplot(time_mat, col = rep(c(4,2,3), length(N)), xlab = 'N', ylab = 'Run Time')

