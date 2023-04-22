ds_resid <- function (yobs, ypred, plotds = TRUE, family = "poisson", phis = NULL, 
          trial.size = 1) 
{
  if (family == "poisson") {
    a <- ppois(yobs - 1, ypred)
    b <- ppois(yobs, ypred)
  }
  if (family %in% c("nbinom", "negbinomial")) {
    a <- pnbinom(yobs - 1, mu = ypred, size = 1/phis)
    b <- pnbinom(yobs, mu = ypred, size = 1/phis)
  }
  if (family == "binomial") {
    a <- pbinom(yobs - 1, trial.size, prob = ypred)
    b <- pbinom(yobs, trial.size, prob = ypred)
  }
  u <- runif(n = length(yobs), min = a, max = b)
  dsres <- qnorm(u)
  group <- as.factor(rep(1:ncol(ypred), each = nrow(ypred)))
  if (plotds) {
    par(mfrow = c(1, 2))
    plot(ypred, dsres, xlab = "Prediction", ylab = "Dunn-Smyth residual", col = group)
    abline(h = 0)
    qqnorm(dsres)
    qqline(dsres)
  }
  return(dsres)
}