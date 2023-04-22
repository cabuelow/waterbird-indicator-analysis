thin = 100
samples = 1000
nChains = 4
set.seed(1)
ptm = proc.time()
m = sampleMcmc(m, samples = samples, thin = thin,
               adaptNf = rep(ceiling(0.4*samples*thin),1),
               transient = ceiling(0.5*samples*thin),
               nChains = nChains, nParallel = nChains,
               initPar = "fixed effects")
computational.time = proc.time() - ptm
filename = file.path(ModelDir, paste("model_", as.character(model), "_",
                                     c("pa","abundance")[modeltype], "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples), ".Rdata", sep = ""))
save(m, file=filename, computational.time)