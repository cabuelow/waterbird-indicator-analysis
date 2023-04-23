get_spprich <- function (hM, Gradient, predY, measure, xlabel = NULL, ylabel = NULL, 
                         index = 1, q = c(0.025, 0.5, 0.975), cicol = rgb(0, 0, 1, 
                                                                          alpha = 0.5), pointcol = "lightgrey", pointsize = 1, 
                         showData = FALSE, jigger = 0, yshow = NA, showPosteriorSupport = TRUE, 
                         main, ...) 
{
  out <- list()
  Pr = NA
  if (is.null(xlabel)) {
    switch(class(hM$X)[1L], matrix = {
      xlabel = colnames(Gradient$XDataNew)[[1]]
    }, list = {
      xlabel = colnames(Gradient$XDataNew[[1]])[[1]]
    })
  }
  switch(class(hM$X)[1L], matrix = {
    xx = Gradient$XDataNew[, 1]
  }, list = {
    if (measure == "Y") {
      xx = Gradient$XDataNew[[index]][, 1]
    } else {
      xx = Gradient$XDataNew[[1]][, 1]
    }
  })
  ngrid = length(xx)
  if (measure == "S") {
    predS = abind::abind(lapply(predY, rowSums), along = 2)
    Pr = mean(predS[ngrid, ] > predS[1, ])
    qpred = apply(predS, c(1), quantile, probs = q, na.rm = TRUE)
    qpred2 = apply(predS, c(1), quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    if (is.null(ylabel)) {
      ylabel = "Summed response"
      if (all(hM$distr[, 1] == 2)) {
        ylabel = "Species richness"
      }
      if (all(hM$distr[, 1] == 3)) {
        ylabel = "Total count"
      }
    }
  }
  if (measure == "Y") {
    tmp = abind::abind(predY, along = 3)
    Pr = mean(tmp[ngrid, index, ] > tmp[1, index, ])
    qpred = apply(tmp, c(1, 2), quantile, probs = q, na.rm = TRUE)
    qpred = qpred[, , index]
    if (is.null(ylabel)) {
      ylabel = hM$spNames[[index]]
    }
  }
  if (measure == "T") {
    if (all(hM$distr[, 1] == 1)) {
      predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)), 
                                                                      hM$nt), ncol = hM$nt))
    }
    else {
      predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a), 
                                                                 hM$nt), ncol = hM$nt))
    }
    predT = abind::abind(predT, along = 3)
    Pr = mean(predT[ngrid, index, ] > predT[1, index, ])
    qpred = apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)
    qpred = qpred[, , index]
    if (is.null(ylabel)) {
      ylabel = hM$trNames[[index]]
    }
  }
  lo = qpred[1, ]
  hi = qpred[3, ]
  me = qpred[2, ]
  lo1 = min(lo, yshow, na.rm = TRUE)
  hi1 = max(hi, yshow, na.rm = TRUE)
  
  datToPlot <- data.frame(x = xx, richness_50 = qpred[2,], richness_97.5 = qpred[3,], richness_2.5 = qpred[1,],
                          richness_75 = qpred2[3,] , richness_25 = qpred2[1,])
  out[[1]] <- datToPlot
  if (showData) {
    switch(class(hM$X)[1L], matrix = {
      XDatacol = which(colnames(Gradient$XDataNew)[[1]] == 
                         colnames(hM$XData))
    }, list = {
      XDatacol = which(colnames(Gradient$XDataNew[[1]])[[1]] == 
                         colnames(hM$XData[[1]]))
    })
    if (measure == "S") {
      pY = rowSums(hM$Y, na.rm = TRUE)
    }
    if (measure == "Y") {
      pY = hM$Y[, index]
    }
    if (measure == "T") {
      if (all(hM$distr[, 1] == 1)) {
        tmp = (exp(hM$Y) %*% hM$Tr)/matrix(rep(rowSums(exp(hM$Y)), 
                                               hM$nt), ncol = hM$nt)
      }
      else {
        tmp = (hM$Y %*% hM$Tr)/matrix(rep(rowSums(hM$Y), 
                                          hM$nt), ncol = hM$nt)
      }
      pY = tmp[, index]
    }
    if (!is.numeric(pY)) 
      pY <- as.numeric(pY)
    switch(class(hM$X)[1L], matrix = {
      pX = hM$XData[, XDatacol]
    }, list = {
      pX = hM$XData[[1]][, XDatacol]
    })
    hi1 <- max(hi1, max(pY, na.rm = TRUE))
    lo1 <- min(lo1, min(pY, na.rm = TRUE))
  }
  
  spp.pred <- data.frame(x = pX, y = pY)
  out[[2]] <- spp.pred
  return(out)
}
