library(spocc)
spocc_download=occ(query='Arabidopsis thaliana', from='gbif', limit=2000);
spocc_download
library(dismo)
wc = getData('worldclim', var='bio', res=5)
library(ENMeval)
library(raster)
occdat <- occ2df(spocc_download)
ext=extent(c(-10, 20, 35, 65))
predictors=crop(wc,ext)
loc=occdat[,c('longitude', 'latitude')]
extr = extract(predictors, loc)
loc = loc[!is.na(extr[,1]),]
eval = ENMevaluate(occ=as.data.frame(loc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
eval@results
which(eval@results$AICc == min(eval@results$AICc))
which(eval@results$avg.test.AUC== max(eval@results$avg.test.AUC))
best=which(eval@results$AICc == min(eval@results$AICc))
best
plot(eval@predictions[[best]])
points(as.data.frame(loc), pch=20, cex =0.1)
est.loc = extract(eval@predictions[[best]], as.data.frame(loc))
est.bg = extract(eval@predictions[[best]], eval@bg.pts)
ev = evaluate(est.loc, est.bg)
thr = threshold(ev)
plot(eval@predictions[[best]] > thr$equal_sens_spec, col = c('lightgrey', 'blue'))
