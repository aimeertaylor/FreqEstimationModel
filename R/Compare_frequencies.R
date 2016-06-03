# Compare frequencies estimated under the model and using counting
load('../FreqResultsStore_counting.RData')

load('../MCMCFreqResults/FreqResultsStore_med_CI.RData')
FreqResultsStore[FreqResultsStore == 0] <- NA

#-----------------------------------------------
# Overall
#-----------------------------------------------
plot(x = FreqResultsStore_counting[,'50%',,,], y = FreqResultsStore[,'50%',,,], pch = 20,
     bty = 'n', cex = 0.5, xlab = 'Frequency estimates, simple counting',
     ylab = 'Posterior median frequency estimates')
segments(y0 = FreqResultsStore[,'2.5%',,,],
         y1 = FreqResultsStore[,'97.5%',,,],
         x0 = FreqResultsStore_counting[,'50%',,,],
         x1 = FreqResultsStore_counting[,'50%',,,],
         adjustcolor('black', alpha = 0.5), lwd = 0.5)


#-----------------------------------------------
# By site
# Note, confidence interval around mutant
#-----------------------------------------------
cols <- rainbow(4, alpha = 0.5); names(cols) <- sites
plot(NULL, ylim = c(0,1), xlim = c(0,1), pch = 20,
     bty = 'n', ylab = 'Frequency estimates, simple counting',
     xlab = 'Posterior median frequency estimates')
legend('top', legend = sites, col = cols, bty = 'n', pch = 20,
       horiz = TRUE, cex = 0.8, x.intersp = 1)

for(site in sites){
points(y = FreqResultsStore_counting[,'50%',,site,],
       x = FreqResultsStore[,'50%',,site,], pch = 20,
       cex = 0.5, col = cols[site])

# segments(y0 = FreqResultsStore[,'2.5%',,site,],
#          y1 = FreqResultsStore[,'97.5%',,site,],
#          x0 = FreqResultsStore_counting[,'50%',,site,],
#          x1 = FreqResultsStore_counting[,'50%',,site,],
#          col = cols[site], lwd = 0.5)
#
# segments(y0 = FreqResultsStore[,'50%',,site,],
#          y1 = FreqResultsStore[,'50%',,site,],
#          x0 = FreqResultsStore_counting[,'2.5%',,site,],
#          x1 = FreqResultsStore_counting[,'97.5%',,site,],
#          col = cols[site], lwd = 0.5)
}
abline(a = 0, b=1)
