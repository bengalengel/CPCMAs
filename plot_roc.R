
library(ROCR)

get_roc_data <- function(fname, cutoff) {
    #cutoff <- -7.26
    dtab <- read.table(fname,header=T)
    labels <- dtab$meas > cutoff
    pred <- prediction(dtab$pred, labels)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
  perf
}

plot_panel <- function(fname, main) {
  cutoff <- c(-7.26, -7.0, -6.5, -6.0)
  perf_a <- get_roc_data(fname, cutoff[1])
  x <- unlist(perf_a@x.values)
  y <- unlist(perf_a@y.values)
  plot(x, y, type='l', col='black', xlab='FPR', ylab='TPR', main=main)

  perf_a <- get_roc_data(fname, cutoff[2])
  x <- unlist(perf_a@x.values)
  y <- unlist(perf_a@y.values)
  lines(x,y, type='l', col='orange')

  perf_a <- get_roc_data(fname, cutoff[3])
  x <- unlist(perf_a@x.values)
  y <- unlist(perf_a@y.values)
  lines(x,y, type='l', col='blue')

  perf_a <- get_roc_data(fname, cutoff[4])
  x <- unlist(perf_a@x.values)
  y <- unlist(perf_a@y.values)
  lines(x,y, type='l', col='purple')

  abline(a=0, b=1, col='gray')
  legend("bottomright", legend=c(-7.26, -7.0, -6.5, -6.0), col=c('black','orange','blue','purple'), pch=19)
}


png(filename='plot_roc.png',width=1200,height=1200)
par(mfrow=c(2,2), pty='s', mex=0.8, lwd=2.0,cex=2.0)

plot_panel('scatter.ABL1.xls', 'ABL1')
plot_panel('scatter.GRB2.xls', 'GRB2')
plot_panel('scatter.NCK1.xls', 'NCK1')
plot_panel('scatter.PTP11.xls', 'PTP11')

dev.off()
