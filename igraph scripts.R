install.packages("igraph")
library(igraph)

### produce degree distribution from single numerical vector

data <- read.csv("degree table.csv")
### grab high affinity subset
High <- data$High.Affinity
## Convert into a frequency table
t <- table(High)
## make the degree distribution (probability density function) plot
plot(t/sum(t), xlab ="Degree", ylab="Frequency")

##revised
plot(t/sum(t), xlim=c(1,4), type='h', axes=FALSE, ann=FALSE)
axis(1, at=1:4, lab=c(1:4))
axis(2)
title(xlab="Degree")
title(ylab="Fraction of Peptides")

