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


> row <- rbind(c(1:4),total,high)
> row
[,1] [,2] [,3] [,4]
1    2    3    4
total    5   11   20   45
high    31   25    7    0
> col <- cbind(c(1:4),total,high)
> col
total high
[1,] 1     5   31
[2,] 2    11   25
[3,] 3    20    7
[4,] 4    45    0
> fisher.test(col)
Error in fisher.test(col) : FEXACT error 7.
LDSTP is too small for this problem.
Try increasing the size of the workspace.
> 2e8
[1] 2e+08
> fisher.test(col,workspace=2e8)

Fisher's Exact Test for Count Data

data:  col 
p-value < 2.2e-16
alternative hypothesis: two.sided 

> fisher.test(col,workspace=2e8, hybrid="True")

Fisher's Exact Test for Count Data

data:  col 
p-value < 2.2e-16
alternative hypothesis: two.sided 

> fisher.test(col,workspace=2e8, hybrid="FALSE")

Fisher's Exact Test for Count Data

data:  col 
p-value < 2.2e-16
alternative hypothesis: two.sided 

> fisher.test(col,workspace=2e8, hybrid="TRUE")

Fisher's Exact Test for Count Data

data:  col 
p-value < 2.2e-16
alternative hypothesis: two.sided 

> fisher.test(row,workspace=2e8, hybrid="TRUE")

Fisher's Exact Test for Count Data

data:  row 
p-value < 2.2e-16
alternative hypothesis: two.sided 

> fisher.test(col, simulate.p.value=TRUE, B=1e6)

Fisher's Exact Test for Count Data with simulated p-value (based on 1e+06 replicates)

data:  col 
p-value = 1e-06
alternative hypothesis: two.sided 

> Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
                +               dimnames = list(income=c("< 15k", "15-25k", "25-40k", "> 40k"),
                                                +                               satisfaction=c("VeryD", "LittleD", "ModerateS", "VeryS")))
> Job
satisfaction
income   VeryD LittleD ModerateS VeryS
< 15k      1       3        10     6
15-25k     2       3        10     7
25-40k     1       6        14    12
> 40k      0       1         9    11
> fisher.test(JOb)
Error in is.data.frame(x) : object 'JOb' not found
> fisher.test(Job)

Fisher's Exact Test for Count Data

data:  Job 
p-value = 0.7827
alternative hypothesis: two.sided 

> fisher.test(Job, simulate.p.value=T, B=1e6)

Fisher's Exact Test for Count Data with simulated p-value (based on 1e+06 replicates)

data:  Job 
p-value = 0.7829
alternative hypothesis: two.sided 

> fisher.test(row, simulate.p.value=TRUE, B=1e6)

Fisher's Exact Test for Count Data with simulated p-value (based on 1e+06 replicates)

data:  row 
p-value = 1e-06
alternative hypothesis: two.sided 

> fisher.test(col, simulate.p.value=TRUE, B=1e7)

	Fisher's Exact Test for Count Data with simulated p-value (based on 1e+07 replicates)

data:  col 
p-value = 1e-07
alternative hypothesis: two.sided 

