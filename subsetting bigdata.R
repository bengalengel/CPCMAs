###subsetting bigtable
bigtable <- read.csv("ST2.csv")
##get all KD estimates for each SH2 domain

SHP2N <- subset(bigtable, SH2.domain=="SHP2N", select="Array.ln.KD")


SHP2Nclean <- SHP2N[!is.na(SHP2N)]


SHP2N <- subset(bigtable, SH2.domain=="SHP2N", select="Array.ln.KD")

### don't want to spend the time figuring out how to average duplicates in datatable used excel
SH2s <- read.csv("lnkds.csv")
SHP2N <- SH2s$SHP2N[!is.na(SH2s$SHP2N)]
NCK1 <- SH2s$NCK1[!is.na(SH2s$NCK1)]
GRB2 <- SH2s$GRB2[!is.na(SH2s$GRB2)]
ABL1 <- SH2s$ABL1[!is.na(SH2s$ABL1)]
## can do many tests with this data type
t.test(SHP2N,NCK1)
t.test(SHP2N,ABL1)
t.test(SHP2N,GRB2)
wilcox.test(ABL1,GRB2)
wilcox.test(ABL1,NCK1)
wilcox.test(ABL1,SHP2N)
t.test(ABL1,SHP2N)
ks.test(ABL1,SHP2N)
ks.test(ABL1,GRB2)
ks.test(ABL1,NCK1)
wilcox.test(ABL1,SHP2N)
ks.test(ABL1,SHP2N, alternative="l")
ks.test(ABL1,SHP2N, alternative="g")
ks.test(ABL1,GRB2, alternative="g")
ks.test(ABL1,NCK1, alternative="g")


### but wait is my data normal??? apparently there is controversy over using this test. Ho: The data is drawn from a normal distribution
## if p<.05 then reject this null
shapiro.test(GRB2)##Reject p = .0376 NORMAL QQ plot looks pretty good though
shapiro.test(NCK1)## DON'T REJECT
shapiro.test(SHP2N)## DON'T REJECT
shapiro.test(ABL1)## DON'T REJECT


##Normal QQ plots for diagnostics, they all look pretty good to me....
qqnorm(GRB2)
qqnorm(NCK1)
qqnorm(SHP2N)
qqnorm(ABL1)

## split plots



### need to rearrange into two columns 1 is KDs and the other is SH2 domain (again in excel)
KDANOVA <- read.csv("KDANOVA.csv")
### create a boxplot for visualizing (note ABL1 is a tight distribution)
boxplot(KDANOVA$KD ~ KDANOVA$SH2)

### run anova MAY REQUIRE BALANCED DATA!
SH2.aov <- aov(KDANOVA$KD ~ KDANOVA$SH2)
summary(SH2.aov)
### 'anova' function OK for unequal sample sizes! result is not significant (can't reject that the means are the same)but what is the power?
anova(SH2.aov)

## was not significant but if it was a multiple testing correction would be needed also see 'pairwise.t.test'


TukeyHSD(SH2.aov)

## Kruskal Wallis Test One Way Anova by Ranks nonparametric
kruskal.test(KDANOVA$KD~KDANOVA$SH2) # where y1 is numeric and A is a factor 
## result is NOT SIGNIFICANT therefore can't reject the null which is:such that the probability that a random observation from one group is greater than 
#a random observation from another group is 0.5. 

### Kolmogorov-smirnoff test on 4C2=6

### anderson darling k-sample tests

install.packages("kSamples")
library(kSamples)
### ad.test tests the null that This function uses the Anderson-Darling criterion to test the hypothesis that k independent samples
## with sample sizes n1,n2,..nk arose from a common unspecified distribution function F (x)

ad.out <- ad.test(NCK1,GRB2,ABL1,SHP2N,method="asymptotic")

## RESULT IS SIGNIFICANT! So two-way tests? YES

ad.outNG <- ad.test(NCK1,GRB2,method="asymptotic")
ad.outNA <- ad.test(NCK1,ABL1,method="asymptotic")
ad.outNS <- ad.test(NCK1,SHP2N,method="asymptotic")
ad.outGA <- ad.test(GRB2,ABL1,method="asymptotic")
ad.outGS <- ad.test(GRB2,SHP2N,method="asymptotic")
ad.outAS <- ad.test(ABL1,SHP2N,method="asymptotic")

## get p-values

pNG <- ad.outNG$ad[1,3]
pNA <- ad.outNA$ad[1,3]
pNS <- ad.outNS$ad[1,3]
pGA <- ad.outGA$ad[1,3]
pGS <- ad.outGS$ad[1,3]
pAS <- ad.outAS$ad[1,3]

### insert into vector

p <- c(pNG,pNA,pNS,pGA,pGS,pAS)

### correct with multiple testing method

p.adjust(p, method="bonferroni")

### SIGNIFICANT IS ABL1-GRB2 AND ABL1-NCK1 with conservative bonferroni. ABL1-SHP2 IS VERY CLOSE! BH you can accept it! 

###############overlaid qqplots with qlines###########

##set the range for the qqplots

ylim=range(ABL1,SHP2N,GRB2,NCK1)

###plot each with lines intersecting the 1st and 3rd quantiles
qqnorm(ABL1, ylim=ylim)
qqline(ABL1)
GRB2pts <- qqnorm(GRB2, plot.it=FALSE)
points(GRB2pts, col=2)
qqline(GRB2, col=2)
SHP2Npts <- qqnorm(SHP2N, plot.it=FALSE)
points(SHP2Npts, col=3)
qqline(SHP2N, col=3)
NCK1pts <- qqnorm(NCK1, plot.it=FALSE)
points(NCK1pts, col=4)
qqline(NCK1, col=4)

##summary stats
> summary(SHP2N)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.3573  2.6440  3.6010  4.2610  6.5830  7.8780 
> sd(SHP2N)
[1] 2.338168
> summary(ABL1)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
1.612   3.111   3.963   4.045   4.715   6.892 
> sd(ABL1)
[1] 1.208477
> summary(NCK1)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
1.143   2.918   4.658   4.702   6.418   8.583 
> sd(NCK1)
[1] 2.079857
> summary(GRB2)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-1.025   2.008   4.247   3.914   5.900   7.991 
> sd(GRB2)
[1] 2.430915







