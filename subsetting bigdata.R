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


### need to rearrange into two columns 1 is KDs and the other is SH2 domain (again in excel)
KDANOVA <- read.csv("KDANOVA.csv")
### create a boxplot for visualizing (note ABL1 is a tight distribution)
boxplot(KDANOVA$KD ~ KDANOVA$SH2)

### run anova MAY REQUIRE BALANCED DATA!
SH2.aov <- aov(KDANOVA$KD ~ KDANOVA$SH2)
summary(SH2.aov)
### 'anova' function OK for unequal sample sizes! result is not significant (can't reject that the means are the same)but what is the power?
anova(SH2.aov)

## was not significant but if it was a multiple testing correction would be needed

TukeyHSD(SH2.aov)

## Kruskal Wallis Test One Way Anova by Ranks nonparametric
kruskal.test(KDANOVA$KD~KDANOVA$SH2) # where y1 is numeric and A is a factor 
## result is not significant therefore can't reject the null which is:such that the probability that a random observation from one group is greater than 
#a random observation from another group is 0.5.

