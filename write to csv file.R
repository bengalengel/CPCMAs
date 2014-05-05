###########################################################################
##Writes the regression coefficients in a csv file ###
 
###########################################################################
 
# Regression Coefficient as .csv , can be extended to add other objects too
# v 1.1
# 10/22/2012
# Anupam Anand
##email:  123@AdotB where 123 means anupam and  A=umd, B=edu
 
## reg_model is the regression model, fname is the name of the csv file you want 
regr_tab <- function(reg_model, fname){
 
  # coefficients in dataframe
regr_tab <- data.frame(summary(reg_model)$coefficients)
 
# grab the coefficients
colnames(regr_tab) <- colnames(summary(reg_model)$coefficients)
# get the p-vals 
regr_tab[ ,4] <- ifelse(regr_tab[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tab[ ,4] < .01, "< 0.01", 
                               round(regr_tab[ ,4], 3)))
 
 
# format the table
summary = format(regr_tab, autoformat = 1)
 
# write it as a csv file 
write.csv(summary, paste(fname,"_model_coeff.csv", sep=''))
}
write.csv(summary,"append.csv",TRUE)
Warning message:
In write.table(summary, "append.csv", TRUE, col.names = NA, sep = ",",  :
  appending column names to file
#the above works but not append=TRUE?
> snrvector=apply(counting[,sort],2,sum)
> snrvector
 [1] 0 0 1 2 3 3 3 3 3 3 3 3
> snrvectortest=data.frame(snrvector)
> snrvectortest
   snrvector
1          0
2          0
3          1
4          2
5          3
6          3
7          3
8          3
9          3
10         3
11         3
12         3
> snrvectortest=data.frame(snrvector)[,2]
Error in `[.data.frame`(data.frame(snrvector), , 2) : 
  undefined columns selected
> snrvectortest=data.frame(snrvector)[,1]
> snrvectortest
 [1] 0 0 1 2 3 3 3 3 3 3 3 3
> > transposesnrvect=t(snrvector)
> transposesnrvect
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,]    0    0    1    2    3    3    3    3    3     3     3     3
> combine=data.frame(regr_tab,transposesnrvect)
> combine
  Estimate Std..Error    t.value Pr...t.. X1 X2 X3 X4 X5 X6 X7 X8 X9 X10
a 17.93446 0.05396034 332.363649  < 0.001  0  0  1  2  3  3  3  3  3   3
d  1.17264 0.18448704   6.356218  < 0.001  0  0  1  2  3  3  3  3  3   3
  X11 X12
a   3   3
d   3   3
> write.csv(combine,"transpose.csv")
> > write.csv(combine,"transpose.csv",TRUE)
Warning message:
In write.table(combine, "transpose.csv", TRUE, col.names = NA, sep = ",",  :
  appending column names to file
> write.csv(combine,"transpose.csv",TRUE)
Warning message:
In write.table(combine, "transpose.csv", TRUE, col.names = NA, sep = ",",  :
  appending column names to file
> write.csv(combine,"transpose.csv",TRUE)
Warning message:
In write.table(combine, "transpose.csv", TRUE, col.names = NA, sep = ",",  :
  appending column names to file
> t4$names
NULL
> w$names
NULL
> names[j]
[1] "RQLNpYIQVDLE"
> identifier=data.frame(names[j],ID[j,i])
> identifier
      names.j. ID.j..i.
1 RQLNpYIQVDLE      B02
> supercombo=data.frame(identifier,combine)
> supercombo
      names.j. ID.j..i. Estimate Std..Error    t.value Pr...t.. X1 X2 X3 X4
a RQLNpYIQVDLE      B02 17.93446 0.05396034 332.363649  < 0.001  0  0  1  2
d RQLNpYIQVDLE      B02  1.17264 0.18448704   6.356218  < 0.001  0  0  1  2
  X5 X6 X7 X8 X9 X10 X11 X12
a  3  3  3  3  3   3   3   3
d  3  3  3  3  3   3   3   3
> write.csv(supercombo,"supercombo.csv",TRUE)
Warning message:
In write.table(supercombo, "supercombo.csv", TRUE, col.names = NA,  :
  appending column names to file
> write.csv(supercombo,"supercombo.csv",TRUE)
Warning message:
In write.table(supercombo, "supercombo.csv", TRUE, col.names = NA,  :
  appending column names to file
No help files found matching ‘CSV files’ using fuzzy matching
No help files found matching ‘CSV files’ using fuzzy matching
> write.csv(supercombo,"supercombo.csv",col.names=NA,TRUE)
Warning messages:
1: In write.csv(supercombo, "supercombo.csv", col.names = NA, TRUE) :
  attempt to set 'col.names' ignored
2: In write.table(supercombo, "supercombo.csv", col.names = NA, TRUE,  :
  appending column names to file
> write.table(supercombo,"supercombo2.csv",col.names=FALSE,append=TRUE)
> write.table(supercombo,"supercombo2.csv",sep=',',col.names=FALSE,append=TRUE)
> write.table(supercombo,"supercombo2.csv",sep=',',col.names=FALSE,append=TRUE)
> 