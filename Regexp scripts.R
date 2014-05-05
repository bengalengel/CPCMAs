##Using reexp in R
##Read in file then change data frame to character vector


peptides = read.table("peptidesnodupe.txt",header=F,sep="\t")
vect=peptides[[1]]
class(vect)
[1] "factor"
vectc=as.character(vect)
class(vectc)
[1] "character"

## an example
> grep("pYVNV",vectc)
 [1]  45  46  47  48  49  50  51  52 461 462 463 464 533 534 535 536

## find occurances of xnx motif
xnx=grepl("pY.N.",vectc)
##write to new file
write.table(xnx,"grb2motif.csv",sep=',',col.names=F,row.names=F)