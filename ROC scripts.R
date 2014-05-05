## to calculate AUC values single
> oneuMsingle=read.table("single1.txt",header=F,sep="\t")
> fiveuMsingle=read.table("single5.txt",header=F,sep="\t")
> tenuMsingle=read.table("single10.txt",header=F,sep="\t")
> predonesingle=prediction(oneuMsingle[1],oneuMsingle[2])
> predfivesingle=prediction(fiveuMsingle[1],fiveuMsingle[2])
> predtensingle=prediction(tenuMsingle[1],tenuMsingle[2])
> AUConesingle=performance(predonesingle,"auc")
> AUCfivesingle=performance(predfivesingle,"auc")
> AUCtensingle=performance(predtensingle,"auc")
> AUConesingle
An object of class "performance"
Slot "x.name":
[1] "None"

Slot "y.name":
[1] "Area under the ROC curve"

Slot "alpha.name":
[1] "none"

Slot "x.values":
list()

Slot "y.values":
[[1]]
[1] 0.5602679


Slot "alpha.values":
list()

> AUCfivesingle
An object of class "performance"
Slot "x.name":
[1] "None"

Slot "y.name":
[1] "Area under the ROC curve"

Slot "alpha.name":
[1] "none"

Slot "x.values":
list()

Slot "y.values":
[[1]]
[1] 0.7729167


Slot "alpha.values":
list()

> AUCtensingle
An object of class "performance"
Slot "x.name":
[1] "None"

Slot "y.name":
[1] "Area under the ROC curve"

Slot "alpha.name":
[1] "none"

Slot "x.values":
list()

Slot "y.values":
[[1]]
[1] 0.6315789


Slot "alpha.values":
list()


##make an ROC plot
oneuM=read.table("1uM.txt",header=F,sep="\t")
predone=prediction(oneuM[1],oneuM[2])
perfone=performance(predone,"tpr","fpr")
plot(perfone,col="black",lwd=3,add=T)
## do the above for each plot

##add a dotted bisecting line
abline(0,1,lwd=1,lty=2)

##add a legend to the bottom right to be adjusted for later
legend("bottomright", c("10 然","5 然","1 然"),lty=c(1,1,1),lwd=c(3,3,3),col=c("green","red","black"))
legend("bottomright", c("10 然","1 然","10 然","1 然"),lty=c(1,1,2,2),lwd=c(2.5,2.5,2.5,2.5),col=c("red","black","red","black"))
