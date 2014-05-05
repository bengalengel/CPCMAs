source("function_miao.R")

t0=read.table("NCK1_data.txt",header=T)
t1=read.table("NCK1_background.txt",header=T)

i=1
##choose the array with replicate i
replicate_fore=t0[,9]
replicate_back=t1[,16]

t=t0[replicate_fore==i,]
t_b=t1[replicate_back==i,]

##select out the signal intensity from multiple PMT and match the names
SI_fore=t[,c(25,23,21,19,17,15)]
SI_back=t_b[,c(15,13,11,9,7,5)]

sat_fore=t[,c(24,22,20,18,16,14)]
sat_back=t_b[,c(14,12,10,8,6,4)]

names(SI_back)=names(SI_fore) 
names(sat_back)=names(sat_fore) 


##combine the foreground and background SI from the same array
SI=rbind(SI_fore,SI_back)
sat=rbind(sat_fore,sat_back)







########### Linear Scatter Uncorrected ###############################
################# lowest/hightest PMT Pairs###############################
par(mfrow=c(2,2))
plot(SI[,6],SI[,5],cex=0.5,xlab="PMT 450",ylab="PMT 550",main="Lowest PMT Pair (original scale)")
points(SI[sat[,5]>=1,2],SI[sat[,5]>=1,1],cex=0.5,col="red")

plot(SI[,2],SI[,1],cex=0.5,xlab="PMT 850",ylab="PMT 950",main="Highest PMT Pair (original scale)")
points(SI[sat[,1]>=1,2],SI[sat[,1]>=1,1],cex=0.5,col="red")




########### Log-log Scatter Uncorrected ###############################
################# lowest/hightest PMT Pairs###############################
plot(SI[,6],SI[,5],cex=0.5,xlab="PMT 450",ylab="PMT 550",main="Lowest PMT Pair (log scale)",log="xy")
points(SI[sat[,5]>=1,2],SI[sat[,5]>=1,1],cex=0.5,col="red")
plot(SI[,2],SI[,1],cex=0.5,xlab="PMT 850",ylab="PMT 950",main="Highest PMT Pair (log scale)",log="xy")
points(SI[sat[,1]>=1,2],SI[sat[,1]>=1,1],cex=0.5,col="red")

dev.copy(pdf,"scatter.pdf")
dev.off()
dev.off()


################# Linear Scatter Uncorrected ALL PMT Pairs#############################


plot(SI[,2],SI[,1],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="Globle Scatter Plot for ALL PMT Pairs",col=1,xlim=c(0,2*10^8),ylim=c(0,2*10^8))
for(i in 1:5){
points(SI[,6],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+1)
}
for(i in 1:4){
points(SI[,5],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+6)
}
for(i in 1:3){
points(SI[,4],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+10)
}

for(i in 1:2){
points(SI[,3],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+13)
}

dev.copy(pdf,"scatter2.pdf")
dev.off()
dev.off()

######## Linear Scatter Uncorrected ALL PMT Pairs (Near Origin)#############


plot(SI[,2],SI[,1],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="Local Scatter Plot for ALL PMT Pairs",col=1,xlim=c(0,10^6),ylim=c(0,10^6))

for(i in 1:5){
points(SI[,6],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+1)
}
for(i in 1:4){
points(SI[,5],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+6)
}
for(i in 1:3){
points(SI[,4],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+10)
}

for(i in 1:2){
points(SI[,3],SI[,i],cex=0.2,xlab="Lower PMT",ylab="Higher PMT",main="ALL PMT Pairs",col=i+13)
}


dev.copy(pdf,"scatter3.pdf")
dev.off()
dev.off()





#################### MA plot Uncorrected ###############################


################# lowest/hightest PMT Pairs###############################
par(mfrow=c(2,2))

plot(0.5*log(SI[,6])+0.5*log(SI[,5]),log(SI[,5])-log(SI[,6]),cex=0.5,xlab="A",ylab="M",main="Lowest PMT Pair")
points(0.5*log(SI[sat[,5]>=1,6])+0.5*log(SI[sat[,5]>=1,5]),log(SI[sat[,5]>=1,5])-log(SI[sat[,5]>=1,6]),cex=0.5,col="red")

plot(0.5*log(SI[,1])+0.5*log(SI[,2]),log(SI[,1])-log(SI[,2]),cex=0.5,xlab="A",ylab="M",main="Highest PMT Pair")
points(0.5*log(SI[sat[,1]>=1,1])+0.5*log(SI[sat[,1]>=1,2]),log(SI[sat[,1]>=1,1])-log(SI[sat[,1]>=1,2]),cex=0.5,col="red")

dev.copy(pdf,"MA_pair.pdf")
dev.off()
dev.off()


################# ALL PMT Pairs###############################

plot(0.5*log(SI[,2])+0.5*log(SI[,1]),log(SI[,1])-log(SI[,2]),cex=0.2,xlab="A",ylab="M",main="MA Plot for ALL PMT Pairs",col=1,xlim=c(10,20),ylim=c(-5,10))
for(i in 1:5){
points(0.5*log(SI[,6])+0.5*log(SI[,i]),log(SI[,i])-log(SI[,6]),cex=0.2,main="ALL PMT Pairs",col=i+1)
}
for(i in 1:4){
points(0.5*log(SI[,5])+0.5*log(SI[,i]),log(SI[,i])-log(SI[,5]),cex=0.2,main="ALL PMT Pairs",col=i+6)
}
for(i in 1:3){
points(0.5*log(SI[,4])+0.5*log(SI[,i]),log(SI[,i])-log(SI[,4]),cex=0.2,main="ALL PMT Pairs",col=i+10)
}

for(i in 1:2){
points(0.5*log(SI[,3])+0.5*log(SI[,i]),log(SI[,i])-log(SI[,3]),cex=0.2,main="ALL PMT Pairs",col=i+13)
}


dev.copy(pdf,"MA_all.pdf")
dev.off()
dev.off()


############################### Offset_corrected MA Plot ######################
preset=as.matrix(SI)/10^6
saturation=(sat>=1)
scan=multiPMT(preset,saturation)
offset_correct=scan$offset_correct*10^6


####corrected MA plot (with saturation points)##########################
plot(0.5*log(offset_correct[,2])+0.5*log(offset_correct[,1]),log(offset_correct[,1])-log(offset_correct[,2]),cex=0.2,xlab="A",ylab="M",main="MA Plot for All PMT Pairs",col=1,xlim=c(10,20),ylim=c(-5,10))
for(i in 1:5){
points(0.5*log(offset_correct[,6])+0.5*log(offset_correct[,i]),log(offset_correct[,i])-log(offset_correct[,6]),cex=0.2,main="ALL PMT Pairs",col=i+1)
}
for(i in 1:4){
points(0.5*log(offset_correct[,5])+0.5*log(offset_correct[,i]),log(offset_correct[,i])-log(offset_correct[,5]),cex=0.2,main="ALL PMT Pairs",col=i+6)
}
for(i in 1:3){
points(0.5*log(offset_correct[,4])+0.5*log(offset_correct[,i]),log(offset_correct[,i])-log(offset_correct[,4]),cex=0.2,main="ALL PMT Pairs",col=i+10)
}

for(i in 1:2){
points(0.5*log(offset_correct[,3])+0.5*log(offset_correct[,i]),log(offset_correct[,i])-log(offset_correct[,3]),cex=0.2,main="ALL PMT Pairs",col=i+13)
}


dev.copy(pdf,"MA_correct1.pdf")
dev.off()
dev.off()





############## Scatter plot after correction and imputation############################

par(mfrow=c(2,2))
plot(offset_correct[,6],offset_correct[,5],cex=0.5,xlab="PMT 450",ylab="PMT 550",main="Lowest PMT Pair (original scale)")
points(offset_correct[sat[,5]>=1,2],offset_correct[sat[,5]>=1,1],cex=0.5,col="red")

plot(offset_correct[,2],offset_correct[,1],cex=0.5,xlab="PMT 850",ylab="PMT 950",main="Highest PMT Pair (original scale)")
points(offset_correct[sat[,1]>=1,2],offset_correct[sat[,1]>=1,1],cex=0.5,col="red")




######################## on log scale############################
plot(offset_correct[,6],offset_correct[,5],cex=0.5,xlab="PMT 450",ylab="PMT 550",main="Lowest PMT Pair (log scale)",log="xy")
points(offset_correct[sat[,5]>=1,2],offset_correct[sat[,5]>=1,1],cex=0.5,col="red")
plot(offset_correct[,2],offset_correct[,1],cex=0.5,xlab="PMT 850",ylab="PMT 950",main="Highest PMT Pair (log scale)",log="xy")
points(offset_correct[sat[,1]>=1,2],offset_correct[sat[,1]>=1,1],cex=0.5,col="red")



dev.copy(pdf,"scatter_pair_corrected.pdf")
dev.off()
dev.off()






################# density #####################
par(mfrow=c(2,2))
plot(density(log(SI[,6])),ylim=c(0,1.5),main="density of raw SI",xlab=expression(log(SI)))
for (i in 1:5){
lines(density(log(SI[,i])))
}



plot(density(log(offset_correct[,6])),ylim=c(0,1.5),main="density of bias corrected SI",xlab=expression(log(SI-a_i)))
for (i in 1:5){
lines(density(log(offset_correct[,i])))
}



normalize=scan$normalize

plot(density(log(normalize[,6])),ylim=c(0,1.5),main="density of scale/bias corrected SI",xlab=expression(log(SI-a_i)/b_i))
for (i in 1:5){
lines(density(log(normalize[,i])))
}

dev.copy(pdf,"density.pdf")
dev.off()
dev.off()


########################## MA plot after bias/scale correction ###############################
plot(0.5*log(normalize[,2])+0.5*log(normalize[,1]),log(normalize[,1])-log(normalize[,2]),cex=0.2,xlab="A",ylab="M",main="MA Plot for ALL PMT Pairs",col=1,xlim=c(-1,10),ylim=c(-5,5))
for(i in 1:5){
points(0.5*log(normalize[,6])+0.5*log(normalize[,i]),log(normalize[,i])-log(normalize[,6]),cex=0.2,main="ALL PMT Pairs",col=i+1)
}
for(i in 1:4){
points(0.5*log(normalize[,5])+0.5*log(normalize[,i]),log(normalize[,i])-log(normalize[,5]),cex=0.2,main="ALL PMT Pairs",col=i+6)
}
for(i in 1:3){
points(0.5*log(normalize[,4])+0.5*log(normalize[,i]),log(normalize[,i])-log(normalize[,4]),cex=0.2,main="ALL PMT Pairs",col=i+10)
}

for(i in 1:2){
points(0.5*log(normalize[,3])+0.5*log(normalize[,i]),log(normalize[,i])-log(normalize[,3]),cex=0.2,main="ALL PMT Pairs",col=i+13)
}

dev.copy(pdf,"MA_offset_bias_corrected.pdf")
dev.off()
dev.off()


###################inter-array correction##############################

##before
t0=read.table("CBD_data.txt",header=T)
replicate=t0[,8]
BSV=t0[,27]
par(mfrow=c(2,2))
boxplot(BSV~replicate,log="y",main="before inter-array normalization")
plot(BSV[replicate==1],BSV[replicate==2],log="xy",cex=0.5,xlab="replicate 1",ylab="replicate 2")
abline(0,1)
plot(BSV[replicate==1],BSV[replicate==3],log="xy",cex=0.5,xlab="replicate 1",ylab="replicate 3")
abline(0,1)
plot(BSV[replicate==2],BSV[replicate==3],log="xy",cex=0.5,xlab="replicate 2",ylab="replicate 3")
abline(0,1)


dev.copy(pdf,"CBD_inter_array.pdf")
dev.off()
dev.off()

######after
replicate=t0[,8]
BSV=t0[,28]
par(mfrow=c(2,2))
boxplot(BSV~replicate,log="y",main="after inter-array normalization")
plot(BSV[replicate==1],BSV[replicate==2],log="xy",cex=0.5,xlab="replicate 1",ylab="replicate 2")
abline(0,1)
plot(BSV[replicate==1],BSV[replicate==3],log="xy",cex=0.5,xlab="replicate 1",ylab="replicate 3")
abline(0,1)
plot(BSV[replicate==2],BSV[replicate==3],log="xy",cex=0.5,xlab="replicate 2",ylab="replicate 3")
abline(0,1)


dev.copy(pdf,"CBD_inter_array_after.pdf")
dev.off()
dev.off()



#####intra-array normalization######################################
block=t0[,1]
block=as.factor(block)
levels(block)=c(1,3,5,7,9,11,2,4,6,8,10,12)
BSV=t0[,28]
boxplot(BSV~block,log="y",main="before intra-array normalization",xlab="Block",ylab="SI")
dev.copy(pdf,"intra_array_before.pdf")
dev.off()
dev.off()

BSV=t0[,29]
boxplot(BSV~block,log="y",main="after intra-array normalization",xlab="Block",ylab="SI")
dev.copy(pdf,"intra_array_after.pdf")
dev.off()
dev.off()







######################### curve fitting #######################
x=read.table("Abl1_data.txt",header=T)
names=c("RQLNpYIQVDLE","FESRpYQQPFED","PPALpYAEPLDS","PSSVpYVPDEWE","NGLNpYIDLDLV",
"SPGEpYVNIEFG","KEEGpYELPYNP","DDPSpYVNVQNL","EDDGpYDVPKPP","DQHDpYDSVASD","AEPQpYEEIPIY")
##corresponding ID
ID=matrix(rep(0,22),nrow=11)
ID[1,]=c("B02","P19")
ID[2,]=c("A23","P16")
ID[3,]=c("A24","P17")
ID[4,]=c("B03","P20")
ID[5,]=c("A21","P14")
ID[6,]=c("A19","P12")
ID[7,]=c("A20","P13")
ID[8,]=c("A12","P09")
ID[9,]=c("A09","P07")
ID[10,]=c("B13","P24")
ID[11,]=c("A16","P10")


######################### original scale (curve) #######################
##The first ID
i=1

plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration (nM)",ylab="normalized SI",main="Curve Fitting (I)")

for (j in 1:11){
if ((j!=5)&(j!=3)){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]
matrix(t1[,5],nrow=3,byrow=T)
y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##points(exp(x0),exp(y0),col=j,cex=0.5,pch=j)
a=coef(reg2)[1]
d=coef(reg2)[2]
curve(exp(llog(log(x),a,d)),from=0.5,to=5000,col=j,add=T,log="x")
}
}

legend("topright",col=c(1:2,4,6:11),lty=1,paste(names,ID[,1])[c(1:2,4,6:11)],cex=0.8)

dev.copy(pdf,"curve_fitting_1.pdf")
dev.off()
dev.off()





######################### the second ID #######################

plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration(nM)",ylab="normalized SI",main="Curve Fitting (II)")

i=2
for (j in 1:11){
if (j!=5){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]
matrix(t1[,5],nrow=3,byrow=T)
y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##points(exp(x0),exp(y0),col=j,cex=0.5,pch=j)
a=coef(reg2)[1]
d=coef(reg2)[2]
curve(exp(llog(log(x),a,d)),from=0.5,to=5000,col=j,add=T,log="x")
}
}


legend("topright",col=c(1:2,4:11),lty=1,paste(names,ID[,2])[c(1:2,4:11)],cex=0.8)

dev.copy(pdf,"curve_fitting_2.pdf")
dev.off()
dev.off()






######################### original scale (points) #######################

i=1


plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration(nM)",ylab="normalized SI",main="Average Points")
con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
for (j in 1:11){
if ((j!=5)&(j!=3)){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]

t1[(t1[,5]==0)|(t1[,2]!=0),4]=NA
matrix=matrix(t1[,4],nrow=3,byrow=T)
y=colMeans(matrix,na.rm=T)

points(con,y,col=j,cex=0.5,pch=j)

}
}

legend("topright",col=c(1:2,4,6:11),pch=c(1:2,4,6:11),paste(names,ID[,1])[c(1:2,4,6:11)],cex=0.8)

dev.copy(pdf,"fitting_points_1.pdf")
dev.off()
dev.off()


##
i=2
plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration(nM)",ylab="normalized SI",main="Average Points")
con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
for (j in 1:11){
if (j!=5){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]

t1[(t1[,5]==0)|(t1[,2]!=0),4]=NA
matrix=matrix(t1[,4],nrow=3,byrow=T)
y=colMeans(matrix,na.rm=T)

points(con,y,col=j,cex=0.5,pch=j)

}
}

legend("topright",col=c(1:4,6:11),pch=c(1:4,6:11),paste(names,ID[,1])[c(1:4,6:11)],cex=0.8)


dev.copy(pdf,"fitting_points (II).pdf")
dev.off()
dev.off()








######################### log scale (curve) #######################


i=1


plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration(nM)",ylab="normalized SI",log="xy",main="fitting curve I (log scale)")

for (j in 1:11){
if ((j!=5)&(j!=3)){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]
matrix(t1[,5],nrow=3,byrow=T)
y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##points(exp(x0),exp(y0),col=j,cex=0.5,pch=j)
a=coef(reg2)[1]
d=coef(reg2)[2]
curve(exp(llog(log(x),a,d)),from=0.5,to=5000,col=j,add=T,log="x")
}
}

legend("bottomright",col=c(1:2,4,6:11),lty=1,paste(names,ID[,1])[c(1:2,4,6:11)],cex=0.8)


dev.copy(pdf,"fitting_curve_log_1.pdf")
dev.off()
dev.off()



plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",log="xy",main="Fitting Curve II (log scale)",xlab="Concentration",yalb="normalized SI")

i=2
for (j in 1:11){
if (j!=5){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]
matrix(t1[,5],nrow=3,byrow=T)
y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##points(exp(x0),exp(y0),col=j,cex=0.5,pch=j)
a=coef(reg2)[1]
d=coef(reg2)[2]
curve(exp(llog(log(x),a,d)),from=0.5,to=5000,col=j,add=T,log="x")
}
}


legend("bottomright",col=c(1:2,4:11),lty=1,paste(names,ID[,2])[c(1:2,4:11)],cex=0.8)


dev.copy(pdf,"fitting_curve_log_2.pdf")
dev.off()
dev.off()







######################### log scale (average points) #######################
i=1


plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",xlab="Concentration(nM)",log="xy",ylab="normalized SI",main="average points I (log scale)")
con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
for (j in 1:11){
if ((j!=5)&(j!=3)){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]

t1[(t1[,5]==0)|(t1[,2]!=0),4]=NA
matrix=matrix(t1[,4],nrow=3,byrow=T)
y=colMeans(matrix,na.rm=T)

points(con,y,col=j,cex=0.5,pch=j)

}
}

legend("bottomright",col=c(1:2,4,6:11),pch=c(1:2,4,6:11),paste(names,ID[,1])[c(1:2,4,6:11)],cex=0.8)

dev.copy(pdf,"fitting_points_log_1.pdf")
dev.off()
dev.off()




######points

i=2
plot(1,1,xlim=c(0.5,5000),ylim=c(10^2,2*10^8),type="n",log="xy",xlab="Concentration (nM)",ylab="normalized SI",main="average points II (log scale)")
con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
for (j in 1:11){
if (j!=5){
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,39)]

t1[(t1[,5]==0)|(t1[,2]!=0),4]=NA
matrix=matrix(t1[,4],nrow=3,byrow=T)
y=colMeans(matrix,na.rm=T)

points(con,y,col=j,cex=0.5,pch=j)

}
}

legend("bottomright",col=c(1:4,6:11),pch=c(1:4,6:11),paste(names,ID[,1])[c(1:4,6:11)],cex=0.8)


dev.copy(pdf,"fitting_points_log_2.pdf")
dev.off()
dev.off()






