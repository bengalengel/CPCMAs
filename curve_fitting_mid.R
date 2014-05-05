##read in the data
x=read.table("Abl1_data.txt",header=T)
y=read.table("GRB2_data.txt",header=T)
z=read.table("NCK1_data.txt",header=T)
w=read.table("SHP2N_data.txt",header=T)


##11 peptides of interest (from excel file)
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



##Remark

##1. j=1,...,11 is the index in the peptide name list; i=1,2 is the index for corresponding id. I've just listed the curve fitting for the first id for all peptides. You can check the second id by changing i to 2.


##2. The main command for curve fitting is
##reg=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=1))
##which is a log-logistic curve fitting after log transformation. 
##This is almost equivalent to weighted least square:
##reg1=nls(exp(y0)~a/(1+exp(d-x0)),start=list(a=10^7,d=1),weight=1/exp(2*y0))

#Note that both the parameter d corresponds to log(kd). You should expect "exp(d)" from this approach is approximately the same as estimated "d" from following approach:
##reg=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))
##You can use the command above to estimate kd and sd for kd directly.##



##3. Any estimates from unconvergent fitting is invalid. Also you might need to tune the start value (a, b) in order to get convergent result. Also you should expect "exp(d)" from approach (1) is similar to "d" from approach(2)


##4. All fitting information, such as p-value, residual standard error, convergence or not, etc, will be printed on screen after typing in summary() command.

 


##########################################################


i=1
#####################j=1 RQLNpYIQVDLE#########################
j=1
arraykd=rep(0,4)
arraykdsd=rep(0,4)

##column 5 is concentration; column 7 is flag; column 13 is protein; column 44 is corrected signal intensity; column 39 is SNR filter

t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]

##Count the number of consecutively increasing concentrations scores above noise (either SNR irregular or SNR cirvular). 

con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
sort=sort.int(con,index=T,decreasing=F)$ix

##Abl1
counting=matrix((t1[,5]==1)&(t1[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

##GRB2
counting=matrix((t2[,5]==1)&(t2[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

##NCK1
counting=matrix((t3[,5]==1)&(t3[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

##SHP2N
counting=matrix((t4[,5]==1)&(t4[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)


y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
##print out all the fitting information
summary(reg1)


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
summary(reg2)

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
summary(reg3)

y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
summary(reg4)


arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]

write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)






#######################j=2 FESRpYQQPFED#################################
j=2
arraykd=rep(0,4)
arraykdsd=rep(0,4)

t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]



write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)








############################# j=3 PPALpYAEPLDS########################
j=3 

arraykd=rep(0,4)
arraykdsd=rep(0,4)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]


y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]

write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)







###########################################################

j=4
arraykd=rep(0,3)
arraykdsd=rep(0,3)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]

y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
##x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
##reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg4))[2,1]


arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg4))[2,2]





write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)







###########################################################

j=5

arraykd=rep(0,3)
arraykdsd=rep(0,3)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


##y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
##x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
##reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]


arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]




write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)







###########################################################

j=6

arraykd=rep(0,4)
arraykdsd=rep(0,4)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]


y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=1))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=1))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=1))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=1))




arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]



write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)




#############################################################


j=7

arraykd=rep(0,4)
arraykdsd=rep(0,4)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]

write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)




###########################################################

j=8

arraykd=rep(0,4)
arraykdsd=rep(0,4)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg3))[2,1]
arraykd[4]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg3))[2,2]
arraykdsd[4]=coef(summary(reg4))[2,2]

write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)



###########################################################

j=9

arraykd=rep(0,3)
arraykdsd=rep(0,3)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
##x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
##reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg4))[2,2]


write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)



#############################################################

j=10

arraykd=rep(0,3)
arraykdsd=rep(0,3)
t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]


y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
##x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
##reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg4))[2,2]



write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)




#############################################################

j=11

arraykd=rep(0,3)
arraykdsd=rep(0,3)

t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]



y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))



y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

##y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
##x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
##reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))


y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))

arraykd[1]=coef(summary(reg1))[2,1]
arraykd[2]=coef(summary(reg2))[2,1]
arraykd[3]=coef(summary(reg4))[2,1]

arraykdsd[1]=coef(summary(reg1))[2,2]
arraykdsd[2]=coef(summary(reg2))[2,2]
arraykdsd[3]=coef(summary(reg4))[2,2]


write.table(arraykd,"kd.txt",row.names=F,col.names=F,append=T)
write.table(arraykdsd,"kdsd.txt",row.names=F,col.names=F,append=T)

#############################



