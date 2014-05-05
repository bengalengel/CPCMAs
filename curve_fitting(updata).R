source("function_miao.R")


x=read.table("Abl1_data.txt",header=T)
y=read.table("GRB2_data.txt",header=T)
z=read.table("NCK1_data.txt",header=T)
w=read.table("SHP2N_data.txt",header=T)


##11 peptides of interest (from excel file); add more if you want see all 144 peptides.
names=c("RQLNpYIQVDLE","FESRpYQQPFED","PPALpYAEPLDS","PSSVpYVPDEWE","NGLNpYIDLDLV",
"SPGEpYVNIEFG","KEEGpYELPYNP","DDPSpYVNVQNL","EDDGpYDVPKPP","DQHDpYDSVASD","AEPQpYEEIPIY")

##corresponding ID; change 11 to 144, 22 to 144*2 and more corresonding ID if you want to see all 144 peptides. 

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





##################### j=1 corresponding RQLNpYIQVDLE, i=1 corresponding ID="B02"#########################
i=1
j=1



t1=x[(x[,4]==names[j])&(x[,8]==ID[j,i]),c(5,7,13,44,46)]
t2=y[(y[,4]==names[j])&(y[,8]==ID[j,i]),c(5,7,13,44,46)]
t3=z[(z[,4]==names[j])&(z[,8]==ID[j,i]),c(5,7,13,49,51)]
t4=w[(w[,4]==names[j])&(w[,8]==ID[j,i]),c(5,7,13,39,41)]
 

con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
sort=sort.int(con,index=T,decreasing=F)$ix


##SHP2N
counting=matrix((t4[,5]==1)&(t4[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

y0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),4])
x0=log1p(t4[(t4[,5]==1)&(t4[,2]==0),1])
##parameterization: log(kd)
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
##parameterization: kd
##reg1=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))

summary(reg1)




##Abl1
counting=matrix((t1[,5]==1)&(t1[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

y0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),4])
x0=log1p(t1[(t1[,5]==1)&(t1[,2]==0),1])
reg2=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
##reg2=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))
summary(reg2)



##GRB2
counting=matrix((t2[,5]==1)&(t2[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

y0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),4])
x0=log1p(t2[(t2[,5]==1)&(t2[,2]==0),1])
reg3=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
##reg3=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))
summary(reg3)




##NCK1
counting=matrix((t3[,5]==1)&(t3[,2]==0),nrow=3,byrow=T)
apply(counting[,sort],2,sum)

y0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),4])
x0=log1p(t3[(t3[,5]==1)&(t3[,2]==0),1])
reg4=nls(y0~a-log(1+exp(d-x0)),start=list(a=7,d=4))
##reg4=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))
summary(reg4)

remove(reg1)
remove(reg2)
remove(reg3)
remove(reg4)

