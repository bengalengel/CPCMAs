source("function_miao.R")

x=read.table("Abl1_data.txt",header=T)
y=read.table("GRB2_data.txt",header=T)
z=read.table("NCK1_data.txt",header=T)
w=read.table("SHP2N_data.txt",header=T)


for (i in 2:143){

t1=x[(x[,4]==list[i,1])&(x[,8]==list[i,2]),c(4,5,7,8,13,44,45,46)]
t2=y[(y[,4]==list[i,1])&(y[,8]==list[i,2]),c(4,5,7,8,13,44,45,46)]
t3=z[(z[,4]==list[i,1])&(z[,8]==list[i,2]),c(4,5,7,8,13,49,50,51)]
t4=w[(w[,4]==list[i,1])&(w[,8]==list[i,2]),c(4,5,7,8,13,39,40,41)]

remove(reg1)
con=c(5000.0,100.0,2500.0,50.0,1000.0,10.0,750.0,5.0,500.0,1.0,250.0,0.5)
sort=sort.int(con,index=T,decreasing=F)$ix

##SHP2N
##snr vector circular
counting=matrix((t4[,7]==1)&(t4[,3]==0),nrow=3,byrow=T)
snr_vect1=t(apply(counting[,sort],2,sum))

##snr vector circular and irregular
counting2=matrix((t4[,8]==1)&(t4[,3]==0),nrow=3,byrow=T)
snr_vect2=t(apply(counting2[,sort],2,sum))

##Combine into table with identifier and snrvectors and write to the appropriate file
name=unique(t4[,1])
ID=unique(t4[,4])
protein=unique(t4[,5])

y0=log1p(t4[(t4[,8]==1)&(t4[,3]==0),6])
x0=log1p(t4[(t4[,8]==1)&(t4[,3]==0),2])

#######parameterization: kd?????
#######reg1=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))
####reg2=nls(y0~a-log(1+d*exp(-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit,weight=1/exp(2*y0))
## This form allows for weighted fitting on Y vs logX and estimation of the logged parameter value

##reg2=nls(exp(y0)~a*(1/(1+exp(d-x0))),start=list(a=max(exp(y0),na.rm=T),d=4),na.action=na.omit,weight=1/exp(2*y0))

##parameterization: log(kd)
##Error handling for Kd estimates

pos_err=tryCatch(nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit), error=function(e) e)

pos_err2=tryCatch(nls(exp(y0)~a*(1/(1+exp(d-x0))),start=list(a=max(exp(y0),na.rm=T),d=4),na.action=na.omit,weight=1/exp(2*y0)), error=function(e) e)


##If no error do fit for parameterization: log(kd)
if(!inherits(pos_err, "error")){
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit)

# coefficients in dataframe
regr_tablog <- data.frame(summary(reg1)$coefficients)

# get the p-vals 
regr_tablog[ ,4] <- ifelse(regr_tablog[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tablog[ ,4] < .01, "< 0.01", 
                               round(regr_tablog[ ,4], 3)))
##convert regr tables into one line
sub1=regr_tablog[1,1:4]
sub2=regr_tablog[2,1:4]
log_result=c(sub1,sub2)
}else{ 
log_result="no convergence log"
}

##If no error do fit for parameterization: linear(kd)
if(!inherits(pos_err2, "error")){
reg2=nls(exp(y0)~a*(1/(1+exp(d-x0))),start=list(a=max(exp(y0),na.rm=T),d=4),na.action=na.omit,weight=1/exp(2*y0))

# coefficients in dataframe
regr_tablinear <- data.frame(summary(reg2)$coefficients)

# get the p-vals 
regr_tablinear[ ,4] <- ifelse(regr_tablinear[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tablinear[ ,4] < .01, "< 0.01", 
                               round(regr_tablinear[ ,4], 3)))
##convert regr tables into one line
sub1=regr_tablinear[1,1:4]
sub2=regr_tablinear[2,1:4]
linear_result=c(sub1,sub2)
}else{ 
linear_result="no convergence linear"
}





out=data.frame(name,ID,protein,snr_vect1,snr_vect2,log_result,linear_result)

write.table(t(colnames(out)),"fits.csv",sep=',',col.names=F,row.names=F,append=T)
write.table(out,"fits.csv",sep=',',col.names=F,row.names=F,append=T)

#clear holders that may not be replaced upon iteration
remove(reg1)
remove(reg2)
remove(regr_tablog)
remove(regr_tablinear)
remove(pos_err)
remove(pos_err2)
remove(log_result)
remove(linear_result)



##Abl1

##snr vector circular
counting=matrix((t1[,7]==1)&(t1[,3]==0),nrow=3,byrow=T)
snr_vect1=t(apply(counting[,sort],2,sum))

##snr vector circular and irregular
counting2=matrix((t1[,8]==1)&(t1[,3]==0),nrow=3,byrow=T)
snr_vect2=t(apply(counting2[,sort],2,sum))

##Combine into table with identifier and snrvector
name=unique(t1[,1])
ID=unique(t1[,4])
protein=unique(t1[,5])

y0=log1p(t1[(t1[,8]==1)&(t1[,3]==0),6])
x0=log1p(t1[(t1[,8]==1)&(t1[,3]==0),2])
#######parameterization: kd
#######reg1=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))

##parameterization: log(kd)
##Error handling
pos_err=tryCatch(nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit), error=function(e) e)

##If no error do fit
if(!inherits(pos_err, "error")){
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit)

# coefficients in dataframe
regr_tab <- data.frame(summary(reg1)$coefficients)
# get the p-vals 
regr_tab[ ,4] <- ifelse(regr_tab[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tab[ ,4] < .01, "< 0.01", 
                               round(regr_tab[ ,4], 3)))

out=data.frame(name,ID,protein,snr_vect1,snr_vect2,regr_tab)
}else{ 
out=data.frame(name,ID,protein,snr_vect1,snr_vect2,"no convergence")
}
#write.table(t(colnames(out)),"fits.csv",sep=',',col.names=F,row.names=F,append=T)
write.table(out,"fits.csv",sep=',',col.names=F,row.names=F,append=T)

#clear holders that may not be replaced upon iteration
remove(reg1)
remove(regr_tab)



##GRB2

##snr vector circular
counting=matrix((t2[,7]==1)&(t2[,3]==0),nrow=3,byrow=T)
snr_vect1=t(apply(counting[,sort],2,sum))

##snr vector circular and irregular
counting2=matrix((t2[,8]==1)&(t2[,3]==0),nrow=3,byrow=T)
snr_vect2=t(apply(counting2[,sort],2,sum))

##Combine into table with identifier and snrvector
name=unique(t2[,1])
ID=unique(t2[,4])
protein=unique(t2[,5])

y0=log1p(t2[(t2[,8]==1)&(t2[,3]==0),6])
x0=log1p(t2[(t2[,8]==1)&(t2[,3]==0),2])
#######parameterization: kd
#######reg1=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))

##parameterization: log(kd)
##Error handling
pos_err=tryCatch(nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit), error=function(e) e)

##If no error do fit
if(!inherits(pos_err, "error")){
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit)

# coefficients in dataframe
regr_tab <- data.frame(summary(reg1)$coefficients)

# get the p-vals 
regr_tab[ ,4] <- ifelse(regr_tab[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tab[ ,4] < .01, "< 0.01", 
                               round(regr_tab[ ,4], 3)))


out=data.frame(name,ID,protein,snr_vect1,snr_vect2,regr_tab)
}else{ 
out=data.frame(name,ID,protein,snr_vect1,snr_vect2,"no convergence")
}
#write.table(t(colnames(out)),"fits.csv",sep=',',col.names=F,row.names=F,append=T)
write.table(out,"fits.csv",sep=',',col.names=F,row.names=F,append=T)

#clear holders that may not be replaced upon iteration
remove(reg1)
remove(regr_tab)



##NCK1

##snr vector circular
counting=matrix((t3[,7]==1)&(t3[,3]==0),nrow=3,byrow=T)
snr_vect1=t(apply(counting[,sort],2,sum))

##snr vector circular and irregular
counting2=matrix((t3[,8]==1)&(t3[,3]==0),nrow=3,byrow=T)
snr_vect2=t(apply(counting2[,sort],2,sum))

##Combine into table with identifier and snrvector
name=unique(t3[,1])
ID=unique(t3[,4])
protein=unique(t3[,5])

y0=log1p(t3[(t3[,8]==1)&(t3[,3]==0),6])
x0=log1p(t3[(t3[,8]==1)&(t3[,3]==0),2])
#######parameterization: kd
#######reg1=nls(y0~a-log(1+d*exp(-x0)),start=list(a=7,d=4))

##parameterization: log(kd)
##Error handling
pos_err=tryCatch(nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit), error=function(e) e)

##If no error do fit
if(!inherits(pos_err, "error")){
reg1=nls(y0~a-log(1+exp(d-x0)),start=list(a=max(y0,na.rm=T),d=4),na.action=na.omit)

# coefficients in dataframe
regr_tab <- data.frame(summary(reg1)$coefficients)

# get the p-vals 
regr_tab[ ,4] <- ifelse(regr_tab[ ,4] < .001, "< 0.001", 
                        ifelse(regr_tab[ ,4] < .01, "< 0.01", 
                               round(regr_tab[ ,4], 3)))

out=data.frame(name,ID,protein,snr_vect1,snr_vect2,regr_tab)
}else{ 
out=data.frame(name,ID,protein,snr_vect1,snr_vect2,"no convergence")
}
#write.table(t(colnames(out)),"fits.csv",sep=',',col.names=F,row.names=F,append=T)
write.table(out,"fits.csv",sep=',',col.names=F,row.names=F,append=T)

#clear holders that may not be replaced upon iteration
remove(reg1)
remove(regr_tab)
}



