
source("function_miao.R")

######################################################################################
######################## CBD Array####################################################
t0=read.table("CBD_data tab delimited.txt",header=T,sep="\t")
t1=read.table("CBD_background tab delimited.txt",header=T, sep="\t")

for (i in 1:3){

##choose the array with replicate i
replicate_fore=t0[,8]
replicate_back=t1[,10]

t=t0[replicate_fore==i,]
t_b=t1[replicate_back==i,]

##select out the signal intensity from multiple PMT Voltages and match the names
SI_fore=t[,c(14,12,10)]
SI_back=t_b[,c(9,7,5)]

sat_fore=t[,c(13,11,9)]
sat_back=t_b[,c(8,6,4)]

names(SI_back)=names(SI_fore) 
names(sat_back)=names(sat_fore) 


##combine the foreground and background SI from the same array
SI=rbind(SI_fore,SI_back)
sat=rbind(sat_fore,sat_back)


##exclude the saturation points (>=1%) and calibrate the SI
preset=as.matrix(SI)/10^6
sat=(sat>=1)
scan=multiPMT(preset,sat)


##estimate both background and foreground SI and overwrite the corresponding columns
t0[replicate_fore==i,25]=scan$estimate[1:(144*12),]*10^6
t1[replicate_back==i,11]=scan$estimate[(144*12+1):(144*12+20*12),]*10^6


##estimate the background mean and background subtracted SI; overwrite the corresponding columns
background=tapply(t1[replicate_back==i,11],t_b[,1],mean)
t0[replicate_fore==i,26]=rep(background,rep(144,12))
t0[replicate_fore==i,27]=t0[replicate_fore==i,25]-t0[replicate_fore==i,26]

##names of the new columns
colnames(t0)[25:27]=c("Estimate_foreground","Estimate_background","Background_subtracted_SI")
colnames(t1)[11]=c("Estimate_SI")

}


write.table(t0,"CBD_data.txt",row.names=F)
write.table(t1,"CBD_background.txt",row.names=F)






#####################################################################################
######################## Abl1 Array##################################################
t0=read.table("Abl1_data tab.txt",header=T,sep="\t")
t1=read.table("Abl1_background.txt",header=T)

for (i in 1:3){

##choose the array with replicate i
replicate_fore=t0[,9]
replicate_back=t1[,14]

t=t0[replicate_fore==i,]
t_b=t1[replicate_back==i,]

##select out the signal intensity from multiple PMTs and match the names
SI_fore=t[,c(23,21,19,17,15)]
SI_back=t_b[,c(13,11,9,7,5)]

sat_fore=t[,c(22,20,18,16,14)]
sat_back=t_b[,c(12,10,8,6,4)]

names(SI_back)=names(SI_fore) 
names(sat_back)=names(sat_fore) 


##combine the foreground and background SI from the same array
SI=rbind(SI_fore,SI_back)
sat=rbind(sat_fore,sat_back)


##exclude the saturation points (>=1%) and calibrate the SI
preset=as.matrix(SI)/10^6
sat=(sat>=1)
scan=multiPMT(preset,sat)


##estimate both background and foreground SI; overwrite the columns
t0[replicate_fore==i,40]=scan$estimate[1:(144*12),]*10^6
t1[replicate_back==i,15]=scan$estimate[(144*12+1):(144*12+20*12),]*10^6


##overwrite the background mean and background subtracted SI
background=tapply(t1[replicate_back==i,15],t_b[,1],mean)
t0[replicate_fore==i,41]=rep(background,rep(144,12))
t0[replicate_fore==i,42]=t0[replicate_fore==i,40]-t0[replicate_fore==i,41]

##names of added columns
colnames(t0)[40:42]=c("Estimate_foreground","Estimate_background","Background_subtracted_SI")
colnames(t1)[15]=c("Estimate_SI")

}


write.table(t0,"Abl1_data.txt",row.names=F)
write.table(t1,"Abl1_background.txt",row.names=F)






##################################################################################
######################## GRB2 Array###############################################
t0=read.table("GRB2_data tab.txt",header=T,sep="\t")
t1=read.table("GRB2_background.txt",header=T)

for (i in 1:3){

##choose the array with replicate i
replicate_fore=t0[,9]
replicate_back=t1[,14]

t=t0[replicate_fore==i,]
t_b=t1[replicate_back==i,]

##select out the signal intensity from multiple PMT and match the names
SI_fore=t[,c(23,21,19,17,15)]
SI_back=t_b[,c(13,11,9,7,5)]

sat_fore=t[,c(22,20,18,16,14)]
sat_back=t_b[,c(12,10,8,6,4)]

names(SI_back)=names(SI_fore) 
names(sat_back)=names(sat_fore) 


##combine the foreground and background SI from the same array
SI=rbind(SI_fore,SI_back)
sat=rbind(sat_fore,sat_back)

    
##exclude the saturation points (>=1%) and calibrate the SI
preset=as.matrix(SI)/10^6
sat=(sat>=1)
scan=multiPMT(preset,sat)


##estimate both background and foreground SI
t0[replicate_fore==i,40]=scan$estimate[1:(144*12),]*10^6
t1[replicate_back==i,15]=scan$estimate[(144*12+1):(144*12+20*12),]*10^6


##overwrite the background mean and background subtracted SI
background=tapply(t1[replicate_back==i,15],t_b[,1],mean)
t0[replicate_fore==i,41]=rep(background,rep(144,12))
t0[replicate_fore==i,42]=t0[replicate_fore==i,40]-t0[replicate_fore==i,41]

##names of added columns
colnames(t0)[40:42]=c("Estimate_foreground","Estimate_background","Background_subtracted_SI")
colnames(t1)[15]=c("Estimate_SI")

}


write.table(t0,"GRB2_data.txt",row.names=F)
write.table(t1,"GRB2_background.txt",row.names=F)







##################################################################################
######################## NCK1 Array###############################################
t0=read.table("NCK1_data tab.txt",header=T,sep="\t")
t1=read.table("NCK1_background.txt",header=T)

for (i in 1:3){

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


##exclude the saturation points (>=1%) and calibrate the SI
preset=as.matrix(SI)/10^6
sat=(sat>=1)
scan=multiPMT(preset,sat)


##estimate both background and foreground SI
t0[replicate_fore==i,45]=scan$estimate[1:(144*12),]*10^6
t1[replicate_back==i,17]=scan$estimate[(144*12+1):(144*12+20*12),]*10^6


##overwrite the background mean and background subtracted SI
background=tapply(t1[replicate_back==i,17],t_b[,1],mean)
t0[replicate_fore==i,46]=rep(background,rep(144,12))
t0[replicate_fore==i,47]=t0[replicate_fore==i,45]-t0[replicate_fore==i,46]

##names of added columns
colnames(t0)[45:47]=c("Estimate_foreground","Estimate_background","Background_subtracted_SI")
colnames(t1)[17]=c("Estimate_SI")

}


write.table(t0,"NCK1_data.txt",row.names=F)
write.table(t1,"NCK1_background.txt",row.names=F)








##################################################################################
########################SHP2N Array###############################################
t0=read.table("SHP2N_data tab.txt",header=T,sep="\t")
t1=read.table("SHP2N_background.txt",header=T)

for (i in 1:3){

##choose the array with replicate i
replicate_fore=t0[,9]
replicate_back=t1[,12]

t=t0[replicate_fore==i,]
t_b=t1[replicate_back==i,]

##select out the signal intensity from multiple PMT and match the names
SI_fore=t[,c(21,19,17,15)]
SI_back=t_b[,c(11,9,7,5)]

sat_fore=t[,c(20,18,16,14)]
sat_back=t_b[,c(10,8,6,4)]

names(SI_back)=names(SI_fore) 
names(sat_back)=names(sat_fore) 


##combine the foreground and background SI from the same array
SI=rbind(SI_fore,SI_back)
sat=rbind(sat_fore,sat_back)


##exclude the saturation points (>1%) and calibrate the SI
preset=as.matrix(SI)/10^6
sat=(sat>=1)
scan=multiPMT(preset,sat)


##estimate both background and foreground SI
t0[replicate_fore==i,35]=scan$estimate[1:(144*12),]*10^6
t1[replicate_back==i,13]=scan$estimate[(144*12+1):(144*12+20*12),]*10^6


##overwrite the background mean and background subtracted SI
background=tapply(t1[replicate_back==i,13],t_b[,1],mean)
t0[replicate_fore==i,36]=rep(background,rep(144,12))
t0[replicate_fore==i,37]=t0[replicate_fore==i,35]-t0[replicate_fore==i,36]

##names of added columns
colnames(t0)[35:37]=c("Estimate_foreground","Estimate_background","Background_subtracted_SI")
colnames(t1)[13]=c("Estimate_SI")

}


write.table(t0,"SHP2N_data.txt",row.names=F)
write.table(t1,"SHP2N_background.txt",row.names=F)



