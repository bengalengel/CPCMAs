source("function_miao.R")


###############################################################################
########################CBD Array###############################################
t0=read.table("CBD_data.txt",header=T)

flag=(t0[,6]==0)
replicate=t0[,8]
SI=t0[,27]
filter=rep(1,dim(t0)[1])

##conduct inter_array normalization
norm=inter_norm(filter,flag,SI,replicate,F)
t0[,28]=norm$inter_normalized
##write the normalized value into the corresponding column
colnames(t0)[28]=c("Inter_normalized SI")

##estimate the intra_array scales
block=t0[,1]
SI=t0[,28]
norm=block_norm(filter,flag,SI,block,F)
block_scale=norm$block_scale


##conduct intra_array normalization
norm=intra_norm(block,SI,block_scale)
t0[,29]=norm$intra_normalized
##write the normalized value into the corresponding column
colnames(t0)[29]=c("Intra_normalized SI")


write.table(t0,"CBD_data.txt",row.names=F)




#################################################################################
########################Abl1 Array###############################################
t0=read.table("Abl1_data.txt",header=T)

block=t0[,1]
filter=t0[,39]
flag=(t0[,7]==0)
replicate=t0[,9]
SI=t0[,42]


##conduct inter_array normalization
norm=inter_norm(filter,flag,SI,replicate,F)
t0[,43]=norm$inter_normalized
colnames(t0)[43]=c("Inter_normalized SI")


##read the intra_array scales from CBD array
SI=t0[,43]
block_scale=read.table("block_scale.txt",header=T)
block_scale=as.matrix(block_scale)
##conduct intra_array normalization
norm=intra_norm(block,SI,block_scale)
t0[,44]=norm$intra_normalized
colnames(t0)[44]=c("Intra_normalized SI")


write.table(t0,"Abl1_data.txt",row.names=F)








#################################################################################
########################GRB2 Array###############################################
t0=read.table("GRB2_data.txt",header=T)

block=t0[,1]
filter=t0[,39]
flag=(t0[,7]==0)
replicate=t0[,9]
SI=t0[,42]


##conduct inter_array normalization
norm=inter_norm(filter,flag,SI,replicate,F)
t0[,43]=norm$inter_normalized
colnames(t0)[43]=c("Inter_normalized SI")


##conduct intra_array normalization
SI=t0[,43]
block_scale=read.table("block_scale.txt",header=T)
block_scale=as.matrix(block_scale)


norm=intra_norm(block,SI,block_scale)
t0[,44]=norm$intra_normalized
colnames(t0)[44]=c("Intra_normalized SI")

##write the result
write.table(t0,"GRB2_data.txt",row.names=F)





#################################################################################
########################NCK1 Array###############################################
t0=read.table("NCK1_data.txt",header=T)

block=t0[,1]
filter=t0[,44]
flag=(t0[,7]==0)
replicate=t0[,9]
SI=t0[,47]


##conduct inter_array normalization
norm=inter_norm(filter,flag,SI,replicate,F)
t0[,48]=norm$inter_normalized
colnames(t0)[48]=c("Inter_normalized SI")


##conduct intra_array normalization
SI=t0[,48]
block_scale=read.table("block_scale.txt",header=T)
block_scale=as.matrix(block_scale)
norm=intra_norm(block,SI,block_scale)
t0[,49]=norm$intra_normalized
colnames(t0)[49]=c("Intra_normalized SI")


##write the result
write.table(t0,"NCK1_data.txt",row.names=F)





#################################################################################
########################SHP2N Array###############################################
t0=read.table("SHP2N_data.txt",header=T)


block=t0[,1]
filter=t0[,34]
flag=(t0[,7]==0)
replicate=t0[,9]
SI=t0[,37]

##conduct inter_array normalization
norm=inter_norm(filter,flag,SI,replicate,F)
t0[,38]=norm$inter_normalized
colnames(t0)[38]=c("Inter_normalized SI")

##conduct intra_array normalization
SI=t0[,38]
block_scale=read.table("block_scale.txt",header=T)
block_scale=as.matrix(block_scale)
norm=intra_norm(block,SI,block_scale)
t0[,39]=norm$intra_normalized
colnames(t0)[39]=c("Intra_normalized SI")

##write the result
write.table(t0,"SHP2N_data.txt",row.names=F)




