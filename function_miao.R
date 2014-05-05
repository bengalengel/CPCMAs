######install R libraries and packages
library(aroma.light)
library(vsn)
library(MASS)
library(drc)




#############################################################################
##Calibrate signal intensity from multiple PMT voltage settings##
##Input:
##preset: An N*K matrix (K>=2) where the columns represent the signal intensities (both foreground and background) of an array from multiple PMT voltages 
##saturation: An N*K matrix where the columns represent if the corresponding spot is saturated at that PMT voltage settings. 
##Output: 
##estimate: An N*1 matrix of estimated signal intensity; the default contrast is the first column in the input matrix preset.
##normalize: An N*k matrix where the columns represent the signal intensities after offset and scale correction
##offset_correct: An N*k matrix where the columns represent the signal intensities after only offset correction
############################################################################

multiPMT=function(preset,saturation){
n=dim(preset)[2]
##set the saturated points as missing data
for (i in 1:n){
preset[saturation[,i]==T,i]=NA
}
normalize=offset_correct=preset
##estimated signal intensity from multiple PMT voltages
estimate=calibrateMultiscan.matrix(preset)

    
for (i in 1:n){    
offset_correct[,i]=preset[,i]-attr(estimate,"modelFit")$a[i]
normalize[,i]=(preset[,i]-attr(estimate,"modelFit")$a[i])/attr(estimate,"modelFit")$b[i]
    offset_correct[saturation[,i]==T,i]=estimate[saturation[,i]==T]*attr(estimate,"modelFit")$b[i]
    normalize[saturation[,i]==T,i]=estimate[saturation[,i]==T]
}    
    

return(list(estimate=estimate,normalize=normalize,offset_correct=offset_correct))
}










#############################################################################
##Conduct inter-array normalization over three replicates##
##input:
##filter: An N*1 matrix where the columns represent SNR score (1 means SNR>=3,0 otherwise).
##flag: An N*1 matrix where the columns represent the flag index from raw data.
##SI: An N*1 matrix where the columns represent background_subtracted SI from all three replicates.
##replicate: An N*1 matrix indexing which replicate it's from.
##trans: F is using robust linear regression; T if using glog transformation.
##Output:
##inter_normalized: An N*1 matrix of signal intensity after inter_array normalization.
##scale: inter_normalization scale for each replicate. 
##A scale.txt file containing three inter_normalization scale is created. Be sure to rename it (in a protein-spcific way) to avoid being rewitten in subsequent analysis.
############################################################################


inter_norm=function(filter,flag,SI,replicate,trans){

##exclude invalid data in the calculation of inter-array scale
non_neg=(SI>=0)
valid=(filter)&(flag)&(non_neg)
result=SI
result[!valid]=NA


##group by replicates
inter_normalized=cbind(result[replicate==1],result[replicate==2],result[replicate==3])
index=which((is.na(inter_normalized[,1])==F)&(is.na(inter_normalized[,2])==F)&(is.na(inter_normalized[,3])==F))


preset=inter_normalized[index,]/10^5
    
if (trans==T){
##calculate the inter-array scale based on variance stablelizing transformation (glog)
model=vsn2(preset)
scale=exp(coef(model)[,,2])[1]/exp(coef(model)[,,2])
}else{
##calculate the inter-array scale based on robust linear regression
reg=rlm(log(c(preset))~as.factor(col(preset)))
scale=c(1,exp(coef(reg))[2:3])
}

write.table(scale,"scale.txt",row.names=F,col.names="inter_array scale")

##conduct inter-array normalization 
for (i in 1:3){
result[replicate==i]=SI[replicate==i]/scale[i]
}

return(list(inter_normalized=result,scale=scale))

}







#############################################################################
##Conduct intra-array normalization over 12 blocks from the same array##
##input:
##block: An N*1 matrix indexing the block.
##SI: An N*1 matrix where the columns represent inter_normalized SI. 
##block_scale: An 12*1 matrix containing intra-array normalization scale from CBD array.
##Output:
##intra_normalized: An N*1 matrix of signal intensity after intra_array normalization.
############################################################################



intra_norm=function(block,SI,block_scale){
    
result=SI

for (i in 1:12){
result[block==i]=SI[block==i]/block_scale[i]
}

return(list(intra_normalized=result))
}







#############################################################################
##Estimate the intra-array normalization scale from three CBD arrays. Inter_array normalization over three replicates should be conducted beforehand.  
##input:
##filter: An N*1 matrix where the columns represent SNR score from CBD arrays (1 means SNR>=3; 0 otherwise).
##flag: An N*1 matrix where the columns represent the flag index from CBD arrays.
##SI: An N*1 matrix where the columns represent inter_array normalized SI from CBD arrays.
##Block: An N*1 matrix indexing the block
##trans: F is using robust linear regression; T if using glog transformation.
##Output:
##block_scale: estimated intra_normalization scale from three CBD arrays. 
##A block_scale.txt file containing 12 intra_normalization scale is created. 
############################################################################


block_norm=function(filter,flag,SI,block,trans){
non_neg=(SI>=0)
valid=(filter)&(flag)&(non_neg)
result=SI
result[!valid]=NA


##group by blocks
index=rep(T,length(SI[block==1]))
block_normalized=NULL
for (i in 1:12){
block_normalized=cbind(block_normalized,result[block==i])
index=index&(is.na(block_normalized[,i])==F)
}
index=which(index)


preset=block_normalized[index,]/10^5

if (trans==T){
model=vsn2(preset)
scale=exp(coef(model)[,,2])[1]/exp(coef(model)[,,2])
}
else{
reg=rlm(log(c(preset))~as.factor(col(preset)))
scale=c(1,exp(coef(reg))[2:12])
}

write.table(scale,"block_scale.txt",row.names=F,col.names="block_scale")


return(list(block_scale=scale))
}




############################################################################
##log-logistic function with d=log(kd) as parameter
############################################################################
llog0=function(x,a,d){
    y=a/(1+exp(d-x))
    return(y)
}



############################################################################
##log-transformed log-logistic function with d=log(kd) as parameter
############################################################################
llog=function(x,a,d){
    y=a-log(1+exp(d-x))
    return(y)
}



