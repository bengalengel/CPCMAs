## read files I can change filenames as needed.

VAVlow = read.table("VAV1 100 1 nM.txt",header=T,sep="\t")
VAVlowbg = read.table("VAV1 100 1 nM BG.txt",header=T,sep="\t")

## For each block calculate average of first row of background file. also calc s.d.

sub1=subset(VAVlowbg, Block==1 & Row==1)
sub2=subset(VAVlowbg, Block==2 & Row==1)
sub3=subset(VAVlowbg, Block==3 & Row==1)
sub4=subset(VAVlowbg, Block==4 & Row==1)
sub5=subset(VAVlowbg, Block==5 & Row==1)
sub6=subset(VAVlowbg, Block==6 & Row==1)
sub7=subset(VAVlowbg, Block==7 & Row==1)
sub8=subset(VAVlowbg, Block==8 & Row==1)
sub9=subset(VAVlowbg, Block==9 & Row==1)
sub10=subset(VAVlowbg, Block==10 & Row==1)
sub11=subset(VAVlowbg, Block==11 & Row==1)
sub12=subset(VAVlowbg, Block==12 & Row==1)



m1=mean(sub1[[35]])
sd1=sd(sub1[[35]])
m2=mean(sub2[[35]])
sd2=sd(sub2[[35]])
m3=mean(sub3[[35]])
sd3=sd(sub3[[35]])
m4=mean(sub4[[35]])
sd4=sd(sub4[[35]])
m5=mean(sub5[[35]])
sd5=sd(sub5[[35]])
m6=mean(sub6[[35]])
sd6=sd(sub6[[35]])
m7=mean(sub7[[35]])
sd7=sd(sub7[[35]])
m8=mean(sub8[[35]])
sd8=sd(sub8[[35]])
m9=mean(sub9[[35]])
sd9=sd(sub9[[35]])
m10=mean(sub10[[35]])
sd10=sd(sub10[[35]])
m11=mean(sub11[[35]])
sd11=sd(sub11[[35]])
m12=mean(sub12[[35]])
sd12=sd(sub12[[35]])






sub1norm=subset(VAVlow, Block==1, select = 35)
bgs1=sub1norm-m1
snr1=(sub1norm-m1/sd1)




