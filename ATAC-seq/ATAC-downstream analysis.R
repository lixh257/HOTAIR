#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("DiffBind")
library(DiffBind)
####制作一个存放数据位置信息的表格

###帮助文档信息
> tamoxifen <- dba(sampleSheet="tamoxifen.csv") %>%
+ dba.blacklist() %>%
+ dba.count() %>%
+ dba.normalize() %>%
+ dba.contrast() %>%
+ dba.analyze()



setwd("my/")
dbObj <- dba(sampleSheet="SampleSheet1.csv")
DBA_object <- dba.count(dbObj)###构建亲和矩阵
DBA_object <- dba.normalize(DBA_object)
consensusObj <- dba.peakset(DBA_object, bRetrieve=TRUE)
res=as.data.frame(consensusObj)
write.table(res,"DBA_object_all_peak_all_sample_new.txt",sep="\t",quote=F)


##wilcox.test
setwd("F:/HOTAIR/HOTAIR-ALL/ATAC-seq/2021-11-10整理/4-motif")
coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
allpeak=read.table("DBA_object_all_peak_all_sample_new.txt",sep="\t",header=T,row.names=1)
##筛选的阈值是：
1.在OE vs Ctrl 里面是显著的，在OE KO vs OE里面得是显著的
###2.然后FC在前两组相反
###3.秩和检验
OE=rownames(coldata)[which(coldata[,1]=="HOTAIR-OE")]
OE_peak=allpeak[,OE]
OE_KO=rownames(coldata)[which(coldata[,1]=="HOTAIR-OE YTHDF3-KO")]
OE_KO_peak=allpeak[,OE_KO]
Ctrl=rownames(coldata)[which(coldata[,1]=="Ctrl")]
Ctrl_peak=allpeak[,Ctrl]
KO=rownames(coldata)[which(coldata[,1]=="YTHDF3-KO")]
KO_peak=allpeak[,KO]

res=matrix(0,ncol=3,nrow=length(allpeak[,1]))
for(i in 1:length(allpeak[,1])){
  p_OE=wilcox.test(as.numeric(OE_peak[i,]),as.numeric(Ctrl_peak[i,]))
  PP_OE=p_OE$p.value
  p_OE_KO=wilcox.test(as.numeric(OE_KO_peak[i,]),as.numeric(OE_peak[i,]))
  PP_OE_KO=p_OE_KO$p.value
  p_KO=wilcox.test(as.numeric(KO_peak[i,]),as.numeric(Ctrl_peak[i,]))
  PP_KO=p_KO$p.value
  res[i,]=c(PP_OE,PP_OE_KO,PP_KO)
}
res_new=res
res_new[,1]=p.adjust(res[,1],method="fdr",n=length(res[,1]))
res_new[,2]=p.adjust(res[,2],method="fdr",n=length(res[,1]))
res_new[,3]=p.adjust(res[,3],method="fdr",n=length(res[,1]))
RES=res_new[intersect(which(res_new[,1]<0.05),which(res_new[,2]<0.05)),]
DEG=intersect(which(res[,1]<0.05),which(res[,2]<0.05))
RES1=allpeak[intersect(which(res[,1]<0.05),which(res[,2]<0.05)),]
Pvalue_RES=res[intersect(which(res[,1]<0.05),which(res[,2]<0.05)),]##P值阈值

##计算FC
OE_mean=apply(OE_peak,1,mean)
OE_KO_mean=apply(OE_KO_peak,1,mean)
Ctrl_mean=apply(Ctrl_peak,1,mean)
KO_mean=apply(KO_peak,1,mean)
FC=matrix(0,ncol=3,nrow=length(allpeak[,1]))
for(i in 1:length(allpeak[,1])){
  FC_OE_ctrl=OE_mean[i]/Ctrl_mean[i]
  FC_OEKO_OE=OE_KO_mean[i]/OE_mean[i]
  FC_KO_Ctrl=KO_mean[i]/Ctrl_mean[i]
  FC[i,]=c(FC_OE_ctrl,FC_OEKO_OE,FC_KO_Ctrl)
}
FC_res=FC[DEG,]
FC_res1=FC_res[which((FC_res[,1]<0.5&FC_res[,2]>2)|(FC_res[,1]>2&FC_res[,2]<0.5)),]##FC阈值
RES2=RES1[which((FC_res[,1]<0.5&FC_res[,2]>2)|(FC_res[,1]>2&FC_res[,2]<0.5)),]
Pvalue_RES1=Pvalue_RES[which((FC_res[,1]<0.5&FC_res[,2]>2)|(FC_res[,1]>2&FC_res[,2]<0.5)),]
###没考虑第三组的差异性不显著这个条件，因为想要保留更多的peak，而且这个不是必要条件
FC_res1[which(FC_res1[,3]==0),3]=1
colnames(FC_res1)=c("OE_vs_Ctrl","OE-KO_vs_OE","KO_vs_Ctrl")
write.table(FC_res1,"585peaks_wilcox_test_p_0.05_FC.txt",sep="\t",quote=F)
write.table(RES2,"585peaks_wilcox_test_p_0.05_peak.txt",sep="\t",quote=F)
write.table(Pvalue_RES1,"585peaks_wilcox_test_p_0.05_Pvalue.txt",sep="\t",quote=F)

#########Plots

FC=read.table("585peaks_wilcox_test_p_0.05_FC.txt",sep="\t",header=T,row.names=1)
peak=read.table("585peaks_wilcox_test_p_0.05_peak.txt",sep="\t",header=T,row.names=1)
Pvalue=read.table("585peaks_wilcox_test_p_0.05_Pvalue.txt",sep="\t",header=T,row.names=1)
YTH=intersect(which(FC[,8]>2|FC[,8]<0.5),which(Pvalue[,3]<0.05))
peak_DEG1=peak
FC_DEG1=FC
clu=rep(0,585)
clu[which(FC[,6]>2)]="Up"
clu[which(FC[,6]<0.5)]="Down"
clu[YTH]="YTHDF3 DEGs"
FC_res=log2(FC_DEG1[,6:8])
FC_res1=cbind(FC_DEG1[,1:5],FC_res,clu)
write.table(FC_res1,"585peak_FC_label_result.txt",sep="\t",row.names = T,col.names = T,quote=F)

peak=read.table("585peaks_wilcox_test_p_0.05_peak.txt",sep="\t",header=T,row.names=1)
peak_res2=c()
for(i in 1:length(FC[,1])){
  peak_res1=peak[which(rownames(peak)==rownames(FC)[i]),6:26]
  peak_res2=rbind(peak_res2,peak_res1)
}
RES=c()
for(i in 1:length(FC[,1])){
  peak_res1=peak[which(rownames(peak)==rownames(FC)[i]),]
  RES=rbind(RES,peak_res1)
}
write.table(RES,"585peak_order_peak_result.txt",sep="\t",quote=F)

####Fig 5 B
coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
annotation_col=data.frame(Condition=factor(coldata[,1]))
annotation_row=data.frame(regulate=factor(FC[,9]))
annotation_colors =list(Condition=c("HOTAIR-OE YTHDF3-KO"="Tan1","HOTAIR-OE"="Salmon",
                                    "YTHDF3-KO"="SkyBlue","Ctrl"="PaleGreen1"),
                        regulate=c("Up"="MistyRose","Down"="Honeydew","YTHDF3 DEGs"="LemonChiffon"))
rownames(peak_res2)=rownames(annotation_row)
colnames(peak_res2)=rownames(annotation_col)
pheatmap(peak_res2, show_colnames= T, show_rownames= F, scale= "row", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=F,
         annotation_col= annotation_col, 
         annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         #cutree_row = 2,
         labels_col=rownames(coldata),width=5,height=5,
         fontsize_row=8,fontsize_col=8,
         cutree_col=4
         
)


##############################Fig 5D
##########
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' 585peaks_OE-Down-92peak.narrowPeak > 585peaks_OE-Down-92peak.homer_peaks.tmp
/md01/lixh/software/Homer/bin/findMotifsGenome.pl 585peaks_OE-Down-92peak.homer_peaks.tmp hg19 585peaks_OE-Down-92peak_motifDir -len 8,10,12
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' 585peaks_OE-UP-482peak.narrowPeak > 585peaks_OE-UP-482peak.homer_peaks.tmp
/md01/lixh/software/Homer/bin/findMotifsGenome.pl 585peaks_OE-UP-482peak.homer_peaks.tmp hg19 585peaks_OE-UP-482peak_motifDir -len 8,10,12
###########这些是linux,得到结果继续下游分析


peak1=read.table("585peak_OE_UP_482peak_knownResults.txt",sep="\t",header=T)
peak1=peak1[which(peak1[,4]<0.01),]
peak2=read.table("585peak_OE_Down_92peak_knownResults.txt",sep="\t",header=T)
peak2=peak2[which(peak2[,4]<0.01),]

res_peak1=peak1[,c(1,4,8,10)]
res_peak11=cbind(rep("UP-482peak",length(res_peak1[,1])),res_peak1,
                 (as.numeric(sub("%", "", res_peak1[,3]))/100)/(as.numeric(sub("%", "", res_peak1[,4]))/100))
res_peak2=peak2[,c(1,4,8,10)]
res_peak22=cbind(rep("Down-96peak",length(res_peak2[,1])),res_peak2,
                 (as.numeric(sub("%", "", res_peak2[,3]))/100)/(as.numeric(sub("%", "", res_peak2[,4]))/100))

RES_585_up_down=matrix(0,ncol=6,nrow=(length(res_peak11[,1])+length(res_peak22[,1])))
colnames(RES_585_up_down)=c("Type","Name","Pvalue","Target","Background","Rate")
RES_585_up_down[1:length(res_peak11[,1]),]=as.matrix(res_peak11)
RES_585_up_down[(length(res_peak11[,1])+1):(length(res_peak22[,1])+length(res_peak11[,1])),]=as.matrix(res_peak22)

RES_585_up_down=as.data.frame(RES_585_up_down)
write.table(RES_585_up_down,"RES_574_up_down_dependent.txt",sep="\t",quote=F)
RES=matrix(0,ncol=6,nrow=21)
colnames(RES)=c("Type","Name","Pvalue","Target","Background","Rate")
RES[1:10,]=as.matrix(res_peak11[1:10,])
RES[11:21,]=as.matrix(res_peak22[1:11,])

RES=as.data.frame(RES)
#RES=RES[order(RES[,3]),]
RES$Name=factor(RES$Name,levels=c("TEAD3(TEA)","TEAD(TEA)","Bach2(bZIP)","Jun-AP1(bZIP)","Fosl2(bZIP)",
                                  "Fra2(bZIP)","AP-1(bZIP)","JunB(bZIP)","BATF(bZIP)","Fra1(bZIP)",
                                  "Atf3(bZIP)","Fos(bZIP)" ))
#library(ggplot2)
ggplot(RES,aes(x=Type,y=Name,col=as.numeric(Rate),size=-log10(as.numeric(Pvalue))))+
  geom_point()+theme_bw()+
  scale_color_gradient(low="lightblue", high="darkblue")+
  scale_size_area(max_size =11)
  #################################################################

##################################################################Fig 5E Down
peak=read.table("6gene_ATAC_peak.txt",sep="\t",header=T,row.names=1)
coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
annotation_col=data.frame(Condition=factor(coldata[,1]))
#annotation_row=data.frame(log2FC=factor(updown))
annotation_colors =list(Condition=c("HOTAIR-OE YTHDF3-KO"="Plum2","HOTAIR-OE"="DarkOliveGreen1",
                                    "YTHDF3-KO"="SandyBrown","Ctrl"="Khaki2"))

colnames(peak)=rownames(annotation_col)
pheatmap(peak, show_colnames= T, show_rownames= T, scale= "row", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=F,
         annotation_col= annotation_col, 
         #annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         #cutree_row = 2,
         labels_col=rownames(coldata),width=5,height=5,
         fontsize_row=8,fontsize_col=8,
         cutree_col=4
         
)
#############################################################################################

###########################################################################################Sup Fig 5D,5F
setwd("F:/HOTAIR/HOTAIR-ALL/ATAC-seq/2021-11-10整理/")
coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
allpeak=read.table("DBA_object_all_peak_all_sample_new.txt",sep="\t",header=T,row.names=1)
OE=rownames(coldata)[which(coldata[,1]=="HOTAIR-OE")]
OE_peak=allpeak[,OE]
OE_KO=rownames(coldata)[which(coldata[,1]=="HOTAIR-OE YTHDF3-KO")]
OE_KO_peak=allpeak[,OE_KO]
Ctrl=rownames(coldata)[which(coldata[,1]=="Ctrl")]
Ctrl_peak=allpeak[,Ctrl]
KO=rownames(coldata)[which(coldata[,1]=="YTHDF3-KO")]
KO_peak=allpeak[,KO]

res=matrix(0,ncol=3,nrow=length(allpeak[,1]))
for(i in 1:length(allpeak[,1])){
  p_OE=wilcox.test(as.numeric(OE_peak[i,]),as.numeric(Ctrl_peak[i,]))
  PP_OE=p_OE$p.value
  p_OE_KO=wilcox.test(as.numeric(OE_KO_peak[i,]),as.numeric(OE_peak[i,]))
  PP_OE_KO=p_OE_KO$p.value
  p_KO=wilcox.test(as.numeric(KO_peak[i,]),as.numeric(Ctrl_peak[i,]))
  PP_KO=p_KO$p.value
  res[i,]=c(PP_OE,PP_OE_KO,PP_KO)
}
res_new=res
RES1=res_new[which(res_new[,1]<0.05),]
RES2=res_new[which(res_new[,2]<0.05),]
RES3=res_new[which(res_new[,3]<0.05),]

##计算FC
OE_mean=apply(OE_peak,1,mean)
OE_KO_mean=apply(OE_KO_peak,1,mean)
Ctrl_mean=apply(Ctrl_peak,1,mean)
KO_mean=apply(KO_peak,1,mean)
FC=matrix(0,ncol=3,nrow=length(allpeak[,1]))
for(i in 1:length(allpeak[,1])){
  FC_OE_ctrl=OE_mean[i]/Ctrl_mean[i]
  FC_OEKO_OE=OE_KO_mean[i]/OE_mean[i]
  FC_KO_Ctrl=KO_mean[i]/Ctrl_mean[i]
  FC[i,]=c(FC_OE_ctrl,FC_OEKO_OE,FC_KO_Ctrl)
}
###提取FC 
FC_res1=FC[which(res_new[,1]<0.05),]
FC_res11=FC_res1[which(FC_res1[,1]<0.5|FC_res1[,1]>2),]
###提取peak set
peak_RES1=allpeak[which(res_new[,1]<0.05),]
peak_RES11=peak_RES1[which(FC_res1[,1]<0.5|FC_res1[,1]>2),]

###输出FC和peak set
colnames(FC_res11)=c("OE_vs_Ctrl","OE-KO_vs_OE","KO_vs_Ctrl")
write.table(FC_res11,"OE_vs_Ctrl_21420peaks_wilcox_p_0.05_FC.txt",sep="\t",quote=F)
write.table(peak_RES11,"OE_vs_Ctrl_21420peaks_wilcox_p_0.05_peakset.txt",sep="\t",quote=F)

peak1=read.table("OE_vs_Ctrl_21420peaks_wilcox_p_0.05_peakset.txt",sep="\t",header=T,row.names=1)
peak1_res=peak1[,6:26]
coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
annotation_col=data.frame(Condition=factor(coldata[,1]))
#annotation_row=data.frame(log2FC=factor(updown))
annotation_colors =list(Condition=c("HOTAIR-OE YTHDF3-KO"="Tan1","HOTAIR-OE"="Salmon",
                                    "YTHDF3-KO"="SkyBlue","Ctrl"="PaleGreen1"))

colnames(peak1_res)=rownames(annotation_col)
list=pheatmap(peak1_res, show_colnames= T, show_rownames= F, scale= "row", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_col= annotation_col, 
         #annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         cutree_row = 2,
         labels_col=rownames(coldata),width=5,height=5,
         fontsize_row=8,fontsize_col=8,
         cutree_col=4
         
)

row_cluster=cutree(list$tree_row,k=2)
newOrder=peak1[list$tree_row$order,]
#newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder),names(row_cluster))]
#colnames(newOrder)[ncol(newOrder)+1]="Cluster"
write.table(newOrder,"row_cluster_order_Peak_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")

###label出YTHDF3/Ctrl差异的peak
FC=read.table("585peaks_wilcox_test_p_0.05_FC.txt",sep="\t",header=T,row.names=1)
peak=read.table("585peaks_wilcox_test_p_0.05_peak.txt",sep="\t",header=T,row.names=1)
Pvalue=read.table("585peaks_wilcox_test_p_0.05_Pvalue.txt",sep="\t",header=T,row.names=1)
YTH=intersect(which(FC[,8]>2|FC[,8]<0.5),which(Pvalue[,3]<0.05))
peak_DEG1=peak
FC_DEG1=FC
clu=rep(0,585)
clu[which(FC[,6]>2)]="Up"
clu[which(FC[,6]<0.5)]="Down"
clu[YTH]="YTHDF3 DEGs"
FC_res=log2(FC_DEG1[,6:8])
FC_res1=cbind(FC_DEG1[,1:5],FC_res,clu)
write.table(FC_res1,"585peak_FC_label_result.txt",sep="\t",row.names = T,col.names = T,quote=F)


FC=read.table("585peak_FC_label_result.txt",sep="\t",header=T,row.names = 1)
FC_res=FC[,6:8]
annotation_row=data.frame(regulate=factor(FC[,9]))
#annotation_row=data.frame(log2FC=factor(updown))
annotation_colors =list(Condition=c("OE.DF3-KO"="Tan1","OE.DF3"="Salmon",
                                    "Ctrl.DF3-KO"="SkyBlue","Ctrl.DF3"="PaleGreen1"),
                        regulate=c("Up"="MistyRose","Down"="Honeydew","YTHDF3 DEGs"="LemonChiffon"))
rownames(FC_res)=rownames(annotation_row)
pheatmap(FC_res, show_colnames= T, show_rownames= F, scale= "none", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=F,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("SteelBlue", "white", "Firebrick3"))(20),
         #cutree_row = 2,
         #width=5,height=5,
         #labels_row = gene8[,1],
         annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         fontsize_row=8,fontsize_col=8,
         #cutree_col=4
         #filename="695gene_616trans_heatmap.pdf"
)##Sup Fig 5E

peak=read.table("585peaks_wilcox_test_p_0.05_peak.txt",sep="\t",header=T,row.names=1)
peak_res2=c()
for(i in 1:length(FC[,1])){
  peak_res1=peak[which(rownames(peak)==rownames(FC)[i]),6:26]
  peak_res2=rbind(peak_res2,peak_res1)
}
RES=c()
for(i in 1:length(FC[,1])){
  peak_res1=peak[which(rownames(peak)==rownames(FC)[i]),]
  RES=rbind(RES,peak_res1)
}
write.table(RES,"585peak_order_peak_result.txt",sep="\t",quote=F)

coldata=read.table("coldata_ATAC.txt",sep="\t",header=T,row.names=1)
annotation_col=data.frame(Condition=factor(coldata[,1]))
annotation_row=data.frame(regulate=factor(FC[,9]))
annotation_colors =list(Condition=c("HOTAIR-OE YTHDF3-KO"="Tan1","HOTAIR-OE"="Salmon",
                                    "YTHDF3-KO"="SkyBlue","Ctrl"="PaleGreen1"),
                        regulate=c("Up"="MistyRose","Down"="Honeydew","YTHDF3 DEGs"="LemonChiffon"))

rownames(peak_res2)=rownames(annotation_row)
colnames(peak_res2)=rownames(annotation_col)
pheatmap(peak_res2, show_colnames= T, show_rownames= F, scale= "row", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=F,
         annotation_col= annotation_col, 
         annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         #cutree_row = 2,
         labels_col=rownames(coldata),width=5,height=5,
         fontsize_row=8,fontsize_col=8,
         cutree_col=4
         
)###Sup Fig 5D
