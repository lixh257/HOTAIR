##Wilocx.test，(P<0.05,|log2FC|>1)
TPM=read.table("22sample_TPM.txt",sep="\t",header=T,row.names=1)#count
colData=read.table("coldata-subgroup.txt",sep="\t",header=T,row.names=1)##样本信息
Ctrl=TPM[,1:6]
Ctrl_mean=apply(Ctrl,1,mean)
OE=TPM[,7:12]
OE_mean=apply(OE,1,mean)
KO=TPM[,13:16]
KO_mean=apply(KO,1,mean)
OE_KO=TPM[,17:22]
OE_KO_mean=apply(OE_KO,1,mean)
name=rownames(TPM)
FC1=log2((OE_mean+0.001)/(Ctrl_mean+0.001))
FC2=log2((OE_KO_mean+0.001)/(OE_mean+0.001))
FC3=log2((KO_mean+0.001)/(Ctrl_mean+0.001))
res1=cbind(FC1,FC2,FC3)
colnames(res1)=c("OE_vs_Ctrl","OE_KO_vs_OE","KO_vs_Ctrl")
write.table(res1,"3_group_log2FC_result_all.txt",sep="\t",quote=F)

pvalue=matrix("NA",ncol=3,nrow=length(TPM[,1]))
rownames(pvalue)=rownames(TPM)
colnames(pvalue)=c("OE_vs_Ctrl","OE_KO_vs_OE","KO_vs_Ctrl")
for(i in 1:length(TPM[,1])){
  a=wilcox.test(as.numeric(OE[i,]),as.numeric(Ctrl[i,]))
  pvalue[i,1]=a$p.value
  b=wilcox.test(as.numeric(OE[i,]),as.numeric(OE_KO[i,]))
  pvalue[i,2]=b$p.value
  c=wilcox.test(as.numeric(KO[i,]),as.numeric(Ctrl[i,]))
  pvalue[i,3]=c$p.value
}
write.table(pvalue,"3_group_pvalue_result_all.txt",sep="\t",quote=F)

###Plots
###Only test HOTAIR-OE/Ctrl
########################################################################Sup Fig5C,E
setwd("F:/HOTAIR/HOTAIR-ALL/overlap-RNA-ATAC/2022-1-17-重新画图/")
FC=read.table("3_group_log2FC_result_all.txt",sep="\t",header=T,row.names = 1)
Pvalue=read.table("3_group_pvalue_result_all.txt",sep="\t",header=T,row.names = 1)
TPM=read.table("22sample_TPM.txt",sep="\t",header=T,row.names = 1)
DEG=intersect(which(abs(FC[,1])>1),which(Pvalue[,1]<0.05))
TPM_DEG=TPM[DEG,]
colData=read.table("coldata-subgroup.txt",sep="\t",header=T,row.names=1)##样本信息
colData<-data.frame(colData,paste(colData$HOTAIR,colData$YTHDF,sep="."))
colnames(colData)[3] <- "HOTAIR_YTHDF"
annotation_col=data.frame(Condition=factor(colData$HOTAIR_YTHDF))
annotation_colors =list(Condition=c("OE.DF3-KO"="Tan1",
                                    "OE.DF3"="Salmon","Ctrl.DF3-KO"="SkyBlue","Ctrl.DF3"="PaleGreen1"),log2FC=c("UP"="Salmon","Down"="LightSkyBlue"))
colnames(TPM_DEG)=rownames(annotation_col)
list=pheatmap(TPM_DEG, show_colnames= T, show_rownames= F, scale= "row", fontsize= 6.5,
              clustering_method ="complete",
              cluster_cols=F,cluster_rows=T,
              annotation_col= annotation_col, 
              #annotation_row= annotation_row,
              annotation_colors= annotation_colors,
              breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
              col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
              cutree_row = 2,
              labels_col=rownames(colData),width=5,height=5,
              fontsize_row=8,fontsize_col=8,
              #cutree_col=2
)
colnames(TPM_DEG)=colnames(TPM)
row_cluster=cutree(list$tree_row,k=2)
newOrder=TPM_DEG[list$tree_row$order,]
write.table(newOrder,"row_cluster_order_TPM_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")


FC_DEG=FC[DEG,]
colnames(FC_DEG)=rownames(annotation_col)
list=pheatmap(FC_DEG, show_colnames= T, show_rownames= F, scale= "row", fontsize= 6.5,
              clustering_method ="complete",
              cluster_cols=F,cluster_rows=T,
              annotation_col= annotation_col, 
              #annotation_row= annotation_row,
              annotation_colors= annotation_colors,
              breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
              col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
              cutree_row = 2,
              labels_col=rownames(colData),width=5,height=5,
              fontsize_row=8,fontsize_col=8,
              #cutree_col=2
)
colnames(FC_DEG_DEG)=colnames(FC)
row_cluster=cutree(list$tree_row,k=2)
newOrder=FC_DEG_DEG[list$tree_row$order,]
write.table(newOrder,"row_cluster_order_FC_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")

###############################################################Fig5 A
YTH=intersect(which(abs(FC[,3])>1),which(Pvalue[,3]<0.05))
YTHDEG=intersect(rownames(FC)[YTH],RNA[,1])
TPM_DEG1=TPM[RNA[,1],]
TPM_DEG1_mean=TPM_mean[RNA[,1],]
FC_DEG1=FC[RNA[,1],]
clu=rep(0,119)
clu[which(FC_DEG1[,1]>1)]="Up"
clu[which(FC_DEG1[,1]<0)]="Down"
clu[which(rownames(FC_DEG1)%in%YTHDEG)]="YTHDF3 DEGs"

FC_res=cbind(FC_DEG1,clu)
write.table(FC_res,"119_FC_result_group_label.txt",sep = "\t",quote=F)
###
FC_DEG2=read.table("119_FC_result_group_label.txt",sep = "\t",header=T,row.names = 1)
RNA=read.table("119gene_symbol_ID_result.txt",sep="\t",header=F,row.names = 1)
symbol=RNA[rownames(FC_DEG2),1]
TPM_DEG2=TPM[rownames(FC_DEG2),]

colData<-data.frame(colData,paste(colData$HOTAIR,colData$YTHDF,sep="."))
colnames(colData)[3] <- "HOTAIR_YTHDF"
annotation_col=data.frame(Condition=factor(colData$HOTAIR_YTHDF))

annotation_row=data.frame(regulate=factor(FC_DEG2[,4]))
#annotation_row=data.frame(log2FC=factor(updown))
annotation_colors =list(Condition=c("OE.DF3-KO"="Tan1","OE.DF3"="Salmon",
                        "Ctrl.DF3-KO"="SkyBlue","Ctrl.DF3"="PaleGreen1"),
                        regulate=c("Up"="MistyRose","Down"="Honeydew","YTHDF3 DEGs"="LemonChiffon"))

colnames(TPM_DEG2)=rownames(annotation_col)
rownames(TPM_DEG2)=rownames(annotation_row)
#library(pheatmap)
list=pheatmap(TPM_DEG2, show_colnames= T, show_rownames= T, scale= "row", fontsize= 6.5,
              clustering_method ="complete",
              cluster_cols=F,cluster_rows=F,
              annotation_col= annotation_col, 
              annotation_row= annotation_row,
              annotation_colors= annotation_colors,
              breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
              col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
              #cutree_row = 2,
              labels_col=rownames(colData),labels_row=symbol,
              #width=5,height=5,
              fontsize_row=5,fontsize_col=5,
              #cutree_col=2
)
res_TPM=cbind(symbol,TPM_DEG2)
write.table(res_TPM,"row_cluster_order_TPM_DEG_119genes.txt",row.names=T,col.names=T,quote=F,sep="\t")

#############################################################################Fig5E UP
gene=c("NLRP2","P2RY10","GRHL2","PAGE5","FBLIM1","AKR1C3")
TPM_res3=TPM_res2[gene,]
for(i in 1:6){
  W=wilcox.test(as.numeric(TPM_res3[i,1:6]),as.numeric(TPM_res3[i,13:16]))
  WP=W$p.value
  print(WP)
}
write.table(TPM_res3,"RNA_ATAC_overlap_6genes_TPM.txt",sep="\t",quote=F)
annotation_col=data.frame(Condition=factor(colData$HOTAIR_YTHDF))
#annotation_row=data.frame(log2FC=factor(updown))
annotation_colors =list(Condition=c("OE.DF3-KO"="Plum2",
                                    "OE.DF3"="DarkOliveGreen1","Ctrl.DF3-KO"="SandyBrown","Ctrl.DF3"="Khaki2"),log2FC=c("UP"="Salmon","Down"="LightSkyBlue"))

colnames(TPM_res3)=rownames(annotation_col)
pheatmap(TPM_res3, show_colnames= T, show_rownames= T, scale= "row", fontsize= 6.5,
         clustering_method ="complete",
         cluster_cols=F,cluster_rows=T,
         annotation_col= annotation_col, 
         #annotation_row= annotation_row,
         annotation_colors= annotation_colors,
         breaks = c(seq(-2,0,length.out=10),seq(0.001,2,length.out=10)),
         col = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         #cutree_row = 2,
         labels_col=rownames(colData),width=5,height=5,
         fontsize_row=8,fontsize_col=8,
         cutree_col=4
         
)
