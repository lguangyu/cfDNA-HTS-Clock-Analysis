# this script is used to analysis and visualize following figures:
# figures 1-2
# figures S1,3-7

# required packages
library(data.table)
library(ggplot2)
library(pheatmap)
library(scales)
library(ComplexHeatmap)
library(tidyverse)
library(jsonlite)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggsankey)
library(dplyr)
library(ggradar)
library(ggalluvial)
library(Metrics)
library(ggsignif)
library(circlize)
library(patchwork)
library(linkET)
library(vegan)
library(ggnewscale)
library(irr)
library(parallel)
library(ggridges)
library(GGally)
library(plotthis)
library(ggcorrplot)
library(cowplot)
library(RColorBrewer)

# set datasets name and colors
dataset_all<-c("gDNA_MSA","gDNA_EPICv2","gDNA_Galaxy","gDNA_Twist","cfDNA_Galaxy","cfDNA_Twist")
dataset_color<-c("#7BA882","#00A79D","#8A92C2","#1A75BC","#BB7113","#BE1D2C")
source_color<-c("#939ca3","#fae8e0")
gender_color<-c("#43b0f1", "#fb6f92")
rep_color<-c("#03658c", "#f29f05")
icc_color<-c("#ef7c8e","#fae8e0","#43b0f1","#a3b7ca")
depth_color<-c("#b2dfeb","#62bed9","#5da7ba","#578e9b","#4e757e")

#### 1. Method-shared CpGs and clock coverage ####
# 1.1 shared CpGs by four methods applied in the paper and HM450K
msa_bed<-fread("data/bed/msa.bed",sep="\t",header=F,data.table = F)
hm450_bed<-fread("data/bed/hm450.bed",sep="\t",header=F,data.table = F)
epic2_bed<-fread("data/bed/epic2.bed",sep="\t",header=F,data.table = F)
galaxy_bed<-fread("data/bed/galaxy.bed",sep="\t",data.table = F)
twist_bed<-fread("data/bed/twist.bed",sep="\t",data.table = F)

colnames(msa_bed)<-c("chr","start","end","probe")
colnames(hm450_bed)<-c("chr","start","end","probe")
colnames(epic2_bed)<-c("chr","start","end","probe")
colnames(galaxy_bed)<-c("chr","start","end")
colnames(twist_bed)<-c("chr","start","end")

msa_bed[["pos"]]<-paste0(msa_bed[,1],"_",msa_bed[,2],"_",msa_bed[,3])
hm450_bed[["pos"]]<-paste0(hm450_bed[,1],"_",hm450_bed[,2],"_",hm450_bed[,3])
epic2_bed[["pos"]]<-paste0(epic2_bed[,1],"_",epic2_bed[,2],"_",epic2_bed[,3])
galaxy_bed[["pos"]]<-paste0(galaxy_bed[,1],"_",galaxy_bed[,2],"_",galaxy_bed[,3])
twist_bed[["pos"]]<-paste0(twist_bed[,1],"_",twist_bed[,2],"_",twist_bed[,3])

tmp_list <- list(
  MSA=msa_bed$pos,
  HM450K=hm450_bed$pos,
  EPICv2=epic2_bed$pos,
  Galaxy=galaxy_bed$pos,
  Twist=twist_bed$pos
)

m <- make_comb_mat(tmp_list)
cs <- comb_size(m)
labels<-names(tmp_list)

pdf("figure/fig_S1a.pdf",width=8,height = 6,onefile = F)
ht = draw(UpSet(m,
                top_annotation = upset_top_annotation(m,ylab="Overlapped CpGs",
                                                      ylim = c(0, max(cs)*1.1),
                                                      annotation_name_rot = 90,
                                                      annotation_name_side = "left",
                                                      height = unit(10, "cm"),
                                                      gp = gpar(
                                                        col = "black",
                                                        fill = "white",
                                                        lwd = 2
                                                      )),
                left_annotation = upset_left_annotation(
                  m,
                  gp = gpar(
                    col = "black",
                    fill = "white",
                    lwd = 2
                  ),width = unit(3, "cm")),
                bg_col = rev(dataset_color[1:5]),set_order = labels))
od = column_order(ht)
decorate_annotation("intersection_size", 
                    {grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + 
                                 unit(2., "pt"), default.units = "native", just = "left", 
                               gp = gpar(fontsize = 10),rot = 90)})
dev.off()

# 1.2 clock coverage
msa_probe<-fread("data/bed/msa.bed",sep="\t",header=F,data.table = F)
epic2_probe<-fread("data/bed/epic2.bed",sep="\t",header=F,data.table = F)
galaxy_probe<-fread("data/bed/galaxy_probe.bed",sep="\t",data.table = F)
twist_probe<-fread("data/bed/twist_probe.bed",sep="\t",data.table = F)

num_cpg_msa<-nrow(msa_bed)
num_cpg_epic2<-nrow(epic2_bed)
num_cpg_galaxy<-nrow(galaxy_bed)
num_cpg_twist<-nrow(twist_bed)

# array clock information
clock_info<-fread("data/clock_info/epigenetic_clocks.txt",sep="\t",header=T,data.table=F)
rownames(clock_info)<-clock_info$Clock
all_clock<-clock_info$Clock
all_clock<-setdiff(all_clock,"magenet")
gp_clock<-c("gpage_10","gpage_30","gpage_71","gpage_a","gpage_b","gpage_c")
pyaging_clock<-setdiff(all_clock[1:53],gp_clock)
clock_cpg<-fromJSON("data/clock_info/all_clock.cpg-dump.json")

# clock coverage by reference
clock_coverage_ref<-data.frame()
for (clock in all_clock){
  if (clock %in% pyaging_clock){tmp_clock_cpg<-clock_cpg[clock][[1]]}else{
    tmp_clock_cpg<-fread(paste0("data/clock_info/",clock,".txt"),sep="\t",header=F,data.table = F)$V1
  }
  num_cpg<-length(tmp_clock_cpg)
  if (clock %in% c("timeseqage","bsclock","komakiage")){
    tmp_clock_pos<-fread(paste0("data/clock_info/",clock,".bed"),sep="\t",header=F,data.table = F)
    num_cpg<-nrow(tmp_clock_pos)
  }
  if ("female" %in% tmp_clock_cpg){num_cpg<-num_cpg-1}
  if ("age" %in% tmp_clock_cpg){num_cpg<-num_cpg-1}
  clock_msa<-intersect(tmp_clock_cpg,msa_probe[,4])
  covered_cpg_msa<-length(clock_msa)
  percent_covered_msa<-100*covered_cpg_msa/num_cpg
  
  clock_epic2<-intersect(tmp_clock_cpg,epic2_probe[,4])
  covered_cpg_epic2<-length(clock_epic2)
  percent_covered_epic2<-100*covered_cpg_epic2/num_cpg
  
  clock_galaxy<-intersect(tmp_clock_cpg,galaxy_probe[,4])
  covered_cpg_galaxy<-length(clock_galaxy)
  if (clock %in% c("timeseqage","bsclock","komakiage")){
    tmp_clock_pos<-fread(paste0("data/clock_info/",clock,"_galaxy.txt"),sep="\t",header=F,data.table = F)
    covered_cpg_galaxy<-nrow(tmp_clock_pos)
  }
  percent_covered_galaxy<-100*covered_cpg_galaxy/num_cpg
  
  clock_twist<-intersect(tmp_clock_cpg,twist_probe[,4])
  covered_cpg_twist<-length(clock_twist)
  if (clock %in% c("timeseqage","bsclock","komakiage")){
    tmp_clock_pos<-fread(paste0("data/clock_info/",clock,"_twist.txt"),sep="\t",header=F,data.table = F)
    covered_cpg_twist<-nrow(tmp_clock_pos)
  }
  percent_covered_twist<-100*covered_cpg_twist/num_cpg
  
  tmp<-c(percent_covered_msa,percent_covered_epic2,percent_covered_galaxy,percent_covered_twist)
  tmp<-t(data.frame(tmp))
  colnames(tmp)<-c("MSA","EPICv2","Galaxy","Twist")
  if (nrow(clock_coverage_ref)==0){clock_coverage_ref<-tmp}else{clock_coverage_ref<-rbind(clock_coverage_ref,tmp)}
}
clock_coverage_ref<-apply(clock_coverage_ref,2,as.numeric)
rownames(clock_coverage_ref)<-all_clock

# clock coverage by observation (53 array clocks)
clock_coverage_obs<-data.frame()
for (i in c(1:6)){
  if (i<3){
    methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_methy_matrix.txt.gz"),sep="\t",data.table=F)
    rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
    methy_matrix_tmp<-methy_matrix_tmp[,-1]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x), na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
  }else{
    methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_depth_matrix.txt.gz"),sep="\t",data.table=F)
    rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
    methy_matrix_tmp<-methy_matrix_tmp[,-1]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x) & x > 4, na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
  }
  clock_coverage_dataset<-data.frame()
  for (clock in all_clock[1:53]){
    if (clock %in% pyaging_clock){tmp_clock_cpg<-clock_cpg[clock][[1]]}else{
      tmp_clock_cpg<-fread(paste0("data/clock_info/",clock,".txt"),sep="\t",header=F,data.table = F)$V1
    }
    num_cpg<-length(tmp_clock_cpg)
    if ("female" %in% clock_cpg){num_cpg<-num_cpg-1}
    if ("age" %in% clock_cpg){num_cpg<-num_cpg-1}
    
    methy_matrix_cpg<-methy_matrix_tmp[rownames(methy_matrix_tmp)%in%tmp_clock_cpg,]
    if (i<3){
      cover_cpg_mean <-apply(methy_matrix_cpg, 2, function(x) sum(!is.na(x), na.rm = TRUE))
    }else{
      cover_cpg_mean <-apply(methy_matrix_cpg, 2, function(x) sum(!is.na(x) & x > 4, na.rm = TRUE))
    }
    cover_cpg_mean<-round(mean(cover_cpg_mean),digits = 0)
    percent_cover<-100*cover_cpg_mean/num_cpg
    
    tmp<-c(clock,percent_cover)
    tmp<-t(data.frame(tmp))
    colnames(tmp)<-c("clock",dataset_all[i])
    if (nrow(clock_coverage_dataset)==0){clock_coverage_dataset=tmp}else{clock_coverage_dataset<-rbind(clock_coverage_dataset,tmp)}
  }
  rownames(clock_coverage_dataset)<-clock_coverage_dataset[,1]
  if (nrow(clock_coverage_obs)==0){clock_coverage_obs<-clock_coverage_dataset
  }else{
    clock_coverage_obs<-cbind(clock_coverage_obs,clock_coverage_dataset[,2])}
}
clock_coverage_obs<-clock_coverage_obs[,-1]
clock_coverage_obs<-apply(clock_coverage_obs,2,as.numeric)
colnames(clock_coverage_obs)<-dataset_all
rownames(clock_coverage_obs)<-all_clock[1:53]

# clock coverage by observation (4 sequencing clocks)
{
clock_coverage_tmp<-data.frame()
clock<-"emseqage"
for (i in c(1:6)){
  if (i<3){
    methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_methy_matrix.txt.gz"),sep="\t",data.table=F)
    rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
    methy_matrix_tmp<-methy_matrix_tmp[,-1]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x), na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
  }else{
    methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_depth_matrix.txt.gz"),sep="\t",data.table=F)
    rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
    methy_matrix_tmp<-methy_matrix_tmp[,-1]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x) & x > 4, na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
  }
  clock_coverage_dataset<-data.frame()
  tmp_clock_cpg<-fread("data/clock_info/emseqage.txt",sep="\t",header=F,data.table = F)
  tmp_clock_cpg<-tmp_clock_cpg$V1
  num_cpg<-length(tmp_clock_cpg)
  methy_matrix_cpg<-methy_matrix_tmp[rownames(methy_matrix_tmp)%in%tmp_clock_cpg,]
  if (i<3){
    cover_cpg_mean <-apply(methy_matrix_cpg, 2, function(x) sum(!is.na(x), na.rm = TRUE))
  }else{
    cover_cpg_mean <-apply(methy_matrix_cpg, 2, function(x) sum(!is.na(x) & x > 4, na.rm = TRUE))
  }
  cover_cpg_mean<-round(mean(cover_cpg_mean),digits = 0)
  percent_cover<-100*cover_cpg_mean/num_cpg
  tmp<-c(clock,percent_cover)
  tmp<-t(data.frame(tmp))
  colnames(tmp)<-c("clock",dataset_all[i])
  if (nrow(clock_coverage_dataset)==0){clock_coverage_dataset=tmp}else{clock_coverage_dataset<-rbind(clock_coverage_dataset,tmp)}
  rownames(clock_coverage_dataset)<-clock_coverage_dataset[,1]
  if (nrow(clock_coverage_tmp)==0){clock_coverage_tmp<-clock_coverage_dataset
  }else{
    clock_coverage_tmp<-cbind(clock_coverage_tmp,clock_coverage_dataset[,2])}
}
clock_cpg_num<-data.frame(clock=c("timeseqage","bsclock","komakiage"),num=c(405,4527,235))
for (clock in c("timeseqage","bsclock","komakiage")){
clock_coverage_dataset<-data.frame()
for (i in c(1:6)){
  if (i<3){
    methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_methy_matrix.txt.gz"),sep="\t",data.table=F)
    rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
    methy_matrix_tmp<-methy_matrix_tmp[,-1]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x), na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
    
    tmp_clock_cpg<-fread(paste0("/mnt/md128/huangwl/ref/Clocks/",clock,"_probe.txt"),sep="\t",header=F,data.table = F)$V7
    num_cpg<-length(tmp_clock_cpg)
    methy_matrix_cpg<-methy_matrix_tmp[rownames(methy_matrix_tmp)%in%tmp_clock_cpg,]
    
  }else{
    methy_matrix_tmp<-fread(paste0("/mnt/md128/huangwl/project/mtAge/data/",dataset_all[i],"/processed/matrix/",dataset_all[i],"_depth_",clock,".matrix"),sep="\t",data.table=F)
    methy_matrix_tmp<-methy_matrix_tmp[,4:ncol(methy_matrix_tmp)]
    num_cpg_mean <-apply(methy_matrix_tmp, 2, function(x) sum(!is.na(x) & x > 4, na.rm = TRUE))
    num_cpg_mean<-round(mean(num_cpg_mean),digits = 0)
  }
  if (i<3){
    cover_cpg_mean <-apply(methy_matrix_cpg, 2, function(x) sum(!is.na(x), na.rm = TRUE))
  }else{
    cover_cpg_mean <-num_cpg_mean
  }
  cover_cpg_mean<-round(mean(cover_cpg_mean),digits = 0)
  num_cpg<-clock_cpg_num[clock_cpg_num$clock==clock,"num"]
  percent_cover<-100*cover_cpg_mean/num_cpg
  tmp<-c(clock,percent_cover)
  tmp<-t(data.frame(tmp))
  colnames(tmp)<-c("clock",dataset_all[i])
  if (nrow(clock_coverage_dataset)==0){clock_coverage_dataset=tmp}else{clock_coverage_dataset<-cbind(clock_coverage_dataset,tmp[,2])}
}
clock_coverage_tmp<-rbind(clock_coverage_tmp,clock_coverage_dataset)
}
rownames(clock_coverage_tmp)<-clock_coverage_tmp[,1]
clock_coverage_tmp<-clock_coverage_tmp[,-1]
clock_coverage_tmp<-apply(clock_coverage_tmp,2,as.numeric)
colnames(clock_coverage_tmp)<-dataset_all
rownames(clock_coverage_tmp)<-c("emseqage","timeseqage","bsclock","komakiage")
}

clock_coverage_obs<-rbind(clock_coverage_obs,clock_coverage_tmp)
write.table(clock_coverage_obs,"./data/clock_info/clock_coverage_obs.txt",sep="\t",row.names=T,col.names=T,quote=F)

# combined and plot
clock_coverage_all<-cbind(clock_coverage_ref,clock_coverage_obs)
clock_coverage_all<-clock_coverage_all[,c("MSA","gDNA_MSA","EPICv2","gDNA_EPICv2",
                                          "Galaxy","gDNA_Galaxy","cfDNA_Galaxy",
                                          "Twist","gDNA_Twist","cfDNA_Twist")]
annotation_col <- data.frame(Type=factor(c("Datasheet","Observation",
                                           "Datasheet","Observation",
                                           "Datasheet","Observation","Observation",
                                           "Datasheet","Observation","Observation"),levels=c("Datasheet","Observation")))
levels(annotation_col$Type)<-c("Reference","Observation")
rownames(annotation_col)<-colnames(clock_coverage_all)
colnames(annotation_col)<-"Coverage Calc. By"
annotation_row<-clock_info[all_clock,c("Generation/Prediction","Category")]
annotation_colors <- list("Coverage Calc. By" = c(Reference = "#d6eadf", Observation = "#eac4d5"),
                          "Category"=c("Hi-Cov-All"="#244CBC","Hi-Repro"="#CB232A",
                                       "Hi-Cov-nonMSA"="#5DB9DD","Other"="#C2C4C6","Sequencing-Clock"="#096377"),
                          "Generation/Prediction"=c(First="#fbf2c4",Second="#F5DE24",Third="#AECE34",
                                                    Fourth="#00743F",Gestational="#F6C3DB",Other="#3D3D3D"))
pdf("figure/fig_S1b.pdf",width=7,height = 8,onefile = F)
pheatmap::pheatmap(clock_coverage_all,
                   color=colorRampPalette(c("#f402f1", "#ffffff","#62cff4"))(100),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   gaps_col = c(2,4,7),
                   number_format = "%.0f",
                   angle_col=45,
                   annotation_names_row=FALSE,
                   annotation_names_col=FALSE,
                   legend = TRUE,
                   border = "grey60",
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   annotation_colors = annotation_colors,
                   main = "Percentage of Covered Clock CpGs")
dev.off()

#### 2. PCA, correlation across samples and consistency across datasets ####
# 2.1 combined meta data and methylation matrix of all datasets
sample_info_all<-data.frame()
methy_matrix_all<-data.frame()
for (i in c(1:6)){
  sample_info_tmp<-fread(paste0("data/meta/",dataset_all[i],"_meta.txt"),sep="\t",data.table=F)
  sample_info_rep1<-sample_info_tmp[,c(2,3,4,8)]
  colnames(sample_info_rep1)[4]<-"Name"
  sample_info_rep2<-sample_info_tmp[,c(2,3,4,9)]
  colnames(sample_info_rep2)[4]<-"Name"
  sample_info_tmp<-rbind(sample_info_rep1,sample_info_rep2)
  sample_info_tmp[["Rep"]]<-c(rep("Rep1",24),rep("Rep2",24))
  sample_info_tmp[["Dataset"]]<-rep(dataset_all[i],nrow(sample_info_tmp))
  if (nrow(sample_info_all)==0){sample_info_all<-sample_info_tmp}else{sample_info_all<-rbind(sample_info_all,sample_info_tmp)}
  methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_methy_matrix.txt.gz"),sep="\t",data.table=F)
  colnames(methy_matrix_tmp)[1]<-"ID"
  if (nrow(methy_matrix_all)==0){methy_matrix_all<-methy_matrix_tmp}else{methy_matrix_all<-merge(methy_matrix_all,methy_matrix_tmp,by="ID")}
}

sample_info_all$Dataset<-factor(sample_info_all$Dataset,levels=dataset_all)
sample_info_all[["Source"]]<-sample_info_all$Dataset
levels(sample_info_all$Source)<-c("gDNA","gDNA","gDNA","gDNA","cfDNA","cfDNA")
rownames(methy_matrix_all)<-methy_matrix_all[,1]
methy_matrix_all<-methy_matrix_all[,-1]

# 2.2 pca of all samples by method-shared CpGs
df_pca<-t(na.omit(methy_matrix_all))
df_pca <- prcomp(df_pca,)
df_pca_sum <- summary(df_pca)
df<- df_pca$x[,1:2]
pc <- df_pca_sum$importance[2,]*100
df <- as.data.frame(df)
df<-cbind(df,sample_info_all)

# color by Dataset
p<-ggplot(data=df,aes(x=PC1,y=PC2))+
  geom_point(aes(fill=Dataset),color="black",alpha=0.7,size=3,shape=21)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(title="Method-shared CpGs (136,437)",x=paste0("PC1 (",pc[1],"%)"),y=paste0("PC2 (",pc[2],"%)"))+
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Dataset),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(name="",values = dataset_color) +
  scale_fill_manual(name="",values = dataset_color)+
  theme_bw()+
  theme(axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22,angle=90),
        axis.text.y=element_text(size=22),
        axis.text.x=element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size=24,hjust=0.5,vjust=0),
        panel.grid=element_blank())
ggsave("figure/fig_1b.pdf",plot=p,width = 7.1,heigh=5)

# color by Rep
p<-ggplot(data=df,aes(x=PC1,y=PC2,color=Rep))+
  geom_point(aes(fill=Rep),color="black",alpha=0.7,size=3,shape=21)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(title="Method-shared CpGs (136,437)",x=paste0("PC1 (",pc[1],"%)"),y=paste0("PC2 (",pc[2],"%)"))+
  stat_ellipse(data=df,
                geom = "polygon",level=0.95,
                linetype = 2,size=0.5,
                aes(fill=Rep),
                alpha=0.2,
                show.legend = T)+
  scale_color_manual(name="",values = rep_color) +
  scale_fill_manual(name="",values = rep_color)+
  theme_bw()+
  theme(axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22,angle=90),
        axis.text.y=element_text(size=22),
        axis.text.x=element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size=24,hjust=0.5,vjust=0),
        panel.grid=element_blank())
ggsave("figure/fig_1c.pdf",plot=p,width = 6,heigh=5)

# 2.3 correlation of all samples by method-shared CpGs
methy_cor<-cor(methy_matrix_all,use = "complete.obs")
annotation_colors <- list(
  Gender = c("M" = "#43b0f1", "F" = "#fb6f92"),
  Source=c("gDNA"=source_color[1],"cfDNA"=source_color[2]),
  Dataset=c("gDNA_MSA"=dataset_color[1],"gDNA_EPICv2"=dataset_color[2],"gDNA_Galaxy"=dataset_color[3],
            "gDNA_Twist"=dataset_color[4],"cfDNA_Galaxy"=dataset_color[5],"cfDNA_Twist"=dataset_color[6])
)
col_anno <- HeatmapAnnotation(
  df = sample_info_all[,c("Gender","Dataset","Source")],
  col = annotation_colors,
  annotation_name_side = "left",
  show_legend = TRUE
)
pdf("figure/fig_S1e.pdf",width = 5.5,heigh=4)
Heatmap(
  matrix = methy_cor,
  column_title = "Method-shared CpGs (136,437)",
  column_title_side = "top",
  name = "Correlation",
  col = colorRamp2(c(0.75, 0.9, 1), c("#0099ff", "#ffffff", "#ff004c")),
  top_annotation = col_anno,
  cluster_rows = TRUE,        
  cluster_columns = TRUE,     
  show_row_names = FALSE,     
  show_column_names = FALSE,
  heatmap_legend_param = list(
    legend_direction = "vertical",  
    legend_width = unit(5, "cm")
  )
)
dev.off()

# 2.4 consistency across datasets by method-shared CpGs
# mean of each CpG across samples
methy_matrix_mean<-data.frame()
for (i in c(1:6)){
  methy_matrix_tmp<-methy_matrix_all[,sample_info_all[sample_info_all$Dataset==dataset_all[i],"Name"]]
  methy_matrix_tmp<-apply(methy_matrix_tmp,1,function(a){return(mean(a,na.rm=TRUE))})
  methy_matrix_tmp<-data.frame(methy_matrix_tmp)
  if(nrow(methy_matrix_mean)==0){methy_matrix_mean=methy_matrix_tmp}else{methy_matrix_mean<-cbind(methy_matrix_mean,methy_matrix_tmp)}
}
colnames(methy_matrix_mean)<-dataset_all
p <- qcorrplot(correlate(methy_matrix_mean), type = "lower",grid_col = NA) +    
  geom_point(shape=21, size=18, fill = NA, color = "black")+  
  geom_point(aes(size=abs(r), fill=r),               
             shape=21,
             color = "black") +    
  scale_size(range = c(18, 18),guide="none") +        
  new_scale("size")+
  scale_fill_gradientn(limits = c(0.9,1),
                       breaks = seq(0.9,1,0.03),
                       colors = rev(brewer.pal(11, "Spectral"))) +
  guides( fill = guide_colorbar(title = "Pearson's r", 
                                keyheight = unit(4, "cm"),
                                keywidth = unit(0.5, "cm"),
                                order = 2)) + 
  labs(title="Method-shared CpGs (136,437)")+
  theme(legend.box.spacing = unit(0, "pt"),
        axis.text = element_text(size=22),
        legend.position = "right",
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 0)
  )
ggsave("figure/fig_1f.pdf",plot=p,width = 8,heigh=6)

# 2.5 pca of all samples by clock-shared CpGs
clock_cpg_all<-c()
for (clock in all_clock[all_clock$Category=="Hi-Cov-All","Clock"]){
  if (clock %in% names(clock_cpg)){
    tmp_clock_cpg<-clock_cpg[clock][[1]]
  }else{
    tmp_clock_cpg<-fread(paste0("data/clock_info/",clock,".txt"),sep="\t",header=F,data.table = F)
    tmp_clock_cpg<-tmp_clock_cpg$V1
  }
  clock_cpg_all<-c(clock_cpg_all,tmp_clock_cpg)
}
clock_cpg_all<-unique(clock_cpg_all)
methy_matrix_clock<-methy_matrix_all[clock_cpg_all,]
df_pca<-t(na.omit(methy_matrix_clock))
df_pca <- prcomp(df_pca,)
df_pca_sum <- summary(df_pca)
df<- df_pca$x[,1:2]
pc <- df_pca_sum$importance[2,]*100
df <- as.data.frame(df)
df<-cbind(df,sample_info_all)

# color by Dataset
p<-ggplot(data=df,aes(x=PC1,y=PC2))+
  geom_point(aes(fill=Dataset),color="black",alpha=0.7,size=3,shape=21)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(title="Clock-shared CpGs (2,727)",x=paste0("PC1 (",pc[1],"%)"),y=paste0("PC2 (",pc[2],"%)"))+
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Dataset),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(name="",values = dataset_color) +
  scale_fill_manual(name="",values = dataset_color)+
  theme_bw()+
  theme(axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22,angle=90),
        axis.text.y=element_text(size=22),
        axis.text.x=element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size=24,hjust=0.5,vjust=0),
        panel.grid=element_blank())
ggsave("figure/fig_S1c.pdf",plot=p,width = 7.1,heigh=5)

# color by Rep
p<-ggplot(data=df,aes(x=PC1,y=PC2,color=Rep))+
  geom_point(aes(fill=Rep),color="black",alpha=0.7,size=3,shape=21)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(title="Clock-shared CpGs (2,727)",x=paste0("PC1 (",pc[1],"%)"),y=paste0("PC2 (",pc[2],"%)"))+
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Rep),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(name="",values = group_color) +
  scale_fill_manual(name="",values = group_color)+
  theme_bw()+
  theme(axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22,angle=90),
        axis.text.y=element_text(size=22),
        axis.text.x=element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size=24,hjust=0.5,vjust=0),
        panel.grid=element_blank())
ggsave("figure/fig_S1d.pdf",plot=p,width = 6,heigh=5)

# 2.6 correlation of all samples by clock-shared CpGs
methy_cor<-cor(methy_matrix_clock,use = "complete.obs")
pdf("figure/fig_S1f.pdf",width = 5.5,heigh=4)
Heatmap(
  matrix = methy_cor,
  column_title = "Clock-shared CpGs (2,727)",
  column_title_side = "top",
  name = "Correlation",
  col = colorRamp2(c(0.8, 0.9, 1), c("#0099ff", "#ffffff", "#ff004c")),
  top_annotation = col_anno,
  cluster_rows = TRUE,        
  cluster_columns = TRUE,     
  show_row_names = FALSE,     
  show_column_names = FALSE,
  heatmap_legend_param = list(
    legend_direction = "vertical",  
    legend_width = unit(5, "cm")
  )
)
dev.off()
  
# 2.7 consistency across datasets by clock-shared CpGs
methy_matrix_mean<-data.frame()
for (i in c(1:6)){
  methy_matrix_tmp<-methy_matrix_clock[,sample_info_all[sample_info_all$Dataset==dataset_all[i],"Name"]]
  methy_matrix_tmp<-apply(methy_matrix_tmp,1,function(a){return(mean(a,na.rm=TRUE))})
  methy_matrix_tmp<-data.frame(methy_matrix_tmp)
  if(nrow(methy_matrix_mean)==0){methy_matrix_mean=methy_matrix_tmp}else{methy_matrix_mean<-cbind(methy_matrix_mean,methy_matrix_tmp)}
}
colnames(methy_matrix_mean)<-dataset_all
methy_matrix_mean<-na.omit(methy_matrix_mean)
p<- qcorrplot(correlate(methy_matrix_mean), type = "lower",grid_col = NA) +    
  geom_point(shape=21, size=18, fill = NA, color = "black")+  
  geom_point(aes(size=abs(r), fill=r),               
             shape=21,
             color = "black") +    
  scale_size(range = c(19, 18),guide="none") +        
  new_scale("size")+
  scale_fill_gradientn(limits = c(0.9,1),
                       breaks = seq(0.9,1,0.03),
                       colors = rev(brewer.pal(11, "Spectral"))) +
  guides( fill = guide_colorbar(title = "Pearson's r", 
                                keyheight = unit(4, "cm"),
                                keywidth = unit(0.5, "cm"),
                                order = 2)) + 
  labs(title="Clock-shared CpGs (2,727)")+
  theme(legend.box.spacing = unit(0, "pt"),
        axis.text = element_text(size=22),
        legend.position = "right",
        legend.text = element_text(size=20),
        legend.title = element_text(size=22),
        plot.title = element_text(size = 24, hjust = 0.5, vjust = 0)
  )
ggsave("figure/fig_S1h.pdf",plot=p,width = 8,heigh=6)

#### 3. ICC calculated by beta values or predicted age ####
# 3.1 ICC of CpGs
calculate_icc<- function(input_data) {
  rep1 <- as.numeric(input_data[1:24])
  rep2 <- as.numeric(input_data[25:48])
  valid_pairs <- !(is.na(rep1) | is.na(rep2))
  rep1 <- rep1[valid_pairs]
  rep2 <- rep2[valid_pairs]
  if (length(rep1) < 2) {
    return(NA)
  }
  tmp_mean<-mean(input_data,na.rm=TRUE)
  tmp_sd<-sd(input_data,na.rm=TRUE)
  icc_input <- cbind(rep1, rep2)
  icc_result <- icc(icc_input, model="twoway", type="agreement", unit="single")
  results <- list(
    Mean=tmp_mean,
    SD=tmp_sd,
    ICC = icc_result$value,
    CI_lower = icc_result$lbound,
    CI_upper = icc_result$ubound,
    F_test = icc_result$Fvalue,
    p_value = icc_result$p.value
  )
  return(results)
}
cl <- makeCluster(8)
clusterExport(cl, c("icc", "calculate_icc"))
cpg_icc<-data.frame()
for (i in c(1:6)){
  methy_matrix_tmp<-fread(paste0("data/matrix/",dataset_all[i],"_methy_matrix.txt.gz"),sep="\t",data.table=F)
  rownames(methy_matrix_tmp)<-methy_matrix_tmp[,1]
  methy_matrix_tmp<-methy_matrix_tmp[,-1]
  icc_values <- parApply(cl, methy_matrix_tmp, 1, calculate_icc)
  result_tmp <- as.data.frame(do.call(rbind, icc_values))
  result_tmp<-apply(result_tmp, 2, as.numeric)
  result_tmp<-data.frame(result_tmp)
  result_tmp[["CpG"]]<-rownames(methy_matrix_tmp)
  result_tmp[["Dataset"]]<-rep(dataset_all[i],nrow(methy_matrix_tmp))
  if (nrow(cpg_icc)==0){cpg_icc<-result_tmp}else{cpg_icc<-rbind(cpg_icc,result_tmp)}
}
stopCluster(cl)
cpg_icc$Group <- cut(
  cpg_icc$ICC,
  breaks = c(-Inf, 0.5, 0.75, 0.9, Inf),
  labels = c("Poor(<0.5)", "Moderate(0.5-0.75)", "Good(0.75-0.9)", "Excellent(≥0.9)"),
  include.lowest = TRUE
)
cpg_icc<-na.omit(cpg_icc)
write.table(cpg_icc,"result/cpg_icc.txt",sep="\t",row.names = F,col.names = T,quote=F)
cpg_icc<-fread("result/cpg_icc.txt",sep="\t",header=T,data.table=F)

# 3.2 ICC of predicted age
cl <- makeCluster(8)
clusterExport(cl, c("icc", "calculate_icc"))
age_icc<-data.frame()
for (i in c(1:6)){
  pred_age_tmp<-fread(paste0("result/",dataset_all[i],"_predicted_age.txt"),sep="\t",data.table=F)
  rownames(pred_age_tmp)<-pred_age_tmp$Name
  pred_age_emseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_emseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_emseqage)<-pred_age_emseqage[,1]
  pred_age_emseqage<-pred_age_emseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["emseqage"]]<-pred_age_emseqage$Predicted_Age
  pred_age_timeseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_timeseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_timeseqage)<-pred_age_timeseqage[,1]
  pred_age_timeseqage<-pred_age_timeseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["timeseqage"]]<-pred_age_timeseqage$Predicted_Age
  pred_age_tmp<-pred_age_tmp[,all_clock[1:55]]
  pred_age_tmp<-t(pred_age_tmp)
  icc_values <- parApply(cl, pred_age_tmp, 1, calculate_icc)
  result_tmp <- as.data.frame(do.call(rbind, icc_values))
  result_tmp<-apply(result_tmp, 2, as.numeric)
  result_tmp<-data.frame(result_tmp)
  result_tmp[["Clock"]]<-rownames(pred_age_tmp)
  result_tmp[["Dataset"]]<-rep(dataset_all[i],nrow(pred_age_tmp))
  if (nrow(age_icc)==0){age_icc<-result_tmp}else{age_icc<-rbind(age_icc,result_tmp)}
}
stopCluster(cl)
age_icc$Group <- cut(
  age_icc$ICC,
  breaks = c(-Inf, 0.5, 0.75, 0.9, Inf),
  labels = c("Poor(<0.5)", "Moderate(0.5-0.75)", "Good(0.75-0.9)", "Excellent(≥0.9)"),
  include.lowest = TRUE
)
age_icc<-na.omit(age_icc)
write.table(age_icc,"result/age_icc.txt",sep="\t",row.names = F,col.names = T,quote=F)

# fig_1g
age_icc$Clock<-factor(age_icc$Clock,levels=rev(unique(age_icc$Clock)))
age_icc$Dataset<-factor(age_icc$Dataset,levels=dataset_all)
age_icc[["Group2"]]<-factor(age_icc$Group,levels=unique(age_icc$Group))
levels(age_icc$Group2)<-rev(c("Poor","Moderate","Good","Excellent"))
p <- ggplot(age_icc, aes(y = Clock, x = ICC)) +
  geom_point(aes(fill=Group2),color="black",pch=21,size = 3,alpha=0.7) +
  geom_errorbarh(aes(color=Group2,xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = c(0.5, 0.75, 0.9), linetype = "dashed", color = "gray60") +
  scale_color_manual(values = icc_color) +
  scale_fill_manual(values = icc_color) +
  labs(y = "", x = "ICC", 
       title = "",
       color = "",
       fill="") +
  theme_bw() +
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         strip.background = element_rect(fill="white"),
         strip.text = element_text(size=22,face="bold"),
         axis.title = element_text(size = 22),
         axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
         axis.text.y = element_text(size = 22, color = "black"),
         axis.ticks = element_line(color = "black"),
         legend.text = element_text(size = 22, color = "black"),
         legend.title = element_text(size = 22, color = "black"),
         legend.position = "bottom",
         plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=1)

ggsave("figure/fig_1g.pdf",plot=p,width = 16.2,heigh=16)

#### 4. Predicted age related analysis ####
# 4.1 consistency
pred_age_all<-data.frame()
for (i in c(1:6)){
  pred_age_tmp<-fread(paste0("result/",dataset_all[i],"_predicted_age.txt"),sep="\t",data.table=F)
  rownames(pred_age_tmp)<-pred_age_tmp$Name
  pred_age_emseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_emseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_emseqage)<-pred_age_emseqage[,1]
  pred_age_emseqage<-pred_age_emseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["emseqage"]]<-pred_age_emseqage$Predicted_Age
  pred_age_timeseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_timeseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_timeseqage)<-pred_age_timeseqage[,1]
  pred_age_timeseqage<-pred_age_timeseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["timeseqage"]]<-pred_age_timeseqage$Predicted_Age
  if (nrow(pred_age_all)==0){pred_age_all<-pred_age_tmp}else{pred_age_all<-rbind(pred_age_all,pred_age_tmp)}
}
rownames(pred_age_all)<-pred_age_all$Name
pred_age_all<-pred_age_all[,-1]

# clock consistency across datasets
all_category<-c("Hi-Cov-All", "Hi-Repro","Hi-Cov-nonMSA")
for (i in c(1:length(all_category))){
  clock_category<-all_category[i]
  methy_matrix_mean<-apply(pred_age_all[,rownames(clock_info[clock_info$Category==clock_category,])],1,function(a){return(mean(a,na.rm=TRUE))})
  names(methy_matrix_mean)<-sample_info_all$Name
  methy_matrix_mean<-data.frame(methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[1],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[2],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[3],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[4],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[5],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[6],"Name"]])
  colnames(methy_matrix_mean)<-dataset_all
  
  p <- qcorrplot(correlate(methy_matrix_mean), type = "lower",grid_col = NA) +    
    geom_point(shape=21, size=18, fill = NA, color = "black")+  
    geom_point(aes(size=abs(r), fill=r),               
               shape=21,
               color = "black") +    
    scale_size(range = c(18, 18), guide = "none") +        
    new_scale("size")+
    scale_fill_gradientn(limits = c(0.9,1),
                         breaks = seq(0.9,1,0.03),
                         colors = rev(brewer.pal(11, "Spectral"))) +
    guides( fill = guide_colorbar(title = "Pearson's r", 
                                  keyheight = unit(4, "cm"),
                                  keywidth = unit(0.5, "cm"),
                                  order = 3)) + 
    labs(title=clock_category)+
    theme(legend.box.spacing = unit(0, "pt"),
          axis.text = element_text(size=22),
          legend.position = "right",
          legend.text = element_text(size=20),
          legend.title = element_text(size=22),
          plot.title = element_text(size = 24, hjust = 0.5, vjust = 0)
    )
  if (i==1){combined_plot<-p}else{combined_plot<-combined_plot+p}
}
combined_plot<-combined_plot+
  plot_layout(nrow = 1, ncol = 3,guides = "collect")
ggsave("figure/fig_1h.pdf",plot=combined_plot,width = 24,heigh=6)

# for each clock
j=0
for (i in c(1:55)){
  clock<-rownames(clock_info)[i]
  methy_matrix_mean<-pred_age_all[,clock]
  names(methy_matrix_mean)<-sample_info_all$Name
  methy_matrix_mean<-data.frame(methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[1],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[2],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[3],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[4],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[5],"Name"]],
                                methy_matrix_mean[sample_info_all[sample_info_all$Dataset==dataset_all[6],"Name"]])
  colnames(methy_matrix_mean)<-dataset_all
  methy_matrix_mean<-na.omit(methy_matrix_mean)
  tmp<-correlate(methy_matrix_mean)
  if (min(tmp$r)>=0.7){
    j=j+1
    p <- qcorrplot(correlate(methy_matrix_mean), type = "lower",grid_col = NA) + 
      geom_point(shape=21, size=18, fill = NA, color = "black")+  
      geom_point(aes(size=abs(r), fill=r),               
                 shape=21,
                 color = "black") +    
      scale_size(range = c(18, 18), guide = "none") +        
      new_scale("size")+
      scale_fill_gradientn(limits = c(0.7,1),
                           breaks = seq(0.7,1,0.1),
                           colors = rev(brewer.pal(11, "Spectral"))) +
      guides( fill = guide_colorbar(title = "Pearson's r", 
                                    keyheight = unit(4, "cm"),
                                    keywidth = unit(0.5, "cm"),
                                    order = 3)) +
      labs(title=clock)+
      theme(legend.box.spacing = unit(0, "pt"),
            axis.text = element_text(size=22),
            legend.position = "right",
            legend.text = element_text(size=20),
            legend.title = element_text(size=22),
            plot.title = element_text(size = 24, hjust = 0.5, vjust = 0)
      )
    if (j==1){combined_plot<-p}else{combined_plot<-combined_plot+p}
  }
}
combined_plot<-combined_plot+
  plot_layout(nrow = 4, ncol = 5,guides = "collect")
ggsave("figure/fig_S5.pdf",plot=combined_plot,width = 40,heigh=24)

# 4.2 RD/Metrics scatter of Hi-Repro clocks
hi_repro_clock<-clock_info[clock_info$Category=="Hi-Repro","Clock"]
for (i in c(1:6)){
  pred_age_tmp<-fread(paste0("result/",dataset_all[i],"_predicted_age.txt"),sep="\t",data.table=F)
 
  # scatter between replicates
  plot_list <- list()
  for (j in c(1:length(hi_repro_clock))) {
    clock<-hi_repro_clock[j]
    rep1<-pred_age_tmp[1:24,clock]
    rep2<-pred_age_tmp[25:48,clock]
    df<-data.frame(rep1,rep2)
    colnames(df)<-c("rep1","rep2")
    df<-na.omit(df)
    cor=round(cor(df$rep1,df$rep2),2)
    mae=round(mae(df$rep1,df$rep2),2)
    medae=round(median(abs(df$rep1-df$rep2)),2)
    rmse=round(rmse(df$rep1,df$rep2),2)
    metrics_label <- paste(
      paste0("MRD = ", mae),
      sep = "\n"
    )
    if (i==1){title=clock}else{title=""}
    if (j==1){y_lab=""}else{y_lab=""}
    if (i==6){x_lab=""}else{x_lab=""}
    p<-ggplot(df, aes(x = rep1, y = rep2)) +
      geom_point(fill=dataset_color[i],color="black",shape = 21,size=3) +
      geom_abline(slope = 1,intercept = 0,linetype="dashed",color="grey")+
      geom_smooth(method = "lm", se = FALSE,color="grey") +
      stat_cor(aes(rep1, rep2),size=7) +
      labs(title=title,x = x_lab, y = y_lab) + 
      annotate("text", x = max(df$rep1), y = (min(df$rep2)+0.1*(max(df$rep2)-min(df$rep2))),
               label = metrics_label, hjust = 1, vjust = 0.5, size = 7,
               color = "black") +
      theme_bw()+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_line(),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 22, color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.title = element_blank(),
            plot.title = element_text(size = 24, hjust = 0.5, vjust = 0, margin = margin(b = 20)))
    plot_list[[j]] <- p
    if (j==1){plot_dataset<-p}else{plot_dataset<-plot_dataset+p}
  }
  plot_dataset<-plot_dataset+plot_layout(nrow = 1, ncol = 6)
  if (i==1){plot_result<-plot_dataset}else(plot_result<-plot_result/plot_dataset)
  combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 6)
  plot_result<-plot_result+plot_layout(heights = c(1, 1, 1, 1, 1, 1))&
    theme(plot.margin = margin(10, 10, 10, 10, "pt"))
  ggsave("figure/fig_S3.pdf", plot_result, width = 30, height = 30)
  
  # scatter between CA and BA
  plot_list <- list()
  for (j in c(1:length(hi_repro_clock))) {
    clock<-hi_repro_clock[j]
    ca<-pred_age_tmp[,"Age"]
    ba<-pred_age_tmp[,clock]
    df<-data.frame(ca,ba,pred_age_tmp$Rep)
    colnames(df)<-c("CA","BA","Rep")
    df<-na.omit(df)
    df1<-df[df$Rep=="Rep1",]
    df2<-df[df$Rep=="Rep2",]
    cor=round(cor(df1$CA,df1$BA),2)
    mae=round(mae(df1$CA,df1$BA),2)
    medae=round(median(abs(df1$CA-df1$BA)),2)
    rmse=round(rmse(df1$CA,df1$BA),2)
    cor2=round(cor(df2$CA,df2$BA),2)
    mae2=round(mae(df2$CA,df2$BA),2)
    medae2=round(median(abs(df2$CA-df2$BA)),2)
    rmse2=round(rmse(df2$CA,df2$BA),2)
    metrics_label<- paste("Rep1/Rep2",
                          paste0("R = ", cor,"/",cor2),
                          paste0("MAE = ", mae,"/",mae2),
                          sep = "\n"
    )
    
    if (i==1){title=clock}else{title=""}
    if (j==1){y_lab=""}else{y_lab=""}
    if (i==6){x_lab=""}else{x_lab=""}
    p<-ggplot(df, aes(x = CA, y = BA)) +
      geom_point(aes(fill=Rep),color="black",shape = 21,size=3) +
      geom_abline(slope = 1,intercept = 0,linetype="dashed",color="grey")+
      geom_smooth(aes(fill=Rep,color=Rep),method = "lm", se = FALSE) +
      scale_color_manual(values=group_color)+
      scale_fill_manual(values=group_color)+
      labs(title=title,x = x_lab, y = y_lab) + 
      annotate("text", x = max(df$CA), y = (min(df$BA)+0.15*(max(df$BA)-min(df$BA))),
               label = metrics_label, hjust = 1, vjust = 0.5, size = 7,
               color = "black") +
      theme_bw()+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_line(),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 22, color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "null",
            plot.title = element_text(size = 24, hjust = 0.5, vjust = 0, margin = margin(b = 20)))
    plot_list[[i]] <- p
    if (j==1){plot_dataset1<-p}else{plot_dataset1<-plot_dataset1+p}
  }
  plot_dataset1<-plot_dataset1+plot_layout(nrow = 1, ncol = 6)
  if (i==1){plot_result1<-plot_dataset1}else(plot_result1<-plot_result1/plot_dataset1)
  plot_result1<-plot_result1+plot_layout(heights = c(1, 1, 1, 1, 1, 1))&
    theme(plot.margin = margin(10, 10, 10, 10, "pt"))
  ggsave("figure/fig_S7.pdf", plot_result1, width = 30, height = 30)
}

# 4.3 boxplot of RD for all clocks
rep_diff_all<-data.frame()
for (i in c(1:6)){
  pred_age_tmp<-fread(paste0("result/",dataset_all[i],"_predicted_age.txt"),sep="\t",data.table=F)
  rownames(pred_age_tmp)<-pred_age_tmp$Name
  pred_age_emseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_emseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_emseqage)<-pred_age_emseqage[,1]
  pred_age_emseqage<-pred_age_emseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["emseqage"]]<-pred_age_emseqage$Predicted_Age
  pred_age_timeseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_timeseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_timeseqage)<-pred_age_timeseqage[,1]
  pred_age_timeseqage<-pred_age_timeseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["timeseqage"]]<-pred_age_timeseqage$Predicted_Age
  pred_age_rep1<-pred_age_tmp[1:24,] %>%
    pivot_longer(cols =clock_info$Clock[1:55], names_to = "clock", values_to = "rep1")
  pred_age_rep2<-pred_age_tmp[25:48,] %>%
    pivot_longer(cols =clock_info$Clock[1:55], names_to = "clock", values_to = "rep2")
  rep_diff_tmp<-data.frame(pred_age_rep1,pred_age_rep2[,"rep2"])
  rep_diff_tmp[["diff"]]<-abs(rep_diff_tmp$rep2-rep_diff_tmp$rep1)
  rep_diff_tmp[["dataset"]]<-rep(dataset_all[i],nrow(rep_diff_tmp))
  if (nrow(rep_diff)==0){rep_diff_all<-rep_diff_tmp}else{rep_diff_all<-rbind(rep_diff_all,rep_diff_tmp)}
}
rep_diff_all<-na.omit(rep_diff_all)

# Hi-Repro
rep_diff<-rep_diff_all[rep_diff_all$clock%in%clock_info[clock_info$Category=="Hi-Repro","Clock"],]
rep_diff$clock<-factor(rep_diff$clock,levels=rev(unique(rep_diff$clock)))
rep_diff$dataset<-factor(rep_diff$dataset,levels=dataset_all)
p<-ggplot(rep_diff, aes(x = diff, y = clock)) +
  geom_boxplot(color="#CB232A",fill="#CB232A",alpha=0.3)+
  labs(title = "Hi-Repro",x = "Replicate Difference",y = NULL)+
  theme_bw() +
  theme(
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "null",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~dataset,nrow=1)
ggsave("figure/fig_S4a.pdf",plot=p,width = 16,heigh=5)

tmp <- rep_diff %>%
  group_by(dataset, clock) %>%
  summarise(
    median_value = median(diff, na.rm = TRUE),
    .groups = "drop"
  )

# Hi-Cov-All
rep_diff<-rep_diff_all[rep_diff_all$clock%in%all_clock[all_clock$Category=="Hi-Cov-All","Clock"],]
rep_diff$clock<-factor(rep_diff$clock,levels=rev(unique(rep_diff$clock)))
rep_diff$dataset<-factor(rep_diff$dataset,levels=dataset_all)
p<-ggplot(rep_diff, aes(x = diff, y = clock)) +
  geom_boxplot(color="#244CBC",fill="#244CBC",alpha=0.3)+
  labs(title = "Hi-Cov-All",x = "Replicate Difference",y = NULL)+
  theme_bw() +
  theme(
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "null",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~dataset,nrow=1)
ggsave("figure/fig_S4b.pdf",plot=p,width = 16,heigh=10)

# Hi-Cov-nonMSA
rep_diff<-rep_diff_all[rep_diff_all$clock%in%all_clock[all_clock$Category=="Hi-Cov-nonMSA","Clock"],]
rep_diff$clock<-factor(rep_diff$clock,levels=rev(unique(rep_diff$clock)))
rep_diff$dataset<-factor(rep_diff$dataset,levels=dataset_all)
p<-ggplot(rep_diff, aes(x = diff, y = clock)) +
  geom_boxplot(color="#5DB9DD",fill="#5DB9DD",alpha=0.3)+
  scale_fill_manual(name="Model",values = rev(model_color))+
  labs(title = "Hi-Cov-nonMSA",x = "Replicate Difference",y = NULL)+
  scale_x_continuous(limits =c(0,80),breaks=c(0,25,50,75)) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "null",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~dataset,nrow=1)
ggsave("figure/fig_S4c.pdf",plot=p,width = 16,heigh=8)

# Sequencing−Clock
rep_diff<-rep_diff_all[rep_diff_all$clock%in%sequence_clock,]
rep_diff$clock<-factor(rep_diff$clock,levels=rev(unique(rep_diff$clock)))
rep_diff$dataset<-factor(rep_diff$dataset,levels=dataset_all)
p<-ggplot(rep_diff, aes(x = diff, y = clock)) +
  geom_boxplot(color="#096377",fill="#096377",alpha=0.3)+
  scale_fill_manual(name="Model",values = rev(model_color))+
  labs(title = "Sequencing−Clock",x = "Replicate Difference",y = NULL)+
  theme_bw() +
  theme(
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "null",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~dataset,nrow=1)
ggsave("figure/fig_S4d.pdf",plot=p,width = 16,heigh=3)

# 4.4 R and MAE plot
calculate_metrics <- function(pred_age,clock_use) {
  metrics <- matrix(nrow = length(clock_use), ncol = 4)
  colnames(metrics) <- c("R","MAE", "MedAE", "RMSE")
  rownames(metrics) <- clock_use
  age_true<-pred_age$Age
  for (i in 1:length(clock_use)) {
    clock<-clock_use[i]
    age_pred <- pred_age[,clock]
    result<-data.frame("true"=age_true,"pred"=age_pred)
    result<-na.omit(result)
    r<-cor(result$true,result$pred)
    mae<-mean(abs(result$true - result$pred))
    medae<-median(abs(result$true - result$pred))
    rmse <- sqrt(mean((result$true - result$pred)^2))
    metrics[i, ] <- c(r, mae, medae, rmse)
  }
  return(as.data.frame(metrics))
}
metrics_all<-data.frame()
for (i in c(1:6)){
  pred_age_tmp<-fread(paste0("result/",dataset_all[i],"_predicted_age.txt"),sep="\t",data.table=F)
  rownames(pred_age_tmp)<-pred_age_tmp$Name
  
  pred_age_emseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_emseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_emseqage)<-pred_age_emseqage[,1]
  pred_age_emseqage<-pred_age_emseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["emseqage"]]<-pred_age_emseqage$Predicted_Age
  pred_age_timeseqage<-fread(paste0("result/",dataset_all[i],"/",dataset_all[i],"_predicted_age_timeseqage.txt"),sep="\t",header=T,data.table=F)
  rownames(pred_age_timeseqage)<-pred_age_timeseqage[,1]
  pred_age_timeseqage<-pred_age_timeseqage[rownames(pred_age_tmp),]
  pred_age_tmp[["timeseqage"]]<-pred_age_timeseqage$Predicted_Age
  
  for (rep in c("Rep1","Rep2")){
    metrics_rep<-calculate_metrics(pred_age_tmp[pred_age_tmp$Rep==rep,],all_clock$Clock[1:55])
    metrics_rep[["Rep"]]<-rep(rep,nrow(metrics_rep))
    metrics_rep[["Dataset"]]<-rep(dataset_all[i],nrow(metrics_rep))
    if (nrow(metrics_all)==0){metrics_all<-metrics_rep}else{metrics_all<-rbind(metrics_all,metrics_rep)}
  }
}
metrics_all$Dataset<-factor(metrics_all$Dataset,levels=unique(metrics_all$Dataset))
metrics_all[["Clock"]]<-rep(all_clock$Clock[1:55],12)

# Hi_Repro
metrics_use<-metrics_all[metrics_all$Clock%in%clock_info[clock_info$Category=="Hi-Repro","Clock"],]
metrics_order<-metrics_use[metrics_use$Dataset=="gDNA_EPICv2",]
metrics_order<-metrics_order[order(metrics_order$MAE,decreasing = T),]

metrics_use$Clock<-factor(metrics_use$Clock,levels=unique(metrics_order$Clock))
metrics_use$Dataset<-factor(metrics_use$Dataset)

p <-ggplot(metrics_use,aes(x = MAE, y = Clock)) +
  geom_segment(aes(y = Clock, yend = Clock, x =0, xend = MAE), linewidth =0.8, color ='grey80') +
  geom_point(aes(fill = R,pch=Rep), color="black",pch=21,size =6) +
  scale_fill_gradient2(low="#3498db",mid="white",high="#871c1c",midpoint = 0.8,
                       limits = c(0.3, 1),breaks = c(0.3,0.5,0.7,0.9),oob = scales::squish)+
  scale_x_continuous(limits =c(0,30),expand =c(0,1)) +
  labs(title ="Hi-Repro", x ="MAE", y =NULL) +
  theme_bw(base_size =18) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=1)

ggsave("figure/fig_2a.pdf", p, width = 18,heigh=5)

# Hi-Cov-All
metrics_use<-metrics_all[metrics_all$Clock%in%clock_info[clock_info$Category=="Hi-Cov-All","Clock"],]
metrics_order<-metrics_use[metrics_use$Dataset=="gDNA_EPICv2",]
metrics_order<-metrics_order[order(metrics_order$MAE,decreasing = T),]
metrics_use$Clock<-factor(metrics_use$Clock,levels=unique(metrics_order$Clock))
metrics_use$Dataset<-factor(metrics_use$Dataset)

p <-ggplot(metrics_use,aes(x = MAE, y = Clock)) +
  geom_segment(aes(y = Clock, yend = Clock, x =0, xend = MAE), linewidth =0.8, color ='grey80') +
  geom_point(aes(fill = R,pch=Rep), color="black",pch=21,size =6) +
  scale_fill_gradient2(low="#3498db",mid="white",high="#871c1c",midpoint = 0.8,
                       limits = c(0.3, 1),breaks = c(0.3,0.5,0.7,0.9),oob = scales::squish)+
  scale_x_continuous(limits =c(0,30),expand =c(0,1)) +
  labs(title ="Hi-Cov-All", x ="MAE", y =NULL) +
  theme_bw(base_size =18) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=1)

ggsave("figure/fig_S6a.pdf", p, width = 18,heigh=10)

# Hi-Cov-All
metrics_use<-metrics_all[metrics_all$Clock%in%clock_info[clock_info$Category=="Hi-Cov-nonMSA","Clock"],]
metrics_order<-metrics_use[metrics_use$Dataset=="gDNA_EPICv2",]
metrics_order<-metrics_order[order(metrics_order$MAE,decreasing = T),]
metrics_use$Clock<-factor(metrics_use$Clock,levels=unique(metrics_order$Clock))
metrics_use$Dataset<-factor(metrics_use$Dataset)

# metrics_use[metrics_use$MAE>30,"MAE"]<-30
p <-ggplot(metrics_use,aes(x = MAE, y = Clock)) +
  geom_segment(aes(y = Clock, yend = Clock, x =0, xend = MAE), linewidth =0.8, color ='grey80') +
  geom_point(aes(fill = R,pch=Rep), color="black",pch=21,size =6) +
  scale_fill_gradient2(low="#3498db",mid="white",high="#871c1c",midpoint = 0.8,
                       limits = c(0.3, 1),breaks = c(0.3,0.5,0.7,0.9),oob = scales::squish)+
  labs(title ="Hi-Cov-nonMSA", x ="MAE", y =NULL) +
  theme_bw(base_size =18) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=1)

ggsave("figure/fig_S6b.pdf", p, width = 18,heigh=8)

# Sequencing-Clock
metrics_use<-metrics_all[metrics_all$Clock%in%c("emseqage","timeseqage"),]
metrics_order<-metrics_use[metrics_use$Dataset=="gDNA_EPICv2",]
metrics_order<-metrics_order[order(metrics_order$MAE,decreasing = T),]

metrics_use$Clock<-factor(metrics_use$Clock,levels=unique(metrics_order$Clock))
metrics_use$Dataset<-factor(metrics_use$Dataset)

p <-ggplot(metrics_use,aes(x = MAE, y = Clock)) +
  geom_segment(aes(y = Clock, yend = Clock, x =0, xend = MAE), linewidth =0.8, color ='grey80') +
  geom_point(aes(fill = R,pch=Rep), color="black",pch=21,size =6) +
  scale_fill_gradient2(low="#3498db",mid="white",high="#871c1c",midpoint = 0.8,
                       limits = c(0.3, 1),breaks = c(0.3,0.5,0.7,0.9),oob = scales::squish)+
  scale_x_continuous(limits =c(0,50),expand =c(0,1)) +
  labs(title ="Sequencing-Clock", x ="MAE", y =NULL) +
  theme_bw(base_size =18) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=1)

ggsave("figure/fig_2f.pdf", p, width = 18,heigh=3)

# 4.5 scatter between ICC and R/MAE
rownames(age_icc)<-paste0(age_icc$Dataset,"_",age_icc$Clock)

plot_data1<-metrics_all[metrics_all$Rep=="Rep1",]
rownames(plot_data1)<-paste0(plot_data1$Dataset,"_",plot_data1$Clock)
plot_data1<-plot_data1[rownames(age_icc),]
plot_data1[["ICC"]]<-age_icc$ICC

plot_data2<-metrics_all[metrics_all$Rep=="Rep2",]
rownames(plot_data2)<-paste0(plot_data2$Dataset,"_",plot_data2$Clock)
plot_data2<-plot_data2[rownames(age_icc),]
plot_data2[["ICC"]]<-age_icc$ICC

plot_data<-rbind(plot_data1,plot_data2)
plot_data<-plot_data[plot_data$Clock%in%clock_info[clock_info$Category=="Hi-Cov-All","Clock"],]
p<-ggplot(plot_data, aes(x = ICC, y = R)) +
  geom_point(aes(fill=Rep),color="black",shape = 21,size=3,alpha = 0.7) +
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  labs(title=category,x = "ICC", y = "R") + 
  geom_smooth(aes(fill=Rep,color=Rep),method = "lm", se = FALSE) +
  stat_cor(
    aes(label = ..r.label..),
    label.x = 0.6,
    label.y = 0.6,
    size = 6,
    color = "black",
    show.legend = FALSE
  ) +
  stat_cor(
    aes(label = ..p.label..),
    label.x = 0.6,
    label.y = 0.52,
    size = 6,
    color = "black",
    show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "null",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=2)
ggsave("figure/fig_2e.pdf", p, width = 9, height = 7)

p<-ggplot(plot_data, aes(x = ICC, y = MAE)) +
  geom_point(aes(fill=Rep),color="black",shape = 21,size=3,alpha = 0.7) +
  geom_smooth(aes(fill=Rep,color=Rep),method = "lm", se = FALSE) +
  stat_cor(aes(ICC, MAE),size=6) +
  scale_color_manual(name="",values=group_color)+
  scale_fill_manual(name="",values=group_color)+
  labs(title=category,x = "ICC", y = "MAE") + 
  theme_bw(base_size =18) +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size=22,face="bold"),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(angle = 90, size = 22, color = "black", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 22, color = "black"),
    legend.position = "bottom",
    plot.title = element_text(size = 24, hjust = 0.5, vjust = 0))+
  facet_wrap(~Dataset,nrow=2)

ggsave("figure/fig_2d.pdf", p, width = 9, height = 8)

#### 5. validation by GSE232346 and GSE245628 ####
metrics1<-fread("result/metrics_GSE232346.txt",sep="\t",header=T,data.table=F)
metrics2<-fread("result/metrics_GSE245628.txt",sep="\t",header=T,data.table=F)
plot_data<-rbind(metrics1,metrics2)
plot_data$clock<-factor(plot_data$clock,levels=c("hannum","horvath2013","skinandblood","dnamphenoage","timeseqage"))
plot_data<-plot_data[plot_data$clock!="timeseqage",]
plot_data[["label_cor"]]<-as.character(round(plot_data$cor,2))
plot_data[["label_mae"]]<-as.character(round(plot_data$mae,2))
colors<-c("#AD3B8F","#46C1BB","#F0545D","#94C73D","#d48a50")
colors<-colors[1:4]
p1 <- ggplot(plot_data, aes(x = clock, y = cor, color = clock)) +  
  geom_point(size = 4) +  
  geom_segment(aes(x = clock, xend = clock, y = 0, yend = cor),
               linewidth = 1.5) +  
  geom_text(aes(y = cor, label = label_cor),size = 8, angle = 0, nudge_y = 0.1, hjust = 0.5, vjust = 0) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5), expand = c(0, 0)) +  
  scale_color_manual(values = colors) +  
  labs(y = "R") +  
  theme_classic2() +   
  theme(legend.position = "none",
        strip.text = element_text(size=22),
        axis.title.y = element_text(size = 22, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 1.5, color = "black"),
        axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
        axis.ticks.length.y = unit(0.2, "cm"))+
  facet_wrap(~platform,nrow=2)

p2 <- ggplot(plot_data, aes(x = clock, y = 0.1, fill = clock)) + 
  geom_tile(width = 0.95) +  
  scale_x_discrete(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +  
  scale_fill_manual(values = colors) +  
  theme_bw() +  
  theme(legend.position = "none", 
        panel.border = element_rect(linewidth = 2),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 22,color = "black",family = "sans",hjust = 1,vjust = 1,angle = 45))

p3 <- ggplot(plot_data, aes(x = clock, y = mae, color = clock)) +  
  geom_point(size = 4) +  
  geom_segment(aes(x = clock, xend = clock, y = 0, yend = mae),
               linewidth = 1.5) +  
  geom_text(aes(y = mae, label = label_mae),size = 8, angle = 0, nudge_y = 5, hjust = 0.5, vjust = 0) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 65, 20), expand = c(0, 0)) +  
  scale_color_manual(values = colors) +  
  labs(y = "MAE") +  
  theme_classic2() +   
  theme(legend.position = "none",
        strip.text = element_text(size=22),
        axis.title.y = element_text(size = 22, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 1.5, color = "black"),
        axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
        axis.ticks.length.y = unit(0.2, "cm"))+
  facet_wrap(~platform,nrow=2)

p12 <- p1 / p2 + plot_layout(heights = c(1, 0.05))
p32 <- p3 / p2 + plot_layout(heights = c(1, 0.05))

ggsave("figure/fig_2g.pdf", plot = p12, width = 6, height = 6)
ggsave("figure/fig_2h.pdf", plot = p32, width = 6, height = 6)
