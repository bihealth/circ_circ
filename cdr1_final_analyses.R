############ install and load all libraries to run the script

library(readxl)
library(dplyr)
library(ggplot2)
library(tmod)
library(tidyr)
library(ggpubr)
library(ggrepel)
library(imputeTS)
library(cowplot)
library(clusterProfiler)



#### read linear gene count data

covariate <-read.table("./data/covariate_file_scn.txt", sep='\t', header=T) %>% tbl_df() %>% mutate(filename=sub("mapping/feature_counts/","", filename))%>% mutate(filename=sub(".all_mates.*","", filename))

linear_all_scn<-readRDS("./data/linear_all_scn.rds")


twin<-readxl::read_xls("./data/elife-10518-supp2-v2.xls") %>% filter(`WGCNA module`=="lightsteelblue1")

genes<-readxl::read_xls("./data/elife-10518-supp2-v2.xls") %>% filter(`WGCNA module`=="lightsteelblue1") %>% arrange(-FPKM)  %>% dplyr::select(Ensembl_ID,Gene_ID)

colnames(genes)<-c("ensembls", "genes")



aln_stats_scn <- read.table("./data/SCN_stats.csv", sep='\t', header=T) %>% tibble::as_tibble() %>% left_join(covariate, by=c("sample_id"="filename")) %>% mutate(sample=label) %>% dplyr::select(sample, key, value)

scn_reads_all <- aln_stats_scn %>% filter(key=="Number of input reads") %>% mutate(mean.v=mean(value)) %>% mutate(norm_factor=mean.v/value)
scn_reads_mu <- aln_stats_scn %>% filter(key=="Uniquely mapped reads number")  %>% mutate(mean.v=mean(value)) %>% mutate(norm_factor=mean.v/value)
scn_reads_um <- aln_stats_scn %>% filter(key=="Number of unmapped reads")

### read circRNA count data 

circ_all_scn<-readRDS("data/circ_all_scn.rds")


# Figure 1 A (circRNA count per time-point)

circ_count_scn<-circ_all_scn %>% group_by(circRNA_ID, time) %>% count() %>% filter(n>1) %>% group_by(time) %>% count()
circ_count_scn1 <- circ_count_scn %>% filter(time=="ZT2") %>% mutate(time="ZT2.1")

circ_count_scn<-rbind(circ_count_scn, circ_count_scn1)

circ_count_scn <- circ_count_scn %>% mutate(replicate=0) %>% dplyr::select(time, replicate, n) %>% mutate(type="unique")

colnames(circ_count_scn) <- c("time", "replicate", "count", "type")


repr_scn_circs <- (circ_all_scn %>% group_by(circRNA_ID, time) %>% count() %>% filter(n>1))$circRNA_ID %>% unique()



scn_Fig.A_data1 <-circ_all_scn %>% filter(circRNA_ID %in% repr_scn_circs) %>% group_by(time, replicate) %>% summarise(count = n()) %>% ungroup() %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")

scn_Fig.A_data2 <- circ_all_scn %>% filter(circRNA_ID %in% repr_scn_circs) %>% group_by(time, replicate) %>% summarise(count = n()) %>% ungroup()

scn_Fig.A_data <- bind_rows(scn_Fig.A_data1, scn_Fig.A_data2) %>% mutate(type="per_replicate")

scn_Fig.A_data<-scn_Fig.A_data %>% rbind(scn_Fig.A_data, circ_count_scn)


scn_Fig.A <- scn_Fig.A_data %>% ggline(x="time", y="count", add = c("mean_sd","jitter"), shape = 5, color="type", palette = "jco", point.size = 3)  + ylab("number of circRNAs per replicate")+ annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22","ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ylab("# of circRNAs")+xlab("zeitgeber(time)") + ggtitle("Number of detected circRNAs") +ylim(0,600)
scn_Fig.A


########### Figure 1 B (overall circRNA expression in the SCN)

sumCDR1reads <- sum((circ_all_scn  %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% dplyr::select(X.junction_reads))$X.junction_reads)


scn_Fig.B <- circ_all_scn %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() %>% mutate(Reads.per.circRNA=replace(Reads.per.circRNA, Reads.per.circRNA<10, 10)) %>% ggplot(aes(x=Reads.per.circRNA))+geom_histogram(color="black", fill="white") + scale_x_log10(breaks=c(0,1,5,10,100,1000,sumCDR1reads))+ scale_y_log10(breaks=c(0,1,5,10,100,1000,5000)) + geom_text(aes(label="Cdr1as", x=sumCDR1reads, y=1 + 2), size=3.5, col="darkred") + geom_segment(aes(x=sumCDR1reads, xend=sumCDR1reads, y=1 + 1.5, yend=1 + 0.25),arrow=arrow(length=unit(2, "mm")), col="darkred") + xlab("aggregate (54 samples) reads") + ylab("# of circRNAs") +ggtitle("Overall circRNA expression") + theme(text = element_text(size=15), plot.title = element_text(size=15, face="bold", hjust = 0.5))
scn_Fig.B




all_w.cdr <- circ_all_scn %>% group_by(sample_id) %>% summarize(read.count=sum(X.junction_reads))  %>% mutate(reads="All circRNAs")

all_wo.cdr <- circ_all_scn %>% filter(circRNA_ID != "chrX:61183248|61186174") %>% group_by(sample_id) %>% summarize(read.count=sum(X.junction_reads)) %>% mutate(reads="Cdr1as (excluded)")



######## Figure 1 C (circular reads per sample)

 all_w.cdr.1 <-all_w.cdr %>% bind_rows(all_wo.cdr) %>% left_join(scn_reads_all, by=c("sample_id"="sample"))  %>% mutate(read.count.norm=read.count * norm_factor) %>% tidyr::separate(sample_id, c("time","replicate")) %>% mutate(sample_id=paste(time, replicate, sep="_")) 
 all_w.cdr.2 <-all_w.cdr %>% bind_rows(all_wo.cdr) %>% left_join(scn_reads_all, by=c("sample_id"="sample"))  %>% mutate(read.count.norm=read.count * norm_factor) %>% tidyr::separate(sample_id, c("time","replicate")) %>% mutate(sample_id=paste(time, replicate, sep="_")) %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")
 
 scn_Fig.C_data <- bind_rows(all_w.cdr.1, all_w.cdr.2)


scn_Fig.C <- scn_Fig.C_data %>% ggline(x="time", y="read.count.norm", add = c("mean_sd","jitter"), shape = 5, color="reads", palette = "jco", point.size = 3) + annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22", "ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of head-to-tail reads")+xlab("zeitgeber(time)") + ggtitle("Circular reads per sample")  + theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))  + ylim(0,1500)
scn_Fig.C

########## Figure 1 D (Cdr1as expression per time-point)

cdr1.as <- circ_all_scn  %>% dplyr::select(circRNA_ID, X.junction_reads, sample_id) %>% spread(sample_id, X.junction_reads) %>% na_replace(0) %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% gather(sample_id, X.junction_reads, ZT10_1:ZT6_9) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(HeadToTail_reads=X.junction_reads * norm_factor)  %>% tidyr::separate(sample_id, c("time","replicate")) %>% mutate(sample_id=paste(time, replicate, sep="_"))


cdr1.as.1 <- cdr1.as %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")

cdr1.as_data <- bind_rows(cdr1.as, cdr1.as.1)

scn_Fig.D <- cdr1.as_data %>% ggline(x="time", y="HeadToTail_reads", add = c("mean_sd", "jitter"), col="darkred", shape = 5, point.size=3) + annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22", "ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ggtitle("Cdr1as expression")  +xlab("zeitgeber(time)")+ylab("norm. # of head-to-tail reads") #+stat_compare_means(comparisons = list(c("ZT14", "ZT10"), c("ZT14", "ZT18"), c("ZT2", "ZT6"), c("ZT6", "ZT10"), c("ZT18", "ZT22")), label = "p.signif", label.y = c(500,450, 350,100, 400))


scn_Fig.D


##################### Figure 2 A and B


# read the tables from Piwecka et al., Sience 2017
# Piwecka et al sequenced poltyA select and total RNA-Seq data of Cdr1as KO mice
#The first (outcommented) block is the polyA data

#cortex<-readxl::read_xlsx("./data/aam8526_Piwecka_SM_table_S4.xlsx", sheet=1, skip=8) %>% tbl_df()
#cereb<-readxl::read_xlsx("./data/aam8526_Piwecka_SM_table_S4.xlsx", sheet=2) %>% tbl_df()
#ofl<-readxl::read_xlsx("./data/aam8526_Piwecka_SM_table_S4.xlsx", sheet=3) %>% tbl_df()
#hippo<-readxl::read_xlsx("./data/aam8526_Piwecka_SM_table_S4.xlsx", sheet=4) %>% tbl_df()


cortex<-readxl::read_xlsx("./data/aam8526_piwecka_sm_table_s3.xlsx", sheet=1, skip=8) %>% tbl_df()
cereb<-readxl::read_xlsx("./data/aam8526_piwecka_sm_table_s3.xlsx", sheet=2) %>% tbl_df()
ofl<-readxl::read_xlsx("./data/aam8526_piwecka_sm_table_s3.xlsx", sheet=3) %>% tbl_df()
hippo<-readxl::read_xlsx("./data/aam8526_piwecka_sm_table_s3.xlsx", sheet=4) %>% tbl_df()

# read tables from Xu et al., Neuron 2022

genes_30min<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=1) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Induced"))$gene
genes_1h<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=2) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Induced"))$gene
genes_1h_r<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=2) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Reduced"))$gene
genes_3h<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=3) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Induced"))$gene
genes_3h_r<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=3) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Reduced"))$gene
genes_6h<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=4) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Induced"))$gene
genes_6h_r<-(readxl::read_xlsx("./data/1-s2.0-S0896627321005705-mmc2.xlsx", sheet=4) %>% tbl_df() %>% mutate(gene=`...1`)  %>% filter(Class=="Reduced"))$gene

mods <- data.frame(ID=c("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7"), Title=c("30 min induced", "1h induced", "1h reduced", "3h induced", "3h reduced", "6h induced", "6h reduced"))

m2g<-list(ID1=genes_30min, ID2=genes_1h, ID3=genes_1h_r, ID4=genes_3h, ID5=genes_3h_r,  ID6=genes_6h, ID7=genes_6h_r)
mymset <- makeTmod(mods, m2g)


#cortex <- cortex %>% mutate(color="grey") %>% mutate(color=ifelse((log2FoldChange < 0) , "darkblue", "darkred")) %>% as.data.frame()
#cereb <- cereb %>% mutate(color="grey") %>% mutate(color=ifelse((log2FoldChange < 0) , "darkblue", "darkred")) %>% as.data.frame()
#hippo <- hippo %>% mutate(color="grey") %>% mutate(color=ifelse((log2FoldChange < 0) , "darkblue", "darkred")) %>% as.data.frame() 
#ofl <- ofl %>% mutate(color="grey") %>% mutate(color=ifelse((log2FoldChange < 0) , "darkblue", "darkred")) %>% as.data.frame() 


cortex <- cortex %>% mutate(color="grey")  %>% as.data.frame()
cereb <- cereb %>% mutate(color="grey")  %>% as.data.frame()
hippo <- hippo %>% mutate(color="grey")  %>% as.data.frame() 
ofl <- ofl %>% mutate(color="grey")  %>% as.data.frame() 

cortex$color[cortex$padj<0.05 & cortex$log2FoldChange < 0]= "blue"
cortex$color[cortex$padj<0.05 & cortex$log2FoldChange > 0]= "red"


hippo$color[hippo$padj<0.05 & hippo$log2FoldChange < 0]= "blue"
hippo$color[hippo$padj<0.05 & hippo$log2FoldChange > 0]= "red"


cereb$color[cereb$padj<0.05 & cereb$log2FoldChange < 0]= "blue"
cereb$color[cereb$padj<0.05 & cereb$log2FoldChange > 0]= "red"


ofl$color[ofl$padj<0.05 & ofl$log2FoldChange < 0]= "blue"
ofl$color[ofl$padj<0.05 & ofl$log2FoldChange > 0]= "red"

############ Figure 2 A

tmod::evidencePlot(cortex$symbol, mset=mymset, m="ID1", gene.labels=T, gene.colors=cortex$color,  main="Xu et al., Neuron 2021 \ninduced 30min after light exposure", gl.cex=0.8)

########### Figure 2 B

tmod::evidencePlot(ofl$symbol, mset=mymset, m="ID1", gene.labels=T, gene.colors=ofl$color,  main="Xu et al., Neuron 2021 \ninduced 30min after light exposure", gl.cex=0.8)


################## supplementary Xu et al., light induced


tmod::evidencePlot(hippo$symbol, mset=mymset, m="ID1", gene.labels=T, gene.colors=hippo$color,  main="30min after light exposure", gl.cex=0.8)
tmod::evidencePlot(cereb$symbol, mset=mymset, m="ID1", gene.labels=T, gene.colors=cereb$color,  main="30min after light exposure", gl.cex=0.8)




####### Supplementary Figure 1A

cdr1_lin<-linear_all_scn %>% filter(Geneid=="ENSMUSG00000090546") %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts=counts*norm_factor) %>% ggline(x="time", y="norm_counts", add = c("mean_sd","jitter"), shape = 5, palette = "jco", point.size = 3, color="darkred")+ggtitle("Cdr1 linear read counts")+ylab("normalized counts")+ annotate("rect", xmin = 3.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.5)+ylim(0,7500)

cdr1_lin

####### Supplementary Figure 1B


cyrano_lin<- linear_all_scn %>% filter(Geneid=="ENSMUSG00000085438") %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts=counts*norm_factor) %>% ggline(x="time", y="norm_counts", add = c("mean_sd","jitter"), shape = 5, palette = "jco", point.size = 3, color="darkred")+ggtitle("Cyrano read counts")+ylab("normalized counts")+ annotate("rect", xmin = 3.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.5)+ylim(0,2500)

cyrano_lin

###### Figure 1 F

twin_histplot<-linear_all_scn %>% left_join(scn_reads_mu, by=c("sample_id"="sample")) %>% mutate(fpkm=counts*10^9/value) %>% mutate(fpkm=fpkm/Length)  %>% right_join(twin, by=c("Geneid"="Ensembl_ID"))  %>% group_by(Geneid) %>% summarize(mean.fpkm=mean(fpkm)) %>% arrange(-mean.fpkm) %>%ggplot(aes(x=log2(mean.fpkm+1))) + geom_histogram(color="black", fill="white")+  theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("# of genes")+ggtitle("Expression of the twin-peak genes")+  geom_vline(xintercept=3.57)+ geom_vline(xintercept=6.41)

twin_histplot

######  Figure 2 C

top50 <-linear_all_scn %>% left_join(scn_reads_mu, by=c("sample_id"="sample")) %>% mutate(fpkm=counts*10^9/value) %>% mutate(fpkm=fpkm/Length)  %>% right_join(twin, by=c("Geneid"="Ensembl_ID"))  %>% group_by(Geneid) %>% summarize(mean.fpkm=mean(fpkm)) %>% arrange(-mean.fpkm) %>% head(50)

top100 <-linear_all_scn %>% left_join(scn_reads_mu, by=c("sample_id"="sample")) %>% mutate(fpkm=counts*10^9/value) %>% mutate(fpkm=fpkm/Length)  %>% right_join(twin, by=c("Geneid"="Ensembl_ID"))  %>% group_by(Geneid) %>% summarize(mean.fpkm=mean(fpkm)) %>% arrange(-mean.fpkm) %>% head(100)

top200 <-linear_all_scn %>% left_join(scn_reads_mu, by=c("sample_id"="sample")) %>% mutate(fpkm=counts*10^9/value) %>% mutate(fpkm=fpkm/Length)  %>% right_join(twin, by=c("Geneid"="Ensembl_ID"))  %>% group_by(Geneid) %>% summarize(mean.fpkm=mean(fpkm)) %>% arrange(-mean.fpkm) %>% head(200)




tmp.mean<-linear_all_scn %>% filter(Geneid %in% top50$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()


total_plot_data.1 <- linear_all_scn %>% filter(Geneid %in% top50$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>%  group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean)  %>% mutate(relative_expr = log2(group.mean/tmp.mean)) 

total_plot_data.2 <- total_plot_data.1 %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")

total_plot_data<- bind_rows(total_plot_data.1, total_plot_data.2)

total_plot <- total_plot_data %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred") + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22", "ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5)) + annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Pemborke et al., 2015")  +ylim(-3,3)

total_plot

#plot_empty <- plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')



######## Figure 2 D

##### load the dark-dark SCN RNA-Seq data from Cheng et al., study

samples<-c("CT0_3xSCN_WT1","CT0_3xSCN_WT2","CT0_3xSCN_WT3","CT0_3xSCN_WT4","CT0_3xSCN_WT5","CT6_3xSCN_WT1","CT6_3xSCN_WT2","CT6_3xSCN_WT3","CT6_3xSCN_WT4","CT6_3xSCN_WT5", "CT12_3xSCN_WT1","CT12_3xSCN_WT2","CT12_3xSCN_WT3","CT12_3xSCN_WT4","CT12_3xSCN_WT5","CT18_3xSCN_WT1","CT18_3xSCN_WT2","CT18_3xSCN_WT3","CT18_3xSCN_WT4","CT18_3xSCN_WT5")



linear_all_scn_pa<-readRDS("./data/linear_all_scn_pa.rds")



aln_stats_scn_pa <- read.csv("./data/SCN_stats_pa.csv", header=T)

scn_reads_all_pa <- aln_stats_scn_pa %>% filter(key=="Number of input reads")%>% mutate(mean.v=mean(value)) %>% mutate(norm_factor=mean.v/value)
scn_reads_mu_pa <- aln_stats_scn_pa %>% filter(key=="Uniquely mapped reads number")



tmp.mean_pa <-linear_all_scn_pa %>% filter(Geneid %in% top50$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()


polyA_plot.data_1 <- linear_all_scn_pa %>% filter(Geneid %in% top50$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean_pa) %>% mutate(relative_expr = log2(group.mean/tmp.mean))

polyA_plot.data_2 <- polyA_plot.data_1 %>% subset(time=="CT0") %>% mutate(time="CT0.1")

polyA_plot.data <- bind_rows(polyA_plot.data_1, polyA_plot.data_2)

polyA_plot <- polyA_plot.data %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred")  + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Cheng et al., 2019") + scale_x_discrete(limits=c("CT0", "CT6","CT12","CT18","CT0.1")) +ylim(-3,3)




###################### Supplementary Figure 2 



tmp.mean_pa_100 <-linear_all_scn_pa %>% filter(Geneid %in% top100$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()


tmp.mean_pa_200 <-linear_all_scn_pa %>% filter(Geneid %in% top200$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()


polyA_plot.data100_1 <- linear_all_scn_pa %>% filter(Geneid %in% top100$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean_pa_100) %>% mutate(relative_expr = log2(group.mean/tmp.mean))

polyA_plot.data100_2 <- polyA_plot.data100_1 %>% subset(time=="CT0") %>% mutate(time="CT0.1")

polyA_plot.data100 <- bind_rows(polyA_plot.data100_1, polyA_plot.data100_2)

polyA_plot100 <- polyA_plot.data100 %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred")  + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Cheng et al., 2019") + scale_x_discrete(limits=c("CT0", "CT6","CT12","CT18","CT0.1")) +ylim(-3,3)

polyA_plot.data200_1 <- linear_all_scn_pa %>% filter(Geneid %in% top200$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all_pa, by=c("sample_id")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean_pa_200) %>% mutate(relative_expr = log2(group.mean/tmp.mean))

polyA_plot.data200_2 <- polyA_plot.data200_1 %>% subset(time=="CT0") %>% mutate(time="CT0.1")

polyA_plot.data200 <- bind_rows(polyA_plot.data200_1, polyA_plot.data200_2)

polyA_plot200 <- polyA_plot.data200 %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred")  + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Cheng et al., 2019") + scale_x_discrete(limits=c("CT0", "CT6","CT12","CT18","CT0.1")) +ylim(-3,3)

tmp.mean_100<-linear_all_scn %>% filter(Geneid %in% top100$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()

tmp.mean_200<-linear_all_scn %>% filter(Geneid %in% top200$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>% group_by(genes) %>% summarize(tmp.mean = mean(norm_counts)) %>% ungroup()

total_plot_data100.1 <- linear_all_scn %>% filter(Geneid %in% top100$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>%  group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean_100)  %>% mutate(relative_expr = log2(group.mean/tmp.mean)) 

total_plot_data100.2 <- total_plot_data100.1 %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")

total_plot_data100<- bind_rows(total_plot_data100.1, total_plot_data100.2)

total_plot100 <- total_plot_data100 %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred") + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22", "ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5)) + annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Pemborke et al., 2015")  +ylim(-3,3)

total_plot_data200.1 <- linear_all_scn %>% filter(Geneid %in% top200$Geneid) %>% left_join(genes, by=c("Geneid"="ensembls")) %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts = counts * norm_factor) %>%  group_by(genes, time) %>% summarize(group.mean=mean(norm_counts)) %>% left_join(tmp.mean_200)  %>% mutate(relative_expr = log2(group.mean/tmp.mean)) 

total_plot_data200.2 <- total_plot_data200.1 %>% subset(time=="ZT2") %>% mutate(time="ZT2.1")

total_plot_data200<- bind_rows(total_plot_data200.1, total_plot_data200.2)

total_plot200 <- total_plot_data200 %>% ggline(x="time", y="relative_expr", group="genes", palette = "jco", point.size = 0.2, color="darkred") + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22", "ZT2.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5)) + annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5)+ggtitle("Pemborke et al., 2015")  +ylim(-3,3)


total_plot100 

total_plot200 

polyA_plot100 

polyA_plot200


############  Supplementary Figure 1  (Cdr1 and Cyrano)


cdr1_lin<-linear_all_scn %>% filter(Geneid=="ENSMUSG00000090546") %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts=counts*norm_factor)

cdr1_lin.2 <- cdr1_lin %>% filter(time=="ZT2") %>% mutate(time="ZT2.1")

cdr1_lin.data<-bind_rows(cdr1_lin, cdr1_lin.2)

cdr1_lin<- cdr1_lin.data %>% ggline(x="time", y="norm_counts", add = c("mean_sd","jitter"), shape = 5, palette = "jco", point.size = 3, color="darkred")+ggtitle("Cdr1 linear read counts")+ylab("normalized counts")+ annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5)+ylim(0,7500) + scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22","ZT2.1")) 


cyrano_lin<- linear_all_scn %>% filter(Geneid=="ENSMUSG00000085438") %>% left_join(scn_reads_all, by=c("sample_id"="sample")) %>% mutate(norm_counts=counts*norm_factor) 

cyrano_lin.2 <- cyrano_lin %>%   filter(time=="ZT2") %>% mutate(time="ZT2.1")

cyrano_lin.data<- bind_rows(cyrano_lin, cyrano_lin.2)

cyrano_lin<-cyrano_lin.data %>% ggline(x="time", y="norm_counts", add = c("mean_sd","jitter"), shape = 5, palette = "jco", point.size = 3, color="darkred")+ggtitle("Cyrano read counts")+ylab("normalized counts")+ annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha=0.5)+ylim(0,2500)+ scale_x_discrete(limits=c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22","ZT2.1")) 

cowplot::plot_grid(cdr1_lin, cyrano_lin, ncol=2, labels=c("A", "B"))


######### Figure 3 A (Cortex, exonic and intronic rates)


f<-readRDS("data/frontal_cortex_exonic_intronic_rates.RDS")

f.1 <- f %>% filter(time=="F7") %>% mutate(time="F7.1")

f<- bind_rows(f,f.1)

f_ie<-f %>% ggline(x="time", y="rate", add = c("mean_sd","jitter"), shape = 5, color="anno", palette = "jco", point.size = 3) + annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("F7" ,"F11", "F15", "F19", "F23", "F3", "F7.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("rate of reads")+xlab("time") + ggtitle("Cortex") + theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.6), legend.text = element_text(size=10)) 

f_ie

######### Figure 3 A (Hippocampus, exonic and intronic rates)


h<-readRDS("data/hippocampus_exonic_intronic_rates.RDS")

h.1 <- h %>% filter(time=="H7") %>% mutate(time="H7.1")

h <- bind_rows(h,h.1)

h_ie<-h %>% ggline(x="time", y="rate", add = c("mean_sd","jitter"), shape = 5, color="anno", palette = "jco", point.size = 3) + annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("H7", "H11", "H15", "H19", "H23", "H3", "H7.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("rate of reads")+xlab("time") + ggtitle("Hippocampus") + theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.6), legend.text = element_text(size=10)) 


######## Figure 3 B and C 

######## Overlap with circbase circRNAs

samples_hippo <- c("H3_1", "H7_1","H11_1","H15_1","H19_1","H23_1", "H3_2", "H7_2","H11_2","H15_2","H19_2","H23_2", "H3_3", "H7_3","H11_3","H15_3","H19_3","H23_3")
samples_cortex<- c("F3_1", "F7_1","F11_1","F15_1","F19_1","F23_1", "F3_2", "F7_2","F11_2","F15_2","F19_2","F23_2", "F3_3", "F7_3","F11_3","F15_3","F19_3","F23_3")


circ_all_hippo <- readRDS("./data/hipppocampus_circ_counts.rds")

rc<-read.table("./data/hippo_mapped_stats.tsv")

colnames(rc)<-c("sample_id", "read_count")

avgLibSize<-mean(rc$read_count)
rc$read_count_norm = avgLibSize/rc$read_count
rc <-tbl_df(rc)

hippo_reproducible_circs <- circ_all_hippo %>% dplyr::select(circRNA_ID, time) %>% group_by(circRNA_ID, time) %>% count() %>% filter(n>1)
hippo_reproducible_circs <- hippo_reproducible_circs$circRNA_ID %>% unique()
circ_all_hippo_repr <- circ_all_hippo %>% filter(circRNA_ID %in% hippo_reproducible_circs)

circbase_hippo <- read.table("./data/circbase_hippo_mm10.bed") %>% tbl_df() %>% dplyr::select(V1:V4, V10, V6)

colnames(circbase_hippo) <- c("chr", "start", "end", "id", "read_count", "strand")
circbase_hippo <- circbase_hippo %>% mutate(circRNA_ID=paste(chr,":",start+1, "|", end, sep=""))

circ_all_cortex <- readRDS("./data/cortex_circ_counts.rds")

#rc_cortex<-read.table("~/projects/2020-05-22_XXXXXXXXXXXXXXX/circs_sge/cortex_total_reads")
rc_cortex<-read.table("./data/cortex_mapped_stats.tsv")

colnames(rc_cortex)<-c("sample_id", "read_count")

avgLibSize_cortex<-mean(rc_cortex$read_count)

rc_cortex$read_count_norm = avgLibSize_cortex/rc_cortex$read_count
rc_cortex <-tbl_df(rc_cortex)



cortex_reproducible_circs <- circ_all_cortex %>% dplyr::select(circRNA_ID, time) %>% group_by(circRNA_ID, time) %>% count() %>% filter(n>1)
cortex_reproducible_circs <- cortex_reproducible_circs$circRNA_ID %>% unique()
circ_all_cortex_repr <- circ_all_cortex %>% filter(circRNA_ID %in% cortex_reproducible_circs)

circbase_cortex <- read.table("./data/circbase_cortex_mm10.bed") %>% tbl_df() %>% dplyr::select(V1:V4, V10, V6)

colnames(circbase_cortex) <- c("chr", "start", "end", "id", "read_count", "strand")
circbase_cortex <- circbase_cortex %>% mutate(circRNA_ID=paste(chr,":",start+1, "|", end, sep=""))



plot.cortex_circular_reads.1  <- circ_all_cortex %>% group_by(sample_id) %>% summarize(AllHeadToTailReads=sum(X.junction_reads)) %>% 
  left_join(rc_cortex, by=c("sample_id")) %>% mutate(AllHeadToTail_norm=AllHeadToTailReads * read_count_norm) %>%
  mutate(time=gsub("\\_.*","",sample_id)) 
  plot.cortex_circular_reads.2  <-  plot.cortex_circular_reads.1 %>% filter(time=="F7") %>% mutate(time="F7.1")
  
  plot.cortex_circular_reads.data <- bind_rows(plot.cortex_circular_reads.1, plot.cortex_circular_reads.2)
  
 plot.cortex_circular_reads <- plot.cortex_circular_reads.data %>% ggline(x="time", y="AllHeadToTail_norm", add = c("mean_sd","jitter"), shape = 5, color="darkred", palette = "jco", point.size = 3) + annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("F7", "F11", "F15", "F19", "F23", "F3", "F7.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of head-to-tail reads")+xlab("time of the day") + ggtitle("CircRNAs in cortex")  +
  theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10)) + ylim(0,30000)  

  
  
  plot.cortex_cdr1.1 <- circ_all_cortex %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% left_join(rc_cortex)  %>% mutate(htt_reads=X.junction_reads * read_count_norm)
  plot.cortex_cdr1.2 <-  plot.cortex_cdr1.1 %>% filter(time=="F7") %>% mutate(time="F7.1")
  
  
  plot.cortex_cdr1.data <- bind_rows(plot.cortex_cdr1.1, plot.cortex_cdr1.2)
  
plot.cortex_cdr1 <- plot.cortex_cdr1.data %>% ggline(x="time", y="htt_reads", add = c("mean_sd", "jitter"), col="darkred", shape = 5, point.size=3) + guides(fill=FALSE) + scale_x_discrete(limits=c("F7", "F11", "F15", "F19", "F23","F3", "F7.1")) + 
  annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5) + ylim(0,2300)+
  theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of head-to-tail reads")+xlab("time of the day") + ggtitle("Cdr1as in cortex")  +
  theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))  


circ_cortex_expr<-circ_all_cortex %>% dplyr::select(circRNA_ID, X.junction_reads, sample_id) %>% tidyr::spread(sample_id, X.junction_reads) %>% 
  imputeTS::na_replace(0) %>% tidyr::gather( sample_id, counts,  F11_1:F7_3) %>% left_join(rc_cortex) %>% mutate(norm.read_count=counts*read_count_norm) %>% mutate(time=gsub("\\_.*","",sample_id)) %>% 
  group_by(circRNA_ID, time) %>% summarize(norm.read_count = mean(norm.read_count)) %>% tidyr::spread(time, norm.read_count) %>% filter(grepl("chr", circRNA_ID))

  
 plot.hippo_circular_reads.1 <- circ_all_hippo %>% group_by(sample_id) %>% summarize(AllHeadToTailReads=sum(X.junction_reads)) %>% 
  left_join(rc, by=c("sample_id")) %>% mutate(AllHeadToTail_norm=AllHeadToTailReads * read_count_norm) %>%
  mutate(time=gsub("\\_.*","",sample_id)) 
  
  plot.hippo_circular_reads.2<-plot.hippo_circular_reads.1 %>% filter(time=="H7") %>% mutate(time="H7.1")
  
  plot.hippo_circular_reads.data <- bind_rows(plot.hippo_circular_reads.1, plot.hippo_circular_reads.2)
  
plot.hippo_circular_reads <- plot.hippo_circular_reads.data %>% ggline(x="time", y="AllHeadToTail_norm", add = c("mean_sd","jitter"), shape = 5, color="darkred", palette = "jco", point.size = 3) + annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("H7", "H11", "H15", "H19", "H23", "H3", "H7.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of head-to-tail reads")+xlab("time of the day") + ggtitle("circRNAs in hippocampus")  + theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))  +ylim(0,30000)

  #ggplot(aes(x=time, y=AllHeadToTail_norm, fill=time)) +  geom_boxplot() + guides(fill=FALSE)

plot.hippo_cdr1.1 <- circ_all_hippo %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% left_join(rc)  %>% mutate(htt_reads=X.junction_reads * read_count_norm)
  
plot.hippo_cdr1.2 <- plot.hippo_cdr1.1 %>% filter(time=="H7") %>% mutate(time="H7.1") 
  
  plot.hippo_cdr1.data <- bind_rows(plot.hippo_cdr1.1, plot.hippo_cdr1.2)
  
plot.hippo_cdr1<- plot.hippo_cdr1.data %>%  ggpubr::ggline(x="time", y="htt_reads", add = c("mean_sd", "jitter"), col="darkred", shape = 5, point.size=3) + guides(fill=FALSE) + scale_x_discrete(limits=c("H7", "H11", "H15", "H19", "H23","H3", "H7.1")) + 
  annotate("rect", xmin = 3.75, xmax = 6.75, ymin = -Inf, ymax = Inf, alpha=0.5)+ylim(0,2000) +
  theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of head-to-tail reads")+xlab("time of the day") + ggtitle("Cdr1as in hippocampus")  +
  theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))

  
circ_hippo_expr<-circ_all_hippo %>% dplyr::select(circRNA_ID, X.junction_reads, sample_id) %>% tidyr::spread(sample_id, X.junction_reads) %>% 
  imputeTS::na_replace(0) %>% tidyr::gather( sample_id, counts,  H11_1:H7_3) %>% left_join(rc) %>% mutate(norm.read_count=counts*read_count_norm) %>% mutate(time=gsub("\\_.*","",sample_id)) %>% 
  group_by(circRNA_ID, time) %>% summarize(norm.read_count = mean(norm.read_count)) %>% tidyr::spread(time, norm.read_count) %>% filter(grepl("chr", circRNA_ID))
  

  
  
 plot.hippo_circular_reads
 plot.cortex_circular_reads
 plot.hippo_cdr1
 plot.cortex_cdr1
  

  
  write.table(circ_all_cortex, "cortex_circRNAs.tsv", sep="\t", quote=F, col.names=T, row.names=F)
  write.table(circ_all_hippo, "hippocampus_circRNAs.tsv", sep="\t", quote=F, col.names=T, row.names=F)
  

##################################### supplement linear vs circular in hippo and cortex
lin_all_hippo_geo<-list()

for (i in samples_hippo){
  sample_lin<-paste("~/projects/2020-05-22_Andranik_Ivanov_Microglia_Analysis/results_2020_05_22/mapping/feature_counts/", i, ".all_mates/out/feature_counts.", i, ".all_mates.feature_counts.jcounts", sep="")
  lin_all_hippo_geo[[i]] <- read.table(sample_lin, header=T, sep='\t', comment.char = "") %>% tbl_df()
  lin_all_hippo_geo[[i]] <- lin_all_hippo_geo[[i]] %>% mutate(sample_id=i) %>% mutate(time=gsub("\\_.*","",i))
    colnames(lin_all_hippo_geo[[i]])[10]= "counts"
    colnames(lin_all_hippo_geo[[i]])[10]= "time"
    colnames(lin_all_hippo_geo[[i]])[11]= "sample_id"
    
  }


lin_all_hippo_geo <- lin_all_hippo_geo %>% bind_rows()
lin_all_hippo_geo <- filter(lin_all_hippo_geo, !is.na(gene_id))


lin_all_hippo<-list()

for (i in samples_hippo){
  sample_lin<-paste("~/projects/2020-05-22_Andranik_Ivanov_Microglia_Analysis/results_2020_05_22/mapping/feature_counts/", i, ".all_mates/out/feature_counts.", i, ".all_mates.feature_counts.jcounts", sep="")
  lin_all_hippo[[i]] <- read.table(sample_lin, header=T, sep='\t', comment.char = "") %>% tbl_df() %>%
    dplyr::select(-SecondaryGenes, -Site1_chr, -Site1_location, -Site1_strand, -Site2_chr, -Site2_location, -Site2_strand)
  colnames(lin_all_hippo[[i]])=c("gene_id", "counts")
  lin_all_hippo[[i]] <- lin_all_hippo[[i]] %>% mutate(sample_id=i) %>% mutate(time=gsub("\\_.*","",i))
}


lin_all_hippo <- lin_all_hippo %>% bind_rows()
lin_all_hippo <- filter(lin_all_hippo, !is.na(gene_id))

lin_all_hippo%>% group_by(sample_id) %>% summarize(linJ.counts=sum(counts)) %>% 
  left_join(rc, by=c("sample_id")) %>% mutate(linJ.counts=linJ.counts * read_count_norm) %>%
  mutate(time=gsub("\\_.*","",sample_id)) %>%
  ggline(x="time", y="linJ.counts", add = c("mean_sd","jitter"), shape = 5, color="darkred", palette = "jco", point.size = 3) + 
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("H11", "H15", "H19", "H23", "H3", "H7")) +
  theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of linear junction reads")+xlab("time of the day") + ggtitle("Linear exon junction reads\n per sample in Hippocampus")  +
  theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))  

lin_hippo_expr<-lin_all_hippo %>% left_join(rc) %>% mutate(norm_counts=counts*read_count_norm) %>% group_by(time, gene_id) %>% summarize(counts_sum=sum(norm_counts)) %>% tidyr::spread(time, counts_sum) 




lin_all_cortex<-list()

for (i in samples_cortex){
  sample_lin<-paste("~/projects/2020-05-22_Andranik_Ivanov_Microglia_Analysis/results_2020_05_22/mapping/feature_counts/", i, ".all_mates/out/feature_counts.", i, ".all_mates.feature_counts.jcounts", sep="")
  lin_all_cortex[[i]] <- read.table(sample_lin, header=T, sep='\t', comment.char = "") %>% tbl_df() %>%
    dplyr::select(-SecondaryGenes, -Site1_chr, -Site1_location, -Site1_strand, -Site2_chr, -Site2_location, -Site2_strand)
  colnames(lin_all_cortex[[i]])=c("gene_id", "counts")
  lin_all_cortex[[i]] <- lin_all_cortex[[i]] %>% mutate(sample_id=i) %>% mutate(time=gsub("\\_.*","",i))
}


lin_all_cortex <- lin_all_cortex %>% bind_rows()
lin_all_cortex <- filter(lin_all_cortex, !is.na(gene_id))

lin_all_cortex%>% group_by(sample_id) %>% summarize(linJ.counts=sum(counts)) %>% 
  left_join(rc_cortex, by=c("sample_id")) %>% mutate(linJ.counts=linJ.counts * read_count_norm) %>%
  mutate(time=gsub("\\_.*","",sample_id)) %>%
  ggline(x="time", y="linJ.counts", add = c("mean_sd","jitter"), shape = 5, color="darkred", palette = "jco", point.size = 3) + 
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha=0.5) + scale_x_discrete(limits=c("F11", "F15", "F19", "F23", "F3", "F7")) +
  theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))+ylab("norm. # of linear junction reads")+xlab("time of the day") + ggtitle("Linear exon junction reads\n per sample in Cortex")  +
  theme(legend.title = element_blank(), plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.position = c(0.25,0.9), legend.text = element_text(size=10))  

lin_cortex_expr<-lin_all_cortex %>% left_join(rc_cortex) %>% mutate(norm_counts=counts*read_count_norm) %>% group_by(time, gene_id) %>% summarize(counts_sum=sum(norm_counts)) %>% tidyr::spread(time, counts_sum) 





circ_hippo_id_conversion<-circ_all_hippo %>% dplyr::select(circRNA_ID, gene_id) %>% mutate(gene_id=sub(",.*", "", gene_id))
circ_cortex_id_conversion<-circ_all_cortex %>% dplyr::select(circRNA_ID, gene_id) %>% mutate(gene_id=sub(",.*", "", gene_id))


colnames(lin_cortex_expr)<-c("gene_id", "F11_l", "F15_l", "F19_l", "F23_l", "F3_l", "F7_l")
cortex_lin_circ_expr <- circ_cortex_expr %>% left_join(circ_cortex_id_conversion %>% unique()) %>% left_join(lin_cortex_expr, by=c("gene_id"))



colnames(lin_hippo_expr)<-c("gene_id", "H11_l", "H15_l", "H19_l", "H23_l", "H3_l", "H7_l")
hippo_lin_circ_expr <- circ_hippo_expr %>% left_join(circ_hippo_id_conversion %>% unique()) %>% left_join(lin_hippo_expr, by=c("gene_id"))



hippo_19_15 <- hippo_lin_circ_expr %>% filter(H15+H19>30) %>% ggplot(aes(x=log2(H19+5)-log2(H15+5), log2(H19_l+5)-log2(H15_l+5)))+geom_point() +xlab("circ LFC") + ylab("linear LFC")+
  theme(plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.text = element_text(size=10)) + ggtitle("Hippocampus: 15:00 vs 19:00") +geom_vline(xintercept = 0) + geom_hline(yintercept =0)+xlim(-2,2)+ylim(-2,2)

cortex_23_19 <- cortex_lin_circ_expr %>% filter(F23+F19>30) %>% ggplot(aes(x=log2(F23+5)-log2(F19+5), log2(F23_l+5)-log2(F19_l+5)))+geom_point() +xlab("circ LFC") + ylab("linear LFC")+
  theme(plot.title = element_text(size=15, face="bold", hjust = 0.5), legend.text = element_text(size=10)) + ggtitle("Cortex: 23:00 vs 19:00") +geom_vline(xintercept = 0) + geom_hline(yintercept =0)+xlim(-2,2)+ylim(-2,2)

  
pdf("supplementary_fig_lin_vs_circs.pdf", height=5, width=12)
cowplot::plot_grid(cortex_23_19, hippo_19_15, labels=c("A","B"), ncol=2)
 dev.off()


#################################### Figure 4 



circ_all_fly =list()
fly_samples<-read.table("~/projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/sra.csv")

#fly_samples<-read.table("~/cluster_projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/sra.csv")
colnames(fly_samples)=c("sample", "pheno", "time", "repl")
fly_samples<-fly_samples %>% tbl_df()
fly_samples<-fly_samples %>% mutate(sample=as.character(sample), pheno=as.character(pheno), time=as.character(time), repl=as.character(repl)) %>% mutate(id=paste(pheno, time, repl, sep="_"))

total_reads<-read.table("~/projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/read.number.csv")

#total_reads<-read.table("~/cluster_projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/read.number.csv")
colnames(total_reads)<-c("sample", "total_reads")
total_reads$total_reads<-total_reads$total_reads/4
fly_samples <- fly_samples %>% left_join(total_reads)
avg.lib.size<-mean(fly_samples$total_reads)

fly_samples <-mutate(fly_samples, norm_factor=total_reads/avg.lib.size)

for (i in fly_samples$sample){
  sample_circs<-paste("~/projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/sea_snap_results_2020_11_10/mapping/ciri/", i, ".all_mates", "/out/ciri.", i, ".all_mates.tsv", sep="")
  circ_all_fly[[i]] <- read.table(sample_circs, header=T, sep='\t', comment.char = "") %>% tbl_df() %>% mutate(sample_id=i) %>% mutate(sample=i) %>% left_join(fly_samples) #%>% mutate(gene_id=gsub(pattern = "\\,",replacement = "", x = gene_id))

}

fly_circRNAs <- circ_all_fly %>% bind_rows()


top100_fly <- (fly_circRNAs %>% group_by(circRNA_ID) %>% summarize(sum_count=sum(X.junction_reads)) %>% arrange(-sum_count) %>% head(100))$circRNA_ID


tmp.mean_fly_100 <-fly_circRNAs %>% filter(circRNA_ID %in% top100_fly)  %>% mutate(counts_norm=X.junction_reads/norm_factor) %>% group_by(circRNA_ID) %>% summarize(tmp.mean = mean(counts_norm)) %>% ungroup()


pdf("jk.pdf")

#fly_circRNAs %>% filter(circRNA_ID %in% top100_fly)  %>% group_by(circRNA_ID, time) %>% summarize(group.mean=mean(X.junction_reads/norm_factor))%>% left_join(tmp.mean_fly_100) %>% mutate(relative_expr = log2((group.mean+5)/(tmp.mean+5))) %>% ggline(x="time", y="relative_expr", group="circRNA_ID", palette = "jco", point.size = 0.2, color="darkred") + scale_x_discrete(limits=c("ZT0", "ZT6", "ZT12", "ZT18"))

dev.off()


fly_circRNAs %>% mutate(id2=paste(pheno, time, sep="_")) %>% ggplot(aes(y=junction_reads_ratio, x=id2))+geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

fly_circRNAs_jr <- fly_circRNAs %>% dplyr::select(circRNA_ID, X.junction_reads, id) 
fly_circRNAs_jr$X.junction_reads <- fly_circRNAs_jr$X.junction_reads %>% tidyr::replace_na( 0)
fly_circRNAs_jr <- fly_circRNAs_jr %>% spread(id, X.junction_reads) %>% mutate_each(funs(replace(., is.na(.), 0)))
fly_circRNAs_jr_norm <- fly_circRNAs_jr %>% gather (id, counts,CS_ZT0_rep1:PER_ZT6_rep2) %>% left_join(fly_samples, by=c("id")) %>% mutate(counts_norm=counts/norm_factor) %>% dplyr::select(circRNA_ID, id, counts_norm) %>% spread( id, counts_norm)


fly_circRNAs_rat <- fly_circRNAs %>% dplyr::select(circRNA_ID, junction_reads_ratio, id) 
fly_circRNAs_rat$junction_reads_ratio <- fly_circRNAs_rat$junction_reads_ratio %>% replace_na( 0)
fly_circRNAs_rat <- fly_circRNAs_rat %>% spread(id, junction_reads_ratio) %>% mutate_each(funs(replace(., is.na(.), 0)))


fly_circRNAs_njr <- fly_circRNAs %>% dplyr::select(circRNA_ID, X.non_junction_reads, id) 
fly_circRNAs_njr$X.non_junction_reads <- fly_circRNAs_njr$X.non_junction_reads %>% replace_na( 0)
fly_circRNAs_njr <- fly_circRNAs_njr %>% spread(id, X.non_junction_reads) %>% mutate_each(funs(replace(., is.na(.), 0)))
fly_circRNAs_njr_norm <- fly_circRNAs_njr %>% gather (id, counts, CS_ZT0_rep1:PER_ZT6_rep2) %>% left_join(fly_samples, by=c("id")) %>% mutate(counts_norm=counts/norm_factor) %>% dplyr::select(circRNA_ID, id, counts_norm) %>% spread( id, counts_norm)




fly_circRNAs_jr_norm_repr <- fly_circRNAs_jr_norm %>% filter((CS_ZT0_rep1>1 & CS_ZT0_rep2>1) | (CS_ZT12_rep1>1 & CS_ZT12_rep2>1)  | (CS_ZT18_rep1>1 & CS_ZT18_rep2>1) | (CS_ZT6_rep1>1 & CS_ZT6_rep2>1) | (PER_ZT0_rep1>1 & PER_ZT0_rep2>1) | (PER_ZT12_rep1>1 & PER_ZT12_rep2>1)  | (PER_ZT18_rep1>1 & PER_ZT18_rep2>1) | (PER_ZT6_rep1>1 & PER_ZT6_rep2>1)) 

#pdf("jk.pdf")

fly_LD <- fly_circRNAs_jr_norm_repr %>% ggplot(aes(log2((CS_ZT0_rep1+CS_ZT0_rep2)/2), log2((CS_ZT18_rep1+CS_ZT18_rep2)/2)))+geom_point()+xlim(0,10)+ylim(0,10)+geom_abline(intercept =0 , slope = 1)+geom_abline(intercept =-1 , slope = 1, linetype = 2)+geom_abline(intercept =1 , slope = 1, linetype = 2)+xlab("log2(ZT0)")+ylab("log2(ZT18)") +ggtitle("CircRNAs upon light induction")+theme( axis.text.y=element_text(size=12), axis.text.x=element_text(size=12))

#dev.off()


fly_circRNAs_jr_repr <- fly_circRNAs_jr %>% filter((CS_ZT0_rep1>1 & CS_ZT0_rep2>1) | (CS_ZT12_rep1>1 & CS_ZT12_rep2>1)  | (CS_ZT18_rep1>1 & CS_ZT18_rep2>1) | (CS_ZT6_rep1>1 & CS_ZT6_rep2>1) | (PER_ZT0_rep1>1 & PER_ZT0_rep2>1) | (PER_ZT12_rep1>1 & PER_ZT12_rep2>1)  | (PER_ZT18_rep1>1 & PER_ZT18_rep2>1) | (PER_ZT6_rep1>1 & PER_ZT6_rep2>1)) 

fly_circRNAs_jr_repr %>% gather(sample, counts, CS_ZT0_rep1:PER_ZT6_rep2) %>% ggplot(aes(sample, counts))+geom_boxplot()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

lin_all_fly<-list()

for (i in fly_samples$sample){
  sample_lin<-paste("~/projects/2019-01-25_Andranik_Ivanov/fly_data_hughes_etal2012/sea_snap_results_2020_11_10/mapping/feature_counts/", i, ".all_mates/out/feature_counts.", i, ".all_mates.feature_counts.jcounts", sep="")
  lin_all_fly[[i]] <- read.table(sample_lin, header=T, sep='\t', comment.char = "") %>% tbl_df() %>%
    dplyr::select(-SecondaryGenes, -Site1_chr, -Site1_location, -Site1_strand, -Site2_chr, -Site2_location, -Site2_strand)
  colnames(lin_all_fly[[i]])=c("gene_id", "counts")
  lin_all_fly[[i]] <- lin_all_fly[[i]] %>% mutate(sample_id=i) 
}


lin_all_fly <- lin_all_fly %>% bind_rows()
lin_all_fly <- filter(lin_all_fly, !is.na(gene_id))

lin_fly_expr<-lin_all_fly %>% left_join(fly_samples, by=c("sample_id"="sample")) %>% mutate(norm_counts=counts * norm_factor) %>% group_by(id, gene_id) %>% summarize(counts_sum=sum(norm_counts)) %>% tidyr::spread(id, counts_sum) %>% replace(is.na(.), 0)

circ_fly_id_conversion<-fly_circRNAs %>% dplyr::select(circRNA_ID, gene_id) %>% mutate(gene_id=sub(",.*", "", gene_id))

lin_fly_expr.2 <- lin_fly_expr

colnames(lin_fly_expr.2)<-c("gene_id", paste(colnames(lin_fly_expr)[-1],"_l", sep=""))


fly_lin_circ_expr <- fly_circRNAs_jr_norm_repr %>% left_join(circ_fly_id_conversion %>% unique()) %>% left_join(lin_fly_expr.2, by=c("gene_id"))



fly_lin_circ_expr_avg <- fly_lin_circ_expr %>% mutate(CS_ZT0=(CS_ZT0_rep1+CS_ZT0_rep2)/2, PER_ZT0=(PER_ZT0_rep1+PER_ZT0_rep2)/2, CS_ZT6=(CS_ZT6_rep1+CS_ZT6_rep2)/2, PER_ZT6=(PER_ZT6_rep1+PER_ZT6_rep2)/2, CS_ZT12=(CS_ZT12_rep1+CS_ZT12_rep2)/2, PER_ZT12=(PER_ZT12_rep1+PER_ZT12_rep2)/2, CS_ZT18=(CS_ZT18_rep1+CS_ZT18_rep2)/2, PER_ZT18=(PER_ZT18_rep1+PER_ZT18_rep2)/2, CS_ZT0_l=(CS_ZT0_rep1_l+CS_ZT0_rep2_l)/2, PER_ZT0_l=(PER_ZT0_rep1_l+PER_ZT0_rep2_l)/2, CS_ZT6_l=(CS_ZT6_rep1_l+CS_ZT6_rep2_l)/2, PER_ZT6_l=(PER_ZT6_rep1_l+PER_ZT6_rep2_l)/2, CS_ZT12_l=(CS_ZT12_rep1_l+CS_ZT12_rep2_l)/2, PER_ZT12_l=(PER_ZT12_rep1_l+PER_ZT12_rep2_l)/2, CS_ZT18_l=(CS_ZT18_rep1_l+CS_ZT18_rep2_l)/2, PER_ZT18_l=(PER_ZT18_rep1_l+PER_ZT18_rep2_l)/2)



sumMBLreads <- sum((fly_circRNAs  %>% filter(circRNA_ID == "chr2R:17275410|17276063") %>% dplyr::select(X.junction_reads))$X.junction_reads)
sumPKAC3reads <- sum((fly_circRNAs  %>% filter(circRNA_ID == "chr3L:15946908|15947579") %>% dplyr::select(X.junction_reads))$X.junction_reads)

fly_hist<-fly_circRNAs %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() %>% mutate(Reads.per.circRNA=replace(Reads.per.circRNA, Reads.per.circRNA<10, 10)) %>% ggplot(aes(x=Reads.per.circRNA))+geom_histogram(color="black", fill="white") + scale_y_log10(breaks=c(0,1,5,10,100,1000,5000)) + geom_text(aes(label="Mbl", x=sumMBLreads, y=1 + 2), size=3.5, col="darkred") + geom_segment(aes(x=sumMBLreads, xend=sumMBLreads, y=1 + 1.5, yend=1 + 0.25),arrow=arrow(length=unit(2, "mm")), col="darkred") + xlab("aggregate (all samples) reads") + ylab("# of circRNAs") +ggtitle("Overall circRNA expression") + theme(text = element_text(size=15), plot.title = element_text(size=15, face="bold", hjust = 0.5)) + geom_text(aes(label="PKA-C3", x=sumPKAC3reads, y=1 + 2), size=3.5, col="darkred") + geom_segment(aes(x=sumPKAC3reads, xend=sumPKAC3reads, y=1 + 1.5, yend=1 + 0.25),arrow=arrow(length=unit(2, "mm")), col="darkred") 


mbl_linSS_plot.1<-lin_fly_expr %>% filter(gene_id=="FBgn0265487") %>% gather(id, counts, CS_ZT0_rep1:PER_ZT6_rep2)  %>% left_join(fly_samples, by=c("id")) %>% mutate(normalized.counts=counts/norm_factor)
mbl_linSS_plot.2<-mbl_linSS_plot.1 %>% filter(time=="ZT0") %>% mutate(time="ZT0.1")
mbl_linSS_data <- bind_rows(mbl_linSS_plot.1, mbl_linSS_plot.2)

mbl_linSS_plot <- mbl_linSS_data %>% ggline(x="time", y="normalized.counts", add = c("mean_sd", "jitter"), col="pheno", shape = 5, size=1.2, point.size=5) + annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf, alpha=0.5)  + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, alpha=0.5)+ scale_x_discrete(limits=c("ZT0", "ZT6", "ZT12", "ZT18","ZT0.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ylim(0,2000)



mbl_HtT_plot.1 <- fly_circRNAs_jr_repr %>% filter(circRNA_ID=="chr2R:17275410|17276063") %>% gather(id, counts, CS_ZT0_rep1:PER_ZT6_rep2) %>% left_join(fly_samples, by=c("id")) %>% mutate(normalized.counts=counts/norm_factor)
mbl_HtT_plot.2 <- mbl_HtT_plot.1 %>% filter(time=="ZT0") %>% mutate(time="ZT0.1")
mbl_HtT_data <- bind_rows(mbl_HtT_plot.1, mbl_HtT_plot.2)


mbl_HtT_plot <- mbl_HtT_data %>% ggline(x="time", y="normalized.counts", add = c("mean_sd", "jitter"), col="pheno", shape = 5, size=1.2, point.size=5) + annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf, alpha=0.5)  + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, alpha=0.5)+ scale_x_discrete(limits=c("ZT0", "ZT6", "ZT12", "ZT18","ZT0.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ylim(0,1000)


mbl_fig<-cowplot::plot_grid(mbl_HtT_plot, mbl_linSS_plot, labels=c("circular junctions", "linear junctions"), ncol=1)



pka_linSS_plot.1 <- lin_fly_expr %>% filter(gene_id=="FBgn0000489") %>% gather(id, counts, CS_ZT0_rep1:PER_ZT6_rep2)  %>% left_join(fly_samples, by=c("id")) %>% mutate(normalized.counts=counts/norm_factor)
pka_linSS_plot.2 <- pka_linSS_plot.1 %>% filter(time=="ZT0") %>% mutate(time="ZT0.1")

pka_linSS_data<- bind_rows(pka_linSS_plot.1, pka_linSS_plot.2)


pka_linSS_plot <- pka_linSS_data %>% ggline(x="time", y="normalized.counts", add = c("mean_sd", "jitter"), col="pheno", shape = 5, size=1.2, point.size=5) + annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf, alpha=0.5)  + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, alpha=0.5)+ scale_x_discrete(limits=c("ZT0", "ZT6", "ZT12", "ZT18", "ZT0.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ylim(0,2000)

pka_HtT_plot.1 <- fly_circRNAs_jr_repr %>% filter(circRNA_ID=="chr3L:15946908|15947579") %>% gather(id, counts, CS_ZT0_rep1:PER_ZT6_rep2) %>% left_join(fly_samples, by=c("id")) %>% mutate(normalized.counts=counts/norm_factor)
pka_HtT_plot.2 <- pka_HtT_plot.1 %>% filter(time=="ZT0") %>% mutate(time="ZT0.1")

pka_HtT_data<- bind_rows(pka_HtT_plot.1, pka_HtT_plot.2)

pka_HtT_plot<- pka_HtT_data %>% ggline(x="time", y="normalized.counts", add = c("mean_sd", "jitter"), col="pheno", shape = 5, size=1.2, point.size=5) + annotate("rect", xmin = 3, xmax = 5, ymin = -Inf, ymax = Inf, alpha=0.5)  + annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, alpha=0.5)+ scale_x_discrete(limits=c("ZT0", "ZT6", "ZT12", "ZT18","ZT0.1")) + theme(text = element_text(size=15), axis.text.y=element_text(size=12), axis.text.x=element_text(size=12), plot.title = element_text(size=15, face="bold", hjust = 0.5))+ylim(0,1000)

 
pka_c3_fig<-cowplot::plot_grid(pka_HtT_plot, pka_linSS_plot, labels=c("circular junctions", "linear junctions"), ncol=1)


fly_low=cowplot::plot_grid(mbl_fig, pka_c3_fig, ncol=2, labels = c("C", "D"), label_x = 0.7, label_y = 1)

pdf("mbl_figure5.pdf", height=13, width=9)

fly_up<-cowplot::plot_grid(fly_hist, fly_LD, ncol=2, labels=c("A","B"))

cowplot::plot_grid(fly_up, fly_low, ncol=1, rel_heights = c(1,1.7))

dev.off()

#add(mbl, jp, mamo, Cadps, Nrg)

######################################## Microglia figure



aln_stats_mg <- read.table("~/projects/2019-01-25_Andranik_Ivanov/sea_snap_mapping/mapping_stats.tsv", sep='\t', header=T)%>% tibble::as_tibble() %>% left_join(covariate_m, by=c("sample_id"="filename")) %>% mutate(sample=label) %>% dplyr::select(sample, key, value)



mg_reads_all <- aln_stats_mg %>% filter(key=="Number of input reads") %>% mutate(mean.v=mean(value)) %>% mutate(norm_factor=mean.v/value)
mg_reads_mu <- aln_stats_mg %>% filter(key=="Uniquely mapped reads number")  %>% mutate(mean.v=mean(value)) %>% mutate(norm_factor=mean.v/value)



covariate_m <-read.table("~/projects/2019-01-25_Andranik_Ivanov/sea_snap_mapping/covariate_file.txt", sep='\t', header=T) %>% tbl_df() %>% mutate(filename=sub("mapping/feature_counts/","", filename))%>% mutate(filename=sub(".all_mates.*","", filename))


cova_h<-covariate_m %>% filter(brain_region=="Hippocampus")
cova_f<-covariate_m %>% filter(brain_region=="Frontal_cortex")




microglia_circ_hippo <-list()


for (i in cova_h$filename){
#  print (i)
  sample <-paste("~/projects/2019-01-25_Andranik_Ivanov/sea_snap_mapping/mapping/ciri/", i, ".all_mates/out/ciri.", i, ".all_mates.tsv", sep="")
  if(file.exists(sample)){
  microglia_circ_hippo[[i]] <- read.table(sample, header=T, sep='\t', comment.char = "") %>% tbl_df() %>% mutate(sample_id=i) %>% left_join(cova_h, by=c("sample_id"="filename")) %>% mutate(sample_id=label)  %>% dplyr::select(-group, -label, -junction_reads_ID) 

  }else{
    print (sample)
    }
}

microglia_circ_hippo <- microglia_circ_hippo %>%  bind_rows() 


microglia_circ_cortex <- microglia_circ_cortex %>%  bind_rows() 




rpc_hippo <-microglia_circ_hippo %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() 
max_rpc_hippo<-max(rpc_hippo$Reads.per.circRNA)

sumCDR1reads_hippo <- (rpc_hippo   %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% dplyr::select(Reads.per.circRNA))$Reads.per.circRNA


hippoM_cdr1 <- microglia_circ_hippo %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() %>% mutate(Reads.per.circRNA=replace(Reads.per.circRNA, Reads.per.circRNA<10, 10)) %>% ggplot(aes(x=Reads.per.circRNA))+geom_histogram(color="black", fill="white")+ scale_x_log10(breaks=c(0,1,5,10,100,1000,5000))+ scale_y_log10(breaks=c(0,1,10,100,1000,5000)) + xlab("aggregate reads") + ylab("# of circRNAs") +ggtitle("Overall circRNA expression") + theme(text = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5)) # + geom_text(aes(label="Cdr1as", x=sumCDR1reads_hippo, y=1 + 2), size=3.5, col="darkred") + geom_segment(aes(x=sumCDR1reads_hippo, xend=sumCDR1reads_hippo, y=1 + 1.5, yend=1 + 0.25),arrow=arrow(length=unit(2, "mm")), col="darkred")





microglia_circ_cortex <-list()


for (i in cova_f$filename){
#  print (i)
  sample <-paste("~/projects/2019-01-25_Andranik_Ivanov/sea_snap_mapping/mapping/ciri/", i, ".all_mates/out/ciri.", i, ".all_mates.tsv", sep="")
  if(file.exists(sample)){
  microglia_circ_cortex[[i]] <- read.table(sample, header=T, sep='\t', comment.char = "") %>% tbl_df() %>% mutate(sample_id=i) %>% left_join(cova_f, by=c("sample_id"="filename")) %>% mutate(sample_id=label)  %>% dplyr::select(-group, -label, -junction_reads_ID) 

  }else{
    print (sample)
    }
}

microglia_circ_cortex <- microglia_circ_cortex %>%  bind_rows() 




rpc_cortex <-microglia_circ_cortex %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() 
max_rpc_cortex<-max(rpc_cortex$Reads.per.circRNA)

sumCDR1reads_cortex <- (rpc_cortex   %>% filter(circRNA_ID == "chrX:61183248|61186174") %>% dplyr::select(Reads.per.circRNA))$Reads.per.circRNA


cortexM_cdr1 <- microglia_circ_cortex %>% group_by(circRNA_ID) %>% summarise(Reads.per.circRNA=sum(X.junction_reads)) %>% ungroup() %>% mutate(Reads.per.circRNA=replace(Reads.per.circRNA, Reads.per.circRNA<10, 10)) %>% ggplot(aes(x=Reads.per.circRNA))+geom_histogram(color="black", fill="white")+ scale_x_log10(breaks=c(0,1,5,10,100,1000,5000))+ scale_y_log10(breaks=c(0,1,10,100,1000,5000))  + xlab("aggregate reads") + ylab("# of circRNAs") +ggtitle("Overall circRNA expression") + theme(text = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5)) #+ geom_text(aes(label="Cdr1as", x=sumCDR1reads_cortex, y=1 + 2), size=3.5, col="darkred") + geom_segment(aes(x=sumCDR1reads_cortex, xend=sumCDR1reads_cortex, y=1 + 1.5, yend=1 + 0.25),arrow=arrow(length=unit(2, "mm")), col="darkred")




microglia_circ <-list()


for (i in covariate_m$filename){
#  print (i)
  sample <-paste("~/projects/2019-01-25_Andranik_Ivanov/sea_snap_mapping/mapping/ciri/", i, ".all_mates/out/ciri.", i, ".all_mates.tsv", sep="")
  if(file.exists(sample)){
  microglia_circ[[i]] <- read.table(sample, header=T, sep='\t', comment.char = "") %>% tbl_df() %>% mutate(sample_id=i) %>% left_join(covariate_m, by=c("sample_id"="filename")) %>% mutate(sample_id=label)  %>% dplyr::select(-group, -label, -junction_reads_ID) 

  }else{
    print (sample)
    }
}

microglia_circ <- microglia_circ %>%  bind_rows() 




cdr1_asM_hippo <-microglia_circ %>% filter(circRNA_ID == "chrX:61183248|61186174")  %>% left_join(mg_reads_all, by=c("sample_id"="sample")) %>% filter(brain_region=="Hippocampus") %>% mutate(norm.counts= X.junction_reads*norm_factor)%>% ggplot(aes(x=time, y=norm.counts))+geom_boxplot(fill=c("white", "darkgrey"))  + theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5), axis.title.y = element_text(size =15), axis.title.x = element_text(size =15)) +ylab("normalized read counts")+ggtitle("Cdr1as") 


cdr1_asM_cortex <-microglia_circ %>% filter(circRNA_ID == "chrX:61183248|61186174")  %>% left_join(mg_reads_all, by=c("sample_id"="sample"))%>% filter(brain_region=="Frontal_cortex") %>% mutate(norm.counts= X.junction_reads*norm_factor)%>% ggplot(aes(x=time, y=norm.counts))+geom_boxplot(fill=c("white", "darkgrey"))  + theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5), axis.title.y = element_text(size =15), axis.title.x = element_text(size =15)) +ylab("normalized read counts")+ggtitle("Cdr1as") 



all_circs_hippoM <- microglia_circ %>% filter(brain_region=="Hippocampus")  %>% group_by(sample_id) %>% summarize(AllHeadToTailReads=sum(X.junction_reads)) %>% left_join(mg_reads_all, by=c("sample_id"="sample"))%>% mutate(norm.counts= AllHeadToTailReads*norm_factor) %>% left_join(covariate_m, by=c("sample_id"="label"))  %>% ggplot(aes(x=time, y=norm.counts, fill=time)) +  geom_boxplot(fill=c("white", "darkgrey")) + guides(fill=FALSE) + theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5), axis.title.y = element_text(size =15), axis.title.x = element_text(size =15)) +ylab("normalized read counts")+ggtitle("All reads") 



all_circs_cortexM <- microglia_circ %>% filter(brain_region=="Frontal_cortex")  %>% group_by(sample_id) %>% summarize(AllHeadToTailReads=sum(X.junction_reads)) %>% left_join(mg_reads_all, by=c("sample_id"="sample"))%>% mutate(norm.counts= AllHeadToTailReads*norm_factor) %>% left_join(covariate_m, by=c("sample_id"="label"))  %>% ggplot(aes(x=time, y=norm.counts, fill=time)) +  geom_boxplot(fill=c("white", "darkgrey")) + guides(fill=FALSE) + theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), plot.title = element_text(size=15, hjust = 0.5), axis.title.y = element_text(size =15), axis.title.x = element_text(size =15)) +ylab("normalized read counts")+ggtitle("All reads") 






upper_panel=cowplot::plot_grid(cortexM_cdr1, all_circs_hippoM, cdr1_asM_cortex, rel_widths = c(1.5, 0.5, 0.5), ncol=3, labels=c("D","E","F"))
lower_panel=cowplot::plot_grid(hippoM_cdr1, all_circs_cortexM, cdr1_asM_hippo, rel_widths = c(1.5, 0.5, 0.5), ncol=3)




pdf("microglia_Figure4.pdf", height=6, width=9)
cowplot::plot_grid(upper_panel, lower_panel, ncol=1, labels = c("Cortex", "Hippocampus"))
dev.off()

lower_new3 <- cowplot::plot_grid(upper_panel, lower_panel, ncol=1)

upper_new3 <- cowplot::plot_grid(plotlist=list(h_ie, f_ie, plot.hippo_circular_reads, plot.cortex_circular_reads, plot.hippo_cdr1, plot.cortex_cdr1), ncol=2, labels=c("A","","B","", "C",""))
 



pdf("Fig3_new.pdf", height=15, width=10)

cowplot::plot_grid(plotlist=list(upper_new3, lower_new3), ncol=1, rel_heights=c(1.2,0.5))

dev.off()



