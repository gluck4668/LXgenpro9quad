# The nine quad of genes and proteins

#---------------------------------------------------------------------

LXgenpro9quad_p1 <- function(gene_data,protein_data,species){

#list all the packages that have been installed
all_packages <- data.frame(installed.packages())

#To judge whether a package was installed. If not, it will be installed.
pack <- data.frame(c("devtools","BiocManager","openxlsx","dplyr","psych",
                     "ggplot2","pheatmap", "ggrepel","pak") )

# To judge whether a package was included in the all_packages: %in%
pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){
  if(pack[i,2]==FALSE)
    install.packages(pack[i,1],update = F,ask = F)
}
rm(i)

# 批量library
packages <- as.character(pack[,1])

for(i in packages){
  library(i, character.only = T)
}
rm(i)

#-----------------
if("tidyverse" %in% all_packages$Package==FALSE)
  pak::pak("tidyverse/tidyverse")
 library(tidyverse)

#-----------------
BiocManager_pack <- data.frame(c("clusterProfiler",
                                 "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"))
# human: "org.Hs.eg.db"
# mouse: "org.Mm.eg.db"
# rat: "org.Rn.eg.db"

BiocManager_pack$type <- BiocManager_pack[,1] %in% all_packages$Package

for (i in 1:nrow(BiocManager_pack)){
  if(BiocManager_pack[i,2]==FALSE)
    BiocManager::install(BiocManager_pack[i,1],update = F,ask = F)
}

# 批量library
Bio_packages <- as.character(BiocManager_pack[,1])
for(i in Bio_packages){
  library(i, character.only = T)
}
rm(i)

#------处理gene data-------------------------------------------------------
#读入数据；
RNA =read.xlsx(gene_data)
head(RNA)
colnames(RNA)[1] <- c("gene_symbol")
table(duplicated(RNA$gene_symbol))  # 查看gene_symbol是否有重复
RNA <- distinct(RNA,gene_symbol,.keep_all = T) # 去重

# hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
RNA_human <- RNA
RNA_human$gene_symbol <-toupper(RNA$gene_symbol)

# mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
# gene_n <- dplyr::filter(RNA,tolower(str_sub(RNA$gene_symbol,1,1)) %in% 0:9) # 数字开头的基因
# gene_L <- dplyr::filter(RNA,tolower(str_sub(RNA$gene_symbol,1,1)) %in% letters) # 字母开头的基因
#gene_L$gene_symbol <- str_to_title(tolower(gene_L$gene_symbol)) # str_to_title()首字母大写
#RNA_animal <-rbind(gene_n,gene_L)

RNA_animal <- RNA

#--------------------------------------------------------------------------
spe <- trimws(species) %>% tolower()

if (spe=="human"){
                 keytypes(org.Hs.eg.db) # to check the key items
                 geneID<-bitr(RNA_human$gene_symbol,
                 fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
                 OrgDb = org.Hs.eg.db)
                }


if (spe=="rat"){
              keytypes(org.Rn.eg.db) # to check the key items
              geneID<-bitr(RNA_animal$gene_symbol,
              fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
              OrgDb = org.Rn.eg.db)
              }


if (spe=="mouse"){
                keytypes(org.Mm.eg.db) # to check the key items
                geneID<-bitr(RNA_animal$gene_symbol,
                fromType = "SYMBOL",toType =c("SYMBOL","ENTREZID","ENSEMBL"),
                OrgDb = org.Mm.eg.db)
               }

RNA$ID <- toupper(RNA$gene_symbol)
geneID$ID <- toupper(geneID$SYMBOL)

table(duplicated(geneID$SYMBOL))  # 查看symbol是否有重复
geneID <- distinct(geneID,SYMBOL,.keep_all = T) # 去重

RNA_data <- dplyr::inner_join(RNA,geneID,"ID")
RNA_data <- RNA_data[,-c(4:6)]
RNA_data <- RNA_data[,c(4,1:3)]
colnames(RNA_data) <- c("ENSEMBL","gene_symbol","log2FC_gene","pvalue_gene")

#-----------处理 protein data----------------------------------------

protein=read.xlsx(protein_data)
head(protein)

colnames(protein) <- c("UNIPROT","log2FC","pvalue")
protein$log2FC <- log2(protein$log2FC)

table(duplicated(protein$UNIPROT))  # 查看UNIPROT是否有重复
protein <- distinct(protein,UNIPROT,.keep_all = T) # 去重

if (spe=="human"){
               keytypes(org.Hs.eg.db) # to check the key items
               proteinID<-bitr(protein$UNIPROT,
               fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
               OrgDb = org.Hs.eg.db)
              }

if (spe=="rat"){
              keytypes(org.Rn.eg.db) # to check the key items
              proteinID<-bitr(protein$UNIPROT,
              fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
              OrgDb = org.Rn.eg.db)
             }


if (spe=="mouse"){
              keytypes(org.Mm.eg.db) # to check the key items
              proteinID <-bitr(protein$UNIPROT,
              fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
              OrgDb = org.Mm.eg.db)
             }

if(length(bitr(protein$UNIPROT,
        fromType = "UNIPROT",toType =c("UNIPROT","SYMBOL","ENTREZID","ENSEMBL"),
        OrgDb = org.Mm.eg.db))==0)
    species="mouse"

colnames(proteinID)[1] <- c("UNIPROT")

table(duplicated(proteinID$UNIPROT))  # 查看UNIPROT是否有重复
proteinID <- distinct(proteinID,UNIPROT,.keep_all = T) # 去重

protein_data <- dplyr::inner_join(protein,proteinID,"UNIPROT")
protein_data <- protein_data[,c(6,1:3)]
colnames(protein_data) <- c("ENSEMBL","UNIPROT","log2FC_protein","pvalue_protein")

#合并gene_data和protein_data两个表格
table(duplicated(RNA_data$ENSEMBL))
RNA_data <- distinct(.data = RNA_data,ENSEMBL,.keep_all = T)

table(duplicated(protein_data$ENSEMBL))
protein_data <- distinct(.data = protein_data,ENSEMBL,.keep_all = T)

data <- dplyr::inner_join(RNA_data,protein_data,"ENSEMBL")

data <- na.omit(data)

colnames(data) <- c("ENSEMBL","gene_symbol","log2FC_RNA","pvalue_gene","UNIPROT","log2FC_Protein","pvalue_protein")

data$log2FC_RNA <- as.numeric(data$log2FC_RNA)
data$log2FC_Protein <- as.numeric(data$log2FC_Protein)


#对数据进行分组；
#生成显著上下调数据标签；
data$part <- case_when(abs(data$log2FC_RNA) >= 1 & abs(data$log2FC_Protein) >= 1 ~ "part1379",
                       abs(data$log2FC_RNA) < 1 & abs(data$log2FC_Protein) > 1 ~ "part28",
                       abs(data$log2FC_RNA) > 1 & abs(data$log2FC_Protein) < 1 ~ "part46",
                       abs(data$log2FC_RNA) < 1 & abs(data$log2FC_Protein) < 1 ~ "part5")
head(data)


sig_cor <- dplyr::filter(data,part=="part1379")

sig_cor$correlation <- case_when(sig_cor$log2FC_RNA >= 1 & sig_cor$log2FC_Protein >= 1 ~ "positive",
                                 sig_cor$log2FC_RNA <= 1 & sig_cor$log2FC_Protein <= 1 ~ "positive",
                                 sig_cor$log2FC_RNA >= 1 & sig_cor$log2FC_Protein <= 1 ~ "negative",
                                 sig_cor$log2FC_RNA <= 1 & sig_cor$log2FC_Protein >= 1 ~ "negative")

postive_cor <- dplyr::filter(sig_cor,correlation=="positive")
postive_cor <- postive_cor[,-7]

negative_cor <- dplyr::filter(sig_cor,correlation=="negative")
negative_cor <-negative_cor[,-7]

if(dir.exists("analysis result")==FALSE)
  dir.create("analysis result")

write.xlsx(data,"analysis result/data_cor.xlsx")
write.xlsx(postive_cor,"analysis result/postive_cor.xlsx")
write.xlsx(negative_cor,"analysis result/negative_cor.xlsx")

#开始尝试绘图；

p0 <-ggplot(data,aes(log2FC_RNA,log2FC_Protein,color=part))

p1 <- p0+geom_point(size = 2.5)+guides(color="none")
p1

#自定义半透明颜色；
#mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
mycolor <- c("#ff4500","#006400","#20b2aa","#da70d6")
p2 <- p1 + scale_colour_manual(name="",values=mycolor)
p2

titile_text <- paste('Genes/Proteins Nine Quadrantal Diagram')

#添加辅助线；
p3 <- p2+geom_hline(yintercept = c(-1,1),
                    linewidth = 0.6,
                    color = "gray30",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             linewidth = 0.6,
             color = "gray30",
             lty = "dashed")+
  labs(title =titile_text) +
  theme(plot.title = element_text(hjust = 0.5))
p3

}


#-------------------------------------------------------------------------

LXgenpro9quad_p2 <- function(xmin,xmax,ymin,ymax){

  data <- read.xlsx("analysis result/data_cor.xlsx")

  p0 <-ggplot(data,aes(log2FC_RNA,log2FC_Protein,color=part))

  p1 <- p0+geom_point(size = 2.5)+guides(color="none")
  p1

  #mycolor <- c("#3B5FCD","#00CE00","#CD3333","#778899")
  mycolor <- c("#ff4500","#006400","#20b2aa","#da70d6")
  p2 <- p1 + scale_colour_manual(name="",values=mycolor)
  p2

  titile_text <- paste('Genes/Proteins Nine Quadrantal Diagram')

  p3 <- p2+geom_hline(yintercept = c(-1,1),
                      linewidth = 0.6,
                      color = "gray30",
                      lty = "dashed")+
    geom_vline(xintercept = c(-1,1),
               linewidth = 0.6,
               color = "gray30",
               lty = "dashed")+
    labs(title =titile_text) +
    theme(plot.title = element_text(hjust = 0.5))
  p3

  xmin <- as.numeric(xmin)
  xmax <- as.numeric(xmax)
  ymin <- as.numeric(ymin)
  ymax <- as.numeric(ymax)

  p4<-p3+
    scale_y_continuous(expand=expansion(add = c(0.01, 0.01)),
                       limits = c(ymin, ymax),
                       breaks = c(ymin,ymin/2,0,ymax/2,ymax),
                       label = c(as.character(ymin),as.character(ymin/2),"0",
                                 as.character(ymax/2),as.character(ymax)))+
    scale_x_continuous(expand=expansion(add = c(0.01, 0.01)),
                       limits = c(xmin, xmax),
                       breaks = c(xmin,xmin/2,0,xmax/2,xmax),
                       label = c(as.character(xmin),as.character(xmin/2),"0",
                                 as.character(xmax/2),as.character(xmax)))+
    theme_bw()+theme(panel.grid = element_blank())

  p4

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",face="bold",colour ="black",size =14),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=18))+
    theme(axis.text.y = element_text(face="bold",color="black",size=18))

  cor  <- WGCNA::cor(data$log2FC_RNA,data$log2FC_Protein,use = "complete.obs",method ="spearman")
  cor <- round(cor,4)

  lab = paste("correlation=",cor,sep="")
  lab

  cor_value <- geom_text(x=xmin+1,y=ymax-0.5,label = lab, size=4,color="black")

  p5 <- p4+mytheme+xytheme+
    annotate("text", x=-1+(xmin+1)/2, y=1+(ymax-1)/2, label= "①",size=12)+
    annotate("text", x=-1+(xmin+1)/2, y=0.1, label= "④", size=12)+
    annotate("text", x=-1+(xmin+1)/2, y=-1+(ymin+1)/2, label= "⑦",size=12)+

    annotate("text", x=0, y=1+(ymax-1)/2, label= "②",size=12)+
    annotate("text", x=0, y=0.1, label= "⑤",size=14)+
    annotate("text", x=0, y=-1+(ymin+1)/2, label= "⑧",size=12)+

    annotate("text", x=1+(xmax-1)/2, y=1+(ymax-1)/2, label= "③",size=12)+
    annotate("text", x=1+(xmax-1)/2, y=0.1, label= "⑥",size=14)+
    annotate("text", x=1+(xmax-1)/2, y=-1+(ymin+1)/2, label= "⑨",size=12)


  ggsave("analysis result/The genes-proteins nine quadrantal diagram.png", p5,width=1200, height =900, dpi=150,units = "px")

  print("Please see the results in the folder of <analysis result>")

  p5

}


