setwd("/Users/saman/Dropbox/my_publications/rs1/data/")
setwd("/Users/saman/Dropbox/systems_biology_projects/pso_cmm/")
setwd("~/git/psoriasis1/")


setRepositories()

if (!require("BiocManager"))
  install.packages("BiocManager")

if (!require("simpleaffy"))
  install.packages("simpleaffy")

if(!require("annotate"))
  BiocManager::install("annotate")

if(!require("Biobase"))
  BiocManager::install("Biobase")

BiocManager::install("affy")
BiocManager::install("GEOquery")
BiocManager::install("limma")
#BiocManager::install("gcrma")
BiocManager::install("simpleaffy")
BiocManager::install("hgu133plus2.db")
BiocManager::install("Vennerable")


library(pheatmap) #
library(affy)  
#library(gplots)
library(limma)
library(GEOquery)  #
#library(gcrma)
library(BiocGenerics )
library(parallel) 
library(reshape2) #
library(plyr) 
library(Biobase)
library(simpleaffy)
library(annotate) #
library(AnnotationDbi) #
library(hgu133plus2cdf) #
library("hgu133plus2.db") #

# Based on workshop ----

#fetch the dataset from GEO:
series <- "GSE30999"
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
#if you do it often, put the following line in a file in home directory : .Renviron
VROOM_CONNECTION_SIZE=262144

#Mac
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE,
       destdir ="/Users/saman/Dropbox/systems_biology_projects/pso_cmm/workshop/data" )
#Ububtu
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE,
               destdir ="/mnt/4tb/Dropbox/systems_biology_projects/pso_cmm/workshop/data/" )

typeof(gset) #list
gset$GSE30999_series_matrix.txt.gz@experimentData
names(gset)
grepl("GPL", attr(gset, "names"))
if (length(gset)>1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1  #only one platform for this dataset
idx



gset[1]  #the name of platform in Annotation line: GPL570
gset <- gset[idx]
typeof(gset) #still a list
gset <- gset[[idx]]
View(gset)
class(gset) #"ExpressionSet"
typeof(gset) # S4
#gr <- c(rep("N", 85), rep("L", 85))
#length(gr)

## extract the new expression matrix based on the workshop ----

ex <- exprs(gset)
dim(gset) #54675 170
dim(ex) #54675 170
max(ex) #up to 20, no need for log
min(ex) #no zeros

## Check if the expression data is log2 transformed and do if not already
#if eg qx[5] being 0.99 quantile were above 100 it means the anti-log of that value would be 2^100 which is crazy
#and it shows the data has not been log2 transformed yet, so we do it!
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
qx
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)

ex[1:5, 1:5]

if (LogC) { ex[which(ex <= 0)] <- NaN  # none was <=0
exprs(gset) <- log2(ex) }

#get sample information containing sample groups
colnames(pData(gset))
colnames(fData(gset))
colnames(exprs(gset))
rownames(exprs(gset))
sum(grepl("_a_",rownames(exprs(gset))))

pData(gset)$data_processing[1]  #"The data was normalized by GCRMA."
pData(gset)
sampleInfo <- pData(gset)
rownames(sampleInfo)
sampleInfo$title
head(sampleInfo)
sampleInfo <- sampleInfo[,c("title", "geo_accession")]
sampleInfo[1:5, ]
sampleInfo$group <- ""


#grep("_NL", sampleInfo$title[1])
for(i in 1:nrow(sampleInfo)) {
  if(grepl("_NL",sampleInfo$title[i])) {
    sampleInfo$group[i] <- "N"
  }
  else if(grepl("_LS", sampleInfo$title[i])) {
    sampleInfo$group[i] <- "L"
  }
}
class(sampleInfo$group)
sampleInfo$title

sampleInfo[1:5, ]

### add geneSymbol column ----
prbs <- rownames(ex)
geneSymbol <- getSYMBOL(prbs,'hgu133plus2.db')
ex <- data.frame(ex)
ex$probe <- prbs
rownames(ex) 
ex$geneSymbol <- geneSymbol
ex[2,"geneSymbol"] #RFC2
ex[2, 172]         #RFC2
ex[2,"probe"]

## replace sample IDs with group names
sampleInfo$group <- paste0(sampleInfo$group, 1:170)
colnames(ex)[1:170] <- sampleInfo$group
# ex[1:5, 1:5]
# dim(ex)
# ex$geneSymbol
# colnames(ex)


### if zeros and large values ----
#ex <- log2(ex + 1)
#exprs(gset) <- ex  #necessary after log2 transformation, exprs is a set and get function
pdf("./genes_acrossSamples.pdf", width=170)
boxplot(ex[,1:170])
dev.off()

### if boxplot expression distribution per sample was not normalized: ----
#ex <- normalizeQuantiles(ex)
#exprs(gset) <- ex

### heatmap ----
# tail(colnames(ex))
# head(colnames(ex))
# ex[1:5,1:5]
# excor <- cor(ex)
# excor[1:5, 1:5]
# 
# ex_t <- as.data.frame(t(ex[, 1:(ncol(ex)-2)]))
# ex_t[1:5,ncol(ex_t)]
# dim(ex_t)
# ex_t$group <- sampleInfo$group
# 
# pdf("./corr_heat.pdf")
# #load gplots library to add e.g. color to pheatmap
# pheatmap(cor(ex_t[,-ncol(ex_t)][,1:2000]), labels_row = ex_t$group, labels_col = ex_t$group)
# dev.off()
#samples of one single tumor from the same person contain more epigenetic variations and mutations than all other cells 
#of the same person!!


### pca ----

pc <- prcomp(ex)
dim(pc$x) #54675   170  : 54675 points in a 170-dimention space

##?? try the whole MI network construction after subtracting means
##to deal with genes with low variation between groups: subtract the mean expression of each gene
#-->> each row is minus the mean expression of the gene in that row:
ex.scale <- t(scale(t(ex), scale = FALSE, center = TRUE)) #scale acts on columns -> mean=0 and divided by SD for each column -> z-score
mean(ex.scale[1,]) #almost 0

##now pca of samples:
dim(pc$rotation)
class(pc$rotation) #matrix, I want to add group names and..., if you add a character, all types are coerced to charachter
#but dataframe can harbour several classes in different columns
pcr <- data.frame(pc$rotation[,1:3], Group = gr) #maximum 3D plot
pcr[1:5,1:4]
dim(pcr)

#ggplot only works with dataframes as input data. you tell ggplot which column is used as the target property of plot





#### RS1
#AffyBatch object
#(pso.affybatch format) is a HG-U133_Plus_2 microarray with 54675 affy IDs. 
#Each affy ID stands for one probe set with which you measure the expression
#of genes. In total,
pso.affybatch<-ReadAffy()
dim(pso.affybatch)  #1164  1164
#pData(pso.affybatch) :show ordering of samples
#image(pso.affybatch[,1])
#hist(pso.affybatch)
#exprs(pso.affybatch)
#plot(exprs(pso.affybatch)[,1], exprs(pso.affybatch)[,2])

pData(pso.affybatch)

#for plotting start
################
expression.log<-log2(exprs(pso.affybatch))
expression<-exprs(pso.affybatch)
# set parameters and draw the plot
palette(c("#a2cbfa","#f5cbcb", "#AABBCC"))
dev.new(width=4+dim(pso.affybatch)[[2]]/5, height=7)
title <- paste ("GSE30999", '/', annotation(pso.affybatch), " All log-transfromed Samples", sep ='')
boxplot(expression.log, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
#for plotting finish

#len_neg_exp<-length(exprs(pso.affybatch)[which(expression <= 0)])

write.table(expression, file="exp.csv3",row.names=FALSE,col.names=T, quote=F, sep="\t")


#normalize 
data.mass <- mas5(pso.affybatch)
dim(data.mass)
data.mass.exp<-exprs(data.mass)
dim(data.mass.exp)


write.table(data.mass.exp, file="data_mass_exp3.csv",row.names=TRUE,col.names=T, quote=F, sep="\t")
data.mass.exp <- read.table("../data_mass_exp3.csv", sep = "\t")
dim(data.mass.exp)
data.mass.exp[1:5,1:5]

#?? log transorm of Normalized Expression value  
#??necessary to log2 transform already normalized expression data?
#?? some papers say mas5 output is already log transformed but it some papers they do log2 on mas5 output
qxm <- as.numeric(quantile(data.mass.exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogCm <- (qxm[5] > 100) ||(qxm[6]-qxm[1] > 50 && qxm[2] > 0) ||
  (qxm[2] > 0 && qxm[2] < 1 && qxm[4] > 1 && qxm[4] < 2)
if (LogCm) 
{ data.mass.exp[which(data.mass.exp<= 0)] <- NaN 
 exprs(data.mass) <- log2(data.mass.exp)}




write.table(exprs(data.mass), file="data_mass_exp_log3.csv",row.names=TRUE,col.names=T, quote=F, sep="\t")


bs <- c()

for (i in 1:10) {
  if (i <= 5) {bs <- paste0(bs, "0")}
  else {bs <- paste0(bs, "1")}
}
  
bs


# group names for all samples in a series
gsms <- paste0("00000000000000000000000000000000000000000000000000000000",
               "000000000000000000000000000001111111111111111111111111111",
               "111111111111111111111111111111111111111111111111111111111")
nchar(gsms)

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

sml

#set group names
sml <- paste("G", sml, sep="")  
sml

# order samples by group
fl <- as.factor(sml)
fl
labels <- c("Non-Lesional","Lesional")
fl

# set parameters and draw the plot
palette(c("#a2cbfa","#f5cbcb", "#AABBCC"))
dev.new(width=4+dim(pso.affybatch)[[2]]/5, height=7)
title <- paste ("GSE30999", '/', annotation(pso.affybatch), " All log-transfromed Normalized Samples", sep ='')
boxplot(exprs(data.mass), boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")



#set up the data and proceed the analysis
data.mass$sample
data.mass$description<-fl

# design matrix gives a matrix representation of different groups of the data, each column, one group, we have two
# groups therefore two columns and 170 rows, each column has 1s and 0s to indicate the presence of a datapoint(sample) 
#in each group. +0 is used so that the intercept column is omitted from the design matrix
design<-model.matrix(~description+0 ,data.mass)
design
colnames(design)<-levels(fl)
fit<-lmFit(data.mass,design)
fit

## look at limma package documentation which explains how linear models are used for t-test
cont.matrix<-makeContrasts(G1-G0, levels=design)
cont.matrix
fit1 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit1, 0.01)

top<-topTable(fit2, adjust="fdr", sort.by="B", number=Inf)   #Mobasheri has the 2000 limit  #?? what is "B" column ?
#sum(top$P.Value<0.01)  :  28770  ,  p1<-top$P.Value[top$P.Value<0.01]
#sum(top$adj.P.Val<0.01) :27018, p2<-top$adj.P.Val[top$adj.P.Val<0.01]
#which((top$P.Value<0.01))
#dim(top[which(top$P.Value<0.01  & top$adj.P.Val<0.01),]) : 27018 


##extract geneSymbol
#first approach:
feature.data <- gset@featureData@data
feature.data[1:5,1:5]
#second approach
rownames(top)
probeList <- rownames(top)
geneSymbol <- getSYMBOL(probeList,'hgu133plus2.db') 
top$genesymbol<-geneSymbol
top$genesymbol

DEG_selected<-top[which(top$P.Value<=0.01  & top$adj.P.Val<=0.01 & (top$logFC>=1 | top$logFC<=-1) & top$genesymbol != 'NA') ,]
rownames(DEG_selected)

DEG_selected$prob_ID<-rownames(DEG_selected)
uniq_genes<-unique(DEG_selected$genesymbol)

####select uniqe genes from their probes with respect to max logFC of each gene similar group 
## ?? To decide which probe to choose among probs of one genes, some probes of the same gene
#have negative and some have positive logFC !!
split<-split(DEG_selected,DEG_selected$genesymbol)
lapp<-lapply(split,function(chunk)
{ if (chunk$logFC>0) { chunk[which.max(chunk$logFC),] }
else{
chunk[which.min(chunk$logFC),]
}})
length(lapp)
DEG_selected_uniq<-do.call(rbind,lapp)
dim(DEG_selected_uniq)
###**** 2231 should be the number of DEG_selected_uniq_orderB 

#select up and down regulated genes
upreg.df <- subset(DEG_selected_uniq, logFC>0 & adj.P.val<0.05)
DEG_selected_uniq_up<-DEG_selected_uniq[which(DEG_selected_uniq$logFC>0) ,]
DEG_selected_uniq_down<-DEG_selected_uniq[which(DEG_selected_uniq$logFC<0) ,]
DEG_selected_uniq_orderB<-DEG_selected_uniq[order(-DEG_selected_uniq$B),][1:2000,]
old_2000 <- readxl::read_excel("/Users/saman/Dropbox/my_publications/rs1/results/Table_S1_top_2000_DEGs.xlsx")
old_2000$B #is already sorted by B
old_2000_lst <- old_2000$genesymbol
old_2000_lst


length(intersect(rownames(DEG_selected_uniq_orderB), old_2000$...1))
#1992   #?? check 3 missing CELL files with newer date of "Added"
length(intersect(rownames(DEG_selected_uniq_orderB), tops)) # 16

#export
#write.table(top, file="DEG_All.csv", sep=",")
write.table(top, file="DEG_All_3.txt", sep="\t")
write.table(expression, file="exp3.txt", sep="\t")
write.table(data.mass.exp, file="mass_exp3.txt", sep="\t")
write.table(DEG_selected_uniq, file="DEG_selected_uniq_3.txt", sep="\t")
write.table(DEG_selected_uniq_orderB, file="DEG_selected_uniq_B3.txt", sep="\t")
write.table(DEG_selected, file="DEG_selected_3.txt", sep="\t")
write.table(DEG_selected_uniq_up, file="DEG_selected_uniq_up_3.txt", sep="\t")
write.table(DEG_selected_uniq_down, file="DEG_selected_uniq_down_3.txt", sep="\t")


#Read files
data_mass_exp_log3_df <- read.table("/Users/saman/Dropbox/my_publications/rs1/data/DataSet/data_mass_exp_log3.csv", sep = "\t")
data_mass_exp_log3_df[1:5, 1:5]
dim(data_mass_exp_log3_df)
colnames(data_mass_exp_log3_df)[1:10]
rownames(data_mass_exp_log3_df)[1:10]
DEG_selected_uniq_orderB <- read.table("/Users/saman/Dropbox/my_publications/rs1/data/DataSet/DEG_selected_uniq_B3.txt")

###################################################
# intersection between DEGs and Top communicatives
###################################################
tops <- c("CCNA2", "CCNB1", "CDC20", "CDK1", "FOXM1", "MKI67", "STAT3", "EGF", "H2AFX", "IL1B",
          "IFNG", "STAT1","CXCL8", "CXCL10", "MAPK14", "MCL1","UBE2N")
tops

#pso_gda <- readxl::read_excel("/Users/saman/Dropbox/my_publications/rs1/
#data/psoriasis_Vulgaris_evidence-based_assiciated_genes_C0263361_disease_gda_summary.xlsx")
general_pso_gda <- readxl::read_excel(
  "/Users/saman/Dropbox/my_publications/rs1/data/psoriasis_C0033860_disease_gda_summary.xlsx")
length(general_pso_gda$Gene)
general_pso_gda[1:30, c("Gene", "Gene_id")]

tops[tops %in% general_pso_gda$Gene]  #11


### for all DEGs
top2000DEGs <- readxl::read_excel(
  "/Users/saman/Dropbox/my_publications/rs1/journal_genes/supplementary_materials/Table_S1.xlsx")
top2000DEGs$genesymbole
length(intersect(general_pso_gda$Gene , top2000DEGs$genesymbole)) #225
length(top2000DEGs$genesymbole[top2000DEGs$genesymbole %in% general_pso_gda$Gene])  #225
length(unique(top2000DEGs$genesymbole)) #2000

length(new2000_df$genesymbol[new2000_df$genesymbol %in% general_pso_gda$Gene])  # without unification 235
length(new2000_df$genesymbol[new2000_df$genesymbol %in% general_pso_gda$Gene])
length(unique(new2000_df$genesymbol)) # 1582

length(new3000_df$genesymbol[new3000_df$genesymbol %in% general_pso_gda$Gene]) #310
length(intersect(general_pso_gda$Gene, new3000_genes)) #225  #?? why intersect and %in% have different results??
length(unique(new3000_df$genesymbol)) #2362

length(intersect(new2000_genes , old_2000_lst)) # 943 
#shows there are several probes for the same gene with high B if probes new list is not unified
#After unification of genes in DEG_selected: 1076

###intersection between three sets
tops_gda <- tops[tops %in% general_pso_gda$Gene]
top2000DEGs$genesymbol[top2000DEGs$genesymbol %in% tops_gda]




########################### 
#### newly GEO2R run result
############################

newtop <- read.delim("./workshop/results/GSE30999_top_table.tsv", sep = "\t")
dim(newtop)
newtop$ID
colnames(newtop)
rownames(newtop)
probeList <- newtop$ID
geneSymbol <- getSYMBOL(probeList,'hgu133plus2.db')
##if you use annotation data from a soft file like "GSE..._family_soft1",multiple genesymbols separated by "\t"
#should be split and put into new rows using cSplit
#when each row is an ID and genesymbol
#,cSplit splits the extra symbols in new columns for the gene and puts NAs for the rest of IDs in that column
class(geneSymbol)
newtop$genesymbol<-geneSymbol
newtop$genesymbol
dim(newtop)

DEG_selected<-newtop[which(newtop$adj.P.Val<=0.01 & (newtop$logFC>=1 | newtop$logFC<=-1) & newtop$genesymbol != 'NA') ,]
rownames(DEG_selected)


dim(DEG_selected) #3995 with genesymbol and 4899 with Gene.symbol
length(DEG_selected$Gene.symbol)
length(DEG_selected$genesymbol)
unique(DEG_selected$Gene.symbol)
unique(DEG_selected$genesymbol)
length(unique(DEG_selected$Gene.symbol)) #3124
length(unique(DEG_selected$genesymbol)) #3125

#select uniqe genes with respect to max logFC of each gene similar group 
split<-split(DEG_selected,DEG_selected$genesymbol)
class(split)
split[[1]]

lapp<-lapply(split,function(chunk)
{ if (mean(chunk$logFC>0)) { chunk[which.max(chunk$B),] }  #change chunk$logFC into $B
  else{
    chunk[which.min(chunk$logFC),]
  }})
length(lapp)  #3320
DEG_selected_uniq<-do.call(rbind,lapp)
dim(DEG_selected_uniq)
DEG_selected_uniq[1:5,]


new3000_df <- DEG_selected[order(DEG_selected$B, decreasing = TRUE),][1:3000,]
dim(new3000_df)
new3000_df$B
new3000_genes <- new3000_df$genesymbol
new3000_ids <- new3000_df$ID

new2000_u_df <- DEG_selected_uniq[order(DEG_selected_uniq$B, decreasing = TRUE),][1:2000,]
dim(new2000_u_df)
new2000_u_genes <- new2000_u_df$genesymbol
length(new2000_u_genes)
length(intersect(DEG_selected_uniq$genesymbol, old_2000_lst)) #1370
length(intersect(new3000_genes, old_2000_lst))


##############################
#Venn Diagrams:
#############################
#install.packages("VennDiagram")
library(VennDiagram)

# move to new plotting page
grid.newpage()


# create pairwise Venn diagram
draw.pairwise.venn(area1=1308, area2=17,cross.area=11,
                   category=c("Mango","Banana"),fill=c("Red","Yellow"))

# create Venn diagram with three sets
pdf("intersections_Venn.pdf", width = 8, height = 8)
draw.triple.venn(area1=1308, area2=2000, area3=17, 
                 n12=225, n23=17, n13=11, n123=11, 
                 category=c("GDA(1308)","Top_DEGs(2000)","Top_results(17)"),
                 col="Red",fill=c("Green","Yellow","Blue"), cat.cex=1, scaled=TRUE, euler.d = TRUE)


install.packages("eulerr")
library(eulerr)
VennDiag <- euler(c("GDA (1308)" = 1308, "Top DEGs (2000)" = 2000,
                    "Top results (17)" = 17, "GDA (1308)&Top DEGs (2000)" = 255,
                    "Top DEGs (2000)&Top results (17)" = 17, 
                    "GDA (1308)&Top results (17)" = 11, "GDA (1308)&Top DEGs (2000)&Top results (17)" = 11))
plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("Green","Yellow","Blue"))

dev.off()




################
#####CMM1 begings
################
# """
# #for WGCNA:
# We do not recommend filtering genes by differential expression. WGCNA is designed to be an unsupervised analysis
# method that clusters genes based on their expression profiles. Filtering genes by differential expression will 
# lead to a set of correlated genes that will essentially form a single (or a few highly correlated) modules. 
# It also completely invalidates the scale-free topology assumption, so choosing soft thresholding power by scale-free topology fit will fail.
# """
#### Using rs1 expression data
##preparing expression (feature) dataframe
#########################################
allProbes <- rownames(data_mass_exp_log3_df)
allProbes[1:5]
all_geneSymbol <- getSYMBOL(allProbes,'hgu133plus2.db')
length(all_geneSymbol)
data_mass_exp_log3_df$all_genesymbol <- all_geneSymbol
data_mass_exp_log3_df$all_genesymbol[1:5]

###select uniqe genes 
split<-split(data_mass_exp_log3_df,data_mass_exp_log3_df$all_genesymbol)
length(split)  #20174
split$A1CF[1:5,1:5]
rownames(split$A2M)
split$A2M[,1]
split$A2M[,ncol(split$A2M)]

which(is.na(split$A1BG), arr.ind = T)
split$A1BG[which(is.na(split$A1BG), arr.ind = T)[,1],]

lapp <- c()
lapp<-lapply(split,function(chunk)
{ chunk[which.max(rowMeans(chunk[,-ncol(chunk)])), ]} 
)
length(lapp)  #20174
data_mass_exp_log2_uniq_df <- NULL
data_mass_exp_log2_uniq_df<-do.call(rbind,lapp)
dim(data_mass_exp_log2_uniq_df)
data_mass_exp_log2_uniq_df[1:5,1:5]
rownames(data_mass_exp_log2_uniq_df)
rownames(split$A2M)
data_mass_exp_log2_uniq_df[1:3, ncol(data_mass_exp_log2_uniq_df)]

#omit the last column which is the gene names from ..._uniq_df and omit rownames as well from features_df
all_uniq_genes <- data_mass_exp_log2_uniq_df$all_genesymbol
data_mass_exp_log2_uniq_df <- data_mass_exp_log2_uniq_df[, -ncol(data_mass_exp_log2_uniq_df)]
tail( colnames(data_mass_exp_log2_uniq_df), n=1)  #the last element of a vector

features_df <- data_mass_exp_log2_uniq_df[, -ncol(data_mass_exp_log2_uniq_df)]
features_df[1:5, 1:5]

rownames(features_df) <- NULL
colnames(features_df) <- NULL
features_df[1:5, 1:5]
dim(features_df)

dim(data_mass_exp_log2_uniq_df)
sum(rownames(data_mass_exp_log2_uniq_df) %in% all_uniq_genes)
colnames(data_mass_exp_log2_uniq_df)

#export
write.csv(data_mass_exp_log2_uniq_df,
          "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/data_mass_exp_log2_uniq.csv")
write.csv(all_uniq_genes,
          "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/gene_names.csv")
write.csv(features_df,"/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/features.csv")
write.csv(colnames(data_mass_exp_log2_uniq_df),
          "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/sample_names.csv")


BiocManager::install("minet")

# read rs1 expression data
data_mass_exp_log2_uniq_df <- 
  read.csv( "/mnt/4tb/Dropbox/systems_biology_projects/pso_cmm/data/data_mass_exp_log2_uniq.csv")
all_uniq_genes <- read.csv( "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/gene_names.csv")
features_df <- read.csv("/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/features.csv", header = FALSE)
sample_names <- read.csv( "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/sample_names.csv")

## Based on rs1, separate D and N data
########################
data_mass_exp_log2_uniq_df[1:5, 1:5]
rownames(data_mass_exp_log2_uniq_df)

N <- data_mass_exp_log2_uniq_df[,grepl("NL", colnames(data_mass_exp_log2_uniq_df))]
D <- data_mass_exp_log2_uniq_df[,grepl("LS", colnames(data_mass_exp_log2_uniq_df))]
dim(N)  # 20174  86
dim(D)  # 20174  86

rownames(N) <- data_mass_exp_log2_uniq_df$X
rownames(D) <- data_mass_exp_log2_uniq_df$X

#highest genes for N
N$means <- rowMeans(N)
sort(rowMeans(N), decreasing = TRUE, index.return=T)$ix #the first index returne = 199
rowMeans(N)[199] #14.08098 , which is the largest mean among genes means

N_sorted <- N[match(sort(rowMeans(N), decreasing=T), N$means),]
rowMeans(N_sorted)
rownames(N_sorted)
#rownames(N_sorted) <- 1:nrow(N_sorted)
#rownames(N_sorted) 
N_800_lst <- rownames(N_sorted)[1:800]
N_800_lst

#highest genes for D
D$means <- rowMeans(D)
dim(D)
rownames(D)
sort(rowMeans(D), decreasing = TRUE, index.return=T)$ix #the first index returned = 16758
rowMeans(D)[16758] #14.12 , which is the largest mean among genes means of D

D_sorted <- D[match(sort(rowMeans(D), decreasing=T), D$means),]
D_sorted[1:5, 1:5]
D_sorted$means  # 14.12 is also the highest
rowMeans(D_sorted) # 14.12 is also the top
#rownames(N_sorted) <- 1:nrow(N_sorted)
#rownames(N_sorted) 
D_800_lst <- rownames(D_sorted)[1:800]
D_800_lst


# intersectional and uniq lists of N and D
N_D_lst <- c(N_800_lst, D_800_lst)
length(N_D_lst) #1600
inter_lst <- intersect(N_800_lst, D_800_lst)
inter_lst #651
sum(N_D_lst %in% inter_lst)
sum(!(N_D_lst %in% inter_lst))

uniq_N_D_exclu_lst <- N_D_lst[!(N_D_lst %in% inter_lst)]
length(uniq_N_D_exclu_lst)  #298

length(unique(old_2000_lst)) # 2000

uniq_N_D_lst <- unique(N_D_lst)
length(uniq_N_D_lst) # 949


length(unique(c(uniq_N_D_lst, old_2000_lst))) #2827
tot_lst <- unique(c(uniq_N_D_lst, old_2000_lst))
tot_lst

# MINET datasets:

#Normal
mi_N <- NULL
mi_N <- N_sorted[tot_lst, -ncol(N_sorted)]

for (c in colnames(mi_N)[1:3]) {
  mi_N[,c] <- round(mi_N[,c], digits = 2)
  #print(round(mi_N[1:5, c], digits = 2))
}

colnames(mi_N)
rownames(mi_N)
class(mi_N)  # data.frame
mi_N[1:5, 1:5]
dim(mi_N) # 2819  85

# the missing values had no gene name
missing_values_df <- as.data.frame(which(is.na(mi_N), arr.ind = TRUE))
missing_values_df[1:20, 1:2]
missing_genes <- unique(missing_values_df$row)
missing_genes
mi_N[missing_genes,1:3]
mi_N[1405,]
dim(missing_values_df)
mi_N <- na.omit(mi_N)



mi_N_t <- as.data.frame(t(mi_N))
dim(mi_N_t)  #85 2819








#Diseased
mi_D <- D_sorted[tot_lst, -ncol(D_sorted)]
mi_D <- t(mi_D)
dim(mi_D)  #85 2223
mi_D

######################
## run minet with discretization none and estimator spearman
library(minet)
#BiocManager::install("Rgraphviz")
#library(Rgraphviz)

net_N <- minet(mi_N_t) #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N ,"graphNEL") )
net_N[1:5, 1:5]
colnames(net_N[1:5, 1:5])
rownames(net_N[1:5,1:5])
dim(net_N) #2819  2819
class(net_N)  # matrix array
net_N != 0
net_N["ACTB", 2]
#test
comb <- t(combn(colnames(net_N[1:5,1:5]), 2))
data.frame(comb, mi=net_N[comb])
net_N[which(net_N !=0, arr.ind = TRUE) ]
net_N_lower <- lower.tri(net_N[1:5,1:5], diag = FALSE)
net_N_lower
net_N_lower <- net_N[1:5,1:5][net_N_lower]
net_N_lower
net_N_lower[which(net_N_lower != 0, arr.ind = TRUE)]


#net_N_lower
net_N_lower <- lower.tri(net_N, diag = FALSE)
net_N_lower <- net_N[net_N_lower]
net_N_lower
net_N_lower <- net_N_lower[which(net_N_lower != 0, arr.ind = TRUE)]
length(net_N_lower[net_N_lower > 0.2]) #975
length(net_N_lower[net_N_lower > 0.1]) #2786
length(net_N_lower[net_N_lower > 0.05]) #31187
hist(net_N_lower[net_N_lower > 0.1])

trace(minet, edit = TRUE)

###########################
#run with discretization equalfreq and estimator mi.empirical
net_N_emp <- minet(mi_N_t, disc = "equalfreq", estimator = "mi.empirical") 
#, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N_emp ,"graphNEL") )
net_N_emp[1:5, 1:5]
colnames(net_N_emp[1:5, 1:5])
rownames(net_N_emp[1:5,1:5])
dim(net_N_emp) #2819  2819
class(net_N_emp)  # matrix array
net_N_emp != 0
net_N_emp["ACTB", 2]
#test
comb <- t(combn(colnames(net_N_emp[1:5,1:5]), 2))
data.frame(comb, mi=net_N_emp[comb])
net_N_emp[which(net_N_emp !=0, arr.ind = TRUE) ]
net_N_emp_lower <- lower.tri(net_N_emp[1:5,1:5], diag = FALSE)
net_N_emp_lower
net_N_emp_lower <- net_N_emp[1:5,1:5][net_N_emp_lower]
net_N_emp_lower
net_N_emp_lower[which(net_N_emp_lower != 0, arr.ind = TRUE)]


#net_N_emp_lower
net_N_emp_lower <- lower.tri(net_N_emp, diag = FALSE)
net_N_emp_lower <- net_N_emp[net_N_emp_lower]
net_N_emp_lower
net_N_emp_lower <- net_N_emp_lower[which(net_N_emp_lower != 0, arr.ind = TRUE)]
length(net_N_emp_lower[net_N_emp_lower > 0.2]) #2817
length(net_N_emp_lower[net_N_emp_lower > 0.1]) #35930
length(net_N_emp_lower[net_N_emp_lower > 0.05]) #447198
hist(net_N_emp_lower[net_N_emp_lower > 0.2])

#########################################################################
##minet with discretization equalfreq and estimator mi.mm
net_N_mm <- minet(mi_N_t, disc = "equalfreq", estimator = "mi.mm") #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N_mm ,"graphNEL") )
net_N_mm[1:5, 1:5]
colnames(net_N_mm[1:5, 1:5])
rownames(net_N_mm[1:5,1:5])
dim(net_N_mm) #2819  2819
class(net_N_mm)  # matrix array
net_N_mm != 0
net_N_mm["ACTB", 2]
#test
comb <- t(combn(colnames(net_N_mm[1:5,1:5]), 2))
data.frame(comb, mi=net_N_mm[comb])
net_N_mm[which(net_N_mm !=0, arr.ind = TRUE) ]
net_N_mm_lower <- lower.tri(net_N_mm[1:5,1:5], diag = FALSE)
net_N_mm_lower
net_N_mm_lower <- net_N_mm[1:5,1:5][net_N_mm_lower]
net_N_mm_lower
net_N_mm_lower[which(net_N_mm_lower != 0, arr.ind = TRUE)]


#net_N_mm_lower
net_N_mm_lower <- lower.tri(net_N_mm, diag = FALSE)
net_N_mm_lower <- net_N_mm[net_N_mm_lower]
net_N_mm_lower
net_N_mm_lower <- net_N_mm_lower[which(net_N_mm_lower != 0, arr.ind = TRUE)]
length(net_N_mm_lower[net_N_mm_lower > 0.2]) #5785
length(net_N_mm_lower[net_N_mm_lower > 0.1]) #107415
length(net_N_mm_lower[net_N_mm_lower > 0.05]) #692748
hist(net_N_mm_lower[net_N_mm_lower > 0.2])




############################################
## Using new expression data based on GEO2R 
#############################################
## separate D and N data
########################
ex[1:5, 1:5]
dim(ex)
colnames(ex)[grep("N", colnames(ex))]
sum(grepl("NA",ex$geneSymbol))
ex$geneSymbol[grepl("NA", ex$geneSymbol)] #super dangerous, gets all the genes with letters "NA" as well
ex$geneSymbol[1:10]
ex$geneSymbol[!is.na(ex$geneSymbol)][1:10]
!is.na(ex$geneSymbol)[1:10]
ex <- ex[!is.na(ex$geneSymbol),]
ex$geneSymbol[1:10]

N <- ex[,grepl("N", colnames(ex))]
D <- ex[,grepl("L", colnames(ex))]
dim(N)  # 41905  85
dim(D)  # 41905  85

N$probe <- ex$probe
D$probe <- ex$probe

N$geneSymbol <- ex$geneSymbol
D$geneSymbol <- ex$geneSymbol

#highest genes for N
colnames(N[1:86])
N$means <- rowMeans(N[1:85])
sort(rowMeans(N[,1:85]), decreasing = TRUE, index.return=T)$ix #the first index of sorted means vector returned = 6819, 
#$x returns only the values when there's no rowname
rowMeans(N[1:85])[6819] #15.61 , which is the largest mean among genes means
N[which(N[,"means"] == 15.61800), "means"]  #it doesn't work if you don't write the full digits of the value 

match(sort(rowMeans(N[,1:85]), decreasing=T), N$means)
match(sort(N$means, decreasing = T), N$means)
N_sorted <- N[match(sort(rowMeans(N[,1:85]), decreasing=T), N$means),]
rowMeans(N_sorted[,1:85])
rownames(N_sorted)
#rownames(N_sorted) <- NULL
class(rownames(N_sorted) ) #character
N_1000_lst <- N_sorted$probe[1:1000]
N_1000_lst

#highest genes for D
dim(D)
colnames(D)
D$means <- rowMeans(D[1:85])
sort(rowMeans(D[,1:85]), decreasing = TRUE, index.return=T)$ix #the first index of sorted means vector returned = 6819, 
#$x returns only the values when there's no rowname
rowMeans(D[1:85])[6819] #15.565 , which is the largest mean among genes means
D[which(round(D[,"means"], digits = 2) == 15.57), "means"]  #it doesn't work if you don't write the full digits of the value 
D["200092_s_at","means"]

match(sort(rowMeans(D[,1:85]), decreasing=T), D$means)
match(sort(D$means, decreasing = T), D$means)
D_sorted <- D[match(sort(rowMeans(D[,1:85]), decreasing=T), D$means),]
rowMeans(D_sorted[,1:85])
rownames(D_sorted)
#rownames(N_sorted) <- NULL
class(rownames(D_sorted) ) #character
D_1000_lst <- D_sorted$probe[1:1000]
D_1000_lst

###########################################
# intersectional and uniq lists of N and D
N_D_lst <- c(N_1000_lst, D_1000_lst)
length(N_D_lst) #1600
inter_lst <- intersect(N_800_lst, D_800_lst)
inter_lst #651
sum(N_D_lst %in% inter_lst)
sum(!(N_D_lst %in% inter_lst))

uniq_N_D_exclu_lst <- N_D_lst[!(N_D_lst %in% inter_lst)]
length(uniq_N_D_exclu_lst)  #298

length(unique(old_2000_lst)) # 2000

uniq_N_D_lst <- unique(N_D_lst)  #top 2000 highly expressed 
length(uniq_N_D_lst) # 1152


tot_lst_5000 <- unique(c(uniq_N_D_lst, new3000_ids))
length(tot_lst_5000) # 4087
##################
## MINET datasets:

#Normal
mi_N <- NULL
colnames(N_sorted)[1:85]
mi_N <- N_sorted[tot_lst_5000,1:85]

# for (c in colnames(mi_N)[1:3]) {
#   mi_N[,c] <- round(mi_N[,c], digits = 2)
#   #print(round(mi_N[1:5, c], digits = 2))
# }

colnames(mi_N)
rownames(mi_N)
class(mi_N)  # data.frame
mi_N[1:5, 1:5]
dim(mi_N) # 2819  85

# the missing values had no gene name
which(is.na(mi_N), arr.ind = TRUE)
missing_values_df <- as.data.frame(which(is.na(mi_N), arr.ind = TRUE))
missing_values_df[1:20, 1:2]
missing_genes_rowidx <- unique(missing_values_df$row)
missing_genes_rowidx
length(missing_genes_rowidx)
mi_N[missing_genes,1:15]
mi_N[966,]
dim(mi_N) #4087 85
dim(mi_N[-missing_genes_rowidx,])  #3687 85
dim(na.omit(mi_N))                 #3687 85

mi_N[-missing_genes_rowidx,][1:10, 1:10]
na.omit(mi_N)[1:10,1:10]

mi_N[-missing_genes_rowidx,][966,]
na.omit(mi_N)[966,]

mi_N <- mi_N[-missing_genes_rowidx,]
dim(mi_N)  #3687 85

mi_N_t <- as.data.frame(t(mi_N))
dim(mi_N_t)  #85 2819








#Diseased
mi_D <- NULL
colnames(D_sorted)[1:85]
mi_D <- D_sorted[tot_lst_5000,1:85]

colnames(mi_D)
rownames(mi_D)
class(mi_D)  # data.frame
mi_D[1:5, 1:5]
dim(mi_D) # 4087  85

# the missing values had no gene name
which(is.na(mi_D), arr.ind = TRUE)
missing_values_df <- as.data.frame(which(is.na(mi_D), arr.ind = TRUE))
missing_values_df[1:20, 1:2]
missing_genes_rowidx <- unique(missing_values_df$row)
missing_genes_rowidx
length(missing_genes_rowidx)
mi_D[79,]
dim(mi_D) #4087 85
dim(mi_D[-missing_genes_rowidx,])  #3681 85
dim(na.omit(mi_D))                 #3681 85

mi_D[-missing_genes_rowidx,][1:10, 1:10]
na.omit(mi_D)[1:10,1:10]

mi_D[-missing_genes_rowidx,][79,]
na.omit(mi_D)[79,]

mi_D <- mi_D[-missing_genes_rowidx,]
dim(mi_D)  #3681 85

mi_D_t <- as.data.frame(t(mi_D))
dim(mi_D_t)  #85 3681

######################################
##run minet with different parameters
#####################################

library(minet)
#BiocManager::install("Rgraphviz")
#library(Rgraphviz)
########
##Normal

####
##run minet with discretization none and estimator spearman
net_N <- minet(mi_N_t) #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 

#plot( as(net_N ,"graphNEL") )
net_N[1:5, 1:5]
colnames(net_N[1:5, 1:5])
rownames(net_N[1:5,1:5])
dim(net_N) 
class(net_N)  # matrix array
net_N != 0
net_N["200092_s_at", 2]
#test
comb <- t(combn(colnames(net_N[1:5,1:5]), 2))
comb
class(comb) #matrix array
as.data.frame(comb, mi=net_N[comb])
net_N[which(net_N !=0, arr.ind = TRUE) ]
net_N_lower <- lower.tri(net_N[1:5,1:5], diag = FALSE)
net_N_lower
net_N_lower <- net_N[1:5,1:5][net_N_lower]
net_N_lower
net_N_lower[which(net_N_lower != 0, arr.ind = TRUE)]


#net_N_lower
net_N_lower <- lower.tri(net_N, diag = FALSE)
net_N_lower <- net_N[net_N_lower]
net_N_lower
net_N_lower <- net_N_lower[which(net_N_lower != 0, arr.ind = TRUE)]
length(net_N_lower[net_N_lower > 0.2]) #930
length(net_N_lower[net_N_lower > 0.1]) #2953
length(net_N_lower[net_N_lower > 0.05]) #22752
hist(net_N_lower[net_N_lower > 0.2])

trace(minet, edit = TRUE)


# N_emp ----
##run with discretization equalfreq and estimator mi.empirical
net_N_emp <- minet(mi_N_t, disc = "equalfreq", estimator = "mi.empirical") 
net_N_emp <- round(net_N_emp, 2)
write.csv(net_N_emp,"./data/adjacency_net_N_emp.csv")
#, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N_emp ,"graphNEL") )
net_N_emp[1:5, 1:5]
colnames(net_N_emp[1:5, 1:5])
rownames(net_N_emp[1:5,1:5])
dim(net_N_emp) #3687  3687
class(net_N_emp)  # matrix array
net_N_emp != 0
net_N_emp["224500_s_at", 2]

# create edgelist
# comb <- t(combn(colnames(net_N_emp[1:5,1:5]), 2))
# data.frame(comb, mi=net_N_emp[comb])
comb_N_emp <- t(combn(colnames(net_N_emp), 2))
head(comb_N_emp)
el_mi_N_emp <- data.frame(comb_N_emp, mi=net_N_emp[comb_N_emp])
length(el_mi_N_emp$mi[el_mi_N_emp$mi > 0.2])  #3372
el_mi_N_emp_2 <- el_mi_N_emp[el_mi_N_emp$mi >0.2,]
el_mi_N_emp_2$mi <- round(el_mi_N_emp_2$mi, 2)

#save mi>0.2 csv
write.csv(el_mi_N_emp_2, "./data/el_mi_N_emp_2.csv")

##create 4.0
el_mi_N_emp_4 <- el_mi_N_emp[el_mi_N_emp$mi >0.4,]
el_mi_N_emp_4$mi <- round(el_mi_N_emp_4$mi, 2)
#save mi>0.4 csv
write.csv(el_mi_N_emp_4, "./data/el_mi_N_emp_4.csv")


## unuique genes of MI > 0.4 
uniq_N_4 <- c(unique(el_mi_N_emp_4$X1), unique(el_mi_N_emp_4$X2))
length(uniq_N_4) #4352
length(unique(uniq_N_4)) #3640  ?? why?! the combination should've omitted 
uniq_N_4 <- unique(uniq_N_4)
length(uniq_N_4) #3640
# 



#adjacency matrix of mi>0.2
net_N_emp_2 <- net_N_emp[colnames(net_N_emp) %in% uniq_N_2, rownames(net_N_emp) %in% uniq_N_2]
net_N_emp_2 <- round(net_N_emp_2, 2)
dim(net_N_emp_2)
net_N_emp_2[1:5, 1:5]
write.csv(net_N_emp_2, "./data/net_N_emp_2.csv")
write.csv(net_N_emp_2[1:200,1:200], "./data/test_net_N_emp_2.csv")


net_N_emp[which(net_N_emp !=0, arr.ind = TRUE) ]
net_N_emp_lower <- lower.tri(net_N_emp[1:5,1:5], diag = FALSE)
net_N_emp_lower
net_N_emp_lower <- net_N_emp[1:5,1:5][net_N_emp_lower]
net_N_emp_lower
net_N_emp_lower[which(net_N_emp_lower != 0, arr.ind = TRUE)]


#net_N_emp_lower
net_N_emp_lower <- lower.tri(net_N_emp, diag = FALSE)
net_N_emp_lower <- net_N_emp[net_N_emp_lower]
net_N_emp_lower
net_N_emp_lower <- net_N_emp_lower[which(net_N_emp_lower != 0, arr.ind = TRUE)]
length(net_N_emp_lower[net_N_emp_lower > 0.2]) #3372
length(net_N_emp_lower[net_N_emp_lower > 0.1]) #52592
length(net_N_emp_lower[net_N_emp_lower > 0.05]) #728373
hist(net_N_emp_lower[net_N_emp_lower > 0.2])

####
##minet with discretization equalfreq and estimator mi.mm
net_N_mm <- minet(mi_N_t, disc = "equalfreq", estimator = "mi.mm") #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N_mm ,"graphNEL") )
net_N_mm[1:5, 1:5]
colnames(net_N_mm[1:5, 1:5])
rownames(net_N_mm[1:5,1:5])
dim(net_N_mm) #2819  2819
class(net_N_mm)  # matrix array
net_N_mm != 0
net_N_mm["ACTB", 2]
#test
comb <- t(combn(colnames(net_N_mm[1:5,1:5]), 2))
data.frame(comb, mi=net_N_mm[comb])
net_N_mm[which(net_N_mm !=0, arr.ind = TRUE) ]
net_N_mm_lower <- lower.tri(net_N_mm[1:5,1:5], diag = FALSE)
net_N_mm_lower
net_N_mm_lower <- net_N_mm[1:5,1:5][net_N_mm_lower]
net_N_mm_lower
net_N_mm_lower[which(net_N_mm_lower != 0, arr.ind = TRUE)]


#net_N_mm_lower
net_N_mm_lower <- lower.tri(net_N_mm, diag = FALSE)
net_N_mm_lower <- net_N_mm[net_N_mm_lower]
net_N_mm_lower
net_N_mm_lower <- net_N_mm_lower[which(net_N_mm_lower != 0, arr.ind = TRUE)]
length(net_N_mm_lower[net_N_mm_lower > 0.2]) #7205
length(net_N_mm_lower[net_N_mm_lower > 0.1]) #161355
length(net_N_mm_lower[net_N_mm_lower > 0.05]) #1126287
hist(net_N_mm_lower[net_N_mm_lower > 0.2])

###########
##Diseased

#run minet with discretization none and estimator spearman
net_D <- minet(mi_D_t) #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_D)) 

#plot( as(net_D ,"graphDEL") )
net_D[1:5, 1:5]
colnames(net_D[1:5, 1:5])
rownames(net_D[1:5,1:5])
dim(net_D) #2819  2819
class(net_D)  # matrix array
net_D != 0
net_D["200092_s_at", 2]
#test
comb <- t(combn(colnames(net_D[1:5,1:5]), 2))
class(comb) #matrix array
as.data.frame(comb, mi=net_D[comb])
net_D[which(net_D !=0, arr.ind = TRUE) ]
net_D_lower <- lower.tri(net_D[1:5,1:5], diag = FALSE)
net_D_lower
net_D_lower <- net_D[1:5,1:5][net_D_lower]
net_D_lower
net_D_lower[which(net_D_lower != 0, arr.ind = TRUE)]


#net_D_lower
net_D_lower <- lower.tri(net_D, diag = FALSE)
net_D_lower <- net_D[net_D_lower]
net_D_lower
net_D_lower <- net_D_lower[which(net_D_lower != 0, arr.ind = TRUE)]
length(net_D_lower[net_D_lower > 0.2]) #1253
length(net_D_lower[net_D_lower > 0.1]) #3255
length(net_D_lower[net_D_lower > 0.05]) #35701
hist(net_D_lower[net_D_lower > 0.2])

trace(minet, edit = TRUE)

# D_emp ----
#run with discretization equalfreq and estimator mi.empirical
net_D_emp <- minet(mi_D_t, disc = "equalfreq", estimator = "mi.empirical") 
net_D_emp <- round(net_D_emp, 2)
write(net_D_emp,"./data/adjacency_net_D_emp.csv")
#, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_D)) 
#plot( as(net_D_emp ,"graphDEL") )
net_D_emp[1:5, 1:5]
colnames(net_D_emp[1:5, 1:5])
rownames(net_D_emp[1:5,1:5])
dim(net_D_emp) #3681  3681
class(net_D_emp)  # matrix array
net_D_emp != 0
net_D_emp["209220_at", 2]

# teset create edgelist
comb <- t(combn(colnames(net_D_emp[1:5,1:5]), 2))
data.frame(comb, mi=net_D_emp[comb])
net_D_emp[which(net_D_emp !=0, arr.ind = TRUE) ]
net_D_emp_lower <- lower.tri(net_D_emp[1:5,1:5], diag = FALSE)
net_D_emp_lower
net_D_emp_lower <- net_D_emp[1:5,1:5][net_D_emp_lower]
net_D_emp_lower
net_D_emp_lower[which(net_D_emp_lower != 0, arr.ind = TRUE)]

#create edge list
comb_D_emp <- t(combn(colnames(net_D_emp), 2))
head(comb_D_emp)
el_mi_D_emp <- data.frame(comb_D_emp, mi=net_D_emp[comb_D_emp])
length(el_mi_D_emp$mi[el_mi_D_emp$mi > 0.2])  #3262
el_mi_D_emp_2 <- el_mi_D_emp[el_mi_D_emp$mi >0.2,]
el_mi_D_emp_2$mi <- round(el_mi_D_emp_2$mi, 2)
#save csv
write.csv(el_mi_D_emp_2, "./data/el_mi_D_emp_2.csv")

##create 4.0
el_mi_D_emp_4 <- el_mi_D_emp[el_mi_D_emp$mi >0.4,]
el_mi_D_emp_4$mi <- round(el_mi_D_emp_4$mi, 2)
#save mi>0.4 csv
write.csv(el_mi_D_emp_4, "./data/el_mi_D_emp_4.csv")


#unuique D genes of MI > 0.4
comb_tot_D <- c(el_mi_D_emp_4$X1, el_mi_D_emp_4$X2)
length(unique(comb_tot_D)) #3597
uniq_D_4 <- c(unique(el_mi_D_emp_4$X1), unique(el_mi_D_emp_4$X2))
class(uniq_D_4)
length(uniq_D_4) #4349
uniq_D_4 <- (unique(uniq_D_4)) 
length(uniq_D_4)#3597  

#adjacency matrix of mi>0.2
net_D_emp_2 <- net_D_emp[colnames(net_D_emp) %in% uniq_D_2, rownames(net_D_emp) %in% uniq_D_2]
dim(net_D_emp_2)
net_D_emp_2 <- round(net_D_emp_2, 2)
net_D_emp_2[1:5,1:5]
write.csv(net_D_emp_2, "./data/net_D_emp_2.csv")






#net_D_emp_lower
net_D_emp_lower <- lower.tri(net_D_emp, diag = FALSE)
net_D_emp_lower <- net_D_emp[net_D_emp_lower]
net_D_emp_lower
net_D_emp_lower <- net_D_emp_lower[which(net_D_emp_lower != 0, arr.ind = TRUE)]
length(net_D_emp_lower[net_D_emp_lower > 0.2]) #3262
length(net_D_emp_lower[net_D_emp_lower > 0.1]) #48266
length(net_D_emp_lower[net_D_emp_lower > 0.05]) #745228
hist(net_D_emp_lower[net_D_emp_lower > 0.2])

####
##minet with discretization equalfreq and estimator mi.mm
# library(future)
# plan(tweak(workers=4, cluster))
net_D_mm <- minet(mi_D_t, disc = "equalfreq", estimator = "mi.mm") #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_D)) 
#plot( as(net_D_mm ,"graphDEL") )
net_D_mm[1:5, 1:5]
colnames(net_D_mm[1:5, 1:5])
rownames(net_D_mm[1:5,1:5])
dim(net_D_mm) #3681 3681
class(net_D_mm)  # matrix array
net_D_mm != 0
net_D_mm["229686_at", 2]
#test
comb <- t(combn(colnames(net_D_mm[1:5,1:5]), 2))
data.frame(comb, mi=net_D_mm[comb])
net_D_mm[which(net_D_mm !=0, arr.ind = TRUE) ]
net_D_mm_lower <- lower.tri(net_D_mm[1:5,1:5], diag = FALSE)
net_D_mm_lower
net_D_mm_lower <- net_D_mm[1:5,1:5][net_D_mm_lower]
net_D_mm_lower
net_D_mm_lower[which(net_D_mm_lower != 0, arr.ind = TRUE)]


#net_D_mm_lower
net_D_mm_lower <- lower.tri(net_D_mm, diag = FALSE)
net_D_mm_lower <- net_D_mm[net_D_mm_lower]
net_D_mm_lower
net_D_mm_lower <- net_D_mm_lower[which(net_D_mm_lower != 0, arr.ind = TRUE)]
length(net_D_mm_lower[net_D_mm_lower > 0.2]) #6285
length(net_D_mm_lower[net_D_mm_lower > 0.1]) #150457
length(net_D_mm_lower[net_D_mm_lower > 0.05]) #1150014
hist(net_D_mm_lower[net_D_mm_lower > 0.2])


##
#******* on your cytoscape install all Dynet... apps, aMatReader and centiscape
save.image("mi.RData")
load("mi.RData")


######################################################################################
##Using rGraph
#############

#install.packages('igraph')
#install.packages('DirectedClustering')
library(igraph)
library(DirectedClustering)

makeGraph <- function(graph){
  vertexList = c()
  for(i in 1:nrow(graph)){
    vertexList <- c(vertexList, graph[i,1], graph[i,2])
  }
  vertexList = unique(vertexList)
  numericalGraph = matrix("", nrow(graph), ncol(graph))
  numericalGraph[,1] = match(graph[,1], vertexList)
  numericalGraph[,2] = match(graph[,2], vertexList)
  numericalGraph[,3] = graph[,3]
  #uniqueNumGraph = numericalGraph[!duplicated(numericalGraph),]
  edgeList = c()
  for(i in 1:nrow(numericalGraph)){
    edgeList <- c(edgeList, numericalGraph[i,1], numericalGraph[i,2])
  }
  #numericalGraph= numericalGraph[,-3]
  g = make_undirected_graph(edgeList)
  E(g)$weight = graph[,3]
  g
}


# READ AND ANALYZE NETS----

path = "/Users/saman/Dropbox/systems_biology_projects/pso_cmm/data/"

#### 1. Normal net
net.normal <- read.csv(paste(path, 'el_mi_N_emp_2.csv', sep = ''), stringsAsFactors = FALSE)
net.normal <- net.normal[,-1]
g <- makeGraph(net.normal)

V(g)$degree <- degree(g)                        # Degree centrality
V(g)$eig <- evcent(g)$vector                    # Eigenvector centrality
V(g)$betweenness <- betweenness(g)              # Vertex betweenness centrality

net.normal.centrality <- data.frame(row.names   = vertexList[V(g)],
                                    degree      = V(g)$degree,
                                    #closeness   = V(g)$closeness,
                                    betweenness = V(g)$betweenness,
                                    eigenvector = V(g)$eig)

write.csv(file = paste(path,'normal_centralities.csv'), net.normal.centrality)

# compute clustering coefficient 
#Get Adjacency 
A<-get.adjacency(g, sparse=FALSE, attr="weight") 
#Compute Barrat et al. (2004) coefficient 
BarratClust<-ClustBCG(A, "undirected")




#### 2. Disease net
net.disease <- read.csv(paste(path, 'el_mi_D_emp_2.csv', sep = ''), stringsAsFactors = FALSE)
net.disease <- net.disease[,-1]
g.d <- makeGraph(net.disease)

V(g.d)$degree <- degree(g.d)                        # Degree centrality
V(g.d)$eig <- evcent(g.d)$vector                    # Eigenvector centrality
V(g.d)$betweenness <- betweenness(g.d)              # Vertex betweenness centrality

net.disease.centrality <- data.frame(row.names   = vertexList[V(g.d)],
                                     degree      = V(g.d)$degree,
                                     betweenness = V(g.d)$betweenness,
                                     eigenvector = V(g.d)$eig)

write.csv(file = paste(path,'disease_centralities.csv'), net.disease.centrality)


# compute clustering coefficient 
B<-get.adjacency(g.d, sparse=FALSE, attr="weight") 
#Compute Barrat et al. (2004) coefficient 
BarratClust<-ClustBCG(B, "undirected")


#? adjacency to edgelist

