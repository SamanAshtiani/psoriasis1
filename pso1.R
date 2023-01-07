setwd("/Users/saman/Dropbox/my_publications/rs1/data/")
setwd("/Users/saman/Dropbox/systems_biology_projects/pso_cmm/")
setwd("~/git/psoriasis1/")



setRepositories()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

BiocManager::install("minet")
BiocManager::install("hgu133plus2.db")
BiocManager::install("Vennerable")
install.packages("")

library(pheatmap) #
library(impute)
library(affy)  
library(ggplot2) #
library(limma) #
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
library(minet)

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

##to deal with genes with low variation between groups: subtract the mean expression of each gene
#-->> each row is minus the mean expression of the gene in that row:
ex.scale <- as.data.frame(scale(t(ex[,1:170]), scale = FALSE, center = TRUE))
#scale acts on columns -> mean=0 and divided by SD for each column -> z-score
mean(ex.scale[,1]) #almost 0 if scale = TRUE
dim(ex.scale) # 170 54675
class(ex.scale)

##now pca of samples:
pc <- mixOmics::pca(ex.scale, ncomp = 5)
dim(pc$variates$X) #170  5  : 170 points in a 5-dimension space
rownames(pc$variates$X)

#plot pca
pcnames <- sprintf("PC%d (%2.0f%%)", 1:ncol(pc$variates$X), unlist(pc$prop_expl_var)*100)

plot( ggplot(data.frame(pc$variates$X, sampleInfo), aes(x=PC1, y=PC2, group=group, color=group, shape= group )) + geom_point()
+ theme_bw() + xlab(pcnames[1]) + ylab(pcnames[2]) +
theme(legend.position=c(0.9, 0.9), legend.justification=c(0,0), legend.background =element_rect(fill="transparent"))
# + geom_label(aes(label=rownames(pc$variates$X)), label.size = 0.01)
)

#find lesional outliers:
pcvar <- as.data.frame(pc$variates$X[,1:2])
pcvar[(pcvar$PC1<0 & pcvar$PC2>-100 & sampleInfo$group == "L"),]

# GEO2R-based DEGs result ----

newtop <- read.delim("~/Dropbox/systems_biology_projects/pso_cmm/workshop/results/GSE30999_top_table.tsv", sep = "\t")
dim(newtop) #54675 6
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
newtop[1:14,1:5]

length(which(newtop$adj.P.Val<=0.01 & (newtop$logFC>=1 | newtop$logFC<=-1) & newtop$genesymbol != 'NA')) # 4098,index 12: NA
length(newtop$adj.P.Val<=0.01 & (newtop$logFC>=1 | newtop$logFC<=-1) & newtop$genesymbol != 'NA') #54675
length(which(newtop$adj.P.Val<=0.01 & (newtop$logFC>=1 | newtop$logFC<=-1) & !is.na(newtop$genesymbol)))
#There are NA e.g. at index 12
class(newtop$genesymbol[12]) # NAs in this case were characters

DEG_selected<-newtop[which(newtop$adj.P.Val<=0.01 & (newtop$logFC>=1 | newtop$logFC<=-1) & !is.na(newtop$genesymbol)) ,]

rownames(DEG_selected)
colnames(DEG_selected)

dim(DEG_selected) #4098 7
unique(DEG_selected$genesymbol)
length(unique(DEG_selected$genesymbol)) #3226

#select uniqe genes with respect to max logFC of each gene similar group 
split<-split(DEG_selected,DEG_selected$genesymbol)
class(split)
split[[1]]
split[1]

lapp<-lapply(split,function(chunk)
{ if (mean(chunk$logFC > 0)) { chunk[which.max(chunk$B),] }  #changed chunk$logFC into $B
  else{
    chunk[which.min(chunk$logFC),]
  }})
length(lapp)  #3226
DEG_selected_uniq<-do.call(rbind,lapp)
#dim(DEG_selected_uniq) #3226  7
#DEG_selected_uniq[1:5,]
#lapp[1]
#lapp[[1]]
#names(lapp)
    
# KYNU with B factor 178 must be selected
lapp$KYNU
split$KYNU[1,]
temp_mtr
do.call("rbind", lapp)

new2000_df <- DEG_selected[order(DEG_selected$B, decreasing = TRUE),][1:2000,]
dim(new2000_df)
new2000_df$B
new2000_genes <- new2000_df$genesymbol
new2000_ids <- new2000_df$ID

new2000_u_df <- DEG_selected_uniq[order(DEG_selected_uniq$B, decreasing = TRUE),][1:2000,]
#dim(new2000_u_df)  
new2000_u_genes <- new2000_u_df$genesymbol
length(new2000_u_genes)
new2000_u_id <- new2000_u_df$ID





# CMM1 begings ----

## Using new expression data based on GEO2R ----
## separate D and N data ----

ex[1:5, 1:5]
dim(ex)
ex <- ex[!is.na(ex$geneSymbol),]

N <- ex[,grepl("N", colnames(ex))]
D <- ex[,grepl("L", colnames(ex))]
dim(N)  # 43145  85
dim(D)  # 43145  85

N$probe <- ex$probe
D$probe <- ex$probe

N$geneSymbol <- ex$geneSymbol
D$geneSymbol <- ex$geneSymbol

#highest genes for N
colnames(N[1:87])
N$means <- rowMeans(N[1:85])
sort(rowMeans(N[,1:85]), decreasing = TRUE, index.return=T)$ix[1:10] #the first index of sorted means vector returned = 7087, 
#$x returns only the values when there's no rowname
rowMeans(N[1:85])[7087] #15.61 , for probeID: 200092_s_at which is the largest mean among genes means
N[which(round(N[,"means"], digits = 1) == 15.6), "means"]  #it doesn't work if you don't write the full digits of the value 


match(sort(N$means, decreasing = TRUE), N$means)
order(N$means, decreasing = TRUE)
N_sorted <- N[match(sort(rowMeans(N[,1:85]), decreasing=T), N$means),]   
N_500_lst <- N_sorted$probe[1:500]
N_500_lst

#highest genes for D
dim(D) # 43145  87
colnames(D)
D$means <- rowMeans(D[1:85])
sort(rowMeans(D[,1:85]), decreasing = TRUE, index.return=T)$ix[1:10] #the first index of sorted means vector returned = 6819, 
#$x returns only the values when there's no rowname
rowMeans(D[1:85])[7087] #15.565  for probeID: 200092_s_at , which is the largest mean among genes means
which(round(D[,"means"], digits = 2) == 15.57)  #it doesn't work if you don't write the full digits of the value 

D_sorted <- D[order(D$means, decreasing = TRUE),]
rowMeans(D_sorted[,1:85])
rownames(D_sorted)
#rownames(N_sorted) <- NULL
class(rownames(D_sorted) ) #character
D_500_lst <- D_sorted$probe[1:500]
D_500_lst


# intersectional and uniq lists of N and D
N_D_lst <- c(N_500_lst, D_500_lst)
length(N_D_lst) #1000

#uniq_N_D_exclu_lst <- N_D_lst[!(N_D_lst %in% inter_lst)]
#length(uniq_N_D_exclu_lst)  


uniq_N_D_lst <- unique(N_D_lst)  #top 1000 highly expressed 
length(uniq_N_D_lst) # 577


total_lst <- unique(c(uniq_N_D_lst, new2000_u_id))
length(total_lst) # 2554

## MINET datasets----

#Normal
mi_N <- NULL
colnames(N_sorted)[1:85]
mi_N <- N_sorted[total_lst,1:85]

# for (c in colnames(mi_N)[1:3]) {
#   mi_N[,c] <- round(mi_N[,c], digits = 2)
#   #print(round(mi_N[1:5, c], digits = 2))
# }

colnames(mi_N)
dim(mi_N)
class(mi_N)  # data.frame
dim(mi_N) # 2554  85

# dummy dataset for imputation
mi_N_imp <- mi_N[521:532, 1:5]
mi_N_imp[1,1] <- NA
grepl("NA", rownames(mi_N_imp))
is.na(rownames(mi_N[521:532, 1:5]))  #in this dataset NAs are Character
mi_N_imp <- mi_N_imp[!grepl("NA", rownames(mi_N_imp)),]
#na.omit(mi_N_imp) #omits rows with any number of missing values on columns specified, default: all
mi_N_imp <- impute.knn(as.matrix(mi_N_imp) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
mi_N_imp <- mi_N_imp$data

#imputation
mi_N_imp <- mi_N
dim(mi_N_imp)
grepl("NA", rownames(mi_N_imp))
is.na(rownames(mi_N))  #in this dataset NAs are Character
mi_N_imp <- mi_N_imp[!grepl("NA", rownames(mi_N_imp)),]
#na.omit(mi_N_imp) #omits rows with any number of missing values on columns specified, default: all
mi_N_imp <- impute.knn(as.matrix(mi_N_imp) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
mi_N_imp <- mi_N_imp$data
sum(is.na(mi_N_imp)) #0
sum(is.na(mi_N)) #6035


mi_N_imp_t <- as.data.frame(t(mi_N_imp))
dim(mi_N_imp_t)  #85 2819


#Diseased
mi_D <- NULL
colnames(D_sorted)[1:85]
mi_D <- D_sorted[total_lst,1:85]

colnames(mi_D)
rownames(mi_D)
class(mi_D)  # data.frame
mi_D[1:5, 1:5]
dim(mi_D) # 2554  85


#imputation
mi_D_imp <- mi_D
dim(mi_D_imp)
grepl("NA", rownames(mi_D_imp))
is.na(rownames(mi_D))  #in this dataset NAs are Character
mi_D_imp <- mi_D_imp[!grepl("NA", rownames(mi_D_imp)),]
#na.omit(mi_N_imp) #omits rows with any number of missing values on columns specified, default: all
mi_D_imp <- impute.knn(as.matrix(mi_D_imp) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
mi_D_imp <- mi_D_imp$data
sum(is.na(mi_D_imp)) #0
sum(is.na(mi_D)) 
dim(mi_D_imp)
as.data.frame(t(mi_D_imp))[1:5, 1:5]
mi_D_imp_t <- as.data.frame(t(mi_D_imp))
dim(mi_D_imp_t)
rownames(mi_D_imp_t)
# Run minet with different parameters ----
#input to minet: data.frame where columns contain variables/features
#and rows contain outcomes/samples.


library(minet)
#BiocManager::install("Rgraphviz")
#library(Rgraphviz)

# Normal ----

####
##run minet with discretization none and estimator spearman
net_N <- minet(mi_N_imp_t) #, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
# trace(minet, edit = TRUE)


# N_emp ----
##run with discretization equalfreq and estimator mi.empirical
net_N_emp <- minet(mi_N_imp_t, disc = "equalfreq", estimator = "mi.empirical") 
net_N_emp <- round(net_N_emp, 2)
write.csv(net_N_emp,"./data/adjacency_net_N_emp.csv")
#, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
#plot( as(net_N_emp ,"graphNEL") )
temp.emp <- net_N_emp[1:5, 1:5]
temp.emp
ind_mx <- which(temp.emp != 0,arr.ind = TRUE)
ind_mx
rownames(ind_mx)
ind_mx[,1]
ind_mx[,2]
temp.emp[ind_mx[,1], ind_mx[,2]]
temp.emp[ind_mx]


class(which(temp.emp == 0,arr.ind = TRUE)) # matrix array
temp.emp[which(temp.emp == 0, arr.ind=TRUE)]
colnames(net_N_emp[1:5, 1:5])
rownames(net_N_emp[1:5,1:5])
net_N_emp[1:5,1:5]
dim(net_N_emp) # 2483 2483
class(net_N_emp)  # matrix array
length( sum(net_N_emp == 0))
net_N_emp[net_N_emp == 0]


# create edgelist
# comb <- t(combn(colnames(net_N_emp[1:5,1:5]), 2))
# data.frame(comb, mi=net_N_emp[comb])
comb_N_emp <- t(combn(colnames(net_N_emp), 2))
head(comb_N_emp)
el_mi_N_emp <- data.frame(comb_N_emp, mi=net_N_emp[comb_N_emp])
length(el_mi_N_emp$mi[el_mi_N_emp$mi > 0.2])  #2267
el_mi_N_emp_2 <- el_mi_N_emp[el_mi_N_emp$mi >0.2,]
el_mi_N_emp_2$mi <- round(el_mi_N_emp_2$mi, 2)

#save mi>0.2 csv
write.csv(el_mi_N_emp_2, "./el_mi_N_emp_2.csv")

##create edge-list >0.1
el_mi_N_emp_1 <- el_mi_N_emp[el_mi_N_emp$mi >0.1,]
el_mi_N_emp_1$mi <- round(el_mi_N_emp_1$mi, 2)
#save mi>0.1 csv
write.csv(el_mi_N_emp_1, "./el_mi_N_emp_1.csv")


## unuique genes of MI > 0.4 
uniq_N_2<- c(unique(el_mi_N_emp_2$X1), unique(el_mi_N_emp_2$X2))
length(uniq_N_2) #3022
length(unique(uniq_N_2)) #2483  
uniq_N_2 <- unique(uniq_N_2)



#adjacency matrix of mi>0.2
net_N_emp_2 <- net_N_emp[colnames(net_N_emp) %in% uniq_N_2, rownames(net_N_emp) %in% uniq_N_2]
net_N_emp_2 <- round(net_N_emp_2, 2)
dim(net_N_emp_2) #2483 2483
net_N_emp_2[1:5, 1:5]
write.csv(net_N_emp_2, "./net_N_emp_2.csv")
#write.csv(net_N_emp_2[1:200,1:200], "./data/test_net_N_emp_2.csv")


# net_N_emp[which(net_N_emp !=0, arr.ind = TRUE) ]
# net_N_emp_lower <- lower.tri(net_N_emp[1:5,1:5], diag = FALSE)
# net_N_emp_lower
# net_N_emp_lower <- net_N_emp[1:5,1:5][net_N_emp_lower]
# net_N_emp_lower
# net_N_emp_lower[which(net_N_emp_lower != 0, arr.ind = TRUE)]


#net_N_emp_lower
net_N_emp_lower2 <- lower.tri(net_N_emp_2, diag = FALSE)
class(net_N_emp_lower2) # marix array
dim(net_N_emp_lower2) #2483 2483
net_N_emp_lower2[1:5, 1:5]
net_N_emp_2 <- net_N_emp_2[net_N_emp_lower2]
length(net_N_emp_2)
#net_N_emp_lower2 <- net_N_emp_lower2[which(net_N_emp_lower2 != 0, arr.ind = TRUE)]
length(net_N_emp_2[net_N_emp_2 > 0.2]) #2267
length(net_N_emp_lower[net_N_emp_lower > 0.1]) #52592
length(net_N_emp_lower[net_N_emp_lower > 0.05]) #728373
hist(net_N_emp_2[net_N_emp_2 > 0.2])

####
##minet with discretization equalfreq and estimator mi.mm
net_N_mm <- minet(mi_N_t, disc = "equalfreq", estimator = "mi.mm") 
#, method = "mrnet", estimator = "spearman", disc = "none", nbins = sqrt(nrow(mi_N)) 
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


# Diseased ----

# N_emp ----
##run with discretization equalfreq and estimator mi.empirical
net_D_emp <- minet(mi_D_imp_t, disc = "equalfreq", estimator = "mi.empirical") 
net_D_emp <- round(net_D_emp, 2)
write.csv(net_D_emp,"./adjacency_net_D_emp.csv")


# create edgelist
# comb <- t(combn(colnames(net_D_emp[1:5,1:5]), 2))
# data.frame(comb, mi=net_D_emp[comb])
comb_D_emp <- t(combn(colnames(net_D_emp), 2))
head(comb_D_emp)
el_mi_D_emp <- data.frame(comb_D_emp, mi=net_D_emp[comb_D_emp])
length(el_mi_D_emp$mi[el_mi_D_emp$mi > 0.2])  #2323
el_mi_D_emp_2 <- el_mi_D_emp[el_mi_D_emp$mi >0.2,]
el_mi_D_emp_2$mi <- round(el_mi_D_emp_2$mi, 2)

#save mi>0.2 csv
write.csv(el_mi_D_emp_2, "./el_mi_D_emp_2.csv")

##create > 0.1
el_mi_D_emp_1 <- el_mi_D_emp[el_mi_D_emp$mi >0.1,]
el_mi_D_emp_1$mi <- round(el_mi_D_emp_1$mi, 2)
#save mi>0.1 csv
write.csv(el_mi_D_emp_1, "./el_mi_D_emp_1.csv")


## Disease unuique genes of MI > 0.2 
uniq_D_2<- c(unique(el_mi_D_emp_2$X1), unique(el_mi_D_emp_2$X2))
length(uniq_D_2) #3165
length(unique(uniq_D_2)) #2554
uniq_D_2 <- unique(uniq_D_2)



#adjacency matrix of mi>0.2
net_D_emp_2 <- net_D_emp[colnames(net_D_emp) %in% uniq_D_2, rownames(net_D_emp) %in% uniq_D_2]
net_D_emp_2 <- round(net_D_emp_2, 2)
dim(net_D_emp_2) #2554 2554
net_D_emp_2[1:5, 1:5]
write.csv(net_D_emp_2, "./net_D_emp_2.csv")
#write.csv(net_D_emp_2[1:200,1:200], "./data/test_net_D_emp_2.csv")



#net_D_emp_lower
net_D_emp_lower2 <- lower.tri(net_D_emp_2, diag = FALSE)
class(net_D_emp_lower2) # marix array
dim(net_D_emp_lower2) # 2554  2554
net_D_emp_lower2[1:5, 1:5]
net_D_emp_2 <- net_D_emp_2[net_D_emp_lower2]
length(net_D_emp_2)
#net_D_emp_lower2 <- net_D_emp_lower2[which(net_D_emp_lower2 != 0, arr.ind = TRUE)]
length(net_D_emp_2[net_D_emp_2 > 0.2]) #2323

hist(net_D_emp_2[net_D_emp_2 > 0.2])


##
#******* on your cytoscape install all Dynet... apps, aMatReader and centiscape
save.image("mi.RData")
load("mi.RData")


# iGraph codes ----
#are in /Users/saman/Dropbox/systems_biology_projects/pso_cmm/pso1_complete.R

# Cytoscape ---- 
# edge-list csv files are in "/Users/saman/git/psoriasis1"
# Cytoscape: input edgelists via "import network from file system"



