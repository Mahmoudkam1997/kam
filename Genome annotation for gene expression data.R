#Differential expression anlysis and gene annotation 
# needed packages
install.packages("BiocManager")#one time
install.packages("matrixTests")#one time
BiocManager::install("genefilter")#one time
BiocManager::install("org.Hs.eg.db")#one time
BiocManager::install("AnnotationDbi")#one time


library(BiocManager)#one time
library(matrixTests)#each time you open your R
library(genefilter)#each time you open your R
library(org.Hs.eg.db)
######################################
#Reading and manipulating the data file
data = read.table("D:/SRP029880.raw_counts.tsv")
data = data[,-11]
dim(data)
data= as.matrix(data)
class(data)
#explore if is there any missing expression value (empty cell)
sum(is.null(data))

#explore the data distribution using the histogram plot 
#if there is a problem with the data
hist(data, col = "orange", main="Histogram")

#scaling the data using log2 transformation to better visulization 
#scaling variability of the data
# we use (+1) to avoid the infinity character when we log zero valus 
#logs according to scaling you want
hist(log2(data+1), col = "orange", main="Histogram")

#filter low count genes which have row mean lower than 1, (indexing)
#>1, [row, data][50,30], reomve the mean lower than 1
data=data[rowMeans(data) > 1,]

####calculating the fold change####

#caculate the logged mean for each group
#log2 fold change = 10/5 = 2 down regulation in gene expression
#log2 is better for subtraction and is accurate for clarrity of difference 
#no dividing
#1 means the mean of rows, 2 means the mean of columns 
tum.mean = apply((log2(data+1))[,1:5], 1, mean)
norm.mean = apply((log2(data+1))[,6:dim(data)[2]], 1, mean)
#calculate the fold change by taking the difference between the two means
#the difference between logged means equl to the fold change insted of using
#division of unlogged data
fold=tum.mean-norm.mean

#view the distribution of the fold change
hist(fold)

#creat a phenotype table as its rows contain a phenotype ethier tumor or
#normal corrseponding to the columns in the expression data; so as we know 
#that the first 50 column in the expression data are tumor so the first 50
#row in the phenotype will be labeld tum and the other 50 norm
phenotype = as.data.frame(factor(rep(c("Case","Ctrl"), c(5, 5))))
colnames(phenotype)= "grouping"

#making the hypothesis testing using T test for each row(gene)
#first mean - second mean / stdv
t=rowttests(data,phenotype$grouping)

#correct the T test p value using the FDR method
# fdr means each p value / total p values
p.adj=p.adjust(t$p.value, method = "fdr")

#save the result in a dataframe contain the fold change and p adjusted valu
result=as.data.frame(cbind(fold , p.adj))

#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
res.deg=result[result$p.adj < 0.05 & abs(result$fold)>2,]

#export the Degs into your current folder for further analysthis
write.csv(as.matrix(res.deg),file="res.degs.csv", quote=F,row.names=T)
##################################
data = read.csv("D:/res.degs.csv", row.names = 1)
#################################
#for symbol and gene name if it doesn't exist
colimn(org.Hs.eg.db)
annots <- select(org.Hs.eg.db, 
                 keys=row.names(data), 
                 columns=("GENENAME"), keytype="SYMBOL")
###################################
#merging 
resultTable <- merge(data, annots, by.x=0, by.y="SYMBOL")
head(resultTable)
######################################
write.csv(as.matrix(resultTable),file="resultTable.csv", quote=F,row.names=T)
#################
