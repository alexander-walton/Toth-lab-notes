## RNA-seq analysis with DESeq2
## Largely based on Stephen Turner, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

#out of date #source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

setwd("~/Box Sync/bionformatics training with Severin/abundance_counts_foundresses")

dat<-read.table("Foundress_count.txt",header = T,quote = "",row.names = 1)

# Convert to matrix
dat <- as.matrix(dat)
head(dat)

#wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {options(width=as.integer(howWide))}
#wideScreen()
#options(width=as.integer(150))


condition<-factor(c("Non", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Foundress", "Foundress", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Non", "Foundress", "Non", "Foundress", "Non"))
condition=relevel(condition,ref = "Non")

location<-factor(c("gal", "gal", "gal",  "keu",  "one",  "tc3",  "gal",  "gal",  "trac", "gal",  "che",  "tc3",  "gal",  "gal",  "sam",  "dry",  "jen",  "gal", "gal",  "trac", "air",  "fil" ))

jh.treat<-factor(c("jh", "jh", "jh", "ace", "jh", "jh",  "ace", "ace", "ace", "jh", "ace", "ace", "ace", "jh", "jh", "jh", "jh", "jh", "jh", "jh", "jh", "ace"))

#date wasp was sampled
culled<-factor(c("2.25.19", "2.25.19", "3.18.19", "3.18.19", "2.25.19", "3.18.19", "2.25.19", "2.25.19", "2.25.19", "2.25.19", "3.18.19", "3.18.19", "3.18.19", "2.25.19", "3.18.19", "3.18.19", "2.25.19", "3.18.19", "2.25.19", "2.25.19", "3.18.19", "3.18.19"))

#mass, pre-overwintering
mass<-factor(c(0.149, 0.122, 0.144, 0.163, 0.142, 0.177, 0.159, 0.127, 0.143, 0.140, 0.111, 0.152, 0.166, 0.171, 0.165, 0.148, 0.208, 0.154, 0.122, 0.143, 0.120, 0.156))

#the hibernaculum the wasp overwintered in
group<-factor(c("R", "H", "Q", "D", "E", "J", "R", "M", "L", "B", "K", "J", "N", "M", "D", "T", "E", "A", "H", "L", "I", "G"))


#what if there are potential co-foundresses, even though we never saw the particular wasp on the nest? What if there non-foundresses from boxes where there are NO nests vs. where a different wasp built a nest in the same box actually means something? Here I'm defining non-foundresses from nest boxes that contained a nest as "co".
co<-factor(c("co", "co", "non", "found", "found", "found", "found", "co", "found", "found", "non", "co", "found", "found", "found", "found", "co", "non", "found", "co", "found", "non"))

#The amount of time since a foundress started a nest may be important. The nests were of different ages when they were culled, so I have binned nest age into young (1-7 days), mid (8-14 days), or old (15-21 days).
age<-factor(c("non", "non", "non", "old", "mid", "old", "young", "non", "mid", "young", "non", "non", "young", "young", "mid", "young", "non", "non", "young", "non", "old", "non"))

##Mass is not suitable as a factor, so I've binned the wasps into size categories: small (0.12-0.13 g), m (0.14-0.159 g), and large (0.16-0.210 g)

size<-factor(c("m", "s", "m", "l", "m", "l", "m", "s", "m", "m", "s", "m", "l", "l", "l", "m", "l", "m", "s", "m", "s", "m"))

coldata.2<-data.frame(row.names=colnames(dat), condition, location)
head(coldata.2)

coldata.3<-data.frame(row.names=colnames(dat), condition, culled)
head(coldata.3)

coldata.4<-data.frame(row.names=colnames(dat), condition, mass)
head(coldata.4)

coldata.5<-data.frame(row.names=colnames(dat), condition, group)
head(coldata.5)

coldata.6<-data.frame(row.names=colnames(dat), condition, jh.treat)
head(coldata.6)

coldata.7<-data.frame(row.names = colnames(dat), condition, size)
head(coldata.7)


coldata.8<-data.frame(row.names = colnames(dat), condition, co, size, group, jh.treat, age, location)
head(coldata.8)

##### DESEq pipeline, first the design and the next step, normalizing to model fitting
#controlling for location, 
#as per http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata.2, design= ~  location + condition)
dds.orig<-dds

dds <- DESeq(dds)

# Plot Dispersions:
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Principal Components Analysis
plotPCA(rld)


# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Get differential expression results
res <- results(dds)
table(res$padj<0.1)

## Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)


idxl<-sort(coldata.2[,2], index.return = TRUE)$ix
test<-sort(coldata.2[,2], index.return = TRUE)
attributes(test)

idx<-sort(as.matrix(condition), index.return = TRUE)$ix
idx
dat["PFUS07749", idx]

dat["PFUS01025", idx]
dat["PFUS06137", idx]

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts["PFUS07749", idx]
round(normalized_counts["PFUS07749", idx])
round(normalized_counts["PFUS07749", order(location)])

normalized_counts["PFUS10952", idx]
round(normalized_counts["PFUS10952", idx])

normalized_counts["PFUS06137", idx]
round(normalized_counts["PFUS06137", idx])

#### DESEq pipeline, first the design and the next step, normalizing to model fitting
#controlling for cull date (which has the same result as just "condition" without controlling), and other factors

#mass isn't a factor, so that's a problem, right?
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata.8, design= ~ size + age)
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)
write.csv(res, file = "p-values for nest age by size.csv")

dat["PFUS03948", idx]
dat["PFUS02161", idx]
dat["PFUS08871", idx]
dat["PFUS09839", idx]
dat["PFUS11390", idx]



dat["PFUS09093", idx]


##
###Rsubread

BiocManager::install("Rsubread")
library(Rsubread)


#Mulitmapping
setwd("~/Box Sync/bionformatics training with Severin/abundance_counts_foundresses_multimapping/counts3")
datm<-read.table("Foundress_multimap_count.txt",header = T,quote = "",row.names = 1)

# Convert to matrix
datm <- as.matrix(datm)
head(datm)

coldata.m<-data.frame(row.names = colnames(datm), condition, co, size, group, jh.treat, age, location)

ddsm <- DESeqDataSetFromMatrix(countData = datm, colData = coldata.m, design= ~location)

ddsm <- DESeq(ddsm)
# Get differential expression results
resm <- results(ddsm)
table(resm$padj<0.05)

## Order by adjusted p-value
resm <- resm[order(resm$padj), ]
head(resm)



###5-26-20
# QuasiSeq


BiocManager::install('edgeR')
BiocManager::install('fields')
BiocManager::install('QuasiSeq')

library(QuasiSeq)
#put in raw reads
setwd("~/Box Sync/bionformatics training with Severin/abundance_counts_foundresses")

uniq<-read.table("Foundress_count.txt",header = T,quote = "",row.names = 1)

#filt to only truely expressed genes  ---what does this mean?
dataIn<-uniq
dataIn2<-dataIn[which((dataIn[,1]+dataIn[,2]+dataIn[,3])>5 | (dataIn[,4]+dataIn[,5]+dataIn[,6])>5),]
dataIn3<-dataIn2[which((sign(dataIn2[,1])+sign(dataIn2[,2])+sign(dataIn2[,3])>1)| (sign(dataIn2[,4])+sign(dataIn2[,5])+sign(dataIn2[,6])>1)),]
dataIn<-as.matrix(dataIn3)

#20705
#set up Models
condition<-factor(c("Non", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Foundress", "Foundress", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Non", "Foundress", "Non", "Foundress", "Non"))
condition=relevel(condition,ref = "Non")
location<-factor(c("gal", "gal", "gal",  "keu",  "one",  "tc3",  "gal",  "gal",  "trac", "gal",  "che",  "tc3",  "gal",  "gal",  "sam",  "dry",  "jen",  "gal", "gal",  "trac", "air",  "fil" ))
age<-factor(c("non", "non", "non", "old", "mid", "old", "young", "non", "mid", "young", "non", "non", "young", "young", "mid", "young", "non", "non", "young", "non", "old", "non"))
size<-factor(c("m", "s", "m", "l", "m", "l", "m", "s", "m", "m", "s", "m", "l", "l", "l", "m", "l", "m", "s", "m", "s", "m"))

n<-length(condition)
mn<-rep(1,length(condition))
design.list<-vector('list',3) #don't know what this is
design.list1<-model.matrix(~location+condition)
design.list2<-model.matrix(~condition)
design.list3<-model.matrix(~location)
log.offset=log(apply(dataIn,2,quantile,.75)) #what is this? 
# calculation
fit2<-QL.fit(dataIn, design.list3,log.offset=log.offset, Model="NegBin")

#look at results
results<-QL.results(fit2)
hist(results$P.values$QLShrink[,1],breaks=seq(0,1,.005))$counts
hist(results$P.values$QLShrink[,2],breaks=seq(0,1,.005))$counts
hist(results$Q.values$QLShrink[,1],breaks=seq(0,1,0.1))$counts
write.csv(dataIn[which(results$Q.values$QLShrink[,1]<0.2),],file="/Users/lintian/Desktop/2590genes_Q200.csv")
write.csv(dataIn[which(results$Q.values$QLShrink[,1]<0.155),],file="/Users/lintian/Desktop/500genes_Q155.csv")



#6-3-20 combining my data and Andrew Legan's data for a broader comparison. I am just focusing on Andrew's samples that were heads, to match my data. Although, he did separately analyze antennae and heads, which I did not.

setwd("~/Box Sync/bionformatics training with Severin/combined data from me and Legan")

datl<-read.table("waltonlegan_heads_count.txt",header = T,quote = "",row.names = 1)
datl <- as.matrix(datl)
head(datl)

conditionl<-factor(c("Gyne", "Gyne", "Gyne", "Gyne", "Gyne", "Gyne", "Male", "Male", "Male", "Male", "Male", "Male", "Non", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Foundress", "Foundress", "Non", "Non", "Foundress", "Foundress", "Foundress", "Foundress", "Non", "Non", "Foundress", "Non", "Foundress", "Non"))
conditionl=relevel(conditionl,ref = "Non")

coldata.l<-data.frame(row.names=colnames(datl), conditionl)
head(coldata.l)
coldata.l

ddsl <- DESeqDataSetFromMatrix(countData = datl, colData = coldata.l, design= ~ conditionl)


ddsl <- DESeq(ddsl)

# Plot Dispersions:
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(ddsl, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rldl <- rlogTransformation(ddsl)
head(assay(rldl))
hist(assay(rldl))

# Sample distance heatmap
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(conditionl))])
sampleDists <- as.matrix(dist(t(assay(rldl))))

library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[conditionl], RowSideColors=mycols[conditionl],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Get differential expression results
resl <- results(ddsl)
table(res$padj<0.1)


## Order by adjusted p-value
resl <- resl[order(res$padj), ]
head(resl)

