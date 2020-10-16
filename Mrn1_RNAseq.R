options(stringsAsFactors=FALSE)

## load in tsv files from Kalliso output:
KO1pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO1pre.tsv", header=TRUE, sep = "\t")
KO2pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO2pre.tsv", header=TRUE, sep = "\t")
KO3pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO3pre.tsv", header=TRUE, sep = "\t")

WT1pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT1pre2.tsv", header=TRUE, sep = "\t")
WT2pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT2pre.tsv", header=TRUE, sep = "\t")
WT3pre <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT3pre.tsv", header=TRUE, sep = "\t")

KO1post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO1post.tsv", header=TRUE, sep = "\t")
KO2post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO2post.tsv", header=TRUE, sep = "\t")
KO3post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/KO3post.tsv", header=TRUE, sep = "\t")

WT1post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT1post.tsv", header=TRUE, sep = "\t")
WT2post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT2post.tsv", header=TRUE, sep = "\t")
WT3post <- read.delim("/mnt/ingolialab/kswain/NIKS010/RNAseq/Kall_out/WT3post.tsv", header=TRUE, sep = "\t")

## make dataframe with pre-diauxic shift counts for two samples:
pre_counts <- KO1pre[,c("target_id", "est_counts")]
colnames(pre_counts) <- c("Yorf", "KO1pre")
head(pre_counts)
## add in counts from other samples
pre_counts$KO2pre <- KO2pre[match(pre_counts$Yorf, KO2pre$target_id), "est_counts"]
pre_counts$KO3pre <- KO3pre[match(pre_counts$Yorf, KO3pre$target_id), "est_counts"]
pre_counts$WT1pre <- WT1pre[match(pre_counts$Yorf, WT1pre$target_id), "est_counts"]
pre_counts$WT2pre <- WT2pre[match(pre_counts$Yorf, WT2pre$target_id), "est_counts"]
pre_counts$WT3pre <- WT3pre[match(pre_counts$Yorf, WT3pre$target_id), "est_counts"]

## round counts 
library(dplyr)
pre_counts <- pre_counts %>% 
  mutate_if(is.numeric, round)
head(pre_counts)

plot(log10(pre_counts$KO1pre), log10(pre_counts$KO2pre))
cor(pre_counts$KO1pre, pre_counts$KO2pre)
cor(pre_counts$KO2pre, pre_counts$KO3pre)

plot(log10(pre_counts$WT1pre), log(pre_counts$WT2pre))
plot(log10(pre_counts$WT2pre), log(pre_counts$WT3pre))
cor(pre_counts$WT1pre, pre_counts$WT2pre)


## run DESeq2 on precount df
library(DESeq2)

rownames(pre_counts) <- pre_counts$Yorf
pre_counts$Yorf = NULL
head(pre_counts)

conditions = data.frame(row.names=c("KO1pre", "KO2pre", "KO3pre", "WT1pre", "WT2pre", "WT3pre"), 
                        geno=factor(c("KO", "KO", "KO", "WT", "WT", "WT"), levels=c("WT", "KO")))
conditions

dds <- DESeqDataSetFromMatrix(countData = pre_counts, colData = conditions, design =  ~ geno)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

resultsNames(dds)
res <- results(dds, name="geno_KO_vs_WT")
res= as.data.frame(res)
head(res)

KOvWT_pre <- res
head(KOvWT_pre)
KOvWT_pre <- cbind(rownames(pre_counts), data.frame(KOvWT_pre, row.names = NULL))
head(KOvWT_pre)
names(KOvWT_pre)[1]<- "Yorf"

## load in sgd dataset for gene names:
if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))
KOvWT_pre$gene <- sgd[match(KOvWT_pre$Yorf, sgd$name), "gene"]
KOvWT_pre$desc <- sgd[match(KOvWT_pre$Yorf, sgd$name), "desc"]


## run DEseq for Post-diauxic shift WT v KO 
post_counts <- KO1post[,c("target_id", "est_counts")]
colnames(post_counts) <- c("Yorf", "KO1post")
head(post_counts)
## add in counts from other samples
post_counts$KO2post <- KO2post[match(post_counts$Yorf, KO2post$target_id), "est_counts"]
post_counts$KO3post <- KO3post[match(post_counts$Yorf, KO3post$target_id), "est_counts"]
post_counts$WT1post <- WT1post[match(post_counts$Yorf, WT1post$target_id), "est_counts"]
post_counts$WT2post <- WT2post[match(post_counts$Yorf, WT2post$target_id), "est_counts"]
post_counts$WT3post <- WT3post[match(post_counts$Yorf, WT3post$target_id), "est_counts"]

post_counts <- post_counts %>% 
  mutate_if(is.numeric, round)
head(post_counts)

cor(post_counts$KO1post, post_counts$KO2post)
cor(post_counts$KO1post, post_counts$KO3post)
cor(post_counts$WT1post, post_counts$WT2post)
cor(post_counts$WT1post, post_counts$WT3post)

rownames(post_counts) <- post_counts$Yorf
post_counts$Yorf = NULL
head(post_counts)

conditions = data.frame(row.names=c("KO1post", "KO2post", "KO3post", "WT1post", "WT2post", "WT3post"), 
                        geno=factor(c("KO", "KO", "KO", "WT", "WT", "WT"), levels=c("WT", "KO")))
conditions

dds <- DESeqDataSetFromMatrix(countData = post_counts, colData = conditions, design =  ~ geno)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

resultsNames(dds)
res <- results(dds, name="geno_KO_vs_WT")
res= as.data.frame(res)
head(res)

KOvWT_post <- res
KOvWT_post <- cbind(rownames(post_counts), data.frame(KOvWT_post, row.names = NULL))
head(KOvWT_post)
names(KOvWT_post)[1]<- "Yorf"

KOvWT_post$gene <- sgd[match(KOvWT_post$Yorf, sgd$name), "gene"]
KOvWT_post$desc <- sgd[match(KOvWT_post$Yorf, sgd$name), "desc"]
write.csv(KOvWT_post, "~/NIKS010/DESeq_out/KOvWT_post.csv")



## redo with WT pre v post shift: 
WT_counts <- pre_counts
WT_counts <- cbind(rownames(WT_counts), data.frame(pre_counts, row.names = NULL))
names(WT_counts)[1] <- "Yorf"
WT_counts$KO1pre = NULL
WT_counts$KO2pre = NULL
WT_counts$KO3pre = NULL
WT_counts$WT1post <- WT1post[match(WT_counts$Yorf, WT1post$target_id), "est_counts"]
WT_counts$WT2post <- WT2post[match(WT_counts$Yorf, WT2post$target_id), "est_counts"]
WT_counts$WT3post <- WT3post[match(WT_counts$Yorf, WT3post$target_id), "est_counts"]
head(WT_counts)

WT_counts <- WT_counts %>% 
  mutate_if(is.numeric, round)
head(WT_counts)

KO_counts <- pre_counts
KO_counts <- cbind(rownames(KO_counts), data.frame(pre_counts, row.names=NULL))
names(KO_counts)[1] <- "Yorf"
KO_counts$WT1pre=NULL
KO_counts$WT2pre=NULL
KO_counts$WT3pre=NULL
KO_counts$KO1post <- KO1post[match(KO_counts$Yorf, KO1post$target_id), "est_counts"]
KO_counts$KO2post <- KO1post[match(KO_counts$Yorf, KO2post$target_id), "est_counts"]
KO_counts$KO3post <- KO1post[match(KO_counts$Yorf, KO3post$target_id), "est_counts"]
head(KO_counts)

KO_counts <- KO_counts %>% 
  mutate_if(is.numeric, round)
head(KO_counts)

## 
rownames(WT_counts) <- WT_counts$Yorf
WT_counts$Yorf = NULL
head(WT_counts)

conditions = data.frame(row.names=c("WT1pre", "WT2pre", "WT3pre", "WT1post", "WT2post", "WT3post"), 
                        time=factor(c("pre", "pre", "pre", "post", "post", "post"), levels=c("pre", "post")))
conditions

dds <- DESeqDataSetFromMatrix(countData = WT_counts, colData = conditions, design =  ~ time)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

resultsNames(dds)
res <- results(dds, name="time_post_vs_pre")
res= as.data.frame(res)
head(res)

WT_prevpost <- res
WT_prevpost <- cbind(rownames(WT_prevpost), data.frame(WT_counts, row.names=NULL))
names(WT_prevpost)[1] <- "Yorf"

WT_prevpost$gene <- sgd[match(WT_prevpost$Yorf, sgd$name), "gene"]
WT_prevpost$desc <- sgd[match(WT_prevpost$Yorf, sgd$name), "desc"]
write.csv(WT_prevpost, "~/NIKS010/DESeq_out/WT_prevpost.csv")
##



## Repeat for KOvWT post
rownames(KO_counts) <- KO_counts$Yorf
KO_counts$Yorf = NULL
head(KO_counts)

conditions = data.frame(row.names=c("KO1pre", "KO2pre", "KO3pre", "KO1post", "KO2post", "KO3post"), 
                        time=factor(c("pre", "pre", "pre", "post", "post", "post"), levels=c("pre", "post")))
conditions

dds <- DESeqDataSetFromMatrix(countData = KO_counts, colData = conditions, design =  ~ time)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

resultsNames(dds)
res <- results(dds, name="time_post_vs_pre")
res= as.data.frame(res)
head(res)

KO_prevpost <- res
KO_prevpost <- cbind(rownames(KO_prevpost), data.frame(KO_counts, row.names=NULL))
names(KO_prevpost)[1] <- "Yorf"
KO_prevpost$gene <- sgd[match(KO_prevpost$Yorf, sgd$name), "gene"]
KO_prevpost$desc <- sgd[match(KO_prevpost$Yorf, sgd$name), "desc"]
write.csv(WT_prevpost, "~/NIKS010/DESeq_out/KO_prevpost.csv")

pre_counts <- KO1pre[,c("target_id", "est_counts")]
colnames(pre_counts) <- c("Yorf", "KO1pre")
head(pre_counts)
## add in counts from other samples
pre_counts$KO2pre <- KO2pre[match(pre_counts$Yorf, KO2pre$target_id), "est_counts"]
pre_counts$KO3pre <- KO3pre[match(pre_counts$Yorf, KO3pre$target_id), "est_counts"]
pre_counts$WT1pre <- WT1pre[match(pre_counts$Yorf, WT1pre$target_id), "est_counts"]
pre_counts$WT2pre <- WT2pre[match(pre_counts$Yorf, WT2pre$target_id), "est_counts"]
pre_counts$WT3pre <- WT3pre[match(pre_counts$Yorf, WT3pre$target_id), "est_counts"]

## round counts 
library(dplyr)
pre_counts <- pre_counts %>% 
  mutate_if(is.numeric, round)
head(pre_counts)



## plot data: 
library(ggplot2)

## plotting replicates: Pre-counts
png(file="~/NIKS010/DESeq_out/WT1vWT2pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(WT1pre), log10(WT2pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep1 counts)", y="Log10(WT rep2 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/WT2vWT3pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(WT2pre), log10(WT3pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep2 counts)", y="Log10(WT rep3 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/WT1vWT3pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(WT1pre), log10(WT3pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep1 counts)", y="Log10(WT rep3 counts)")
dev.off()


png(file="~/NIKS010/DESeq_out/KO1vKO2pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(KO1pre), log10(KO2pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep1 counts)", y="Log10(Δmrn1 rep2 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/KO2vKO3pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(KO2pre), log10(KO3pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep2 counts)", y="Log10(Δmrn1 rep3 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/KO1vKO3pre.png", width=1000,height=1000,res=144)
ggplot(pre_counts, aes(log10(KO1pre), log10(KO3pre))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene pre-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep1 counts)", y="Log10(Δmrn1 rep3 counts)")
dev.off()

## plotting replicates: post-counts 
head(post_counts)
png(file="~/NIKS010/DESeq_out/WT1vWT2post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(WT1post), log10(WT2post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep1 counts)", y="Log10(WT rep2 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/WT1vWT3post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(WT1post), log10(WT3post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep1 counts)", y="Log10(WT rep3 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/WT3vWT2post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(WT3post), log10(WT2post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("WT counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(WT rep3 counts)", y="Log10(WT rep2 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/KO1vKO2post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(KO1post), log10(KO2post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep1 counts)", y="Log10(Δmrn1 rep2 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/KO1vKO3post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(KO1post), log10(KO3post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep1 counts)", y="Log10(Δmrn1 rep3 counts)")
dev.off()

png(file="~/NIKS010/DESeq_out/KO3vKO2post.png", width=1000,height=1000,res=144)
ggplot(post_counts, aes(log10(KO3post), log10(KO2post))) + geom_point(size=0.5, col="black") + theme_minimal() + ggtitle("Δmrn1 counts per gene post-diauxic shift") +
  theme(axis.title = element_text(size=14), plot.title = element_text(size = 16)) +  labs(x="Log10(Δmrn1 rep3 counts)", y="Log10(Δmrn1 rep2 counts)")
dev.off()

## plotting KOvWT data, Pre- and post-diauxic shift: 

## plot KOvWT pre count data:
KOvWT_pre <- mutate(KOvWT_pre, sig=ifelse(KOvWT_pre$padj < 0.05, "P-adj < 0.05", "P-adj > 0.05"))
plot_KOvWTpre <- ggplot(KOvWT_pre, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(size=0.5, aes(col=sig))  + xlim(-5, 5) + ylim(0, 111) +
  scale_color_manual(values=c("cyan2",  "grey80")) +
  theme_minimal() + ggtitle("Pre-diauxic shift: Δmrn1 v WT") +
  theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16), legend.position = c(0.85, 0.90),
        legend.box.background = element_rect(color="grey", size=0.5), legend.title=element_blank())
plot_KOvWTpre + geom_text_repel(data=filter(KOvWT_pre, padj < 0.05 & log2FoldChange < -0.5 | padj < 0.05 & log2FoldChange > 0.5), aes(label=gene), size=2.5, segment.alpha=0.15)

## plot KOvWT post count data: 
KOvWT_post <- mutate(KOvWT_post, sig=ifelse(KOvWT_post$padj < 0.05, "P-adj < 0.05", "P-adj > 0.05"))
plot_KOvWTpost <- ggplot(KOvWT_post, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(size=0.5, aes(col=sig))  + xlim(-5, 5) + ylim(0, 111) +
  scale_color_manual(values=c("hotpink2",  "grey80")) +
  theme_minimal() + ggtitle("Post-diauxic shift: Δmrn1 v WT") +
  theme(axis.title = element_text(size = 16), plot.title = element_text(size = 16), legend.position = c(0.85, 0.90),
        legend.box.background = element_rect(color="grey", size=0.5), legend.title=element_blank())
plot_KOvWTpost + geom_text_repel(data=filter(KOvWT_post, padj < 0.05 & log2FoldChange < -0.75 | padj < 0.05 & log2FoldChange > 0.75), aes(label=gene), size=2.5, segment.alpha=0.15)


