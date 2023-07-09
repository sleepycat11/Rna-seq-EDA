library(here)

count_data <- read.delim(here("data", "GSM7166269_raw_count_data.txt"), sep='\t',header=T)
head(count_data)

sample_data <- read.delim(here("data", "GSM7166269_sample_data.txt"), sep='\t',header=T)
head(sample_data)

##processing count data
rownames(count_data) <- count_data$ENSG
count_data = count_data[2:7]
head(count_data)

##processing sample data
rownames(sample_data) <- colnames(count_data)
sample_data = sample_data[2]
head(sample_data)

##analysing data

##histogram
hist(count_data$Non_expanded_Tregs_1,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Non-expanded T_reg 1", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 1 Non-expanded T_reg cells")
hist(count_data$Non_expanded_Tregs_2,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Non-expanded T_reg 2", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 2 Non-expanded T_reg cells")
hist(count_data$Non_expanded_Tregs_3,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Non-expanded T_reg 3", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 3 Non-expanded T_reg cells")
hist(count_data$IMP_Tregs_1,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Expanded T_reg 1", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 1 Expanded T_reg cells")
hist(count_data$IMP_Tregs_2,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Expanded T_reg 2", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 2 Expanded T_reg cells")
hist(count_data$IMP_Tregs_3,xlim=c(0,5500),breaks=1000,prob=T,
     xlab="Expanded T_reg 3", cex.lab = 1.4, col = rgb(0.6,0,1,0.5),
     main="Histrogram for sample 3 Expanded T_reg cells")

##density line plot
plot(density(count_data$Non_expanded_Tregs_1),xlim=c(0,5000),
     main="Density plot for sample 1",col='blue')
lines(density(count_data$IMP_Tregs_1),xlim=c(0,5000),col='red')
legend(3000, 0.15, legend=c("non-expanded", "expanded"), fill = c("blue","red"))

plot(density(count_data$Non_expanded_Tregs_2),xlim=c(0,5000),col='blue',
     main="Density plot for sample 2")
lines(density(count_data$IMP_Tregs_2),xlim=c(0,5000),col='red')
legend(3000, 0.15, legend=c("non-expanded", "expanded"), fill = c("blue","red"))


plot(density(count_data$Non_expanded_Tregs_3),xlim=c(0,5000),col='blue',
     main="Density plot for sample 3")
lines(density(count_data$IMP_Tregs_3),xlim=c(0,5000),col='red')
legend(3000, 0.15, legend=c("non-expanded", "expanded"), fill = c("blue","red"))


##boxplot
boxplot(count_data$Non_expanded_Tregs_1,count_data$IMP_Tregs_1)
boxplot(count_data$Non_expanded_Tregs_2,count_data$IMP_Tregs_2)
boxplot(count_data$Non_expanded_Tregs_3,count_data$IMP_Tregs_3)

##loading DESeq2 library
library("DESeq2")

##making DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design = ~expansion)
dds

##finding row indexes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10

##filtering out rows with less than 10 reads
dds <- dds[keep,]
dds

##set factor level
dds$expansion <- relevel(dds$expansion, ref= "no")

##run DESeq analysis
prdds <- DESeq(dds)

##results
res <- results(prdds)
res

summary(res)
resultsNames(prdds)

##list of genes with 10% FDR cut-off and LFC threshold=1
list_set <- subset(res, padj <= 0.1)
list_set <- subset(list_set, abs(log2FoldChange)>=1)

list_set <- list_set[ order(list_set$log2FoldChange, decreasing = TRUE), ]
list_set

names <- rownames(list_set) ##final list
names
write.csv2(names, file = "genes.csv", row.names = F) ##saving list

##making MA plot
plotMA(res, main="MA plot")

##making volcano plot
library("ggplot2")

df = data.frame(res)

df$diffexpressed <- "no" ##adding a column
df$diffexpressed[ df$log2FoldChange > 1 & df$pvalue < 0.05] <- "upregulated"
df$diffexpressed[ df$log2FoldChange < -1 & df$pvalue < 0.05] <- "downregulated"

df$diffexpressed

ggplot(data=df, aes(x=log2FoldChange,y= -log10(pvalue),col=diffexpressed)) +
  geom_vline(xintercept = c(-1,1),col='black',linetype='dashed') +
  geom_hline(yintercept = 0.05,col='black',linetype='dashed') +
  geom_point() +
  scale_color_manual(values=c("#00AFBB","gray","#bb0c00"),
                     labels=c("Down-regulated","non-significant","Upregulated"))

