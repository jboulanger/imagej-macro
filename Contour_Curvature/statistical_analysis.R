 PCA plot for morphological quantifications in brain organoids
filename <- ''
#how to open a matrix (csv) file d55 <- read.csv(filename, row.names=1)
matrix(c(97428.33333, 1139.534808, 94.2894174, 0, 90350.33333, 1097.37621, 94.42556463, 0, 95945, 1125.084258, 95.26941112, 0, 151497.5, 1382.628484, 99.57739267, 0.1), nrow = 4, ncol = 4)
data.matrix <- matrix(c(97428.33333, 1139.534808, 94.2894174, 0, 90350.33333, 1097.37621, 94.42556463, 0, 95945, 1125.084258, 95.26941112, 0, 151497.5, 1382.628484, 99.57739267, 0.1), nrow = 4, ncol = 4)
rownames(data.matrix) <- paste(c("Area", "Perimeter", "Circularity", "Inflection points"), sep="")
colnames(data.matrix) <- paste(c("B1C1", "B1C2", "B1C3", "B1C4"), sep="")
head(data.matrix)
pca <- prcomp(t(data.matrix), scale=TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x),
X=pca$x[,1],
Y=pca$x[,2])
pca.data
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
geom_text() +
xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
theme_bw() +
ggtitle("My PCA Graph")
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores)
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_4_genes <- names(gene_score_ranked[1:4])
top_4_genes

#How to make subclusters in the PCA data
library(stringr)
pca.data$batch = rownames(pca.data) #to make a subcolumns called batch with the rownames
batches = sub(".C*", "", pca.data$batch)
batches = substr(pca.data.play$batch, 1,2) #to take the first two digits
batches = str_extract(pca.data$batch, "[^_]+")#to take everything before the first underscore
pca.data$batches = batches
single = sub("^[^_]*_", "", pca.data$batch) #to extract everything after the first underscore
pca.data.SW318

#How to make a ggplot with dots and colors
ggplot(data=pca.data, aes(x=X, y=Y, color=pca.data$batches, label=single)) +
 geom_point(size=6) + geom_text(vjust = -1, hjust= 0.4, size = 4) +
 xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
 ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
 theme_bw()
#Aesthetics
pca.data$batches <- factor(pca.data$batches, levels = c("Kit", "Fast", "noMG", "Diss")) #to change the order of the samples
geom_point(size=3) + geom_text(vjust = -0.2, hjust= 0.2, size = 2.5) + #vjust highness of the label, hjust right or left
print(i+ labs(colour = "condition")) #to change the clusters label

#How to make a correlation Plot
library(devtools)
eig.val <- get_eigenvalue(pca)
eig.val
res.var <- get_pca_var(pca)
res.var$coord
res.var$contrib
res.var$cos2
library(corrplot)
corrplot(res.var$cos2, is.corr=FALSE)
corrplot(res.var$cos2, is.corr=FALSE, tl.col="black", tl.srt=45)