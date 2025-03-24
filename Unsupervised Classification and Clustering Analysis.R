# =============================================================================
# Unsupervised Classification and Clustering Analysis
# =============================================================================
# Developed by: Mohammad Hassan Darabi
# Bachelor’s Degree in Data Analytics
# University of Campania Luigi Vanvitelli, Caserta, Italy
# Academic Year: 2022-2023 (Second Semester)
#
# This script includes comprehensive unsupervised learning techniques 
# and hierarchical clustering methods applied to various datasets.
# =============================================================================

# --------------------
# Libraries Required
# --------------------
library(cluster)    
library(factoextra)
library(ggplot2)
library(BBmisc)
library(mlr)
library(plsgenomics)

# --------------------
# Global Settings
# --------------------
set.seed(123)

# ============================================================
# SECTION 1 - Hierarchical Clustering on Iris Dataset
# ============================================================

# Inspecting Iris dataset structure
dim(iris)
head(iris)
tail(iris)

# Selecting variables Petal Length and Petal Width
data <- iris[, 3:4]

# ------------------------------------------------------------
# (1.1) Euclidean Distance - Complete Linkage
# ------------------------------------------------------------

dd_euclidean <- dist(data, method = "euclidean")
clusters_complete <- hclust(dd_euclidean, method = "complete")

# Plot dendrogram
plot(clusters_complete, main = "Euclidean - Complete Linkage")

# Cut dendrogram into 3 clusters and visualize
cut_complete <- cutree(clusters_complete, k = 3)
table(cut_complete, iris$Species)
rect.hclust(clusters_complete, k = 3, border = 2:6)

# Visualization of Clusters
ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = Species), size = 4) +
  ggtitle("Real Iris Categories")

ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = factor(cut_complete)), size = 4) +
  ggtitle("Clustering Results (Euclidean-Complete)")

fviz_cluster(list(data = data, cluster = cut_complete))

# ------------------------------------------------------------
# (1.2) Minkowski Distance (p=3) - Average Linkage
# ------------------------------------------------------------

dd_minkowski <- dist(data, method = "minkowski", p = 3)
clusters_average <- hclust(dd_minkowski, method = "average")

# Plot dendrogram
plot(clusters_average, main = "Minkowski (p=3) - Average Linkage")

# Cut dendrogram into 3 clusters and visualize
cut_average <- cutree(clusters_average, k = 3)
table(cut_average, iris$Species)
rect.hclust(clusters_average, k = 3, border = 2:6)

# Visualization of Clusters
ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = Species), size = 4) +
  ggtitle("Real Iris Categories")

ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = factor(cut_average)), size = 4) +
  ggtitle("Clustering Results (Minkowski-Average)")

# ------------------------------------------------------------
# (1.3) Manhattan Distance - Single Linkage
# ------------------------------------------------------------

dd_manhattan <- dist(data, method = "manhattan")
clusters_single <- hclust(dd_manhattan, method = "single")

# Plot dendrogram
plot(clusters_single, main = "Manhattan - Single Linkage")

# Cut dendrogram into 3 clusters and visualize
cut_single <- cutree(clusters_single, k = 3)
table(cut_single, iris$Species)
rect.hclust(clusters_single, k = 3, border = 2:6)

# Visualization of Clusters
ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = Species), size = 4) +
  ggtitle("Real Iris Categories")

ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = factor(cut_single)), size = 4) +
  ggtitle("Clustering Results (Manhattan-Single)")

# ============================================================
# SECTION 2 - Additional Distance Metrics (Reference)
# ============================================================

# Note: For future use/reference only, not executed here.
# library(multivariate)
# dcov(x, y, index = 1.0)
# dcor(x, y, index = 1.0)  # correlation distance
# DCOR(x, y, index = 1.0)    

# ============================================================
# SECTION 3 - Data Normalization
# ============================================================

# Check for missing data
any(is.na(iris))

# Normalization (range 0 to 1)
normalized_data <- normalize(iris, method = "range", range = c(0,1))

# ============================================================
# SECTION 4 - Genomic Data Clustering (Leukemia Dataset)
# ============================================================

# Install and load plsgenomics package
# (installation commented if already installed)
# install.packages("plsgenomics")
library(plsgenomics)

# Load Leukemia dataset
data(leukemia)
# Inspect Leukemia data
class(leukemia)
dim(leukemia$X)
head(leukemia$X)

# Manhattan distance with average linkage
dd_leukemia <- dist(leukemia$X, method = "manhattan")
clusters_leukemia <- hclust(dd_leukemia, method = "average")

# Plot dendrogram
plot(clusters_leukemia, main = "Leukemia Data - Manhattan-Average")

# Cut dendrogram into 2 clusters and compare with true labels
cut_leukemia <- cutree(clusters_leukemia, k = 2)
table(cut_leukemia, leukemia$Y)
rect.hclust(clusters_leukemia, k = 2, border = 2:3)

# Visualization using first two dimensions
fviz_cluster(list(data = leukemia$X, cluster = cut_leukemia))

# ============================================================
# SECTION 5 - Divisive Hierarchical Clustering (DIANA method)
# ============================================================

# DIANA clustering on Iris (Petal Length and Width)
clusters_diana <- diana(data)

# Print and visualize results
print(clusters_diana)
plot(clusters_diana, main = "Divisive Clustering - Iris Data")

# Divisive coefficient (higher means clearer structure)
clusters_diana$dc

# =============================================================================
# SECTION 6 - Determining the Optimal Number of Clusters
# =============================================================================

library(factoextra)
library(gridExtra)

# Using Iris dataset (Petal Length and Width)
data <- iris[, 3:4]

# Elbow method, Silhouette, and Gap statistic
p1 <- fviz_nbclust(data, FUN = hcut, method = "wss", k.max = 10) +
  ggtitle("(A) Elbow method")

p2 <- fviz_nbclust(data, FUN = hcut, method = "silhouette", k.max = 10) +
  ggtitle("(B) Silhouette method")

p3 <- fviz_nbclust(data, FUN = hcut, method = "gap_stat", k.max = 10) +
  ggtitle("(C) Gap statistic")

# Display plots side by side
grid.arrange(p1, p2, p3, nrow = 1)

# =============================================================================
# SECTION 7 - Heatmap Visualization (Leukemia Data)
# =============================================================================

heatmap(leukemia$X,
        scale = "column",
        xlab = "Genes",
        ylab = "Patients",
        main = "Leukemia Heatmap")

# Check specific group condition (optional)
leukemia$Y == 2

# =============================================================================
# SECTION 8 - Circular Dendrogram Visualization
# =============================================================================

library(dendextend)
library(circlize)

# Dendrogram using Iris dataset
dd_euclidean <- dist(data, method = "euclidean")
clusters_avg <- hclust(dd_euclidean, method = "ave")

dend <- as.dendrogram(clusters_avg)
dend <- color_branches(dend, k = 3)

circlize_dendrogram(dend)

# =============================================================================
# SECTION 9 - Interactive Heatmap Visualization
# =============================================================================

library(gplots)
library(d3heatmap)

# Interactive heatmap using Iris dataset (first 4 variables)
gplots::heatmap.2(as.matrix(iris[,1:4]), scale = "row")

d3heatmap(as.matrix(iris[,1:4]),
          dendrogram = "row",
          Rowv = dend,
          colors = "Greens",
          width = 900,
          height = 700,
          show_grid = TRUE)

# =============================================================================
# SECTION 10 - Multiple ggplot Visualizations (Boxplots)
# =============================================================================

library(ggpubr)

# Boxplot for Iris variables by Species
p1 <- ggplot(iris, aes(x=Species, y=Sepal.Length, fill=Species)) + 
  geom_boxplot()

p2 <- ggplot(iris, aes(x=Species, y=Sepal.Width, fill=Species)) + 
  geom_boxplot()

p3 <- ggplot(iris, aes(x=Species, y=Petal.Length, fill=Species)) + 
  geom_boxplot()

p4 <- ggplot(iris, aes(x=Species, y=Petal.Width, fill=Species)) + 
  geom_boxplot()

# Combine boxplots into a single visualization
ggarrange(p1, p2, p3, p4 + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

# =============================================================================
# SECTION 11 - Comparing Clustering Algorithms (Optional Reference)
# =============================================================================

# This section is for reference and exploration purposes only.
# Uncomment and run if needed.

# library(dendextend)
# 
# hclust_methods <- c("single", "complete", "average")
# iris_dendlist <- dendlist()
# 
# for(i in seq_along(hclust_methods)){
#   hc_iris <- hclust(dist(iris[,1:4]), method = hclust_methods[i])   
#   iris_dendlist <- dendlist(iris_dendlist, as.dendrogram(hc_iris))
# }
# 
# names(iris_dendlist) <- hclust_methods
# iris_dendlist

# =============================================================================
# SECTION 12 - Cluster Validation (Reference)
# =============================================================================

# Reference notes on validation using `fpc` package
# Uncomment if needed.

# library(fpc)
# cluster.stats(d, fit1$cluster, fit2$cluster)

# d: Distance matrix
# fit1$cluster, fit2$cluster: Integer vectors from two clusterings

# =============================================================================
# SECTION 13 - COVID-19 Data Clustering (Italy, 2020)
# =============================================================================

# Read data directly (without crawler)
COVID_ITA_REG <- read.csv("COVID_ITA_REG.csv")
View(COVID_ITA_REG)

# Inspect regions
COVID_ITA_REG$denominazione_regione[1:21]

data <- COVID_ITA_REG
data$denominazione_regione <- as.factor(data$denominazione_regione)

# Rename regions for consistency
levels(data$denominazione_regione)[levels(data$denominazione_regione) == "P.A. Trento"] <- "Trento"
levels(data$denominazione_regione)[levels(data$denominazione_regione) == "P.A. Bolzano"] <- "Bolzano"

# Calculate dimensions and prepare matrix
days <- nrow(data) / 21
nreg <- length(levels(data$denominazione_regione))

contagiati <- matrix(NA, nreg, days)

# Populate matrix with cases per region
for(i in 1:nreg) {
  contagiati[i, ] <- data$totale_casi[seq(i, 567, by = nreg)]
}

# Convert matrix to dataframe
contagiati <- as.data.frame(contagiati)
colnames(contagiati) <- paste0("Day_", 1:days)
rownames(contagiati) <- levels(data$denominazione_regione)

# Perform clustering on COVID-19 cases
dd_covid <- dist(contagiati, method = "euclidean")
clusters_covid <- hclust(dd_covid, method = "ave")

# Dendrogram plot
plot(clusters_covid, main = "COVID-19 Clustering (Italian Regions)")

# Cut into 2 clusters
cut_covid <- cutree(clusters_covid, k = 2)
rect.hclust(clusters_covid, k = 2, border = 2:3)

# Evaluate the optimal number of clusters
p1 <- fviz_nbclust(contagiati, FUN = hcut, method = "wss", k.max = 10) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(contagiati, FUN = hcut, method = "silhouette", k.max = 10) +
  ggtitle("(B) Silhouette method")
p3 <- fviz_nbclust(contagiati, FUN = hcut, method = "gap_stat", k.max = 10) +
  ggtitle("(C) Gap statistic")

grid.arrange(p1, p2, p3, nrow = 1)

# =============================================================================
# SECTION 14 - K-means Clustering
# =============================================================================

# ----- (14.1) Iris Dataset (2 Variables) -----

data <- iris[, 3:4]

# Standardize data
data_scaled <- as.data.frame(scale(data))

# Apply K-means clustering (3 centers)
kk <- kmeans(data_scaled, centers = 3, iter.max = 10, nstart = 1)

# Inspect results
kk$centers
kk$cluster
kk$iter
kk$withinss
kk$tot.withinss
kk$betweenss
kk$totss

# Visualization: real vs. predicted clusters
ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = factor(Species)), size = 4) +
  ggtitle("Real Iris Categories")

ggplot(iris, aes(Petal.Length, Petal.Width)) +
  geom_point(aes(colour = factor(kk$cluster)), size = 4) +
  ggtitle("K-means Clustering Results")

fviz_cluster(kk, data_scaled, ellipse.type = "norm")

# ----- (14.2) Genomic Data (Leukemia) -----

library(plsgenomics)
data(leukemia)

scaled_leukemia <- as.data.frame(scale(leukemia$X))

# K-means clustering (5 clusters)
leukemia_kk <- kmeans(leukemia$X, centers = 5, iter.max = 50, nstart = 1)

leukemia_kk$cluster
leukemia_kk$iter
leukemia_kk$tot.withinss

fviz_cluster(leukemia_kk, leukemia$X)

# Determine optimal number of clusters
p1 <- fviz_nbclust(scaled_leukemia, kmeans, method = "wss", k.max = 10) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(scaled_leukemia, kmeans, method = "silhouette", k.max = 10) +
  ggtitle("(B) Silhouette method")
p3 <- fviz_nbclust(scaled_leukemia, kmeans, method = "gap_stat", k.max = 10) +
  ggtitle("(C) Gap statistic")

grid.arrange(p1, p2, p3, nrow = 1)

# K-means with 2 clusters for visualization
leukemia_kk_2 <- kmeans(leukemia$X, centers = 2, iter.max = 50, nstart = 1)
fviz_cluster(leukemia_kk_2, leukemia$X)

# =============================================================================
# SECTION 15 - PAM and CLARA Clustering (Leukemia Data)
# =============================================================================

# ----- (15.1) PAM -----

leukemia_pam <- pam(leukemia$X, 2, metric = "euclidean", stand = FALSE)

leukemia_pam$clusinfo
leukemia_pam$id.med
leukemia_pam$clustering

fviz_cluster(leukemia_pam, leukemia$X)

# ----- (15.2) CLARA -----

leukemia_clara <- clara(leukemia$X, 2, metric = "euclidean", stand = FALSE, samples = 5)
fviz_cluster(leukemia_clara, leukemia$X)

# =============================================================================
# SECTION 16 - Cluster Validity Analysis (Internal and Stability Measures)
# =============================================================================

library(clValid)

data("mouse")
dim(mouse)

express <- mouse[, c("M1", "M2", "M3", "NC1", "NC2", "NC3")]
rownames(express) <- mouse$ID

# ----- (16.1) Internal Validation -----

intern <- clValid(express, 2:6, 
                  clMethods = c("hierarchical", "kmeans", "diana", "pam"),
                  validation = "internal")

summary(intern)

# Plot internal validation results
op <- par(no.readonly=TRUE)
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(intern, legend=FALSE)
plot(nClusters(intern), measures(intern,"Dunn")[,,1], type="n", axes=F, xlab="", ylab="")
legend("center", clusterMethods(intern), col=1:9, lty=1:9, pch=paste(1:9))
par(op)

# ----- (16.2) Stability Measures -----

stab <- clValid(express, 2:6, 
                clMethods = c("hierarchical","kmeans","diana","pam"),
                validation = "stability")

summary(stab)

# Stability validation plots
op <- par(no.readonly=TRUE)
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(stab, measure=c("APN","AD","ADM"), legend=FALSE)
plot(nClusters(stab), measures(stab,"APN")[,,1], type="n", axes=F, xlab="", ylab="")
legend("center", clusterMethods(stab), col=1:9, lty=1:9, pch=paste(1:9))
par(op)

# ----- (16.3) Rank Aggregation -----

library(RankAggreg)

result <- clValid(express, 4:6, 
                  clMethods=c("hierarchical","kmeans","pam"),
                  validation=c("internal","stability"))

res <- getRanksWeights(result)

print(res$ranks[,1:3], quote=FALSE)

# Rank aggregation analysis
if(require("RankAggreg")) {
  CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
  print(CEWS)
  plot(CEWS)
}

# =============================================================================
# SECTION 17 - Fuzzy K-means (BES Data Analysis)
# =============================================================================

library(fclust)
library(FactoMineR)
library(sjPlot)
library(dplyr)
library(readxl)

# Load BES dataset
BES_2015 <- read_excel("BES_2015.xlsx")
data <- data.frame(BES_2015)

rownames(data) <- data[,1]
data <- data[,-1]

data2015 <- select(data, contains("2015"))

# Descriptive statistics
sjt.corr(data2015, corr.method="pearson")

# Determine optimal number of clusters

# Modified Partition Coefficient (MPC)
prove1 <- numeric(9)
for (i in 2:10) {
  fclus_2015 <- FKM(data2015, i)
  prove1[i-1] <- MPC(fclus_2015$U)
}

# Partition Coefficient (PC)
prove2 <- numeric(9)
for (i in 2:10) {
  fclus_2015 <- FKM(data2015, i)
  prove2[i-1] <- PC(fclus_2015$U)
}

# SIL.F measure
prove3 <- numeric(9)
for (i in 2:10) {
  fclus_2015 <- FKM(data2015, i)
  prove3[i-1] <- SIL.F(data2015, fclus_2015$U)
}

# Combine results and visualize
provedata <- data.frame(MPC=prove1, PC=prove2, SIL_F=prove3)

x <- 2:10
plot(x, provedata$MPC, type="l", lty=1, lwd=2, ylim=c(0, max(provedata)+0.5),
     xlim=c(1,10), main="Optimal Number of Groups", 
     xlab="Number of Groups", ylab="Index Value")
lines(x, provedata$PC, type="l", lty=2, lwd=2)
lines(x, provedata$SIL_F, type="l", lty=3, lwd=2)
legend("topright", legend=c("MPC", "PC", "SIL.F"), lty=1:3, lwd=2)

# =============================================================================
# SECTION 18 - Fuzzy Analysis (FANNY Algorithm)
# =============================================================================

fannyx <- fanny(data2015, 2)

summary(fannyx)
plot(fannyx)

fanni <- fanny(data2015, 2, diss=FALSE, memb.exp=2, metric="euclidean",
               stand=FALSE, maxit=500)

fviz_cluster(fanni, repel=TRUE)

fanni$membership
# =============================================================================
# SECTION 19 - Principal Component Analysis (PCA)
# =============================================================================

library(FactoMineR)

# PCA on BES 2015 data
pc <- PCA(data2015, ncp = 2, scale.unit = TRUE)

# Eigenvalues and variance explained
pc$eig
pc$var

# =============================================================================
# SECTION 20 - Fuzzy Cluster Plot (Italian Regions)
# =============================================================================

library(ggplot2)
library(ggrepel)

# PCA Coordinates
coo <- pc$ind$coord
dataf <- data.frame(x1 = coo[,1], x2 = coo[,2])

# Fuzzy cluster visualization
fuzzyplot <- ggplot(dataf, aes(x = x1, y = x2, colour = fanni$membership[,1])) +
  geom_point(size = 3) +
  scale_colour_gradient("Membership Degree (Group 1)", low = "green", high = "red") +
  ggtitle("Fuzzy Clustering of Italian Regions") +
  geom_text_repel(aes(label = rownames(data2015)), size = 4)

fuzzyplot

# =============================================================================
# SECTION 21 - Biclustering with 'biclust' Package
# =============================================================================

library(biclust)

# Load example yeast data
data(BicatYeast)
dim(BicatYeast)

# Cheng and Church Algorithm (Biclustering)
XCC <- biclust(BicatYeast, method = BCCC(), delta = 0.01, alpha = 1.5, number = 5)

# Visualizations of biclusters
parallelCoordinates(x = BicatYeast, bicResult = XCC, number = 1)
parallelCoordinates(x = BicatYeast, bicResult = XCC, number = 2)
parallelCoordinates(x = BicatYeast, bicResult = XCC, number = 3)
parallelCoordinates(x = BicatYeast, bicResult = XCC, number = 4)

drawHeatmap(x = BicatYeast, result = XCC, bicluster = 4)

# =============================================================================
# SECTION 22 - Co-clustering (BlockCluster Package)
# =============================================================================

library(blockcluster)

set.seed(3101)

# Default co-clustering strategy
default_strategy <- coclusterStrategy()
summary(default_strategy)

# Modify co-clustering strategy
new_strategy <- coclusterStrategy(nbtry = 5, nbxem = 10, algo = "XCEMStrategy")

# Co-clustering on binary data
data("binarydata", package = "blockcluster")
coclust_binary <- cocluster(binarydata, datatype = "binary", nbcocluster = c(2, 3))
summary(coclust_binary)
plot(coclust_binary, asp = 0)
plot(coclust_binary, type = "distribution")

# Co-clustering on contingency data
data("contingencydataunknown", package = "blockcluster")
coclust_contingency <- cocluster(contingencydataunknown, datatype = "contingency", nbcocluster = c(2, 3))
summary(coclust_contingency)
plot(coclust_contingency)
plot(coclust_contingency, type = "distribution")

# Co-clustering on continuous data
data("gaussiandata", package = "blockcluster")
coclust_gaussian <- cocluster(gaussiandata, datatype = "continuous", nbcocluster = c(2, 3))
summary(coclust_gaussian)
plot(coclust_gaussian, asp = 0)
plot(coclust_gaussian, type = "distribution")

# Image segmentation (example with 'snake.jpg')
library(ReadImages)
# x <- read.jpeg("snake.jpg") 
# coclust_image <- cocluster(x, datatype = "continuous", nbcocluster = c(3, 3))
# plot(coclust_image)

# Document clustering (example with 'textdatamedcran.dat')
# x <- read.csv("textdatamedcran.dat", header = FALSE)
# x1 <- t(x[, 2:9276])
# coclust_text <- cocluster(x1, datatype = "contingency", nbcocluster = c(2, 2))
# predicted_class <- coclust_text["colclass"]
# original_class <- x[, 1]
# misclassified <- min(sum(abs(predicted_class - original_class)), 2431 - sum(abs(predicted_class - original_class)))

# =============================================================================
# SECTION 23 - BES Data Exercise (Hierarchical and K-means Clustering)
# =============================================================================

library(readxl)
library(cluster)
library(factoextra)

# Load BES 2015 dataset
BES_2015 <- read_excel("BES_2015.xlsx")
data <- as.data.frame(BES_2015)
rownames(data) <- data[, 1]
data <- data[, -1]

# Select 2015 data
data2015 <- dplyr::select(data, contains("2015"))

# Hierarchical clustering (average linkage)
dist_matrix <- dist(data2015, method = "euclidean")
hc_clusters <- hclust(dist_matrix, method = "ave")

# Dendrogram visualization
plot(hc_clusters, main = "Hierarchical Clustering (BES 2015)")

# Determine optimal clusters
p1 <- fviz_nbclust(data2015, FUN = hcut, method = "wss", k.max = 10) +
  ggtitle("(A) Elbow method")
p2 <- fviz_nbclust(data2015, FUN = hcut, method = "silhouette", k.max = 10) +
  ggtitle("(B) Silhouette method")
p3 <- fviz_nbclust(data2015, FUN = hcut, method = "gap_stat", k.max = 10) +
  ggtitle("(C) Gap statistic")

grid.arrange(p1, p2, p3, nrow = 1)

# Cutting dendrogram into 2 clusters
cut_clusters <- cutree(hc_clusters, k = 2)
rect.hclust(hc_clusters, k = 2, border = 2:6)

# Heatmap visualization
heatmap(as.matrix(data2015), scale = "column", xlab = "BES Indicators", ylab = "Regions",
        main = "Heatmap of BES Indicators")

# ----- K-means Clustering (BES Data) -----

kmeans_BES <- kmeans(data2015, centers = 2, iter.max = 10, nstart = 1)

# Cluster inspection
kmeans_BES$cluster
kmeans_BES$iter
kmeans_BES$tot.withinss

# K-means clustering visualization
fviz_cluster(kmeans_BES, data2015, ellipse.type = "norm")

# PCA visualization for K-means interpretation
pc_BES <- PCA(data2015, ncp = 2, scale.unit = TRUE)

