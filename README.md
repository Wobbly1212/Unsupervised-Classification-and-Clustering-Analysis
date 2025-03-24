# 📊 Unsupervised Classification and Clustering Analysis

### Developed by: Mohammad Hassan Darabi  
🎓 **Bachelor’s Degree in Data Analytics**  
🏛️ **University of Campania Luigi Vanvitelli, Caserta, Italy**  
🗓️ **Academic Year: 2022-2023 (Second Semester)**

---

## 👨‍🏫 Instructor
**Prof. Fabrizio Maturo**  
Faculty of Economics, Universitas Mercatorum, Rome, Italy  
📧 Contact: [fabrizio.maturo@unimercatorum.it](mailto:fabrizio.maturo@unimercatorum.it)

---

## 🚀 Project Overview

This project provides an extensive exploration and practical comparison of various unsupervised learning techniques, particularly clustering algorithms. The primary goal is to demonstrate clustering patterns, determine the optimal number of clusters, and effectively visualize the results across different datasets.

---

## 📁 Project Structure

- 🌳 **Hierarchical Clustering** (Agglomerative and Divisive methods)
- 🎯 **K-means Clustering**
- 🔄 **Fuzzy Clustering**
- 📊 **PAM and CLARA** (Partitioning methods)
- 🧩 **Biclustering and Co-clustering**
- ✅ **Cluster Validation and PCA Analysis**

---

## 📚 Datasets Used

1. 🌸 **Iris Dataset**
   - **Variables:** Sepal Length, Sepal Width, Petal Length, Petal Width
   - **Purpose:** Introductory clustering and visualization.

2. 🧬 **Leukemia Genomic Data** (`plsgenomics` package)
   - **Variables:** Patient gene expression profiles
   - **Purpose:** High-dimensional biological data clustering.

3. 🐭 **Mouse Gene Expression Data** (`clValid` package)
   - **Variables:** Gene expression levels
   - **Purpose:** Validation of clustering methods (internal and stability measures).

4. 🍞 **Yeast Microarray Data** (`BicatYeast` dataset)
   - **Variables:** Gene expression across multiple experiments
   - **Purpose:** Biclustering method demonstration.

5. 🦠 **COVID-19 Italian Regional Data**
   - **Variables:** Total COVID-19 cases across Italian regions
   - **Purpose:** Spatio-temporal epidemiological clustering.

6. 🇮🇹 **BES 2015 Data** (Well-being indicators for Italian regions)
   - **Variables:** Socioeconomic and quality-of-life measures
   - **Purpose:** Advanced clustering methods and PCA.

---

## ⚙️ Methods and Techniques Explained

### 1. Hierarchical Clustering
- **Agglomerative Methods:** Complete, Single, Average linkage
- **Divisive Method:** DIANA

### 2. Partitioning Methods
- **K-means:** Identification of optimal clusters (Elbow, Silhouette, Gap statistics)
- **PAM (Partitioning Around Medoids)** and **CLARA (Clustering Large Applications)**

### 3. Fuzzy Clustering
- **Fuzzy K-means and FANNY:** Analyzing membership degrees
- Optimal cluster selection using MPC, PC, SIL.F measures

### 4. Biclustering and Co-clustering
- **Biclustering:** Cheng and Church algorithm for coherent clusters
- **Co-clustering:** Blockcluster method for simultaneous clustering

### 5. Cluster Validation and Rank Aggregation
- **Internal Validation:** Dunn index, Connectivity, Silhouette width
- **Stability Validation:** APN, AD, ADM
- **Rank Aggregation:** Consolidating multiple validation criteria

### 6. Principal Component Analysis (PCA)
- Dimensionality reduction for clearer cluster visualization

---

## 📈 Visualization Techniques

- 🌲 Dendrograms (standard and circular)
- 🔥 Heatmaps (standard and interactive)
- 📊 PCA scatter plots
- 📏 Parallel coordinate plots
- 📌 Multi-panel plots (boxplots, clustering outcomes)

---

## 📌 Dependencies

- `cluster`, `factoextra`, `ggplot2`, `ggpubr`
- `BBmisc`, `mlr`
- `plsgenomics`, `clValid`, `fclust`, `FactoMineR`
- `dendextend`, `circlize`, `gplots`, `d3heatmap`, `ggrepel`
- `biclust`, `blockcluster`, `RankAggreg`
- `gridExtra`, `sjPlot`, `e1071`, `readxl`, `dplyr`

Install dependencies using:
```r
install.packages(c("cluster", "factoextra", "ggplot2", "ggpubr", "BBmisc", "mlr",
                   "plsgenomics", "clValid", "fclust", "FactoMineR", "dendextend",
                   "circlize", "gplots", "d3heatmap", "ggrepel", "biclust", "blockcluster",
                   "RankAggreg", "gridExtra", "sjPlot", "e1071", "readxl", "dplyr"))
```

---

## 🖥️ How to Use
- Clone or download the repository.
- Install the required dependencies listed above.
- Execute scripts sequentially or independently, based on the chosen method or dataset.
- Consult script comments for detailed guidance.

---

## 🎓 Conclusion
This project offers an in-depth illustration of multiple clustering methodologies, validation strategies, and visualization techniques, serving as a valuable resource for educational and practical data analytics applications.

