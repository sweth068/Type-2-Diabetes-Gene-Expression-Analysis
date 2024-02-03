library(GEOquery)

my_id <- "GDS3782"
dat <- getGEO(my_id)

#Extracting the table containing gene expression data
K <- dat@dataTable@table

#Checking the dimensions of the table
dim(K)




##Removing the 1st column in the dataset, making the IDENTIFIER column as the row names and 
##filtering out the duplicate genes and putting it into another dataframe.

#Assigning dataframe K to x
x <- K

#Removing the 1st column in dataframe x
x <- x[, -1]

#Removing row names
rownames(x) <- NULL

#Creating a new dataframe 'unique_df' by removing duplicated rows based on the 'IDENTIFIER' column
unique_df <- x[!duplicated(x$IDENTIFIER),]

#Setting row names of 'unique_df' as the unique values in the 'IDENTIFIER' column
rownames(unique_df) <- unique_df$IDENTIFIER

unique_df$IDENTIFIER <- NULL

head(unique_df)




#Creating a vector with new column names for control and diabetic samples

new_column_names <- c("GSM524151_Control","GSM524152_Control","GSM524153_Control","GSM524154_Control","GSM524155_Control","GSM524156_Control","GSM524157_Control","GSM524158_Control","GSM524159_Control","GSM524160_Control","GSM524161_Diab","GSM524162_Diab","GSM524163_Diab","GSM524164_Diab","GSM524165_Diab","GSM524166_Diab","GSM524167_Diab","GSM524168_Diab","GSM524169_Diab","GSM524170_Diab")  

#Assigning the new column names to the columns of the dataframe 'unique_df'
colnames(unique_df) <- new_column_names




###Normalizing the data
#Performing log2 normalization on the dataframe 'unique_df'

norm_data <- log2(unique_df)
#head(norm_data)

#Calculating pairwise correlation coefficients between columns with complete observations
data.cor <- cor(norm_data, use = 'pairwise.complete.obs')

#Creating a layout for the plot with 5 rows and 2 columns
layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = TRUE))

#Setting outer margins to a smaller value
par(oma = c(3, 0, 3, 0)) 

#Defining a custom color palette for the correlation plot
cx <- rev(colorRampPalette(c("red", "white", "blue"))(25))

#Creating a sequence of values for the legend
leg <- seq(min(data.cor, na.rm = TRUE), max(data.cor, na.rm = TRUE), length.out = 10)

#Creating the correlation plot
image(data.cor, main = "Correlation plot of control & type2 diabetic samples",
      axes = FALSE, col = cx)

#Adding x and y axes with labels
axis(1, at = seq(0, 1, length.out = ncol(data.cor)), labels = dimnames(data.cor)[[2]],
     cex.axis = 0.9, las = 2)
axis(2, at = seq(0, 1, length.out = ncol(data.cor)), labels = dimnames(data.cor)[[2]],
     cex.axis = 0.9, las = 2)

#Creating the color legend with a smaller plot region
#Adjusting the plot region as needed
par(plt = c(0.1, 0.5, 0.1, 0.5))

#Creating a color legend for the correlation values
image(as.matrix(leg), col = cx, axes = FALSE, xlab = 'Pearson correlation')

#Adding x-axis to the color legend with correlation values
tmp <- round(leg, 2)
axis(1, at = seq(0, 1, length.out = length(leg)), labels = tmp, cex.axis = 1)




##Outlier Assessment
#Calculating the mean expression level for each gene (row) and storing it in the variable 'test'
test <- rowMeans(norm_data)

#Calculating quantile values (25%, 50%, and 75%) of the mean expression levels 
quantile_value <- quantile(test, probs = c(.25, .5, .75),na.rm = TRUE)

#Printing the calculated quantile values
print(quantile_value)

#Identifying the row indices of genes with expression levels above the 25% quantile value
test.rows <- which(test > quantile_value[1])

#Creating a new dataframe 'data' by selecting only the rows with expression levels above the 25% quantile
data <- norm_data[test.rows, ]




##Plotting an Average correlation plot to find outliers

#Calculating the correlation matrix for the 'data' dataframe
data.cor <- cor(data, use = 'pairwise.complete.obs')

#Calculating the average correlation for each row in 'data.cor' and storing it in 'dat_new.avg'
dat_new.avg <- apply(data.cor, 1, mean)

#Setting plot parameters for margins
par(oma = c(3, 0.1, 0.1, 0.1))

#Creating a plot with specified parameters
plot(c(1, length(dat_new.avg)), range(dat_new.avg), type = "n", xlab = "", ylab = "Avg r",
     main = "Avg correlation of Control/Type2 Diabetic samples", axes = FALSE)

#Adding points to the plot using 'dat_new.avg' with blue color, white border, and circle shape
points(dat_new.avg, bg = "blue", col = 1, pch = 21, cex = 1.25)

#Adding x-axis labels with gene names
axis(1, at = c(1:length(dat_new.avg)), labels = dimnames(data)[[2]], las = 2, cex.lab = 0.4, cex.axis = 0.6)

#Adding y-axis
axis(2)

#Adding vertical lines at regular intervals
abline(v = seq(0.5, 62.5, 1), col = "grey")




##Finding outliers using hierarchical clustering

#Transposing the 'data' dataframe to create a new dataframe 'dat_new'
dat_new <- t(data)

#Calculating the Euclidean distance between rows of 'dat_new' and creating a distance matrix
dat_new.dist <- dist(dat_new, method = "euclidean")

#Performing hierarchical clustering on the distance matrix using the single linkage method
dat_new.clust <- hclust(dat_new.dist, method = "single")

#Plotting the hierarchical clustering dendrogram
#Labels are set as the names of data points, and label size is adjusted to 0.75
plot(dat_new.clust, labels = names(dat_new), cex = 0.75)




##From the above 2 visualizations - we can see that GSM524165_Diab and GSM524170_Diab are outlier 
##samples and needs to be removed before further analysis

#Creating a new dataframe 'data_no_outlier' by excluding specific columns (GSM524165_Diab, GSM524170_Diab) from the 'data' dataframe
data_no_outlier <- subset(data, select = -c(GSM524165_Diab, GSM524170_Diab))

#Displaying the dimensions (number of rows and columns) of the 'data_no_outlier' dataframe
dim(data_no_outlier)

#Displaying the column names of the 'data_no_outlier' dataframe
colnames(data_no_outlier)




##Differential Testing
#Creating a vector 'labels' to assign group labels (control and diabetic)
labels <- c(rep("control", 10), rep("diabetic", 8))

#Conducting an independent two-sample t-test on the expression levels of the first gene for both classes
result_for_1_gene <- t.test(data_no_outlier[1, labels == 'control'], data_no_outlier[1, labels == 'diabetic'], var.equal = FALSE)

result_for_1_gene

# Conducting an independent two-sample t-test for all genes between the 'control' and 'diabetic' groups
result_for_all_genes <- t.test(data_no_outlier[c(1:10), labels == 'control'], data_no_outlier[c(11:18), labels == 'diabetic'], var.equal = FALSE)

result_for_all_genes 




###To find gene with p-values under threshold 
#Define a function 'my.ttest' to perform t-tests for a gene and return the p-value
my.ttest <- function(v, labels) {
  levels <- unique(labels)
  v1 <- v[labels == levels[1]]
  v2 <- v[labels == levels[2]]
  pval <- t.test(v1, v2, var.equal = FALSE)$p.value
  pval
}

#Apply the 'my.ttest' function to calculate p-values for all genes and store them in 'allpvalues'
allpvalues <- apply(data_no_outlier, 1, my.ttest, labels)

#Adjust p-values for multiplicity using the Benjamini-Hochberg (BH) method
padj <- p.adjust(allpvalues, method = "BH")

#Extract genes with p-values less than the threshold value 0.05
sig.padj <- padj[padj < 0.05]
sig.padj.inds <- which(padj < 0.05)

#Plot a histogram to visualize the distribution of p-values for significant genes
hist(sig.padj, col = "darkmagenta", xlab = "p-values", main = "Genes with p-values under threshold", cex.main = 0.9)

#Grab all the rows of genes with p-values under the threshold from 'data_no_outlier'
select.genes <- data_no_outlier[c(sig.padj.inds), ]




###Dimensionality reduction methods
#Checking the class of the 'select.genes' object
class(select.genes)

#Unlisting the 'select.genes' object to convert it into a vector
select.genes <- unlist(select.genes)

#Subsetting the 'data_no_outlier' dataframe with the previously determined significant genes
data_pvalue <- data_no_outlier[c(select.genes),]

#Displaying the column names of the 'data_pvalue' dataframe
colnames(data_pvalue)

#Rename the first 10 columns of 'data_pvalue' to "c" (control) and columns 11 to 18 to "s" (sample)
colnames(data_pvalue)[1:10] <- "c"
colnames(data_pvalue)[11:18] <- "s"

#Perform hierarchical cluster analysis on the transposed data
dat.hca <- hclust(dist(t(data_pvalue), "man"), method = "median")

#Plot the dendrogram of the hierarchical cluster analysis
plot(dat.hca, main = "HCA of control and sample data with selected genes")




###classification method to classify the samples into their respective classes

#Assigning the column names of the dataframe into a variable called clas
clas <- colnames(data_no_outlier)

#Assigning colors based on sample groups (Control, Diabetes, others)
group_colors <- ifelse(grepl("Control$", clas), "black", ifelse(grepl("Diab$", clas), "red", "blue"))

#Performing Principal Component Analysis (PCA) on the transposed data
data_no_outliers_pca <- prcomp(t(data_no_outlier), cor = FALSE)

#Extracting the loadings of the first three principal components
data_no_outliers_loadings <- data_no_outliers_pca$x[, 1:3]

#Generating a scatter plot of the first two principal components with group colors
plot(
  data_no_outliers_loadings[, 1],
  data_no_outliers_loadings[, 2],
  col = group_colors,
  pch = 16,
  cex = 1.5,
  xlab = 'PC1',
  ylab = 'PC2',
  main = 'PCA plot of data_no_outliers'
)

#Adding a legend to the plot
legend("topright", legend = c("Control", "Diabetes"), col = c("black", "red"), pch = 16)




###Gene information

#Calculating ranks for genes based on adjusted p-values (padj)
gene_ranks <- rank(padj)

#Selecting the top 5 genes with the highest ranks
top_5_ranks <- head(sort(gene_ranks, decreasing = TRUE), 5)

#Printing the names of the top 5 genes
print(top_5_ranks)


#Calculating ranks for genes based on adjusted p-values (padj)
gene_ranks <- rank(padj)

#Selecting the top 5 genes with the lowest ranks
top_5_ranks <- head(sort(gene_ranks, decreasing = FALSE), 5)

#Printing the names of the top 5 genes with the lowest ranks
print(top_5_ranks)