check_expression <- function(row) {           
	samples <- 0
	for (col in 1:length(row)) {
		#Count samples for each gene above threshold
		if (row[col] > log2(15)) {
			samples <- samples + 1
		}
	}
	#For 20% of values above log2(15)
	return ((samples / length(row)) > 0.2)
}

check_variance <- function(row,median_sd) {
	df <- length(row) - 1
	testchi = df*(sd(row)/(median_sd))^2
	chiupper = qchisq((0.99)/2, df, lower.tail = FALSE)
	return (testchi > chiupper)
}

separate <- function(clustdata, filterdata, genenum, num) {
	newmatrix <- matrix(,genenum,0)
	for (x in 1:numsamples) {
		if (clustdata[x] == num) {
			newmatrix <- cbind(newmatrix,filterdata[x])
		}
	}
	return (newmatrix)
}

get_subtypes <- function() {
	colors <- c()
	for (x in 1:numsamples) {
		if (annomatrix[x] == "C3") {
			colors <- c(colors,"red")
		}
		else {
			colors <- c(colors,"blue")
		}
	}
	return (colors)
}

run_t_test <- function() {
	for (x in 1:numgenes) {
		t.test(cluster1[x,],cluster2[x,])
	}
}

read_gedata <- read.csv("/projectnb/bf528/users/group4/project1/code/combat_csv",sep=",")
gedata <- read_gedata[,-1]
rownames(gedata) <- read_gedata[,1]
read_annomatrix <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")	
annomatrix = read_annomatrix$cit.coloncancermolecularsubtype
filter1 <- apply(gedata,1,function(x) check_expression(x))
express <- gedata[filter1,]
median_sd <- median(apply(gedata,1,function(x) sd(x)))
filter2 <- apply(express,1,function(x) check_variance(x,median_sd))
variance <-express[filter2,]
filter3 <- apply(variance,1,function(x) (sd(x)/mean(x)) > 0.186)
covar <- variance[filter3,] #Gene expression data, fully-filtered
write.csv(covar,"filtered_gedata.csv")
numgenes = nrow(covar)
numsamples = ncol(covar)


distdata <- dist(t(covar))		
dendrogram <- hclust(distdata,method="complete")
clusters <- cutree(dendrogram, k=2)
cluster1 <- separate(clusters,covar,numgenes,1)
cluster2 <- separate(clusters,covar,numgenes,2)
heatmap.2(as.matrix(covar),ColSideColors = get_subtypes(),xlab="Patient tumor samples",ylab="Microarray probesets",labRow=c(""),labCol=c(""),tracecol="white")tstatistics <- unlist(lapply(1:numgenes, function(x) t.test(cluster1[x,],cluster2[x,])$statistic))
pvalues <- unlist(lapply(1:numgenes, function(x) t.test(cluster1[x,],cluster2[x,])$p.value))
diffexpress <- data.frame(ID = rownames(covar), t_statistic = tstatistics, p_value = pvalues, adjusted_p_value = p.adjust(pvalues,"fdr"))
write.csv(diffexpress,"tttest.csv")
sigfilter <- diffexpress$adjusted_p_value < 0.05
siggenes <- diffexpress[sigfilter,]
sortedgenes <- siggenes[order(siggenes$adjusted_p_value),]
negfilter <- sortedgenes$t_statistic<0
negatives <-  sortedgenes[negfilter,]
sortednegs <-  negatives[order(negatives$adjusted_p_value),]


#For biologist
varfilter <- apply(gedata,1,function(x) check_variance(x,median_sd))
varonly <- gedata[varfilter,]
write.csv(varonly,"var_filtered.csv")
var_numgenes = nrow(varonly)
var_distdata <- dist(t(varonly))		
var_dendrogram <- hclust(var_distdata,method="complete")
var_clusters <- cutree(var_dendrogram, k=2)
var_cluster1 <- separate(var_clusters,varonly,var_numgenes,1)
var_cluster2 <- separate(var_clusters,varonly,var_numgenes,2)
var_tstatistics <- unlist(lapply(1:var_numgenes, function(x) t.test(var_cluster1[x,],var_cluster2[x,])$statistic))
var_pvalues <- unlist(lapply(1:var_numgenes, function(x) t.test(var_cluster1[x,],var_cluster2[x,])$p.value))
var_diffexpress <- data.frame(ID = rownames(varonly), t_statistic = var_tstatistics, p_value = var_pvalues, adjusted_p_value = p.adjust(var_pvalues,"fdr"))
write.csv(var_diffexpress,"var_ttest.csv")
