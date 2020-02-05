---
title: "programmer_jia"
output:
  pdf_document: default
  html_document: default
---

##install and library packages

```{r}
install.packages("BiocManager")
install.packages("BiocManager")
library("BiocManager", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
BiocManager::install('affy')
BiocManager::install('affyPLM')
BiocManager::install('sva')
BiocManager::install('AnnotationDbi')
BiocManager::install('hgu133plus2.db')
```

```{r}
library("affy", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
library("affyPLM", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
library("sva", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
library("AnnotationDbi", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
library("hgu133plus2.db", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
```


##read array and normalize 

```{r}
dir_cels='/projectnb/bf528/users/group4/project1/samples'
affy_data = ReadAffy(celfile.path=dir_cels)
View(affy_data)
norm_arry <- rma(affy_data)
```

##get information from norm_arry(this is just a general look)

```{r}
getClass(affy_data)
```

##compute Relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) scores of the microarray samples
```{r}
my_plm<- fitPLM(affy_data,normalize = T, background = T)

#boxplot_rle=Mbox(my_plm,main="RLE") #this is a boxplot cause i want to take a look at distribution

sta_rle=RLE(my_plm,type="stats")
#transfer from atom vector to matrix
dim(sta_rle)
dim(sta_rle) <- c(2,134)
my_med=sta_rle[1,]
#plot histgram
hist_rle=hist(my_med)
```


```{r}
sta_nuse=NUSE(my_plm,type="stats")
dim(sta_nuse)
dim(sta_nuse) <- c(2,134)
my_med_nuse=sta_nuse[1,]
hist_nuse=hist(my_med_nuse)
```

##correct for batch effects
```{r}
library(readr)
proj_metadata=read_csv("/project/bf528/project_1/doc/proj_metadata.csv")
norm_matrix=exprs(norm_arry)
head.matrix(norm_matrix)
write.csv(norm_matrix,'norm.csv')
norm <- read_csv("norm.csv")
head(norm)
```

```{r}
model <- model.matrix(~normalizationcombatmod, data = proj_metadata)
combat_norm_matrix <- ComBat(dat = norm_matrix, batch = proj_metadata$normalizationcombatbatch, mod = model)
write.csv(combat_norm_matrix,'combat_norm.csv')
library(readr)
combat_norm <- read_csv("combat_norm.csv")
transed_data=t(combat_norm_matrix)
transed_scale=scale(transed_data, center=T,scale=T)
scaled_data=t(transed_scale)
prcomp_data=prcomp(scaled_data,scale=F,center = F)
prcomp_data$rotation
pca_plot=plot(prcomp_data$rotation[,1:2])
```
