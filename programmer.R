#install and library packages
> install.packages("BiocManager")
> library("BiocManager", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
> BiocManager::install('affy')
> BiocManager::install('affyPLM')
> BiocManager::install('sva')
> BiocManager::install('AnnotationDbi')
> BiocManager::install('hgu133plus2.db')
> library("affy", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
> library("affyPLM", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
> library("sva", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
> library("AnnotationDbi", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")
> library("hgu133plus2.db", lib.loc="/share/pkg.7/r/3.6.0/install/lib64/R/library")

#read array and normalize 
> dir_cels='/projectnb/bf528/users/group4/project1/samples'
> affy_data = ReadAffy(celfile.path=dir_cels)
> View(affy_data)
> norm_arry <- rma(affy_data)

#read information of norm_arry
> getClass(affy_data)
AffyBatch object
size of arrays=1164x1164 features (78 kb)
cdf=HG-U133_Plus_2 (54675 affyids)
number of samples=134
number of genes=54675
annotation=hgu133plus2


> my_plm<- fitPLM(affy_data,normalize = T, background = T)


> boxplot_rle=Mbox(my_plm,main="RLE")

> sta_rle=RLE(my_plm,type="stats")



#transfer from atom vector to matrix
> dim(sta_rle)
[1]   2 134
> dim(sta_rle) <- c(2,134)
> my_med=sta_rle[1,]
#plot histgram
> hist_rle=hist(my_med)

> hist_rle=hist(my_med)
> sta_nuse=NUSE(my_plm,type="stats")
> dim(sta_nuse)
[1]   2 134
> dim(sta_nuse) <- c(2,134)
> my_med_nuse=sta_nuse[1,]
> hist_nuse=hist(my_med_nuse)

#correct for batch effects
> library(readr)> proj_metadata 
<- read_csv("/project/bf528/project_1/doc/proj_metadata.csv")

#get error in combat: no method for coercing this S4 class to a vector
#so I transfer s4 to data frame.
> new_norm=as.data.frame(norm_arry)
########
> norm_matrix=exprs(norm_arry)
> head.matrix(norm_matrix)
> write.csv(norm_matrix,'norm.csv')
> norm <- read_csv("norm.csv")
> View(norm)

> model <- model.matrix(~normalizationcombatmod, data = proj_metadata)
> combat_norm_matrix <- ComBat(dat = norm_matrix, batch = proj_metadata$normalizationcombatbatch, mod = model)
> write.csv(combat_norm_matrix,'combat_norm.csv')
> library(readr)
> combat_norm <- read_csv("combat_norm.csv")

> transed_data=t(combat_norm_matrix)
> transed_scale=scale(transed_data, center=T,scale=T)
> scaled_data=t(transed_scale)
> prcomp_data=prcomp(scaled_data,scale=F,center = F)
> prcomp_data$rotation
> pca_plot=plot(prcomp_data$rotation[,1:2])
