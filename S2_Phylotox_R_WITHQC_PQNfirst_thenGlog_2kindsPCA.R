# install BiocManager if not present
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# install structToolbox and dependencies
#BiocManager::install("structToolbox")
#BiocManager::install("mixOmics")
#a
#BiocManager::install(c('pmp', 'ropls', 'BiocFileCache'))
#install.packages(c('cowplot', 'openxlsx'))
#a
suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  library(BiocFileCache)
  
  # CRAN libraries
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  library(openxlsx)
  library(ggthemes)
  library(RColorBrewer)
  library(ropls)
  library(dplyr)
  library(mixOmics)
})

setwd("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2022.03\\4tp4conc")
MATLABresult = read.xlsx("Combined_withQC_4tp_for4tp4conc_all.xlsx",sheet = 1,colNames= FALSE,rowNames = FALSE);
MassList = MATLABresult$X1;
MassList = MassList[-1];
Mass_list = as.character(MassList);
PhylotoxData = MTBLS79
NumSamples = as.numeric(dim(MATLABresult)[2])-1
NumMetabolites = as.numeric(dim(MATLABresult)[1])-1
nSamples = 1
RawMatrix = as.matrix(MATLABresult)
Type_names = vector();
while (nSamples <= NumSamples) { Type_names[nSamples] = RawMatrix[1,nSamples+1];
nSamples = nSamples +1;

}
noQC = read.xlsx("Combined_withQC_4tp_for4tp4conc_all.xlsx",sheet = 1,colNames= TRUE,rowNames = TRUE);
noQC = as.matrix(noQC)
colnames(noQC) = Type_names
rownames(noQC) = Mass_list
Type = Type_names
OPLS_DA_GROUPS = factor(Type_names)
#Type = vector()
#for (s in 1:NumSamples) { if(Type_names[s] == 'E3'|Type_names[s] =='E2'|Type_names[s] =='E1') {Type[s] = 'E'}}
#for (s in 1:NumSamples) { if(Type_names[s] == 'C3'|Type_names[s] =='C2'|Type_names[s] =='C1') {Type[s] = 'C'}}
#for (s in 1:NumSamples) { if(Type_names[s] == 'C0') {Type[s] = 'C0'}}%%
Batch = vector()
for (s in 1:NumSamples) {Batch[s] = '1'}

Replicates = vector()
NumQC = 0
for (s in 1:NumSamples) {if (Type_names[s]=='QCs') {NumQC = NumQC +1}}
NumE33 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E33') {NumE33 = NumE33 +1}}
NumE32 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E32') {NumE32 = NumE32 +1}}
NumE31 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E31') {NumE31 = NumE31 +1}}
NumE23 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E23') {NumE23 = NumE23 +1}}
NumE22 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E22') {NumE22 = NumE22 +1}}
NumE21 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E21') {NumE21 = NumE21 +1}}
NumE13 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E13') {NumE13 = NumE13 +1}}
NumE12 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E12') {NumE12 = NumE12 +1}}
NumE11 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='E11') {NumE11 = NumE11 +1}}
NumCT3 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='CT3') {NumCT3 = NumCT3 +1}}
NumCT2 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='CT2') {NumCT2 = NumCT2 +1}}
NumCT1 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='CT1') {NumCT1 = NumCT1 +1}}
NumCT0 = 0
for (s in 1:NumSamples) {if (Type_names[s]=='CT0') {NumCT0 = NumCT0 +1}}


for (s in 1:NumQC) {Replicates[s] = as.character(s)}
for (s in 1:NumE33) {Replicates[s+NumQC] = as.character(s)}
for (s in 1:NumE32) {Replicates[s+NumE33+NumQC] = as.character(s)}
for (s in 1:NumE31) {Replicates[s+NumE33+NumE32+NumQC] = as.character(s)}
for (s in 1:NumE23) {Replicates[s+NumE33+NumE32+NumE31+NumQC] = as.character(s)}
for (s in 1:NumE22) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumQC] = as.character(s)}
for (s in 1:NumE21) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumQC] = as.character(s)}
for (s in 1:NumE13) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumQC] = as.character(s)}
for (s in 1:NumE12) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumQC] = as.character(s)}
for (s in 1:NumE11) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumE12+NumQC] = as.character(s)}
for (s in 1:NumCT3) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumE12+NumE11+NumQC] = as.character(s)}
for (s in 1:NumCT2) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumE12+NumE11+NumCT3+NumQC] = as.character(s)}
for (s in 1:NumCT1) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumE12+NumE11+NumCT3+NumCT2+NumQC] = as.character(s)}
for (s in 1:NumCT0) {Replicates[s+NumE33+NumE32+NumE31+NumE23+NumE22+NumE21+NumE13+NumE12+NumE11+NumCT3+NumCT2+NumCT1+NumQC] = as.character(s)}
PhylotoxData@colData@listData[["Sample_Rep"]] = Replicates
PhylotoxData@colData@listData[["Batch"]] = Batch
PhylotoxData@colData@listData[["Class"]] = Type
PhylotoxData@colData@nrows = as.integer(NumSamples)
PhylotoxData@elementMetadata@nrows = as.integer(NumMetabolites)
PhylotoxData@colData@listData[["Class2"]] = Type_names
PhylotoxData@assays@data@listData[[1]] = noQC
PhylotoxData@NAMES= Mass_list
PhylotoxData@colData@rownames = Type_names
DE = as.DatasetExperiment(PhylotoxData)
DE$sample_meta$run_order = 1:nrow(DE)
Type=as.character(DE$sample_meta$Class)
Type[Type != 'QCs'] = 'Sample'
DE$sample_meta$Type = factor(Type)
DE$sample_meta$Batch = factor(DE$sample_meta$Batch)
DE$sample_meta$Class = factor(DE$sample_meta$Class)





M5 = pqn_norm(qc_label='QCs',factor_name='Type') 
  
M5 = model_apply(M5,DE)

PQN_result = M5@normalised@value@assays@data@listData[[1]]
PQN_result = t(PQN_result)
colnames(PQN_result) = Type_names
write.csv(PQN_result,file = "PQN_predictedresult.csv")

M6 = knn_impute(neighbours=5,by='samples') +
  glog_transform(qc_label='QCs',factor_name='Type')

M6 = model_apply(M6,predicted(M5))

GlogedData = predicted(M6)
glog_result = GlogedData@assays@data@listData[[1]]
glog_result = t(glog_result)
colnames(glog_result) = Type_names
write.csv(glog_result,file = "glog_predictedresult.csv")

# PCA
M7  = mean_centre() + PCA(number_components = 2)

# apply model sequence to data
M7 = model_apply(M7,predicted(M6))

# plot pca scores
C = pca_scores_plot(factor_name=c('Sample_Rep','Class'),ellipse='none')
chart_plot(C,M7[2]) + coord_fixed() +guides(colour=FALSE) + scale_shape_manual(values = 10:24)#change to the corresponding number

PhylotoxPCA = M7@models[[2]]@scores@value@assays@data@listData[[1]]
write.csv(PhylotoxPCA,file = "WithQCs_pca.csv")

PhylotoxData2 = MTBLS79

nRemove = 1
while(nRemove <= NumQC){Replicates = Replicates[-1];
Batch = Batch[-1];
Type = Type[-1];
Type_names = Type_names[-1];
M6@models[[2]]@transformed@value@assays@data@listData[[1]]= M6@models[[2]]@transformed@value@assays@data@listData[[1]][-1,];
nRemove = nRemove + 1;
}

PhylotoxData2@colData@listData[["Sample_Rep"]] = Replicates
PhylotoxData2@colData@listData[["Batch"]] = Batch
PhylotoxData2@colData@listData[["Class"]] = Type
PhylotoxData2@colData@nrows = as.integer(NumSamples-NumQC)
PhylotoxData2@elementMetadata@nrows = as.integer(NumMetabolites)
PhylotoxData2@colData@listData[["Class2"]] = Type_names
PhylotoxData2@assays@data@listData[[1]] = t(M6@models[[2]]@transformed@value@assays@data@listData[[1]])
PhylotoxData2@NAMES= Mass_list
PhylotoxData2@colData@rownames = Type_names
DE = as.DatasetExperiment(PhylotoxData2)
DE$sample_meta$run_order = 1:nrow(DE)
Type=as.character(DE$sample_meta$Class)
Type[Type != 'QCs'] = 'Sample'
DE$sample_meta$Type = factor(Type)
DE$sample_meta$Batch = factor(DE$sample_meta$Batch)
DE$sample_meta$Class = factor(DE$sample_meta$Class)



# PCA
M7  = mean_centre() + PCA(number_components = 2)

# apply model sequence to data
M7 = model_apply(M7,DE)

# plot pca scores
C = pca_scores_plot(factor_name=c('Sample_Rep','Class2'),ellipse='none')
chart_plot(C,M7[2]) + coord_fixed() +guides(colour=FALSE) + scale_shape_manual(values = 11:24)#change to the corresponding number

PhylotoxPCA = M7@models[[2]]@scores@value@assays@data@listData[[1]]
write.csv(PhylotoxPCA,file = "WithQCsbutremoved_pca.csv")

gloged = read.csv("glog_predictedresult.csv",header = TRUE,row.names = 1)

gloged = t(gloged)

plsda.datatm <- plsda(gloged,OPLS_DA_GROUPS,ncomp = 3)
jpeg(filename = "plada_result.jpeg",width = 5, height = 5, units = "in",res =300 )

plotIndiv(plsda.datatm,ind.names = FALSE,legend = TRUE, ellipse = FALSE,title = "PLS-DA - Result")
dev.off()


background <- background.predict(plsda.datatm, comp.predicted = 2, dist = "max.dist")
jpeg(filename = "plsda_max.jpeg",width = 5, height = 5, units = "in",res =300)

 plotIndiv(plsda.datatm, comp = 1:2, ind.names = FALSE, title = "Maximum Distance",
          legend = TRUE, background = background, ellipse = FALSE)
 dev.off()
#opls(gloged,OPLS_DA_GROUPS,predI = 1)