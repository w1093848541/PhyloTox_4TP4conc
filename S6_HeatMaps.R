#BiocManager::install("ComplexHeatmap")
#BiocManager::install("circlize")
#BiocManager::install("grid")
suppressPackageStartupMessages({library (ggplot2);
library(openxlsx);
library(ggthemes)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
  library(grid)
})
setwd("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2022.03\\4tp4conc")
path = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2022.03\\4tp4conc"
fileNames = dir(path, pattern = "*sigpeaks_zscored.xlsx")
NumHm = length(fileNames)
NumHm = as.numeric(NumHm)
i = 1

for (i in 1:NumHm) {

data <- read.xlsx(fileNames[i],sheet = 1,colNames= FALSE,rowNames = FALSE);
data <- as.matrix(data)
names <- c("Quality Controls","EC50 Time Point 3","EC50 Time Point 2","EC50 Time Point 1","EC25 Time Point 3","EC25 Time Point 2","EC25 Time Point 1","EC10 Time Point 3","EC10 Time Point 2","EC10 Time Point 1","Control Time Point 3","Control Time Point 2","Control Time Point 1","Control Time Point 0");
colnames(data)= names;
idfile = fileNames[i];
idfile = substr(idfile,1,nchar(idfile)-12);
idfile = paste(idfile,"afterGLOG_IDs.xlsx",collapse = "",sep = "");
ids = read.xlsx(idfile,sheet = 1,colNames= FALSE,rowNames = FALSE);
ids = t(ids);
rownames(data)=ids;
col_fun <- colorRamp2(c(-2,0,2),c("#5296cc","#fdedf6","#f064af"))

hmname = fileNames[i];
hmname = substr(idfile,1,nchar(idfile)-20);
hmname = paste(hmname,"_heatmap.jpeg",collapse = "",sep = "")
jpeg(filename = hmname,width = 6, height = 6, units = "in",res =300)

heatmap1 = Heatmap(data,
        #color
        col = col_fun,
        #height of the dendrogram
        rect_gp = gpar(col = "white", lwd = 1),
        
        row_dend_width = unit(2,"cm"),
        #size of labels
        row_names_gp = gpar(fontsize = 7 , fontface = "italic"),
        column_names_gp = gpar(fontsize= 7 ),
        #note in the cells
        #cell_fun = function(j,i,x,y,width,height,fill) {
          #grid.text(sprintf("%.1f",data[i,j]),x,y,gp = gpar(fontsize = 5))},
        )
print(heatmap1)
dev.off()



}

