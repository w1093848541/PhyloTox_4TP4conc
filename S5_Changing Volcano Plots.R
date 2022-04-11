
library (ggplot2);
library(openxlsx);
library("ggthemes")
library("RColorBrewer")
setwd("C:\\Users\\10938\\Desktop\\Jacobson Lab\\2022.03\\4tp4conc")

path = "C:\\Users\\10938\\Desktop\\Jacobson Lab\\2022.03\\4tp4conc"
fileNames = dir(path, pattern = "*_qvalues.xlsx")
NumHm = length(fileNames)
NumHm = as.numeric(NumHm)
i = 1

for (i in 1:NumHm){

T1data = read.xlsx(fileNames[i],sheet = 1,colNames= FALSE,rowNames = FALSE);



data = T1data[,3:5];

colnames(data) = c("q","log2FoldChange","label")



n = nrow(data)
for (q in 1:n) {if (data[q,1] > 0.05) {data[q,3] = ""}}


q1 <- ggplot(data,aes(log2FoldChange,-log10(q))) +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  geom_vline(xintercept = c(-1.2,1.2),linetype = "dashed")+
  geom_point(aes(size = -log10(q),color = -log10(q))) +
  theme_bw() +
theme(panel.grid = element_blank(),
      legend.position = c(0.99,0.99),
      legend.justification = c(0.99,0.99))  +
  scale_colour_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  guides(col = guide_colourbar(title = "-log10 q-value"),
         size = "none")+
         #size = guide_legend(title = "-log10 q-value"))
 xlab("log2FoldChange(Exposed/Control)")+
   ylab ("-log10(FDR q-value)")+
  geom_text(aes(label = label,color = -log10(q)),size = 3, vjust = 1.5, hjust = 1)

vpname = fileNames[i];
vpname = substr(vpname,1,nchar(vpname)-13);
vpname = paste(vpname,"_volcano.jpeg",collapse = "",sep = "")


ggsave(q1,filename = vpname,width = 16, height = 12)

}


