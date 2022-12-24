#### GO analysis and visualization

# Load all the required libraries
library(tidyverse)
library(clusterProfilter)
library(org.Mm.eg.db)
library(ggnewscale)
library(ggupset)
library(DOSE)
library(enrichplot)
library(ggridges)
library(pathview)
library(writexl) 
library(tidyverse)
library(ggplot2)

# Set working directory 
setwd("/home/jie/Documents/jie21/R BIo")

# Load data 
data1 <-read_csv("geneList_Jie.csv")  # First Column is gene symbol, Second column is log2FC.  
head(data1)

# Clean data, convert gene symbol to ENTREZID
Symbol <- data1[,1]

gene.df <- bitr(Symbol, fromType = "SYMBOL",  # Convert gene symbol to ENTREZID
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db) 
head(gene.df)

data2 <- merge(data1, gene.df, by.x="GeneName", by.y="SYMBOL")
str(data2)

test3 <- data2[,-1] # Remove GeneName column
head(test3)
geneList <- test3[,2]

names(geneList) <- as.character(test3[,1]) # Set first column ENTREZID as vector and character
geneList <- sort(geneList, decreasing = TRUE)

gene <- names(geneList)[abs(geneList)>1.5]
head(gene)

length(gene)
length(geneList)

# Gene Ontology Analysis - Cellular Component 
ego_CC <- enrichGO(gene       = gene,
                universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
str(ego_CC)

# Visualize the top items in Cellular Component
barplot(ego_CC,showCategory = 30)
dotplot(ego_CC,showCategory=30)
upsetplot(ego_CC)

# Export the data for detail analysis to excel
result_CC <- ego_CC@result
str(result_CC)
write_xlsx(result_CC, "XX_result_CC.xlsx") 


# Gene Ontology Analysis - Biological Process
ego_BP <- enrichGO(gene       = gene,
                universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_BP)

# Visualize the top items in Biological Process
barplot(ego_BP,showCategory = 30)
dotplot(ego_BP,showCategory=30)
upsetplot(ego_BP)

# Export the data for detail analysis to excel
result_BP <- ego_BP@result
str(result_BP)
write_xlsx(result_BP, "XX_result_BP.xlsx") 


# Gene Ontology Analysis - Molecular Function
ego_MF <- enrichGO(gene       = gene,
                universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_MF)

# Visualize the top items in Molecular Function
barplot(ego_MF,showCategory = 30)
dotplot(ego_MF,showCategory=30)
upsetplot(ego_MF)

# Export the data for detail analysis to excel
result_MF <- ego_MF@result
str(result_MF)
write_xlsx(result_MF, "XX_result_MF.xlsx") 


# With this excel file exported from above (CC, BP, and MF), will then check in details on individual GO terms, and create a new spreadsheet called "GO_summary" for selected the ones for visualization 
GO_summary <- read_csv("GO_summary.csv")
GO_summary <- GO_summary[GO_summary$pvalue <0.05, ]


# Avoid order change, set column "Description" as factor 
GO_summary$Description <-factor(GO_summary$Description, levels = GO_summary$Description)

# Use ggplot2 to visualize the data
# Dual x axis with both gene count number and p value  
offval <- 5
gg <- ggplot(data = GO_summary, aes(y=Description, x= `-logP_value`)) + 
    geom_bar( stat = "identity",fill= "steelblue", width = 0.4) + 
    geom_line(aes(x=Count * 0.2), group=1, color = "darkred") +
    scale_x_continuous(name="-Log P value", sec.axis = sec_axis(~.*offval, name = "Count")) + 
    scale_y_discrete(limits = rev(levels(summary$Description))) + 
    theme(aspect.ratio = 1/2,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
    ylab("GO Terms")
    
    