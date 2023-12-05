library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(writexl)
library(readxl)
library(stats)

setwd("<enter the path for the working directory")

################################################################################
#################################    GSEA     ##################################
################################################################################

########################### Plot after GSEA #####################################
par(las=2) # make label text perpendicular to axis
par(mar=c(5,20,8,0.5)) # increase y-axis margin. 12 20

####Define Theme
##Plot GSEA results
#text size
t=3.5
m=15
l=12
z=2
mytheme<- theme(legend.position = "bottom",
                legend.text = element_text(size=l),
                text = element_text(family="Calibri", colour = "black"),
                axis.text.y=element_blank(),
                axis.text.x=element_text(size=10, color="black"),
                axis.title.x = element_text(size=12, color="black"),
                axis.ticks=element_blank(),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                plot.title = element_text(face="bold", size=m, hjust = 0.5, family="Calibri", margin = margin(b=z)),
                plot.subtitle = element_text(face="bold", size=m, hjust = 0.5, family="Calibri", margin =margin(t=z, b=z+2)),
                # plot.margin = unit(c(0.2, 0,0.9, 0), units = "cm"))
                plot.margin = unit(c(0.2, 0,0.5, 0), units = "cm"))

adapt <- function (x) {gsub("interferon", "IFN", gsub("reactive oxigen species pathway", "ROS pathway", gsub("epithelial mesenchymal transition", "EMT", gsub("_"," ",gsub("reactome_", "", tolower(x))))))}

# DESeq2 FDR
# Find which pathways are enriched in the model we built 

UP <- fread("<enter the path for gsea_report_for_na_pos.tsv and file name>")
DW <- fread("<enter the path for gsea_report_for_na_neg.tsv and file name>")
UP$class="UP"
DW$class="DW"

UP = UP %>% 
  dplyr::select(NAME, SIZE, NES, ES, `FDR q-val`, `FWER p-val`,class)
UP <- UP[1:30,]

DW = DW %>%
  dplyr::select(NAME, SIZE, NES, ES, `FDR q-val`, `FWER p-val`,class) 
DW <- DW[1:5,]

gsea_report = rbindlist(list(UP,DW))

setorder(gsea_report, NES)
gsea_report$NAME = factor(gsea_report$NAME, gsea_report$NAME)
val<-ceiling(max(abs(gsea_report$NES))/2)
gsea_report$hjust <- ifelse(gsea_report$NES > 0, -val, val)

#Alpha trasparency based on FDR
gsea_report$alpha = ifelse(gsea_report$`FDR q-val` < 0.05, 1, 0.5)
gsea_report$class<-factor(gsea_report$class, levels = c("UP", "DW"))
palette = c("firebrick1","steelblue3")

#Barplot
plot <- ggplot(gsea_report, aes(NAME, NES,
                                label = adapt(NAME))) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = class), show.legend = T) +
  ylim(c(-2*val,2*val)) +
  scale_fill_manual(values = palette, labels=c("Pathway overexpression", "Pathway downregulation"), name=NULL) + 
  theme_bw() +
  geom_text(aes(family="Calibri", y=hjust), size=t) +
  coord_flip() +
  xlab("<enter x ax lable>") + ylab("Normalized Enrichment Score") +
  mytheme +
  geom_hline(yintercept = 0) +
  ggtitle(expression(bold(paste(bolditalic("<enter gene name>"), "<enter other parts of the title>"))))

ggsave(plot = plot, filename = "<enter path to save the plot and plot name>.png", device = "png", width = 40, height = 40, units = "cm")
