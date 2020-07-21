library("R.utils")
library("lumi")
library("limma")
library("lumiMouseAll.db")
library("annotate")
library("tidyr")
library("ggplot2")
library("cowplot")
library("grDevices")
library("RColorBrewer")

load("function.r")

#############################
####    Preprocessing    ####
#############################
#if the data file is compressed, you have to uncompress the file
gunzip("../data/Sample_Probe_Profile.txt.gz")

#Sample_Probe level illumina data file
filename <- "../data/Sample_Probe_Profile.txt"

#sample information file
sampleInfo <- "../data/Sample_Table.txt"

#sample names
sample_names <- "../data/Sample_metadata.txt"

# Reading the raw data
x.lumi <- lumiR.batch(filename, sampleInfoFile=sampleInfo)

#Change sample ID to sample name ([age][dpp/cddp]_[replocate index])
meta <- read.table(sample_names, sep="\t", header=T, stringsAsFactors = FALSE)
colnames(x.lumi@QC$sampleSummary) <- meta$Sample.Name
x.lumi@phenoData@data$sampleID <- meta$Sample.Name
x.lumi@phenoData@data$Sentrix.Barcode <- meta$Sample.Name

#background correction, variance-stabilizing trans-formation (VST), and quantile normalization
lumi.N.Q <- lumiExpresso(x.lumi, QC.evaluation=TRUE)

#####################################################
####    Vidualization of Preprocessing Result    ####
#####################################################
#boxplot
p1 <- qcBoxplot(x.lumi)
p2 <- qcBoxplot(lumi.N.Q)

#extract only legend
legend <- get_legend(p1 + 
    guides(color = guide_legend(nrow = 6))
)

#add title and remove legend
p1n <- p1 + ggtitle("Raw data") + theme(legend.position="none")
p2n <- p2 + ggtitle("Preprocessed data") + theme(legend.position="none")

#Integration of plots
p1n2nl <- plot_grid(p1n, p2n, legend, ncol = 1, scale=0.9, rel_heights=c(2, 2, 0.5)) 
ggsave(p1n2nl, file="qc_boxplot.pdf", width=14, height=14)


##################################################
####    Differentially Expression Analysis    ####
##################################################
#Output Normalized data matrix
dataMatrix <- exprs(lumi.N.Q)

#Remove low quality probes
presentCount <- detectionCall(x.lumi, Th = 0.01)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)

#linear model fitting
label <- lumi.N.Q$Sample.Group 
sampleType <- label
design <- model.matrix(~ 0+factor(sampleType))
colnames(design) <- c("cdpp11", "dpp11", "cdpp13", "dpp13", "cdpp14", "dpp14", "cdpp16", "dpp16","cdpp21","dpp21","dpp7","cdpp9","dpp9")
fit <- lmFit(selDataMatrix, design)

#Differentiall expression analysis
degs <- lumiDEG(fit = fit, control = "dpp9", treatment = "cdpp9")
write.table(degs, file='9.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

degs <- lumiDEG(fit = fit, control = "dpp11", treatment = "cdpp11")
write.table(degs, file='11.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

degs <- lumiDEG(fit = fit, control = "dpp13", treatment = "cdpp13")
write.table(degs, file='13.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

degs <- lumiDEG(fit = fit, control = "dpp14", treatment = "cdpp14")
write.table(degs, file='14.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

degs <- lumiDEG(fit = fit, control = "dpp16", treatment = "cdpp16")
write.table(degs, file='16.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

degs <- lumiDEG(fit = fit, control = "dpp21", treatment = "cdpp21")
write.table(degs, file='21.5dpp_DEGs.txt',quote=F, sep="\t",row.names=F)

