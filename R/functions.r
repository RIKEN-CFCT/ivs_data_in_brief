#--------------------------------------------
# qcBoxplot
#--------------------------------------------
#' Draw a boxplot from a lumiR.batch object.
#' 
#' Draw a boxplot of all samples of a lumiR.batch object.
#' 
#' @param data a lumiR.batch object.
#' @param  title title of boxplot.
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr tidyr %>%
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_y_continuous labs scale_color_manual theme element_blank element_rect element_text
#' @importFrom graphics plot
#' 
#' @return a ggplot object
#' 
#' @keywords boxplot ggplot
#' @export 

qcBoxplot <- function(data, title){
    exp_matrix <- exprs(data)
    df <- data.frame(exp_matrix)
    colnames(df) <- meta$Sample.Name

    mycols1 <- colorRampPalette(brewer.pal(12, "Paired"))(ncol(df)) # color for motif enrichment

    df2 <- df %>% 
        tidyr::pivot_longer(
            cols = "7.5dpp_rep1":"14.5cdpp_rep3",
            names_to = "sample_names",
            values_to = "value")

    scaleFUN <- function(x) sprintf("%.1f", x)

    ## Boxplot of log 10 intensity
    g <- ggplot(df2, aes(x = sample_names, y=log(value,10)))
    g <- g + geom_boxplot(aes(color = sample_names), show.legend = TRUE)
    g <- g + scale_y_continuous(labels=scaleFUN)
    g <- g + labs(x = "Sample", y = "Intensity (log 10)")
    g <- g + scale_color_manual(values = mycols1, limits=c(
        "7.5dpp_rep1", "7.5dpp_rep2", "7.5dpp_rep3",
        "9.5dpp_rep1", "9.5dpp_rep2", "9.5dpp_rep3",
        "9.5cdpp_rep1", "9.5cdpp_rep2", "9.5cdpp_rep3",
        "11.5dpp_rep1", "11.5dpp_rep2", "11.5dpp_rep3",
        "11.5cdpp_rep1", "11.5cdpp_rep2", "11.5cdpp_rep3",
        "13.5dpp_rep1", "13.5dpp_rep2", "13.5dpp_rep3",
        "13.5cdpp_rep1", "13.5cdpp_rep2", "13.5cdpp_rep3",
        "14.5dpp_rep1", "14.5dpp_rep2", "14.5dpp_rep3",
        "14.5cdpp_rep1", "14.5cdpp_rep2", "14.5cdpp_rep3",
        "16.5dpp_rep1", "16.5dpp_rep2", "16.5dpp_rep3",
        "16.5cdpp_rep1", "16.5cdpp_rep2", "16.5cdpp_rep3",
        "21.5dpp_rep1", "21.5dpp_rep2", "21.5dpp_rep3",
        "21.5cdpp_rep1", "21.5cdpp_rep2", "21.5cdpp_rep3"
        ))
    g <- g + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "gray1", fill = NA),
            axis.text = element_text(size = 18),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 18),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.background = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
            )
    plot(g)
    return(g)
}

#--------------------------------------------
# lumiDEG
#--------------------------------------------
#' Identify differentially expressed genes.
#' 
#' From the output expression matrix, this script identifies differentially expressed genes based on  p-values of empirical Bayes moderated t-statistics test
#' 
#' @param data an MArrayLM object or a list object produced by the function lm.series or equivalent. Must contain components coefficients and stdev.unscaled.
#' @param control sample name of control sample of the comparisn
#' @param treatment name of treatment sample of the comparisn
#' @param cutoff cutoff fdr for the differentially expression analysis. defoult is 0.05.
#' 
#' @importFrom limma contrasts.fit eBayes topTable
#' @importFrom annotte getSYMBOL getEG 
#' 
#' @return a data.frame of differentially expression analysis
#' 
#' @keywords seurat 10x
#' @export 

lumiDEG <- function(fit, control, treatment, cutoff = 0.05){
    #empirical Bayes moderated t-statics test
    cont.matrix <- eval(parse(text = paste0("makeContrasts(", treatment, "-" , control, ",levels=design)")))
    fit2<- contrasts.fit(fit, cont.matrix)
    fit3 <- eBayes(fit2)

    #annltation
    geneSymbol <- getSYMBOL(probeList, 'lumiMouseAll.db')
    entrezID <- getEG (probeList, 'lumiMouseAll.db')
    geneName <- sapply(lookUp(probeList, 'lumiMouseAll.db', 'GENENAME'), function(x) x[1])
    fit3$genes <- data.frame(ID=probeList, geneSymbol=geneSymbol, entrezID=entrezID, geneName=geneName, stringsAsFactors=FALSE)

    #Extract DEGs
    totalnum <- nrow(fit3$genes) 
    tt <- topTable(fit3, coef=1, adjust='fdr', p.value = cutoff, sort.by="P", number=totalnum)
    return(tt)
}

