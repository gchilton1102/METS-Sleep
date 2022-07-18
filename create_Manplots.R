## Creating Manhattan plots of the PrediXcan results for every tissue in each population (META and AFR) and each phenotype (chronotype and duration).

#load libraries
library(data.table)
library(dplyr)
library(qqman)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
# reads in the gene_anno file and removes chrM and chrY cause they are not numeric 
gene_anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ') 
gene_anno <- gene_anno[gene_anno$chr != "M"]
gene_anno <- gene_anno[gene_anno$chr != "Y"]
gene_anno$chr <- as.numeric(gene_anno$chr)
#creates column (chr_pos) with chromosome number:start position
gene_anno <- gene_anno %>% mutate(chr_pos = ""%&% gene_anno$chr %&%":"%&% 
                                    gene_anno$start)

Total <- data.frame() #create starting data frame
#list of brain tissues
tissues <- list("Hypothalamus", "Hippocampus", "Amygdala", 
                "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
                "Cerebellar_Hemisphere", "Cerebellum", "Cortex", 
                "Frontal_Cortex_BA9", "Hippocampus", "Hypothalamus", 
                "Nucleus_accumbens_basal_ganglia", "Putamen_basal_ganglia", 
                "Spinal_cord_cervical_c-1", "Substantia_nigra")

#loop creating new data frame 
plots <- list('META','AFR')
phenos <- list('duration','chronotype')
for (plot in plots) {
  Total <- data.frame()
  for (pheno in phenos) {
    Total <- data.frame()
    for (x in tissues) {
      path <- paste("/home/grace/mets-gwas/Brain_Tissues_PrediXcan/PUKBB_", 
                    plot,"_",pheno,"_Brain_",x,".csv", sep="") #finds file with specified tissue
      #opens file with the PrediXcan data from the specified tissue
      Tissue <- fread(path, header= T, sep=",") 
      Tissue <- mutate(Tissue, gene=strtrim(gene, 15)) %>% #gets rid of the decimal point
        mutate(Tissue, tissue= x) #adds tissue column
      Total <- rbind(Total,Tissue) #adds this data frame to the Total data frame
      Total <- arrange(Total,pvalue) #sorts by pval
    }
    newPath <- paste(pheno,'_',plot,'_Total.csv', sep='')
    fwrite(Total, newPath, append= FALSE, quote = 'auto', sep=',')
    
    gwas_file <- Total %>% 
      select(-c(pred_perf_pval,pred_perf_r2,pred_perf_qval)) 
    #edits the ensemble ids so there are no decimal points
    gwas_file <- mutate(gwas_file, gene = strtrim(gene, 15))
    #joins gwas with gene anno based on the ensemble ids
    gwas_file <- left_join(gwas_file,gene_anno, by=c("gene"="gene_id")) 
    #gets rid of NA data
    gwas_file <- na.omit(gwas_file)
    
    #adds column of -log10log(pvalue)
    gwas_file$logP <- as.numeric(-log10(gwas_file$pvalue))
    #adds column with transcription start position of each gene
    gwas_file <- gwas_file %>%  group_by(chr) %>% summarise(chr_len=max(start)) %>% 
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>% 
      left_join(gwas_file, ., by=c("chr"="chr")) %>% 
      arrange(chr, start) %>% mutate(BPcum=start+tot)
    #creates data frame with the middle position of each chromosome in the "center" column
    axisdf <- gwas_file %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
    
    title <- paste(pheno, plot, 'Total', sep=' ')
    pathTitle <- paste("/home/grace/mets-gwas/tissue_total_manplots/",pheno,'_',plot,'_Total.png',sep='')
    library("wesanderson")
    pal <- c(wes_palette("GrandBudapest1", type="continuous"),
             wes_palette("GrandBudapest2", type="continuous"), 
             wes_palette("Moonrise1", type="continuous"),  
             wes_palette("Moonrise2", type="continuous"))
    ggplot(gwas_file,aes(x = BPcum, y = logP)) + 
      geom_point(alpha=0.5, size = 0.5, aes(color = tissue, shape= tissue), show.legend=TRUE) + 
      ggtitle(title) + 
      scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) + 
      scale_color_manual(values = rep(pal)) + 
      scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13)) +
      labs(x = "Chromosome",y = "-log10(p)") + 
      geom_hline(yintercept = -log10(3.82e-6), color = "grey40", linetype = "dashed") + 
      theme(axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5), legend.text= element_text(size=6))
    
    ggsave(pathTitle, plot=last_plot(), width = 10)
    
  }
}
