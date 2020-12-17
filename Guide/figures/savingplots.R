source("./PlottingFunctions.R")

TPM_df <- read_csv("./example_data/Ccera_rsem_tpm.csv")
TPM_df <- TPM_df %>% mutate_if(is.character,as.factor)

TPM_df <- df_clean(TPM_df, class="class", toxin_family="toxin_family", colors=toxin_colors)

TPM_df2 <- t(cmultRepl(t(TPM_df[,4:12]),output = "p-counts"))
TPM_df2 <- cbind(TPM_df[,1:3], TPM_df2)
rownames(TPM_df2) <- TPM_df2$gene_id

FullBarplot(TPM_df2,"CLP2057",class="class")
FullPie(TPM_df2,"CLP2057",class="class")
ToxinBarplot(TPM_df2,"CLP2057",class="class",toxin_family="toxin_family", colors=toxin_colors)
ToxinPie(TPM_df2,"CLP2057",class="class",toxin_family="toxin_family", colors=toxin_colors)
FancyFigure(TPM_df2,"CLP2057",class="class",toxin_family="toxin_family", colors=toxin_colors)

png("figures/TransCompPlot.png", width =480, height=480)
TransCompPlot(TPM_df2, "CLP2057", "CLP2065")
invisible(dev.off())

samples <- read_csv("./example_data/Ccera_samples.csv")
samples <- as.data.frame(samples) %>% mutate_if(is.character,as.factor)
rownames(samples) <- samples$ID
# Log-transforming data
lnToxins <- TPM_df2 %>% filter(class == "Toxin") %>% 
  mutate_if(is.numeric,log) %>% 
  arrange(desc(Average))
rownames(lnToxins) <- lnToxins$gene_id
png("figures/heatmap_toxins.png", width=480, height=720)
pheatmap(lnToxins[,4:11], cluster_rows=F, show_rownames=F,cluster_cols=T, 
         annotation_col=samples[,2:5], annotation_legend=T)
dev.off()

TPM_family_df <- TPM_df2 %>% 
  filter(class=="Toxin") %>%
  group_by(toxin_family,class) %>% 
  summarize_if(is.numeric,sum) %>%
  arrange(desc(Average)) %>% ungroup()
TPM_family_df <- as.data.frame(TPM_family_df)
rownames(TPM_family_df) <- TPM_family_df$toxin_family
# Log transform
lnClasses <- as.data.frame(TPM_family_df) %>%
  mutate_if(is.numeric,log) %>%
  arrange(desc(Average))
rownames(lnClasses) <- lnClasses$toxin_family
png("figures/heatmap_classes.png", width=480, height=720)
pheatmap(lnClasses[,3:10], cluster_rows=F, show_rownames=T,cluster_cols=T, 
         annotation_col=samples[,2:5], annotation_legend=T)
dev.off()

library(ape)
Tree<-read.tree(file="./example_data/Ccera-Phylogeny.nwk")
library(phytools)
png("figures/phyloheatmap_toxins.png", width=720, height=480)
phylo.heatmap(Tree, t(lnToxins[,4:11]), fsize=c(1,0.8,1), standardize=F, labels=F,
              split=c(0.3,0.7), ylim=c(-0.25,1.25), grid=T,
              colors=colorRampPalette(rev(brewer.pal(n = 7,name="RdYlBu")))(100))
dev.off()


PCA <- prcomp(as.data.frame(t(TPM_df2[,4:11])), center=TRUE, scale=TRUE)
PCA_df<-data.frame(samples,PCA$x)
contribs<-round((PCA$sdev^2/sum(PCA$sdev^2))*100,2)
png("figures/pca_toxins.png", width=480, height=480)
ggscatter(PCA_df,"PC1","PC2",color="NontoxPhylo",fill="NontoxPhylo",size=8,
          ellipse = T, ellipse.level = 0.95, ellipse.type = "norm",
          xlab = paste0("PC1 (",contribs[1],"%)"), ylab=paste0("PC2 (",contribs[2],"%)"))
dev.off()
