#!/usr/bin/env Rscript
##### INSTALL PACKAGES #####
# Package names
packages <- c("argparse")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
   install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


##### SETUP ARGPARSE #####
parser <- ArgumentParser()
parser$add_argument("-i","--input", default="example_data/Balte-SB0022CVR.genes.results", 
#parser$add_argument("-i","--input", default="example_data/Ccera_rsem_tpm.csv", 
                    help="Expression quantification results (rsem or kallisto) [default: \"%(default)s\"]")
parser$add_argument("-s","--sample", default="Balte-SB0022CVR",
#parser$add_argument("-s","--sample", default="CLP2057",                    
                    help="Sample name and prefix to use [default: \"%(default)s\"]")
parser$add_argument("--iforgot", action="store_true", default=FALSE,
                    help="Did you label the beginning of toxin sequences with \">toxin\" (case-insensitive). [default off]")
parser$add_argument("--noplot", action="store_false", default=TRUE,
                    help="Don't make plots, only format dataframe [default off]")
parser$add_argument("--noformat", action="store_false", default=TRUE,
                    help="Don't format input data, skip to making plots [default off]")
args <- parser$parse_args()


##### SOURCE PLOTTING FUNCTIONS #####
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "PlottingFunctions.R")
print(paste("Sourcing",other.name))
source(other.name)
#source("PlottingFunctions.R")


##### READ DATA #####
rsem<-normalizePath(args$input)
if(grepl("\\.csv",rsem)){
   rsem<-read_csv(rsem)
}else{
   rsem<-read_delim(rsem,"\t")  
}
rsem<-rsem %>% mutate_if(is.character,as.factor)

##### RENAME/REFORMAT COLUMNS #####
if(args$noformat){
   
   names(rsem)[grep("tpm",names(rsem),ignore.case = T)]=paste0(args$sample,"_tpm")
   names(rsem)[grep("count",names(rsem),ignore.case = T)]=paste0(args$sample,"_count")
   rsem=rsem[,-grep("fpkm",names(rsem),ignore.case = T)]
   
   gene_id <- names(rsem)[1]
   
   ##### CREATE "CLASS" COLUMN #####
   if(args$iforgot){
      toxin_pattern<-paste0("^toxin|toxin|", paste(names(toxin_colors),collapse="|"))
      rsem<-rsem %>% rowwise() %>% mutate(class=ifelse(grepl(toxin_pattern,!!as.name(gene_id),ignore.case = T),"Toxin",
                                          ifelse(grepl("\\W",!!as.name(gene_id)),"Nontoxin","Uncharacterized")))
   }else{
      rsem<-rsem %>% rowwise() %>%mutate(class=ifelse(grepl("^toxin",!!as.name(gene_id),ignore.case = T),"Toxin",
                                          ifelse(grepl("[\\W_]",!!as.name(gene_id)),"Nontoxin","Uncharacterized")))
   }
   
   
   ##### SPLIT TOXINS/OTHER #####
   Toxins<-as.data.frame(rsem %>% filter(class=="Toxin") %>% mutate(toxin_family=""))
   Others<-as.data.frame(rsem %>% filter(class!="Toxin") %>% mutate(toxin_family=class))
   
   
   ##### IDENTIFY BEST TOXIN FAMILY (split & amatch) #####
   tox_names<-names(toxin_colors)
   r0<-vector()
   for(i in 1:nrow(Toxins)){
      x<-as.character(Toxins[[gene_id]][i])
      x<-unlist(strsplit(x, '[\\W_|-]'))
      x<-x[x != ""]
      r1<-tox_names[amatch(toupper(x), toupper(tox_names))]
      r1<-r1[complete.cases(r1)]
      if(length(r1)==0){
         Toxins$toxin_family[i]<-gsub(".*_","",Toxins[[gene_id]][i])
      }else{
         Toxins$toxin_family[i]<-r1[1]
      }
   }
   
   Toxins$toxin_family<-gsub("LAO","LAAO",Toxins$toxin_family)
   Toxins$toxin_family<-gsub("SNACLEC","CTL",Toxins$toxin_family)
   
   
   ##### CHECK FOR IMPOSTER NONTOXINS POSING AS TOXINS #####
   p=0
   imposters<-vector()
   tmp<-Toxins
   while(p<0.05){
      t<-grubbs.test(log(table(tmp$toxin_family)))
      p=t$p.value
      if(p<0.05){
         n=names(t$p.value)
         if(n %in% tox_names){
            tmp<-tmp %>% filter(toxin_family!=n)
         }
         else{
            imposters<-c(imposters,n)
            tmp<-tmp %>% filter(toxin_family!=n)
         }
      }else{}
   }
   
   Toxins<-Toxins %>% mutate(class=ifelse(toxin_family %in% imposters,"Nontoxin",class),
                      toxin_family=ifelse(toxin_family %in% imposters,"Nontoxin",toxin_family))
   imposters<-Toxins %>% filter(class!="Toxin")
   Toxins<-Toxins %>% filter(class=="Toxin") %>% arrange(desc(!!as.name(paste0(args$sample,"_tpm"))))
   
   Others<-rbind(Others,imposters)
   Others<-Others %>% arrange(class,desc(!!as.name(paste0(args$sample,"_tpm"))))
   
   
   ##### ADD ANY NEW TOXINS TO TOXIN COLORS #####
   new_toxs<-sort(setdiff(unique(Toxins$toxin_family), tox_names))
   j=1
   for(i in new_toxs){
      if(j<7){
         names(toxin_colors)[which(names(toxin_colors)==paste0("zOpenColor",j))]<-i
         j=j+1
      }else{
         new_col<-sample(colors(),1)
         names(new_col)<-i
         c(toxin_colors,new_col)
         j=j+1
      }
   }
   
   rsem<-rbind(Toxins,Others)
   rsem<-rsem %>% mutate_if(is.character,as.factor)
   
   
   ##### WRITE DATA #####
   write_csv(rsem,paste0(args$sample,"_expression_format.csv"))
}

##### MAKE PLOTS #####
if(args$noplot){
   
   if(args$noformat){
      column=paste0(args$sample,"_tpm")
      
      rsem[[column]] <- rsem[[column]]+1
      
      pdf(paste0(args$sample,"_transcriptome.pdf"),width = 16, height=10)
      FancyFigure(df=rsem, id=column, class="class", toxin_family="toxin_family",colors=toxin_colors)
      dev.off()
      
      svg(paste0(args$sample,"_transcriptome.svg"),width = 12, height=7)
      FancyFigure(df=rsem, id=column, class="class", toxin_family="toxin_family",colors=toxin_colors)
      dev.off()
      
      unlink("Rplots.pdf")
   }
   else{
      column=args$sample
      
      rsem<-df_clean(rsem, class="class", toxin_family="toxin_family",colors=toxin_colors)
      
      rsem[[column]] <- rsem[[column]]+1
      
      pdf(paste0(column,"_transcriptome.pdf"),width = 16, height=10)
      FancyFigure(df=rsem, id=column, class="class", toxin_family="toxin_family",colors=toxin_colors)
      dev.off()
      
      svg(paste0(column,"_transcriptome.svg"),width = 12, height=7)
      FancyFigure(df=rsem, id=column, class="class", toxin_family="toxin_family",colors=toxin_colors)
      dev.off()
      
      unlink("Rplots.pdf")
   }
}else{unlink("Rplots.pdf")}