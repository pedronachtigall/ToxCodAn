# Venom Transcriptome Plotting Functions & Colors
# Functions by: Rhett M. Rautsaw

# # CRAN Packages
packages<-c("compositions","dplyr", "ggpubr", "outliers", "patchwork", "pheatmap", "RColorBrewer", "readr", "stringdist", "zCompositions")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
 install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Toxin Colors
toxin_colors<-c("#d30b94","3FTx",
                "#201923","BPP",
                "#772b9d","CRISP",
                "#0ec434","CTL",
                "#a4a4a4","Ficolin",
                "#ffffff","FusedToxin",
                "#29bdab","HYAL",
                "#ffdab9","KUN",#e68f66
                "#ffc413","LAAO",
                "#fcff5d","Lacta",
                "#f47a22","MYO",
                "#fffac8","NGF",#ffcba5
                "#5d4c86","NUC",
                "#7cd676","PDE",
                "#f52727","PLA2",
                "#991919","PLA2_neuro",
                "#c3a5b4","PLB",
                "#3998f5","SVMMP",
                "#277da7","SVMP",
                "#277da7","SVMPI",
                "#3750db","SVMPII",
                "#2f2aa0","SVMPIII",
                "#228c68","SVSP",
                "#aaffc3","VEGF",#946aa2
                "#aa6e28","Vespryn",#c56133
                "#cccccc","VF",
                "#f07cab","Waprin",
                "#8B5A2B","zOpenColor1",
                "#235b54","zOpenColor2",
                "#37294f","zOpenColor3",
                "#96341c","zOpenColor4",
                "#632819","zOpenColor5",
                "#b732cc","zOpenColor6",
                "#8ad8e8","zOpenColor7")

toxin_colors_df<-matrix(toxin_colors,nrow=length(toxin_colors)/2,ncol=2,byrow=T)
toxin_colors_df<-as.data.frame(cbind(toxin_colors_df,rep(1,nrow(toxin_colors_df))))
toxin_colors_df$V2 <- factor(toxin_colors_df$V2, levels = toxin_colors_df$V2)
toxin_colors_df$V3 <- as.numeric(toxin_colors_df$V3)
#ggbarplot(toxin_colors_df, "V2","V3", fill=toxin_colors_df$V1, width = 1, xlab="",ylab="", main="Toxin Colors") + rotate_x_text(angle = 45)
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
names(toxin_colors)<-toxin_colors_df$V2

ToxinNon_colors <- colorRampPalette(c('black','gray80'))

## Cleaning Dataframe
df_clean<-function(df=TPM_df2,class="class",toxin_family="toxin_family",colors=toxin_colors){
  df2 <- df %>% mutate(class = gsub("^t.*","Toxin",!!as.name(class), ignore.case = T)) %>%
                mutate(class = gsub("^n.*","Nontoxin",!!as.name(class), ignore.case = T)) %>%
                mutate(class = gsub("^u.*","Uncharacterized",!!as.name(class), ignore.case = T))
  Toxins<-df2 %>% filter(!!as.name(class)=="Toxin")
  Others<-df2 %>% filter(!!as.name(class)!="Toxin")
  tox_names<-names(colors)
  Toxins<-Toxins %>% rowwise() %>% 
    mutate(toxin_family = ifelse(
      is.na(tox_names[amatch(toupper(!!as.name(toxin_family)), toupper(tox_names))]),
      as.character(!!as.name(toxin_family)),
      tox_names[amatch(toupper(!!as.name(toxin_family)), toupper(tox_names))]))
  new_toxs<-sort(setdiff(unique(Toxins$toxin_family), tox_names))
  j=1
  for(i in new_toxs){
    if(j<7){
      names(colors)[which(names(colors)==paste0("zOpenColor",j))]<-i
      j=j+1
    }else{
      new_col<-sample(colors(),1)
      names(new_col)<-i
      c(toxin_colors,new_col)
      j=j+1
    }
  }
  df2<-rbind(Toxins,Others)
  df2 <- as.data.frame(df2) %>% mutate_if(is.character,as.factor)
  assign("toxin_colors",colors)
  return(df2)
}

## Full transcriptome barplot
FullBarplot<-function(df=TPM_df2,id="Average",class="class",print=TRUE){
  tmp_df <- df %>% mutate_if(is.character,as.factor) %>%
    mutate(lnTPM=log(!!as.name(id))) %>%
    arrange(-lnTPM)
  transcript<-names(tmp_df)[1]
  
  pal<-ToxinNon_colors(nlevels(tmp_df[[class]]))
  names(pal)[1]<-"Toxin"
  names(pal)[2]<-"Nontoxin"
  
  A<-ggbarplot(tmp_df, transcript, "lnTPM", sort.val = "desc", sort.by.groups = F,color=NA,
               size=0, width=1, fill=class, palette = pal, legend="right",
               ylab="ln(TPM)") + #, title=paste(id,"Venom Gland Transcriptome")) + 
    rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")
  
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Expression piecharts
FullPie<-function(df=TPM_df2,id="Average",class="class",print=TRUE){
  tmp_sum<-df %>% mutate_if(is.character,as.factor) %>%
    group_by(!!as.name(class)) %>% summarize(per=sum(!!as.name(id))) %>%
    ungroup() %>% mutate(per=round((per/sum(per))*100,2)) %>%
    arrange(per) %>%
    mutate(labels=paste0(!!as.name(class),"\n",per,"%"))
  
  tmp_sum[[class]]<-factor(tmp_sum[[class]],levels=tmp_sum[[class]])
  # tmp_sum[[class]]<-relevel(tmp_sum[[class]],"Nontoxin")
  # tmp_sum[[class]]<-relevel(tmp_sum[[class]],"Toxin")
  
  pal<-ToxinNon_colors(nlevels(df[[class]]))
  names(pal)[1]<-"Toxin"
  names(pal)[2]<-"Nontoxin"
  
  A<-ggpie(tmp_sum,"per",fill = class, palette=pal,
           label="labels",lab.pos = "in", lab.font = c(5,"white","bold"))+
    theme_transparent() +  rremove("legend")
  
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Toxin barplot
ToxinBarplot<-function(df=TPM_df2,id="Average",class="class",toxin_family="toxin_family",colors=toxin_colors,print=TRUE){
  tmp_df <- df %>% mutate_if(is.character,as.factor) %>%
    filter(!!as.name(class)=="Toxin") %>%
    mutate(lnTPM=log(!!as.name(id))) %>%
    arrange(-lnTPM)
  transcript<-names(tmp_df)[1]
  labels <- tmp_df[[transcript]]
  labels<-sub("^.*[^*]","",labels)
  labels<-sub("...","*",labels)
  
  tmp_df[[1]]<-factor(tmp_df[[1]],levels=tmp_df[[1]])
  
  A<-ggbarplot(tmp_df, transcript, "lnTPM", #sort.val = "desc", sort.by.groups = F,
               size=0.5, width=1, fill=toxin_family, palette = colors,legend="right",
               ylab="ln(TPM)", label=labels) +  #title=paste(id,"Venom Gland Transcriptome")
    rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Expression piecharts
ToxinPie<-function(df=TPM_df2,id="Average",class="class",toxin_family="toxin_family",colors=toxin_colors,print=TRUE){
  tmp_sum<-df %>% mutate_if(is.character,as.factor) %>%
    filter(!!as.name(class)=="Toxin") %>% 
    group_by(!!as.name(toxin_family)) %>% 
    summarize(per=sum(!!as.name(id))) %>%
    ungroup() %>% 
    mutate(per=round((per/sum(per))*100,2)) %>% 
    arrange(per) %>%
    mutate(labels=ifelse(per<1,"",paste0(!!as.name(toxin_family),"\n",per,"%")))
  
  tmp_sum[[toxin_family]]<-factor(tmp_sum[[toxin_family]],levels=tmp_sum[[toxin_family]])
  
  A<-ggpie(tmp_sum,"per",fill = toxin_family, palette=colors,
           label="labels",lab.pos = "out") + theme(
             panel.background = element_rect(fill = "transparent"), # bg of the panel
             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
           ) + rremove("legend")
  
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}



## FancyFigure
FancyFigure <- function(df=TPM_df2,id="Average",class="class",toxin_family="toxin_family",colors=toxin_colors){
  A <- FullBarplot(df,id,class,print=FALSE)+rremove("legend")
  B <- FullPie(df,id,class,print=FALSE)
  C <- ToxinBarplot(df,id,class,toxin_family,colors,print=FALSE)
  D <- ToxinPie(df,id,class,toxin_family,colors,print=FALSE)
  
  E<-A + inset_element(B, left = 0, bottom = 0.4, right = 0.6, top = 1.1)
  E$patches$layout$widths  <- 1
  E$patches$layout$heights <- 1
  
  F<-(C+plot_spacer()+plot_layout(widths=c(1,1.5)))/E
  
  G<-F+inset_element(D, left = 0.4, bottom = 0.4, right = 1, top = 2)+plot_annotation(title = id)
  G$patches$layout$widths  <- 1
  G$patches$layout$heights <- 1
  
  #F<-C + D + plot_layout(widths = c(1, 1.25))
  #G<-F/E
  print(G)
}

## Pairwise scatterplot of clr transformed expression data
TransCompPlot<-function(df=TPM_df2,id1="CLP2057", id2="CLP2065",class="class",toxin_family="toxin_family",colors=toxin_colors,print=TRUE){
  #df2<-df %>% mutate_if(is.numeric,clr)
  df2<-df %>% mutate_if(is.numeric, function(x) as.numeric(clr(x)))
  Nontoxins<-df2 %>% filter(!!as.name(class)=="Nontoxin")
  Toxins<-df2 %>% filter(!!as.name(class)=="Toxin")
  
  minX<-min(df2[,c(id1,id2)])
  maxX<-max(df2[,c(id1,id2)])
  
  # Linear Regression
  y=Nontoxins[[id2]]
  x=Nontoxins[[id1]]
  model <- lm(y~x)
  
  # Add predictions 
  new.x=seq(minX,maxX,by = 0.05)
  new.y <- predict(model, data.frame(x=new.x), interval = "prediction")
  pred.int<-as.data.frame(cbind(new.x,new.y))
  colnames(pred.int)[1]<-c(id1)
  
  # Plot Nontoxins with regression line + prediction intervals
  p <- ggscatter(Nontoxins, x=id1, y=id2, color="#8080801A",  add="reg.line", add.params = list(color="black"),fullrange = T)+
    stat_cor(aes(label = paste("Nontoxins:",nrow(Nontoxins), ..rr.label.., ..p.label.., sep = "~` `~")))+
    geom_line(data=pred.int,aes(y = lwr), color = "black", linetype = "dashed")+
    geom_line(data=pred.int,aes(y = upr), color = "black", linetype = "dashed")
  
  # Add Toxins
  p <-ggscatter(Toxins, x=id1, y=id2, color=toxin_family,  palette=toxin_colors, legend="right",ggp=p)+
    stat_cor(data=Toxins, aes(label = paste("Toxins:",nrow(Toxins), ..rr.label.., ..p.label.., sep = "~` `~")), label.y.npc = 0.93)
  
  # Print plot
  if(print==TRUE){
    print(p)
  }
  if(print==FALSE){
    return(p)
  }
}
