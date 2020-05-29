# Venom Transcriptome Plotting Functions & Colors
# Functions by: Rhett M. Rautsaw

library(RColorBrewer)
library(ggpubr)
library(patchwork)
library(compositions)
library(dplyr)

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
                "#8B5A2B","PDI",
                "#f52727","PLA2",
                "#991919","PLA2_neuro",
                "#c3a5b4","PLB",
                "#3998f5","SVMMP",
                "#277da7","SVMPI",
                "#3750db","SVMPII",
                "#2f2aa0","SVMPIII",
                "#228c68","SVSP",
                "#aaffc3","VEGF",#946aa2
                "#aa6e28","Vespryn",#c56133
                "#cccccc","VF",
                "#f07cab","Waprin")
#                "#235b54","zOpenColor1",
#                "#37294f","zOpenColor2",
#                "#96341c","zOpenColor3",
#                "#632819","zOpenColor4",
#                "#b732cc","zOpenColor5",
#                "#8ad8e8","zOpenColor6")

toxin_colors_df<-matrix(toxin_colors,nrow=length(toxin_colors)/2,ncol=2,byrow=T)
toxin_colors_df<-as.data.frame(cbind(toxin_colors_df,rep(1,nrow(toxin_colors_df))))
toxin_colors_df$V2 <- factor(toxin_colors_df$V2, levels = toxin_colors_df$V2)
toxin_colors_df$V3 <- as.numeric(toxin_colors_df$V3)
#ggbarplot(toxin_colors_df, "V2","V3", fill=toxin_colors_df$V1, width = 1, xlab="",ylab="", main="toxin colors") + rotate_x_text(angle = 45)
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
toxin_colors<-palette(as.vector(toxin_colors_df$V1))
names(toxin_colors)<-toxin_colors_df$V2

ToxinNon_colors <- colorRampPalette(c('black','gray80'))

## Full transcriptome barplot
FullTranscriptomePlot<-function(df=TPM_df2,id="Average",class="class",print=TRUE){
  df[[id]]<-log(df[[id]])
  transcript<-names(df)[1]
  
  df[[class]]<-factor(df[[class]],levels=c("Toxin","Nontoxin",levels(df[[class]])[grep("Toxin|Nontoxin",levels(df[[class]]),invert=T)]))
  
  A<-ggbarplot(df, transcript, id, sort.val = "desc", sort.by.groups = F,
                size=0, width=1, fill=class, palette = ToxinNon_colors(nlevels(df[[class]])),legend="right",
                ylab="ln(TPM)", title=paste(id,"Venom Gland Transcriptome")) + 
    rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")

  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Toxin barplot
ToxinBarplot<-function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors,print=TRUE){
  df[[id]]<-log(df[[id]])
  df<-df[df[class]=="Toxin",]
  #df<-subset(df,df[[class]]=="Toxin")
  transcript<-names(df)[1]
  labels <- df[[transcript]]
  labels<-sub("^.*[^*]","",labels)
  labels<-sub("...","*",labels)
  
  A<-ggbarplot(df, transcript, id, sort.val = "desc", sort.by.groups = F,
               size=0.5, width=1, fill=toxin_class, palette = colors,legend="right",
               ylab="ln(TPM)", title=paste(id,"Venom Gland Transcriptome"), label=labels) + 
               rremove("xlab")+rremove("x.axis")+rremove("x.text")+rremove("x.ticks")
  if(print==TRUE){
    print(A)
  }
  if(print==FALSE){
    return(A)
  }
}

## Expression piecharts
ExpressionPie<-function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors,print=TRUE){
  tmp_sum<-aggregate(df[[id]],by=list(class=df[[class]]),FUN=sum)
  tmp_sum$label <- paste0(tmp_sum$class," ",round((tmp_sum$x/sum(tmp_sum$x))*100,2),"%")
  
  tmp_sum[[class]]<-factor(tmp_sum[[class]],levels=c("Toxin","Nontoxin",levels(tmp_sum[[class]])[grep("Toxin|Nontoxin",levels(tmp_sum[[class]]),invert=T)]))
  
  A<-ggpie(tmp_sum,"x",fill = "class", palette=ToxinNon_colors(nlevels(as.factor(df[[class]]))), label="label",lab.pos = "in", lab.font = "white")+rremove("legend")

  df2<-df[df[class]=="Toxin",]
  #df2<-subset(df,df[[class]]=="Toxin")
  tmp_sum<-aggregate(df2[[id]],by=list(class=df2[[toxin_class]]),FUN=sum)
  tmp_sum$label <- paste0(tmp_sum$class," ",round((tmp_sum$x/sum(tmp_sum$x))*100,2),"%")
  tmp_sum<-tmp_sum[order(tmp_sum$x),]
  tmp_sum$class <- factor(tmp_sum$class, levels = tmp_sum$class)
  B<-ggpie(tmp_sum,"x",fill = "class", palette=colors, label="label",lab.pos = "in", lab.font = "white")+rremove("legend")
  
  C<-A + B
  if(print==TRUE){
    print(C)
  }
  if(print==FALSE){
    return(C)
  }
}

## FancyFigure
FancyFigure <- function(df=TPM_df2,id="Average",class="class",toxin_class="toxin_class",colors=toxin_colors){
  A <- FullTranscriptomePlot(df,id,class,print=FALSE)
  B <- ToxinBarplot(df,id,class,toxin_class,colors,print=FALSE)
  C <- ExpressionPie(df,id,class,toxin_class,colors,print=FALSE)
  
  D<-(B + C) / A
  print(D)
}

## Pairwise scatterplot of clr transformed expression data
TransCompPlot<-function(df=TPM_df2,id1="CLP2057", id2="CLP2065",class="class",toxin_class="toxin_class",colors=toxin_colors,print=TRUE){
  df2<-df %>% mutate_if(is.numeric,clr)
  Nontoxins<-df2[df2[class]=="Nontoxin",]
  Toxins<-df2[df2[class]=="Toxin",]
  
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
  p <-ggscatter(Toxins, x=id1, y=id2, color=toxin_class,  palette=toxin_colors, legend="right",ggp=p)+
    stat_cor(data=Toxins, aes(label = paste("Toxins:",nrow(Toxins), ..rr.label.., ..p.label.., sep = "~` `~")), label.y.npc = 0.93)
  
  # Print plot
  if(print==TRUE){
    print(p)
  }
  if(print==FALSE){
    return(p)
  }
}
