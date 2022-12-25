## col1 -> rowname  
egg_us<-function(data){
  rownames(data)<-data[,1]
  data<-data[,-1]
  data<-as.matrix(data)
  return(data)
}

###fread -> list
freadlist<-function(site){
  data<-fread(site,data.table = FALSE)
  out<-list()
  for (i in 1:(dim(data)[2])) {
    temp<-which(data[,i]=="")
    if(length(temp)==0){
      out[[colnames(data)[i]]] <- data[,i]
    }
    else{
      out[[colnames(data)[i]]] <- data[-which(data[,i]==""),i]
    } 
    out[[colnames(data)[i]]]<-na.omit(out[[colnames(data)[i]]])
  }
  return(out)
}

##immune infiltration
Im.inf<-function(data,array=FALSE){
  library(immunedeconv)
  res_quantiseq = deconvolute(data, "quantiseq", tumor = TRUE,arrays=array)        
  res_epic = deconvolute(data,"epic", tumor = TRUE)
  res_xcell = deconvolute(data, "xcell", tumor = TRUE,arrays=array)
  res_mcp_counter<-deconvolute_mcp_counter(data,feature_types = "HUGO_symbols")
  res_timer<-deconvolute_timer(data,indications=rep('lihc',dim(data)[2]))
  quantiseq<-as.data.frame(res_quantiseq)
  rownames(quantiseq)<-quantiseq[,1]
  quantiseq<-quantiseq[,-1]
  epic<-as.data.frame(res_epic)
  rownames(epic)<-epic[,1]
  epic<-epic[,-1]
  xcell<-as.data.frame(res_xcell)
  rownames(xcell)<-xcell[,1]
  xcell<-xcell[,-1]
  mcp_counter<-as.data.frame(res_mcp_counter)
  timer<-as.data.frame(res_timer)
  out<-list(timer=timer,mcp_counter=mcp_counter,xcell=xcell,epic=epic,quantiseq=quantiseq)
  return(out)
}
##limma
DEG_limma <- function(exp,label,log=FALSE){
  data<-exp[,-1]
  rownames(data)<- exp[,1]
  if(log==FALSE){
    data<-log2(data+1)
  }
  aaa<-label[,]
  design <-model.matrix(~0+factor(aaa))#把group设置成一个model matrix#
  colnames(design)=levels(factor(aaa))
  rownames(design)=colnames(data)
  fit <- lmFit(data,design)
  cont.matrix<-makeContrasts(paste0(c("C2","C1"),collapse = "-"),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ###
  ##
  allDEG <- topTable(fit2,n = Inf,p.value = 0.05)
  ##allDEG <- topTable(fit2,n = Inf)
  DEG<-cbind(rownames(allDEG),allDEG[,1]>0)
  #return(allDEG)
  gene_low <- DEG[which(DEG[,2]==TRUE),1]
  gene_high <- DEG[which(DEG[,2]==FALSE),1]
  return(list(gene_low=gene_low,gene_high=gene_high))
}

########violin
library(rlang)

transvar <- function(xglue,x){
  #xglue <- englue("{{ x }}")
  if(grepl("\"", xglue)){
    return(x)
  }
  else{
    if(xglue %in% c("NULL","NA")){
      return(x)
    }
    else{
      return(xglue)
    }
  }
}




library(rlang)
plot.violin <- function(data,x,y,fill=NA,color=NA,
                        split=NA,
                        free="free_x",
                        group=NA,
                        label="p.signif",
                        sizeP=5,
                        hide.ns=T,
                        scale = "width", #"area", "count", "width"
                        palette.fill=NULL,
                        palette.color=NULL,
                        alpha.v=0.5,
                        
                        add=c("violinplot","boxplot"), #"shadow","point"
                        alpha.b=1,
                        width.b=0.3,
                        color.b=T,
                        fill.b=T,
                        lwd.b=1,
                        notch.b=F,
                        outlier.shape = 19,
                        theme=theme_bw(),
                        font_xy_title = 14,
                        font_xy_text = 14,
                        font_legend_title = 14,
                        font_legend_text = 14,
                        font_title=15,
                        legend_title = " ",
                        hjust_title=0,
                        x_label=NULL,
                        y_label=NULL,
                        title=NULL,
                        shadow.col=c("#E64B35FF","#4DBBD5FF"),
                        shadow.alpha=0.15,
                        
                        shape.p=21,
                        size.p=2,
                        alpha.p=1,
                        fill.p=T,
                        color.p=T,
                        
                        angle.x=0,
                        hjust.x=0.5,
                        vjust.x=1
){
  ##string
  x = transvar(englue("{{ x }}"), x)
  y = transvar(englue("{{ y }}"), y)
  fill = transvar(englue("{{ fill }}"), fill)
  color = transvar(englue("{{ color }}"), color)
  ##binary
  bifil = fill %in% colnames(data)
  bicol = color %in% colnames(data)
  ##base_plot
  p1 <- ggplot( data,aes_string(x, y) )
  
  ##violinplot
  if( c("violinplot") %in% add ){
    ###color
    if( is.null(palette.fill) & bifil){
      palette.fill = ggsci::pal_nejm("default")(length(unique(data[,fill])))
    }
    if( is.null(palette.color) & bicol){
      palette.color=ggsci::pal_nejm("default")(length(unique(data[,color])))
    }
    #
    #
    if(bifil){
      if(bicol){
        p1<-p1 + geom_violin(scale=scale,aes_string(fill= fill,color=color),
                             position = position_dodge(1),
                             alpha=alpha.v)+
          scale_fill_manual(values=palette.fill)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1 + geom_violin(scale=scale,aes_string(fill= fill),color=color,
                             position = position_dodge(1),
                             alpha=alpha.v)+
          scale_fill_manual(values=palette.fill)
      }
    }
    else{
      if(bicol){
        p1<-p1 + geom_violin(scale=scale,aes_string(color=color),fill=fill,
                             position = position_dodge(1),
                             alpha=alpha.v)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1 + geom_violin(scale=scale,color=color,fill=fill,
                             position = position_dodge(1),
                             alpha=alpha.v)
      }
    }
  }
  
  
  #### boxplot
  if(c("boxplot") %in% add){
    #
    fill.b = transvar(englue("{{ fill.b }}"), fill.b)
    color.b = transvar(englue("{{ color.b }}"), color.b)
    if(fill.b=="T"){
      fill.b <- fill
    }
    else{
      if(fill.b %in% colnames(data)){
        fill.b <- fill.b
      }
    }
    if(color.b=="T"){
      color.b <- color
    }
    else{
      if(color.b %in% colnames(data)){
        color.b <- color.b
      }
    }
    #
    
    if(length(palette.fill)==0 & (fill.b %in% colnames(data))){
      palette.fill = ggsci::pal_nejm("default")(length(unique(data[,fill.b])))
    }
    if(length(palette.color)==0 & (color.b %in% colnames(data))){
      palette.color=ggsci::pal_nejm("default")(length(unique(data[,color.b])))
    }  
    #
    if(fill.b %in% colnames(data)){
      if(color.b %in% colnames(data)){
        p1<-p1+geom_boxplot(width=width.b,position = position_dodge(1),
                            aes_string(fill=fill.b,color=color.b),
                            alpha=alpha.b,lwd=lwd.b,outlier.shape = outlier.shape,notch = notch.b)+
          scale_fill_manual(values=palette.fill)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1+geom_boxplot(width=width.b,position = position_dodge(1),
                            aes_string(fill=fill.b),color=color.b,
                            alpha=alpha.b,lwd=lwd.b,outlier.shape = outlier.shape,notch = notch.b) + 
          scale_fill_manual(values=palette.fill)
      }
    }
    else{
      if(color.b %in% colnames(data)){
        p1<-p1+geom_boxplot(width=width.b,position = position_dodge(1),
                            aes_string(color=color.b),fill=fill.b,
                            alpha=alpha.b,lwd=lwd.b,outlier.shape = outlier.shape,notch = notch.b)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1+geom_boxplot(width=width.b,position = position_dodge(1),
                            fill=fill.b,color=color.b,
                            alpha=alpha.b,lwd=lwd.b,outlier.shape = outlier.shape,notch = notch.b)
      }
    }
  }
  
  
  
  #### point
  if(c("point") %in% add){
    #
    fill.p = transvar(englue("{{ fill.p }}"), fill.p)
    color.p = transvar(englue("{{ color.p }}"), color.p)
    if(fill.p=="T"){
      fill.p <- fill
    }
    else{
      if(fill.p %in% colnames(data)){
        fill.p <- fill.p
      }
    }
    if(color.p=="T"){
      color.p <- color
    }
    else{
      if(color.p %in% colnames(data)){
        color.p <- color.p
      }
    }
    #
    
    if(length(palette.fill)==0 & (fill.p %in% colnames(data))){
      palette.fill = ggsci::pal_nejm("default")(length(unique(data[,fill.p])))
    }
    if(length(palette.color)==0 & (color.p %in% colnames(data))){
      palette.color=ggsci::pal_nejm("default")(length(unique(data[,color.p])))
    }  
    #
    if(fill.p %in% colnames(data)){
      if(color.p %in% colnames(data)){
        p1<-p1+geom_point(aes_string(fill=fill.p,color=color.p),shape = shape.p, size=size.p, 
                          position = position_jitterdodge(dodge.width=1), alpha=alpha.p)+
          scale_fill_manual(values=palette.fill)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1+geom_point(aes_string(fill=fill.p),color=color.p,
                          shape = shape.p, size=size.p, 
                          position = position_jitterdodge(dodge.width=1), alpha=alpha.p) + 
          scale_fill_manual(values=palette.fill)
      }
    }
    else{
      if(color.p %in% colnames(data)){
        p1<-p1+geom_point(aes_string(color=color.p),fill=fill.p,
                          shape = shape.p, size=size.p, 
                          position = position_jitterdodge(dodge.width=1), alpha=alpha.p)+
          scale_colour_manual(values=palette.color)
      }
      else{
        p1<-p1+geom_point(fill=fill.p,color=color.p,
                          shape = shape.p, size=size.p, 
                          position = position_jitterdodge(dodge.width=1), alpha=alpha.p)
      }
    }
  }
  
  
  
  ###shadow
  ##shading
  if( c("shadow") %in% add ){
    len <- length(unique(data[,x]))
    if(len==2){
      shading_C1 <- data.frame(
        min = c(-Inf),
        max=  c(1.5)
      )
    }
    else{
      shading_C1 <- data.frame(
        min = c(-Inf, seq(from=2.5, to=1*(len-1)+0.5, by=2) ),
        max=  c(seq(from=1.5, to=1*(len)+0.5, by=2))
      )
    }
    shading_C2<-data.frame(
      min = c(seq(from=1.5, to=1*(len-1)+0.5, by=2)),
      max=  c(seq(from=2.5, to=1*(len)+0.5, by=2))
    )
    if(max(shading_C1$max) > max(shading_C2$max)){
      shading_C1$max[length(shading_C1$max)] <- Inf
    }
    else{
      shading_C2$max[length(shading_C2$max)] <- Inf
    }
    
    p1<- p1+ geom_rect(data = shading_C1, inherit.aes = FALSE,
                       aes(xmin = min, xmax = max, ymin = -Inf, 
                           ymax = Inf),fill = shadow.col[1],alpha = shadow.alpha,show.legend=F)+
      geom_rect(data = shading_C2, inherit.aes = FALSE,
                aes(xmin = min, xmax = max, ymin = -Inf, 
                    ymax = Inf),fill=shadow.col[2],alpha = shadow.alpha,show.legend=F)
  }
  
  ##pvalue
  symnum.args<-list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1,1),
                    symbols = c("****", "***", "**", "*", "+","ns"))
  bista = englue("{{ group }}") %in% colnames(data)
  if(bista){
    p1 <- p1 + stat_compare_means(aes(group = {{ group }}),label = label,vjust=0.8,hide.ns=hide.ns,size=sizeP)
  }
  
  
  ###theme and font
  p1 <- p1+theme+ggtitle(title)+
    theme(axis.title.x =element_text(size=font_xy_title,color = "black"), 
          axis.title.y=element_text(size=font_xy_title,color = "black"),
          axis.text.x=element_text(size=font_xy_text,color = "black",angle = angle.x,hjust=hjust.x, vjust=vjust.x),
          axis.text.y=element_text(size=font_xy_text,color = "black"),
          legend.title=element_text(size=font_legend_title,color = "black"),legend.text = element_text(size=font_legend_text,color = "black"),
          plot.title = element_text(hjust=hjust_title, color="black", size=font_title))
  if(length(x_label)!=0){
    p1 <- p1 + labs(x=x_label)
  }
  if(length(y_label)!=0){
    p1 <- p1 + labs(y=y_label)
  }
  
  return(p1)
}


## dataset surplot function
Surplot <- function(data_demo,type,title){
  if(type=="OS"){
    surv_TTP<-surv_fit(Surv(OS.time, OS) ~ subtype, data =data_demo)
    xlab = "Overall Survival (days)"
  }
  else if(type=="PFI"){
    surv_TTP<-surv_fit(Surv(PFI.time, PFI) ~ subtype, data =data_demo)
    xlab = "Progression Free Interval (days)"
  }
  else{
    surv_TTP<-surv_fit(Surv(RFS.time, RFS) ~ subtype, data =data_demo)
    xlab = "Relapse Free Survival (days)"
  }
  TitleNames <- "Subtype"  
  title <- title  
  legend_label <- c(paste0("C1 n=",table(data_demo$subtype)[1]),
                    paste0("C2 n=",table(data_demo$subtype)[2]))
  ggsurvplot(surv_TTP,
             pval = TRUE, 
             pval.size = 6.5,  
             pval.coord=c(3,5),
             conf.int = F, 
             fun = "pct", 
             palette =c("#CD1818","#0F2C67"),
             risk.table = F, 
             xlab = xlab,  
             font.tickslab = 16, 
             font.x = 18,font.y = 18,  
             font.subtitle = 18, 
             font.legend = 15,  
             legend = c(0.8,0.9),  
             legend.labs = legend_label,
             censor.size = 0,  
             size = 1.8,  
             axes.offset = T, 
             legend.title = TitleNames,
             main = "Survival curves",
             submain = title)
}
###
xlsx.list<-function(dir){
  temp<-excel_sheets(dir)
  out<-list()
  for (i in 1:length(temp)) {
    out[[temp[i]]]<-read_excel(dir,sheet=i)
  }
  names(out)<-temp
  return(out)
}
###
imm.select<-function(inf.all,name,data){
  out<-list()
  for (i in 1:length(names(data))) {
    a<-intersect(data[[i]][,1,drop=TRUE],colnames(inf.all[[name]]))
    #data[[i]][,1,drop=TRUE]
    temp<-inf.all[[name]][,c("subtype",a)]
    temp<-del_immune(temp)  
    if(length(temp)!=0){
      colnames(temp)[-1]<-data[[i]][match(colnames(temp)[-1],data[[i]][,1,drop=TRUE]),2,drop=TRUE]
      out[[i]]<-imm.inf(temp)
    }
  }
  names(out)<-names(data)
  out1<-lapply(out,function(x) x$part)
  out2<-lapply(out,function(x) x$out)
  out1<-do.call('rbind',lapply(out1,function(x) x))
  out2<-do.call('rbind',lapply(out2,function(x) x))
  out3<-list(part=out,total=out1,mydata=out2)
  return(out3)
}
###
del_immune<-function(data){
  if(length(grep("epic",colnames(data)))!=0){
    data<-data[,-grep("epic",colnames(data))]
  }
  a<-as.data.frame(apply(data[,-1],2,function(x){
    sum(round(x,10)==0)
  }))
  temp<-which(a>(length(data[,1])*0.2))
  if(length(temp)!=0){
    abc<-data[,-(temp+1)]
  }
  else{
    abc<-data
  }
  return(abc)
}
####
imm.inf<-function(data){
  library(dplyr)
  mydata<-data
  mydata <-reshape2::melt(mydata,id="subtype")
  ##
  pvalue<-compare_means(value~subtype,group.by="variable",method = "t.test",data=mydata)$p
  aaa<-mydata %>% group_by(variable,subtype) %>% summarise(mean = mean(value)) 
  temp<-reshape2::dcast(aaa, variable  ~ subtype,value.var = 'mean' ) %>% 
    mutate(df=C1-C2,pvalue=pvalue,status=((df>0)-0.5)*2,value=status*((pvalue<0.05)+1-1),
           p_dir=pvalue*status)
  temp<-list(part=temp,out=mydata)
  return(temp)
}
###替换函数
replace_egg<-function(data,col=NULL,to1,to2=NULL,from1=0,from2=1){
  if(is.null(to2)){
    data<-matrix(data,ncol = 1)
  }
  if(is.null(to2)){
    data[which(data[,col]==from1),col]<-to1
  }
  else{
    data[which(data[,col]==from1),col]<-to1
    data[which(data[,col]==from2),col]<-to2
  }
  return(data)
}

###GSEA Kegg
gse_kegg<-function(data,label,log=FALSE){
  set.seed(123)
  if(log){
    data1<-data[,-1]
  }
  else{
    data1<-log2(data[,-1]+1) 
  }
  mean1<-apply(data1[,which(label=="C1")], 1, mean)
  mean2<-apply(data1[,which(label=="C2")], 1, mean)
  mean3<-mean1-mean2
  genelist<-data.frame(data[,1],mean3)
  genelist<-as.data.frame(genelist)
  colnames(genelist)<-c("id","logFC")
  geneList_NRVC.sort <- arrange(genelist, desc(logFC))
  geneList_NRVC<-geneList_NRVC.sort$logFC
  names(geneList_NRVC)<-geneList_NRVC.sort[,1]
  out<-gseKEGG(geneList_NRVC,use_internal_data =T)
  return(out)
}
###
gsea_egg<-function(data,label,geneset,log=FALSE){
  if(log){
    data1<-data[,-1]
  }
  else{
    data1<-log2(data[,-1]+1) 
  }
  mean1<-apply(data1[,which(label=="C1")], 1, mean)
  mean2<-apply(data1[,which(label=="C2")], 1, mean)
  mean3<-mean1-mean2
  genelist<-data.frame(data[,1],mean3)
  genelist<-as.data.frame(genelist)
  colnames(genelist)<-c("id","logFC")
  geneList_NRVC.sort <- arrange(genelist, desc(logFC))
  geneList_NRVC<-geneList_NRVC.sort$logFC
  names(geneList_NRVC)<-geneList_NRVC.sort[,1]
  egmt <- GSEA(geneList_NRVC, TERM2GENE=geneset,verbose=F,pvalueCutoff =1.1)
  return(egmt)
}
###cox
cox_s_m<-function(sur,label){
  subtype<-ifelse(label=="C1","1","0")
  subtype<-as.character(subtype)
  cli<-dplyr::mutate(sur[,],subtype)
  #return(cli)
  if(length(cli$grade)==0){
    single_name<-c("stage","age","gender","subtype")
    cox_m<-coxph(Surv(OS.time,OS) ~ stage+age+gender+subtype,data=cli)
  }
  else{
    single_name<-c("stage","age","gender","grade","subtype")
    cox_m<-coxph(Surv(OS.time,OS) ~ stage+age+gender+grade+subtype,data=cli)
  }
  a<-summary(cox_m)
  multi_cox<-cbind(as.data.frame(a$conf.int),as.data.frame(a$coefficients[,5,drop=FALSE]))
  multi_cox1<-multi_cox
  multi_cox<-round(multi_cox,3)
  multi_cox<-data.frame(multi_cox,other=paste0(multi_cox[,1],"(",multi_cox[,3],"-",multi_cox[,4],")"))
  cox_s<-lapply(single_name,function(name){
    x<-as.formula(paste('Surv(OS.time,OS)~',name))
    cox<-coxph(x,data=cli)
    a<-summary(cox)
    single_cox<-cbind(as.data.frame(a$conf.int),as.data.frame(a$coefficients[,5,drop=FALSE]))
    return(single_cox)
  })
  cox_s<-do.call(rbind,cox_s)
  cox_s1<-cox_s
  cox_s<-round(cox_s,3)
  cox_s<-data.frame(cox_s,other=paste0(cox_s[,1],"(",cox_s[,3],"-",cox_s[,4],")"))
  cox_result<-list(single=cox_s,multi=multi_cox,single1=cox_s1,multi1=multi_cox1)
  return(cox_result)
}
###
APM_24IMMUNE<-function(data,label,name){
  TIS<-apply(data[1:9,],2,mean)
  IIS<-apply(data[-which(rownames(data) %in% c("APM","TFH","tgd")),],2,mean)
  APM<-data["APM",]
  temp<-data.frame(label,TIS,IIS,APM)
  temp<-reshape2::melt(temp)
  temp<-data.frame(temp,Dataset=name)
  return(temp)
}
###IPS
IPS<-function(data,dir){
  library(ggplot2)
  library(grid)
  library(gridExtra)
  data<-egg_us(data)
  ## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
  ## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
  if(name=="TCGA" || name=="ICGC"){
    gene_expression<-as.data.frame(log2(data+1))
  }
  else{
    gene_expression<-as.data.frame(data)
  }
  sample_names<-names(gene_expression)
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  # For different 
  IPSG<-read.table(dir,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
  unique_ips_genes<-as.vector(unique(IPSG$NAME))
  IPS<-NULL
  MHC<-NULL
  CP<-NULL
  EC<-NULL
  SC<-NULL
  AZ<-NULL
  # Gene names in expression file
  GVEC<-row.names(gene_expression)
  # Genes names in IPS genes file
  VEC<-as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind<-which(is.na(match(VEC,GVEC)))
  # List genes missing or differently named
  MISSING_GENES<-VEC[ind]
  dat<-IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }
  for (i in 1:length(sample_names)) {	
    GE<-gene_expression[[i]]
    mGE<-mean(GE)
    sGE<-sd(GE)
    Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
    W1<-IPSG$WEIGHT
    WEIGHT<-NULL
    MIG<-NULL
    k<-1
    for (gen in unique_ips_genes) {
      MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k<-k+1
    }
    WG<-MIG*WEIGHT
    MHC[i]<-mean(WG[1:10])
    CP[i]<-mean(WG[11:20])
    EC[i]<-mean(WG[21:24])
    SC[i]<-mean(WG[25:26])
    AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
    IPS[i]<-ipsmap(AZ[i])
  }
  DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
  #write.table(DF,file="IPS.txt",row.names=FALSE, quote=FALSE,sep="\t")
  return(DF)
}
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

## Assign colors 
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

##
trans_score1 <- function(data,ratio,label,gene_set,len){
  rownames(data)<-data[,1]
  data<-data[,-1]
  data<-as.matrix(data)
  result<- gsva(data,gene_set,method='ssgsea')
  result1 <-result %>% pheatmap:::scale_rows()
  high<-result1[c(1),,drop=FALSE]
  low<-result1[c(2),,drop=FALSE]
  out<-data.frame(C1_score=t(high),C2_score=t(low),C1_C2=t(high-low),ratio=ratio)
  colnames(out)<-c("C1_score","C2_score","C1score_C2score","Ratio")
  out$class <- ifelse(label=="C1" & ((out$C1_score-out$C2_score)>0) =="TRUE","C1",
                      ifelse(label=="C2" & ((out$C1_score-out$C2_score)>0) == "FALSE","C2","Missubtype"))
  return(out)
}
##
point_section <- function(data, x, y, colour,
                          colour_value = c("#CD1818", "#0F2C67", "yellow"),
                          annotate="",
                          sizeP = 2,
                          sizeAnno= 5,
                          font_xy_label = 14,
                          font_xy_text = 14,
                          font_legend_title = 14,
                          font_legend_text = 14,
                          legend_title = " ",
                          x_label=NULL,
                          y_label=NULL,
                          slope = 1,
                          intercept = 0,
                          col_section_up = "#E64B35FF",
                          col_section_down = "#4DBBD5FF",
                          xintercept=0,
                          yintercept=0,
                          xlinetype="dashed",
                          ylinetype="dashed",
                          ablinetype="solid",
                          abline=T,
                          xlwd=0.8,
                          ylwd=0.8,
                          ablwd=1,
                          cor=F,
                          cor_method="spearman",
                          cor_colour="red",
                          xprobs=0.01,
                          yprobs=0.001
){
  requireNamespace("tidyr",quietly = TRUE)
  requireNamespace("rlang",quietly = TRUE)
  requireNamespace("dplyr",quietly = TRUE)
  requireNamespace("ggstatsplot",quietly = TRUE)
  if(cor){
    p1 <- ggscatterstats(data,{{x}},{{y}},marginal=F,results.subtitle=F,
                         point.args = list(size = sizeP))+
      stat_cor(method = cor_method,size=5,color = cor_colour)
  }
  else{
    p1 <- ggplot(data = data,aes(x= {{x}},y= {{y}} ))
  }
  
  p1 <- p1+
    geom_point(size=sizeP, aes(colour=factor({{colour}})))+
    geom_hline(yintercept = yintercept,lwd=xlwd,linetype=xlinetype)+
    geom_vline(xintercept = xintercept,lwd=ylwd,linetype=ylinetype)+
    geom_section(slope=slope, intercept=intercept, above=F,fill=col_section_up,alpha=0.15)+
    geom_section(slope=slope, intercept=intercept, above=T,fill=col_section_down,alpha=0.15)+
    theme_bw()+
    scale_colour_manual(values=colour_value)+
    font("xy", size = font_xy_label)+
    font("xy.text", size = font_xy_text,color = "black")+
    font("legend.text", size = font_legend_title)+
    font("legend.title", size = font_legend_text)+
    labs(colour=legend_title,x=x_label,y=y_label)+
    annotate('text',
             x=quantile(data[,englue("{{ x }}")], probs = xprobs),
             y=quantile(data[,englue("{{ y }}")], probs = yprobs),
             label=annotate,
             size=5,color='black',hjust = "left")
  if(!abline | slope %in% c(Inf,-Inf)){
    p1
  }
  else{
    p1+geom_abline(intercept=intercept,slope=slope,lwd=ablwd,linetype=ablinetype)
  }
}
##
geom_section <- function (mapping = NULL, data = NULL, ..., slope, intercept, above,
                          na.rm = FALSE, show.legend = NA) {
  
  buildPoly <- function(slope, intercept, above, xr, yr){
    # By Joran Elias, @joran https://stackoverflow.com/a/6809174/1870254
    #Find where the line crosses the plot edges
    yCross <- (yr - intercept) / slope
    xCross <- (slope * xr) + intercept
    
    #Build polygon by cases
    if (above & (slope >= 0)){
      rs <- data.frame(x=-Inf,y=Inf)
      if (xCross[1] < yr[1]){
        rs <- rbind(rs,c(-Inf,-Inf),c(yCross[1],-Inf))
      }
      else{
        rs <- rbind(rs,c(-Inf,xCross[1]))
      }
      if (xCross[2] < yr[2]){
        rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,Inf))
      }
      else{
        rs <- rbind(rs,c(yCross[2],Inf))
      }
    }
    if (!above & (slope >= 0)){
      rs <- data.frame(x= Inf,y= -Inf)
      if (xCross[1] > yr[1]){
        rs <- rbind(rs,c(-Inf,-Inf),c(-Inf,xCross[1]))
      }
      else{
        rs <- rbind(rs,c(yCross[1],-Inf))
      }
      if (xCross[2] > yr[2]){
        rs <- rbind(rs,c(yCross[2],Inf),c(Inf,Inf))
      }
      else{
        rs <- rbind(rs,c(Inf,xCross[2]))
      }
    }
    if (above & (slope < 0)){
      rs <- data.frame(x=Inf,y=Inf)
      if (xCross[1] < yr[2]){
        rs <- rbind(rs,c(-Inf,Inf),c(-Inf,xCross[1]))
      }
      else{
        rs <- rbind(rs,c(yCross[2],Inf))
      }
      if (xCross[2] < yr[1]){
        rs <- rbind(rs,c(yCross[1],-Inf),c(Inf,-Inf))
      }
      else{
        rs <- rbind(rs,c(Inf,xCross[2]))
      }
    }
    if (!above & (slope < 0)){
      rs <- data.frame(x= -Inf,y= -Inf)
      if (xCross[1] > yr[2]){
        rs <- rbind(rs,c(-Inf,Inf),c(yCross[2],Inf))
      }
      else{
        rs <- rbind(rs,c(-Inf,xCross[1]))
      }
      if (xCross[2] > yr[1]){
        rs <- rbind(rs,c(Inf,xCross[2]),c(Inf,-Inf))
      }
      else{
        rs <- rbind(rs,c(yCross[1],-Inf))
      }
    }
    return(rs)
  }
  GeomSection <- ggplot2::ggproto("GeomSection", GeomPolygon,
                                  default_aes = list(fill="blue", size=0, alpha=0.2, colour=NA, linetype="dashed"),
                                  required_aes = c("slope", "intercept", "above"),
                                  draw_panel = function(data, panel_params, coord) {
                                    ranges <- coord$backtransform_range(panel_params)
                                    data$group <- seq_len(nrow(data))
                                    data <- data %>% group_by_all %>% do(buildPoly(.$slope, .$intercept, .$above, ranges$x, ranges$y)) %>% unnest
                                    GeomPolygon$draw_panel(data, panel_params, coord)
                                  }
  )
  
  if (missing(mapping) && missing(slope) && missing(intercept) && missing(above)) {
    slope <- 1
    intercept <- 0
    above <- TRUE
  }
  if (!missing(slope) || !missing(intercept)|| !missing(above)) {
    if (missing(slope))
      slope <- 1
    if (missing(intercept))
      intercept <- 0
    if (missing(above))
      above <- TRUE
    data <- data.frame(intercept = intercept, slope = slope, above=above)
    mapping <- aes(intercept = intercept, slope = slope, above=above)
    show.legend <- FALSE
  }
  layer(data = data, mapping = mapping, stat = StatIdentity,
        geom = GeomSection, position = PositionIdentity, show.legend = show.legend,
        inherit.aes = FALSE, params = list(na.rm = na.rm, ...))
}


####
Gene_surplot <- function(data, pdata, gene, time, event, cutoff_method="best",
                         palette =c("#CD1818","#0F2C67"),
                         xlab="Overall Survival (days)",
                         count=T,
                         type.sur="type1",
                         P_HR=F,
                         ...){
  #
  exp <- data[match(gene,data[,1]),-1]
  exp <- as.numeric(exp)
  #
  if(mean(exp)==0){
    return(NULL)
  }
  ##sur data
  pdata <- pdata[,c(time,event)]
  sur <- data.frame(gene=exp,pdata)
  colnames(sur)[1] <- gene
  
  ##
  sur <- structure(sur, class = c(class(sur),cutoff_method))
  #return(sur)
  twosurdata <- sur_cutoff(sur,time,event)
  
  if(P_HR){   
    twosurdata[,3]<-factor(twosurdata[,3],levels = c("low","high"))
    FORMULA <-as.formula(paste("Surv(",time,",",event,") ~ ",gene))
    sdf <- survdiff(formula=FORMULA,data=twosurdata)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    res.cox <- coxph(FORMULA, data = twosurdata)
    HR<-(summary(res.cox))[["coefficients"]][2]
    out <- cbind(p.val,HR) 
    rownames(out) <- gene
    return(out)
  }

  #return(twosurdata)
  #return(twosurdata)
  FORMULA <-as.formula(paste("Surv(",time,",",event,") ~ ",gene))
  surv_TTP<-surv_fit(FORMULA, data =twosurdata)
  ###
  xlab = xlab
  title <- gene
  if(count){
    legend_label <- c(paste0("High n=",table(twosurdata[[gene]])[1]),
                      paste0("Low n=",table(twosurdata[[gene]])[2]))
  }
  else{
    legend_label <- c('High',"Low")
  }
  
  surv_TTP <- structure(surv_TTP, class = c(class(surv_TTP),type.sur))
  p1 <- surplot(surv_TTP,palette,title,xlab,legend_label,...)
  return(p1)
}
sur_cutoff <- function(sur,...){
  UseMethod("sur_cutoff")
}
##
sur_cutoff.best <- function(sur,time,event, ...){
  res.cut1 <- surv_cutpoint(sur, time=time , event =event, variables = colnames(sur)[1])
  res.cat <- surv_categorize(res.cut1)
  return(res.cat)
}
##
sur_cutoff.median <- function(sur,  ...){
  cutoff <- median(sur[,1])
  out <- data.frame(sur[,-1],ifelse(sur[,1]>cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}
##
sur_cutoff.mean <- function(sur,  ...){
  cutoff <- mean(sur[,1])
  out <- data.frame(sur[,-1],ifelse(sur[,1]>cutoff,"high","low"))
  colnames(out)[3] <- colnames(sur)[1]
  return(out)
}
####
surplot <- function(surv_TTP,...){
  UseMethod("surplot")
}
surplot.type1 <- function(surv_TTP,palette,title,xlab,legend_label,legend.coord=c(0.8,0.9),...){
  ggsurvplot(surv_TTP,
             pval = TRUE, 
             pval.size = 6.5,  
             pval.coord=c(3,5),
             conf.int = F, 
             fun = "pct", 
             palette =palette,
             risk.table = F, 
             xlab = xlab,  
             font.tickslab = 16, 
             font.x = 18,font.y = 18,  
             font.subtitle = 18, 
             font.legend = 15,  
             legend =legend.coord,  
             legend.labs = legend_label,
             censor.size = 0,  
             size = 1.8,  
             axes.offset = T,
             main = "Survival curves",
             submain = title)
}

##
#####ligand–receptor
network_ratio<-function(name,data,predict,all_int,ciber,ciber_cbind){
  ciber1<-sapply(ciber_cbind,function(x){
    apply(ciber[[name]][,x,drop=FALSE],1,sum)
  })
  print("11111")
  TCGA<-as.data.frame(t(egg_us(data)))
  TCGA<-data.frame(ciber1,TCGA)
  test<-apply(TCGA,2,function(x){
    low<-(x<=quantile(x,c(1/3)))
    median<-((x>quantile(x,c(1/3))) & (x<=quantile(x,c(2/3))))
    high<-(x>quantile(x,c(2/3)))
    x[low]<-"low"
    x[median]<-"median"
    x[high]<-"high"
    return(x)
  })
  print("22222")
  test<-as.data.frame(test) 
  all_int<-all_int[(all_int[,1,drop=TRUE] %in% colnames(test)) & 
                     (all_int[,2,drop=TRUE] %in% colnames(test)),]  
  label<-predict$Subtype
  high<-test[which(label=="C1"),]
  low<-test[which(label=="C2"),]
  high_ratio<-ratio_concordance(all_int,high)
  print("33333")
  low_ratio<-ratio_concordance(all_int,low)
  out<-list(high=high_ratio,low=low_ratio)
  return(out)
}

ratio_concordance<-function(all_int,data_class){
  temp3<-apply(all_int, 1, function(x){
    temp1<-data_class[,as.character(x)]
    hh<-sum((temp1[,1]=="high") & (temp1[,2]=="high"))
    hl<-sum((temp1[,1]=="high") & (temp1[,2]=="low"))
    lh<-sum((temp1[,1]=="low") & (temp1[,2]=="high"))
    ll<-sum((temp1[,1]=="low") & (temp1[,2]=="low"))
    concordance<-(hh+ll)/(hl+lh+1)
    col1<-sum((temp1[,1]=="high") | (temp1[,1]=="median"))/length(temp1[,1])
    col2<-sum((temp1[,2]=="high") | (temp1[,2]=="median"))/length(temp1[,1])
    ratio<-data.frame(From.ratio=col1,To.ratio=col2,concordance=concordance)
    return(ratio)
  })
  temp3<-do.call(rbind,lapply(temp3, function(x) x))
  node<-data.frame(all_int,temp3)
  return(node)
}
###immunotherapy surplot
imm_surplot <- function(data,dataset,type,title,class=c("dataset","type")){
  if(class=="dataset"){
    data_demo <- data[which(data$dataset %in% dataset),]
  }
  else{
    data_demo <- data[which(data$type %in% type),]
  }
  surv_TTP<-surv_fit(Surv(os_time, os_event) ~ subtype, data =data_demo)
  xlab = "Overall Survival (days)"
  TitleNames <- "Subtype"  
  title <- title  
  legend_label <- c(paste0("C1 n=",table(data_demo$subtype)[1]),
                    paste0("C2 n=",table(data_demo$subtype)[2]))
  ggsurvplot(surv_TTP,
             pval = TRUE, 
             pval.size = 6.5,  
             pval.coord=c(3,5),
             conf.int = F,
             fun = "pct", # 
             palette =c("#CD1818","#0F2C67"),
             risk.table = F, 
             xlab = xlab,  
             font.tickslab = 16, 
             font.x = 18,font.y = 18,  
             font.subtitle = 18, 
             font.legend = 15, 
             legend = c(0.8,0.9), 
             legend.labs = legend_label,
             censor.size = 0, 
             size = 1.8, 
             axes.offset = T, 
             legend.title = TitleNames,
             main = "Survival curves",
             submain = title)
}


##
df_pre<-function(data,reverse=F){
  if(reverse){
    data$out <- paste0("->",data$target)
    data$inner <- paste0(data$source," (",data$dataset,")")
    data$direction <- "reverse"
  }
  else{
    data$out <- paste0(data$source,"->")
    data$inner <- paste0(data$target," (",data$dataset,")")
    data$direction <- "no_reverse"
  }
  return(data)
}
##
LR_plot <- function(df){
  ##out level
  df$out <- factor(df$out,levels = unique(df$out))
  ##inner level
  cells.level <- levels(cellchat@idents$C1)
  sources.use=1:26
  targets.use=1:26
  sources.use <- cells.level[sources.use]
  targets.use <- cells.level[targets.use]
  source.target <- paste0(rep(targets.use,each=2),
                          rep(c(" (C1)"," (C2)"),time=length(targets.use)))
  inner_level<-unique(df$inner)[na.omit(match(source.target,unique(df$inner)))]
  df$inner <- factor(df$inner,inner_level)
  
  index<-unique(df[order(df$out,df$inner),c("dataset","out","inner","direction")])
  #color_level<-index[,1]
  #color_index<-match(color_level, c("C1","C2"))
  #color.text	= c("black","#CD1818")
  color_level <- index$direction
  color_index <- match(color_level, c("reverse","no_reverse"))
  color.text	= c("#CD1818","black")
  
  n.colors=10
  color.use <- tryCatch({
    RColorBrewer::brewer.pal(n = n.colors, name ="Spectral")
  }, error = function(e) {
    (scales::viridis_pal(option = "Spectral", direction = -1))(n.colors)
  })
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  ##dashed line
  dashed <- match(unique(index$out),index$out)[-1]-1+0.5
  ##shading
  shading_C1<-data.frame(
    min = which(index$dataset=="C1")-0.5,
    max= which(index$dataset=="C1")+0.5
  )
  shading_C2<-data.frame(
    min = c(which(index$dataset=="C2")-0.5),
    max=c(which(index$dataset=="C2")+0.5)
  )
  if(min(shading_C1$min) < min(shading_C2$min)){
    shading_C1$min[1] <- -Inf
  }else{
    shading_C2$min[1] <- -Inf
  }
  if(max(shading_C1$max) > max(shading_C2$max)){
    shading_C1$max[length(shading_C1$max)] <- Inf
  }else{
    shading_C2$max[length(shading_C2$max)] <- Inf
  }
  ###
  p1<-ggplot(df, aes(x = interaction(inner,out), y = interaction(interaction_name_2,pathway_name), 
                     color = prob, size = pval)) + geom_point(pch = 16,size=4.5) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5), axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")+
    scale_radius(range = c(min(df$pval), max(df$pval)), 
                 breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                          sort(unique(df$pval))], name = "p-value")+
    scale_colour_gradientn(colors = colorRampPalette(color.use[10:1])(99), 
                           na.value = "white", limits = c(quantile(df$prob,0, na.rm = T),
                                                          quantile(df$prob, 1, na.rm = T)), 
                           breaks = c(quantile(df$prob, 0, na.rm = T), 
                                      quantile(df$prob,1, na.rm = T)), labels = c(round(min(df$prob),3), round(max(df$prob),3))) + 
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))+
    guides(x="axis_nested",y="axis_nested")+
    geom_vline(xintercept = seq(1.5, length(unique(df[,c("out","inner")])[,1]) - 
                                  0.5, 1), lwd = 0.1, colour = "grey95") +
    geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                  0.5, 1), lwd = 0.1, colour = "grey95")+
    
    geom_vline(xintercept = dashed, linetype = "dashed",color = "black", size = 0.3)+
    theme(axis.text.x = element_text(colour = color.text[color_index]),
          ggh4x.axis.nestline.x = element_line(size = 0.8),
          ggh4x.axis.nestline.y = element_line(size = 0.8),
          ggh4x.axis.nesttext.x = element_text(angle = 0,hjust = 0.5,colour = color.text[2:1][color_index][match(unique(index$out),index$out)] ),
    )+
    geom_rect(data = shading_C1, inherit.aes = FALSE,
              aes(xmin = min, xmax = max, ymin = -Inf, 
                  ymax = Inf),fill = "#E64B35FF",alpha = 0.15,show.legend=F)+
    geom_rect(data = shading_C2, inherit.aes = FALSE,
              aes(xmin = min, xmax = max, ymin = -Inf, 
                  ymax = Inf),fill="#4DBBD5FF",alpha = 0.15,show.legend=F)
  
  return(p1)
}
##
plot_point<-function(cellchat,pairLR.use=NULL,sources.use,targets.use,comparison=c(2,1),
                     name="high",max.dataset=NULL,net.up=NULL,DEG_up=FALSE,prob=FALSE,
                     color.text=c("#CD1818","#0F2C67")){
  if(DEG_up){
    comparison=c(1,2)
  }
  if(length(comparison)<=1){
    color.text=color.text
  }
  else{
    color.text=c("#CD1818","#0F2C67")
  }
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = sources.use, 
                          targets.use = targets.use, comparison = comparison,  
                          angle.x = 45, remove.isolate = T,
                          title.name = name,thresh=0.01,
                          max.dataset = max.dataset,
                          color.text=color.text)
  if(prob){
    if(max.dataset==1){
      b<-gg2$data
      b1<-b[which(b$dataset=="C1"),]
    }
    if(max.dataset==2){
      b<-gg2$data
      b1<-b[which(b$dataset=="C2"),]
    }
  }
  else{
    b<-gg2$data
    b1<-left_join(b,net.up[,c(1:4,13)],by=c("source","target","ligand","receptor"))
    b1<-na.omit(b1)
    #b1$source.target<-as.factor(as.character(b1$source.target))
  }
  gg2$data<-b1
  if(DEG_up){
    out<-plot_point(cellchat,pairLR.use,sources.use,
                    targets.use,comparison=max.dataset,name,max.dataset,net.up,DEG_up = F)
    b2<-left_join(out$data,b1[,c(1:4,16)],by=c("source","target","ligand","receptor"))
    b2<-na.omit(b2)
    b1<-b2
    gg2$data<-b1
  }
  color<-gg2[["plot_env"]][["xtick.color"]][match(b1$dataset[match(levels(as.factor(as.character(b1$source.target))),b1$source.target)],
                                                  names(gg2[["plot_env"]][["xtick.color"]]))]
  p1<-gg2 + theme(axis.text.x = element_text(colour = color))
  return(p1)
}


ES_sur_best <- function(ES,data,name,cell,cluster){
  a<-ES[[cell]][,cluster,drop=F]
  sur<-data
  sub1<-data.frame(a,sur)
  colnames(sub1)[1]<-"gene"
  res.cut1 <- surv_cutpoint(sub1, time = "OS.time", event = "OS",
                            variables = c("gene"))
  cut1<-as.numeric(summary(res.cut1)[1])
  low<-which(a<cut1)
  high<-which(a>=cut1)
  #return(low)
  subtype<-data.frame()
  subtype[low,1]<-"low"
  subtype[high,1]<-"high"
  colnames(subtype)<-"subtype"
  sur<-data[,2:3]
  sub_sur<-data.frame(subtype,sur)
  return(sub_sur)
}
####