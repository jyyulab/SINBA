#' @import Biobase limma igraph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot
#' @importFrom reshape melt
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf rgb
#' @importFrom graphics abline arrows axis barplot boxplot hist image layout legend lines mtext par points polygon rect segments strheight stripchart strwidth text
#' @importFrom utils read.delim write.table
#' @importFrom ROCR prediction performance
#' @importFrom  httr POST content
#' @importFrom jsonlite fromJSON
#' @importFrom pbapply pblapply
#' @importFrom ggpubr ggboxplot
#' @importFrom dplyr left_join
#' @importFrom ghql GraphqlClient
#' @importFrom GA ga
#' @importFrom openxlsx write.xlsx read.xlsx
#######################dependent packages##################
#library(plyr)
library(dplyr)
library(ROCR)
library(RColorBrewer)
library(httr)
library(jsonlite)
library(reshape)
library(pbapply)
library(ggpubr)
library(limma)
library(ghql)
library(GA)
library(openxlsx)
#######################Keep internal check##################
check_para <- function(para_name,envir){
  if(base::exists(para_name,envir=envir)==FALSE){message(sprintf('%s missing !',para_name));return(0)}
  if(is.null(base::get(para_name,envir=envir))==TRUE){message(sprintf('%s is NULL !',para_name));return(0)}
  return(1)
}
check_option <- function(para_name,option_list,envir){
  if(!base::get(para_name,envir=envir) %in% option_list){
    message(sprintf('Only accept %s set at: %s !',para_name,base::paste(option_list,collapse=';')));return(0)
  }
  return(1)
}
clean_charVector <- function(x){
  x1 <- names(x)
  x <- as.character(x);
  x[which(x=='')] <- 'NULL';
  x[which(is.null(x)==TRUE)] <- 'NULL'
  x[which(is.na(x)==TRUE)] <- 'NA'
  names(x) <- x1
  x
}
check_SINBA.par <- function(SINBA.par=NULL,step='pre-load'){
  if(class(SINBA.par)!='list'){message('Invalid SINBA.par !');return(FALSE)}
  if(step=='pre-load'){
    n1 <- names(SINBA.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.PLOT','tf.network.file','sig.network.file'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s SINBA.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  if(step %in% c('exp-load','exp-QC','get_AC','combo-AC','ms-tab')){
    n1 <- names(SINBA.par)
    n2 <- setdiff(c('main.dir','project.name','out.dir','out.dir.QC','out.dir.DATA','out.dir.PLOT','cal.eset','tf.network.file','sig.network.file'),n1)
    if(length(n2)>0){
      message(sprintf('Miss %s SINBA.par, please check and re-try !',paste(n2,collapse=';')))
      return(FALSE)
    }
  }
  return(TRUE)
}
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ', '),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}
###########################################################
####helper function####
#' @title get_net2target_list
#' @param net_dat data.frame, must contain four columns with column names "source" (driver) and "target" (target genes), "MI" (mutual information) and "spearman" (spearman correlation coefficient).
#' @return Return a list. The names of the list elements are drivers.
#' Each element is a data frame, contains three columns. "target", target gene names;"MI", mutual information; "spearman", spearman correlation coefficient.
#' @noRd
#' @examples
#' out_net_file<-"path_to_network_file/consensus_network_ncol_.txt"
#' net_dat      <- read.delim(file=out_net_file,stringsAsFactors = FALSE)
#' target_list  <- get_net2target_list(net_dat)
get_net2target_list <- function(net_dat=NULL) {
  all_input_para <- c('net_dat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  all_source <- base::unique(net_dat$source)
  all_target <- lapply(all_source, function(x) {
    n1 <- net_dat[which(net_dat$source == x), base::intersect(c('target', 'MI', 'spearman'),colnames(net_dat))]
    if(class(n1)=='character') n1 <- data.frame('target'=n1,'MI'=1,'spearman'=1,stringsAsFactors=F)
    n1 <- unique(n1)
    if(length(unique(n1$target))!=length(n1$target)){ ## multiple
      t1 <- table(n1$target)
      w1 <- names(which(t1==1)); w2 <- names(which(t1>1))
      w21 <- n1[which(n1$target %in% w1),]
      w22 <- do.call(rbind,lapply(w2,function(x){
        x1 <- n1[which(n1$target==x),]
        x1 <- x1[which.max(x1$MI),]
      }))
      n1 <- rbind(w21,w22)
    }
    rownames(n1) <- n1$target
    return(n1)
  })
  names(all_target) <- all_source
  return(all_target)
}
#' @title combine_names
#' @param x seed gene name which is included in target list
#' @param target_list \code{target_list} is a list contain all the regulons
#' @param partner_name partner genes for combination
#' @return Return a comination name list. seed and partner genes are seperated by colon.
#' @export
#'
combine_names<-function(x,target_list,partner_name){
  if(is.null(partner_name)){
    partner_name<-setdiff(names(target_list),x)
  }else{
    partner_name<-setdiff(partner_name,x)
  }
  y<-partner_name
  name<-paste(x,y,sep = ":")
  return(name)}
#'
#' @title l2l_ftest
#' @param gsc1 a list of genesets, each geneset contains multiple gene names.
#' @param gsc2 a list of genesets, each geneset contains multiple gene names.
#' @param N mumeric, total unique gene number in gsc1 and gsc2.
#' @export
#' @return Return a table of all fisher exact testing results.
#'
l2l_ftest<-function(gsc1,gsc2,N){
  l_ftest<-function(gsc,l,N){
    rs<-data.frame(signature=paste(names(gsc)),n.overlap=NA,n1=NA,n2=NA,N=N,pval.FET=NA,OR.FET=NA,o.genes=NA)
    for(i in 1:nrow(rs)){
      N1<-unique(unlist(l))
      o.genes<-intersect(gsc[[i]],N1)
      rs$o.genes[i]<-paste(o.genes,collapse = ";")
      rs$n.overlap[i]<-length(o.genes)
      rs$n1[i]<-length(N1)
      rs$n2[i]<-length(gsc[[i]])
      rs$pval.FET[i]<-phyper(rs$n.overlap[i],rs$n1[i],rs$N[i]-rs$n1[i],rs$n2[i],lower.tail=FALSE)
      rs$OR.FET[i]=length(o.genes)*N/length(N1)/length(gsc[[i]])
    }
    rs$FDR.BH.FET<-p.adjust(rs$pval.FET,method='BH')
    return(rs)
  }
  #all.rs<-data.frame(matrix(ncol = 9,nrow = 0))
  #colnames(all.rs)<-c("signature","n.overlap","n1","n2","N","pval.FET","FDR.BH.FET","OR.FET","o.genes")
  all.rs<-sapply(gsc1,FUN = l_ftest,gsc=gsc2,N=N,simplify = F)
  return(all.rs)
}
#'
#' @concept activity calculation
#'
es <- function(z, es.method = "mean") {
  if (es.method == "maxmean") {
    n <- base::length(z)
    m1 <- ifelse(sum(z > 0) > 0, sum(z[z > 0]) / n, 0)
    m2 <- ifelse(sum(z < 0) > 0, sum(z[z < 0]) / n, 0)
    if (m1 > -m2)
      es <- m1
    else
      es <- m2
  }
  else if (es.method == 'absmean') {
    es <- base::mean(abs(z),na.rm=TRUE)
  }
  else if (es.method == 'mean') {
    es <- base::mean(z,na.rm=TRUE)
  }
  else if (es.method == 'median') {
    es <- stats::median(z,na.rm=TRUE)
  }
  else if (es.method == 'max') {
    es <- base::max(z,na.rm=TRUE)
  }
  else if (es.method == 'min') {
    es <- base::min(z,na.rm=TRUE)
  }
  return(es)
}
#'
#' @concept activity calculation
#'
do.std <- function(x) {
  x <- x[!is.na(x)]
  (x - base::mean(x,na.rm=TRUE)) / sd(x,na.rm=TRUE)
}
#'
#' Calculate the mean of logged values
#' Calculate mean of logged values in non-log space (return answer in log-space)
#' @param x A vector of values
#' @param ... Other arguments (not used)
#' @return Returns the mean in log-space
#' @concept helper
#' @examples
#' ExpMean(x = c(1, 2, 3))
#' @concept DE
#' @noRd
ExpMean <- function(x, ...) {
  if (inherits(x = x, what = 'AnyMatrix')) {
    return(apply(X = x, FUN = function(i) {log(x = mean(x = exp(x = i) - 1) + 1)}, MARGIN = 1))
  } else {
    return(log(x = mean(x = exp(x = x) - 1) + 1))
  }
}
#'
# internal function to calculate AUC values
#' @importFrom pbapply pblapply
#' @concept deltaROC
#' @noRd
AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
  myAUC <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(DifferentialAUC(
        x = as.numeric(x = data1[x, ]),
        y = as.numeric(x = data2[x, ])
      ))
    }
  ))
  myAUC[is.na(x = myAUC)] <- 0
  iterate.fxn <- ifelse(test = print.bar, yes = pblapply, no = lapply)
  avg_diff <- unlist(x = iterate.fxn(
    X = mygenes,
    FUN = function(x) {
      return(
        base::mean(
          x = as.numeric(x = data1[x, ])
        ) - base::mean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}
#'
# internal function to calculate AUC values
#' @importFrom ROCR prediction performance
#' @concept deltaROC
#' @noRd
DifferentialAUC <- function(x, y) {
  prediction.use <- ROCR::prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- ROCR::performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}
#'
#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using
#' the Student's t-test
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats t.test
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#'
#' @examples
#' cells.1<-colnames(sinba_ac_mat)[1:20]
#' cells.2<-colnames(sinba_ac_mat)[21:40]
#' data.use<-sinba_ac_mat
#' DiffTTest(data.use = data.use,cells.1 = cells.1,cells.2 = cells.2)
#' @keywords internal
#' @noRd
DiffTTest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose,
  ...
) {
  my.sapply <- ifelse(
    test = verbose && future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  p_val <- unlist(
    x = my.sapply(
      X = 1:nrow(data.use),
      FUN = function(x) {
        t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(p_val,row.names = rownames(x = data.use))
  return(to.return)
}
#'
#' Differential expression using Wilcoxon Rank Sum
#'
#' Identifies differentially expressed genes between two groups of cells using a Wilcoxon Rank Sum test. Makes use of limma::rankSumTestWithCorrelation for a more efficient implementation of the wilcoxon test.
#' @param data.use Data matrix to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param verbose Print a progress bar
#' @param ... Extra parameters passed to wilcox.test
#' @return Returns a p-value ranked matrix of putative differentially expressed features.
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @concept DE
#' @noRd
#' @examples
#' cells.1<-colnames(sinba_ac_mat)[1:20]
#' cells.2<-colnames(sinba_ac_mat)[21:40]
#' data.use<-sinba_ac_mat
#' WilcoxDETest(data.use = data.use,cells.1 = cells.1,cells.2 = cells.2)
#'
WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose,
  ...
) {
  my.sapply <- ifelse(
    test = verbose && future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(x = group.info), drop = FALSE]
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
    }
  )
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}
#'
#' Differential expression using limma
#'
#' Identifies differentially expressed genes between two groups of cells using limma method.
#' @param data.use Data matrix to test
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param verbose Print a progress bar
#' @param ... Extra parameters passed to limma
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' features
#'
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#'
#' @examples
#' cells.1<-colnames(sinba_ac_mat)[1:20]
#' cells.2<-colnames(sinba_ac_mat)[21:40]
#' data.use<-sinba_ac_mat
#' LimmaDETest(data.use = data.use,cells.1 = cells.1,cells.2 = cells.2)
#' @concept DE
#' @noRd
#'
LimmaDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose,
  random_effect,
  ...
){
  my.sapply <- ifelse(
    test = verbose && future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  use_samples<-c(cells.1,cells.2)
  data.use <- data.use[, use_samples, drop = FALSE]
  design.mat <-as.data.frame(matrix(NA, nrow = base::length(use_samples), ncol = 1))
  rownames(design.mat) <-use_samples
  colnames(design.mat) <- 'group'
  design.mat[base::intersect(cells.1, use_samples), 'group'] <- 'G1'
  design.mat[base::intersect(cells.2, use_samples), 'group'] <- 'G0'
  #  design <- model.matrix( ~ group + 0, design.mat)
  group <- factor(design.mat$group)
  design <- model.matrix(~0+group);
  colnames(design) <- levels(group);
  rownames(design) <- colnames(data.use)

  if(is.null(random_effect)==TRUE){
    fit <- limma::lmFit(data.use,design)
  }else{
    random_effect <- random_effect[colnames(data.use)]
    corfit <- limma::duplicateCorrelation(data.use,design,block=random_effect)
    fit <- limma::lmFit(data.use,design,block=random_effect,correlation=corfit$consensus)
  }
  contrasts <- limma::makeContrasts(G1-G0,levels=design)
  fit2 <- limma::contrasts.fit(fit,contrasts=contrasts)
  fit2 <- limma::eBayes(fit2,trend=TRUE)
  #summary(decideTests(fit2, method="global"))
  ##
  tT <- limma::topTable(fit2, adjust.method = "fdr", number = Inf,coef=1)
  if(nrow(tT)==1){
    rownames(tT) <- rownames(data.use)
  }
  tT <- base::cbind(ID=rownames(tT),tT,stringsAsFactors=FALSE)
  tT <- tT[rownames(data.use),,drop=FALSE]

  exp_G1 <- base::rowMeans(data.use[,cells.1,drop=FALSE]);
  exp_G0 <- base::rowMeans(data.use[,cells.2,drop=FALSE]);
  w1 <- which(tT$P.Value<=0);
  if(base::length(w1)>0) tT$P.Value[w1] <- .Machine$double.xmin;
  logfc_sign<-ifelse(tT$logFC==0,1,sign(tT$logFC))#JJ modify
  z_val <- sapply(tT$P.Value*logfc_sign,function(x)combinePvalVector(x,twosided = TRUE)[1])#JJ modify
  #z_val <- sapply(tT$P.Value*sign(tT$logFC),function(x)combinePvalVector(x,twosided = TRUE)[1]) #original
  if(is.null(random_effect)==TRUE){
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1)
  }else{
    tT <- base::cbind(tT,'Z-statistics'=z_val,'Ave.G0'=exp_G0,'Ave.G1'=exp_G1,
                      'Ave.G0_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G0'],
                      'Ave.G1_RemoveRandomEffect'=fit@.Data[[1]][rownames(tT),'G1'])
  }
  #if(is.null(G0_name)==FALSE) colnames(tT) <- gsub('Ave.G0',paste0('Ave.',G0_name),colnames(tT))
  #if(is.null(G1_name)==FALSE) colnames(tT) <- gsub('Ave.G1',paste0('Ave.',G1_name),colnames(tT))
  tT <- tT[order(tT$P.Value, decreasing = FALSE), ]

  return(tT)
}
#'
#' @concept DE
#' @noRd
#'
PerformDE <- function(
  object,
  cells.1,
  cells.2,
  features,
  test.use,
  verbose,
  random_effect,
  ...
) {
  if (!test.use %in% c("wilcox","t.test","limma")) {
    stop("please select one of the methods: 'wilcox','t.test','limma'")
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = WilcoxDETest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    't.test' = DiffTTest(
      data.use = object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    'limma' = LimmaDETest(
      data.use=object[features, c(cells.1, cells.2), drop = FALSE],
      cells.1=cells.1,
      cells.2=cells.2,
      verbose = TRUE,
      random_effect=random_effect,
      ...
    ),
    stop("Unknown test: ", test.use)
  )
  return(de.results)
}
#1.
#' @title get_label_manual
#' @concept SINBA_GSEA
#' @noRd
get_label_manual <- function(x){
  x1 <- sapply(x,function(x2){
    x3 <- unlist(strsplit(as.character(x2),""))
    x4 <- length(x3)%/%3 ## add number
    if(x4>0){
      pp <- length(x3)-seq(1,x4)*3; x3[pp] <- paste0(x3[pp],','); paste(x3,collapse="")
    }else{
      x2
    }
  })
  unlist(x1)
}
#2.
#' @title z2col
#' @concept SINBA_GSEA
#' @noRd
z2col <- function(x,n_len=60,sig_thre=0.01,col_min_thre=0.01,col_max_thre=3,
                  blue_col=RColorBrewer::brewer.pal(9,'Set1')[2],
                  red_col=RColorBrewer::brewer.pal(9,'Set1')[1]){
  ## create vector for z-score, can change sig threshold
  x[which(is.na(x)==TRUE)] <- 0
  x[which(x==Inf)]<-  max(x[which(x!=Inf)])+1
  x[which(x==-Inf)]<- min(x[which(x!=-Inf)])-1
  if(col_min_thre<0) col_min_thre<-0.01
  if(col_max_thre<0) col_max_thre<-3
  #c1 <- brewer.pal(9,'Set1')
  c2 <- colorRampPalette(c(blue_col,'white',red_col))(n_len)
  r1 <- 1.05*max(abs(x)) ## -r1~r1
  if(r1 < col_max_thre){
    r1 <- col_max_thre
  }
  if(col_min_thre>r1){
    r2 <- seq(-r1,r1,length.out=n_len+1)
  }else{
    r21 <- seq(-r1,-col_min_thre,length.out=n_len/2)
    r22 <- seq(col_min_thre,r1,length.out=n_len/2)
    r2 <- c(r21,r22)
  }
  x1 <- cut(x,r2)
  names(c2) <- levels(x1)
  x2 <- c2[x1]
  x2[which(abs(x)<sig_thre)] <- 'white'
  x2
}
#3.
#' @title get_transparent
#' @concept SINBA_GSEA
#' @noRd
#'
get_transparent <- function(x,alpha=0.1){
  rgb(t(col2rgb(x)/255),alpha=alpha)
}
#4.
#' @title get_z2p
#' @concept SINBA_GSEA
#' @noRd
#'
get_z2p <- function(x,use_star=FALSE){
  x[which(is.na(x)==TRUE)] <- 0
  #if(is.na(x[1])==TRUE) return('NA')
  x <- abs(x)
  x[which(is.na(x)==TRUE)] <- 0 ##
  if(max(x)<5){
    use_pv <- 1-pnorm(x)
    use_pv<-use_pv*2 #JJ modified
    use_p <- format(use_pv,digits=2,scientific = TRUE)
  }else{
    low_p <- .Machine$double.xmin
    low_z <- sapply(10^(-(1:(1+-log10(low_p)))),combinePvalVector) #function combinePvalVector
    use_pv <- sapply(x,function(x1){
      low_z[2,which(low_z[1,]>=x1)[1]]}
    )
    use_pv<-use_pv*2
    use_p <- format(use_pv, digits=3,scientific = TRUE)
    use_p[which(use_p=='NA')] <- '<1e-308'
    use_p <- as.character(use_p)
  }
  x_star <- rep('',length.out=length(use_pv))
  x_star[which(use_pv<0.05)] <-'*'
  x_star[which(use_pv<0.01)] <-'**'
  x_star[which(use_pv<0.001)] <-'***'
  if(use_star==TRUE) use_p<-paste0(use_p,x_star)
  return(use_p)
}
#5.
#' @title Combine P Values Using Fisher's Method or Stouffer's Method
#' @description \code{combinePvalVector} is a function to combine multiple comparison's P values using Fisher's method or Stouffer's method.
#'
#' @param pvals a vector of numerics, the P values from multiple comparison need to be combined.
#' @param method character, users can choose between "Stouffer" and "Fisher". Default is "Stouffer".
#' @param signed logical, if TRUE, will give a sign to the P value to indicate the direction of testing.
#' Default is TRUE.
#' @param twosided logical, if TRUE, P value is calculated in a one-tailed test.
#' If FALSE, P value is calculated in a two-tailed test, and it falls within the range 0 to 0.5.
#' Default is TRUE.
#' @return Return a vector contains the "Z-statistics" and "P.Value".
#' @examples
#' combinePvalVector(c(0.1,1e-3,1e-5))
#' combinePvalVector(c(0.1,1e-3,-1e-5))
#' @concept SINBA_GSEA
#' @noRd
#'
combinePvalVector <-
  function(pvals,
           method = 'Stouffer',
           signed = TRUE,
           twosided = TRUE) {
    #
    all_input_para <- c('pvals','method','signed','twosided')
    check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
    if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
    check_res <- c(check_option('signed',c(TRUE,FALSE),envir=environment()),
                   check_option('twosided',c(TRUE,FALSE),envir=environment()),
                   check_option('method',c('Stouffer','Fisher'),envir=environment()))
    if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
    #
    #remove NA pvalues
    pvals <- pvals[!is.na(pvals) & !is.null(pvals)]
    pvals[which(abs(pvals)<=0)] <- .Machine$double.xmin
    if (sum(is.na(pvals)) >= 1) {
      stat <- NA
      pval <- NA
    } else{
      if (twosided & (sum(pvals > 1 | pvals < -1) >= 1))
        stop('pvalues must between 0 and 1!\n')
      if (!twosided & (sum(pvals > 0.5 | pvals < -0.5) >= 1))
        stop('One-sided pvalues must between 0 and 0.5!\n')

      if (!signed) {
        pvals <- abs(pvals)
      }

      signs <- sign(pvals)
      signs[signs == 0] <- 1

      if (grepl('Fisher', method, ignore.case = TRUE)) {
        if (twosided & signed) {
          neg.pvals <- pos.pvals <- abs(pvals) / 2
          pos.pvals[signs < 0] <- 1 - pos.pvals[signs < 0]
          neg.pvals[signs > 0] <- 1 - neg.pvals[signs > 0]
        } else{
          neg.pvals <- pos.pvals <- abs(pvals)
        }
        pvals <-
          c(1, -1) * c(
            pchisq(
              -2 * sum(log(as.numeric(pos.pvals))),
              df = 2 * base::length(pvals),
              lower.tail = FALSE
            ) / 2,
            pchisq(
              -2 * sum(log(as.numeric(neg.pvals))),
              df = 2 * base::length(pvals),
              lower.tail = FALSE
            ) / 2
          )
        pval <- base::min(abs(pvals))[1]
        #if two pvals are equal, pick up the first one
        stat <-
          sign(pvals[abs(pvals) == pval])[1] * qnorm(pval, lower.tail = F)[1]
        pval <- 2 * pval
      }
      else if (grepl('Stou', method, ignore.case = TRUE)) {
        if (twosided) {
          zs <- signs * qnorm(abs(pvals) / 2, lower.tail = FALSE)
          stat <- sum(zs) / sqrt(base::length(zs))
          pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
        }
        else{
          zs <- signs * qnorm(abs(pvals), lower.tail = FALSE)
          stat <- sum(zs) / sqrt(base::length(zs))
          pval <- pnorm(abs(stat), lower.tail = FALSE)
        }
      }
      else{
        stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
      }
    }
    return(c(`Z-statistics` = stat, `P.Value` = pval))
  }
#6.
#' @title get_ES
#' @description  get enrichment score
#' @concept SINBA_GSEA
#' @noRd
#'
get_ES <- function(rank_profile=NULL,use_genes=NULL,weighted.score.type=1){
  gene.list <- names(rank_profile)
  correl.vector <- rank_profile
  tag.indicator <- sign(match(gene.list, use_genes, nomatch=0))# notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(use_genes)
  Nm <-  N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (isTRUE(max.ES > - min.ES)){
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
#######################SINBA prepare###############################
#'@title SINBA.analysis.dir.create
#'
#'@description  \code{SINBA.analysis.dir.create} is used to help users create an organized working directory for SINBA analysis.However, it is not essential for the analysis.It creates a hierarchcial working directory and returns a list contains this directory information.
#'
#' This function needs users to define the main working directory and the project's name.
#' It creates a main working directory with a subdirectory of the project.
#' It also automatically creates three subfolders (QC, DATA and PLOT) within the project folder.
#'
#' QC/, storing Quality Control related plots; DATA/, saving data in RData format; PLOT/, storing the plot files generated during the analysis
#'
#' @param project_main_dir character, name or absolute path of the main working directory.
#' @param project_name character, name of the project folder.
#' @param tf.network.file character, the path of the TF network file.
#' @param sig.network.file character, the path of the SIG network file.
#'
#' @examples
#'
#' \dontrun{
#' # Creating a main working directory under the current working directory by folder name
#' SINBA.par <- SINBA.analysis.dir.create("MyMainDir","MyProject")
#' }
#' @export
SINBA.analysis.dir.create <- function(project_main_dir=NULL,
                                      project_name=NULL,
                                      tf.network.file=NULL,
                                      sig.network.file=NULL){
  if(exists('SINBA.par')==TRUE){message('SINBA.par is occupied in the current session,please manually run: rm(SINBA.par) and re-try, otherwise will not change !');
    return(SINBA.par)}
  if(is.null(project_main_dir)==TRUE){message('project_main_dir required, please input and re-try!');return(FALSE)}
  if(is.null(project_name)==TRUE){message('project_name required, please input and re-try!');return(FALSE)}

  SINBA.par <- list()
  SINBA.par$main.dir <- project_main_dir
  SINBA.par$project.name <- project_name
  SINBA.par$out.dir <- sprintf('%s/%s/',SINBA.par$main.dir,SINBA.par$project.name)
  # create output directory
  if (!dir.exists(SINBA.par$out.dir)) {
    dir.create(SINBA.par$out.dir, recursive = TRUE)
  }
  SINBA.par$out.dir.QC <- paste0(SINBA.par$out.dir, '/QC/')
  if (!dir.exists(SINBA.par$out.dir.QC)) {
    dir.create(SINBA.par$out.dir.QC, recursive = TRUE) ## directory for QC
  }
  SINBA.par$out.dir.DATA <- paste0(SINBA.par$out.dir, '/DATA/')
  if (!dir.exists(SINBA.par$out.dir.DATA)) {
    dir.create(SINBA.par$out.dir.DATA, recursive = TRUE) ## directory for DATA
  }
  SINBA.par$out.dir.PLOT <- paste0(SINBA.par$out.dir, '/PLOT/')
  if (!dir.exists(SINBA.par$out.dir.PLOT)) {
    dir.create(SINBA.par$out.dir.PLOT, recursive = TRUE) ## directory for Result Plots
  }
  SINBA.par$tf.network.file <- ''
  SINBA.par$sig.network.file <- ''
  if(is.null(tf.network.file)==FALSE){
    SINBA.par$tf.network.file  <- tf.network.file
  }else{
    message('tf.network.file is required, please input and re-try!');return(FALSE)
  }
  if(is.null(tf.network.file)==FALSE){
    SINBA.par$sig.network.file <- sig.network.file
  }else{
    message('sig.network.file is required, please input and re-try!');return(FALSE)
  }
  if(file.exists(SINBA.par$tf.network.file)){
    message(sprintf('TF network file found in %s',SINBA.par$tf.network.file))
  }else{
    message(sprintf('TF network file not found in %s, please check and re-try !',SINBA.par$tf.network.file))
    return(FALSE)
  }
  if(file.exists(SINBA.par$sig.network.file)){
    message(sprintf('SIG network file found in %s',SINBA.par$sig.network.file))
  }else{
    message(sprintf('SIG network file not found in %s, please check and re-try ',SINBA.par$sig.network.file))
    return(FALSE)
  }
  message(sprintf('Analysis space created, please check %s',SINBA.par$out.dir))
  return(SINBA.par)
}
#' @title SINBA.saveRData
#'
#'@description  \code{SINBA.saveRData} is a function to save complicated list object generated by certain steps of SINBA pipeline (e.g. load network files, 'pre-load', load expression files, 'exp-load', calculate driver pairs activity as 'combo-AC', predict drug pairs table 'ms-tab').This function is not essential, but it is highly suggested for easier pipeline step checkout and reference.
#'
#' Assigning the \code{step} name to save the RData for easier reference.
#' Calling \code{SINBA.loadRData} to load the corresponding step RData, users can avoid repeating the former steps.
#'
#' @param SINBA.par list, stores all related datasets from driver selection and driver combination prediction pipeline step.
#' @param step character, name of the pipeline step decided by user for easier reference.
#'
#' @examples
#' \dontrun{
#' SINBA.par <- list()
#' SINBA.par$out.dir.DATA <- system.file('demo1','/DATA/',package = "SINBA")
#' SINBA.loadRData(SINBA.par=SINBA.par,step='pre-load')
#' SINBA.saveRData(SINBA.par=SINBA.par,step='pre-load')
#' }
#' @export
SINBA.saveRData <- function(SINBA.par=NULL,step='pre-load'){
  if(is.null(SINBA.par)==FALSE){
    check_SINBA.par(SINBA.par = SINBA.par,step=step)
    save(SINBA.par,file=sprintf('%s/SINBA.par.Step.%s.RData',SINBA.par$out.dir.DATA,step))
    message(sprintf('Successful save to %s',sprintf('%s/SINBA.par.Step.%s.RData',SINBA.par$out.dir.DATA,step)))
  }
}
#' @title SINBA.loadRData
#' @description  \code{SINBA.loadRData} is a function reloads RData saved by \code{SINBA.saveRData} function. It prevents user from repeating former pipeline steps.
#'
#' @param SINBA.par list, stores all related datasets from driver selection and driver combination prediction pipeline step.
#' @param step character, name of the pipeline step decided by user for easier reference.It should be previously assigned by user when calling \code{SINBA.saveRData} function.
#'
#' @examples
#' \dontrun{
#' SINBA.par <- list()
#' SINBA.par$out.dir.DATA <- system.file('demo1','/DATA/',package = "SINBA")
#' SINBA.loadRData(SINBA.par=SINBA.par,step='pre-load')
#' SINBA.saveRData(SINBA.par=SINBA.par,step='pre-load')
#' }
#' @export
SINBA.loadRData <- function(SINBA.par=NULL,step='pre-load'){
  if(is.null(SINBA.par)==FALSE){
    check_SINBA.par(SINBA.par = SINBA.par,step=step)
    load(file=sprintf('%s/SINBA.par.Step.%s.RData',SINBA.par$out.dir.DATA,step),.GlobalEnv)
    message(sprintf('Successful load from %s',sprintf('%s/SINBA.par.Step.%s.RData',SINBA.par$out.dir.DATA,step)))
  }
}
###########################################################
#######################SINBA analysis######################
#' @title get.single.network
#' @description  Read SJARACNe Network Result and Return it as List Object
#'
#' \code{get.single.network} reads SJARACNe network construction result and returns a list object including 1.network data frame, 2.driver-to-target list and 3.igraph object wrapped inside.
#'
#' In the demo, "consensus_network_ncol_.txt" file will be read and convert into a list object.
#' This list contains three elements, \code{network_data}, \code{target_list} and \code{igraph_obj}.
#' \code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed.
#' \code{target_list} is a driver-to-target list object. Please check details in \code{get_net2target_list}.
#' \code{igraph_obj} is an igraph object used to save this directed and weighted network.
#' Each edge of the network has two attributes, \code{weight} and \code{sign}.
#' \code{weight} is the "MI (mutual information)" value and \code{sign} is the sign of the spearman
#' correlation coefficient (1, positive regulation; -1, negative regulation).
#'
#' @param network_file character, the path for storing network file. For the output of SJAracne, the name of the network file will be "consensus_network_ncol_.txt" under the output directory.
#'
#' @return Return a list containing three elements, \code{network_dat}, \code{target_list} and \code{igraph_obj}.
#'
#' @examples
#' out_net_file<-"path_to_network_file/consensus_network_ncol_.txt"
#' network  <- get.single.network(network_file=out_net_file)
#'
#' \dontrun{
#' }
#' @export
get.single.network <- function(network_file=NULL){
  all_input_para <- c('network_file')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  net_dat      <- read.delim(file=network_file,stringsAsFactors = FALSE)
  target_list  <- get_net2target_list(net_dat)
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  if('MI' %in% colnames(net_dat)) igraph_obj   <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  if('spearman' %in% colnames(net_dat)) igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}
#' @title get.combined.network
#' @description This function will merge two networks into one network prepared for SINBA activity calculation.
#' @param network1 network object, which is generated by \code{get.single.network} contaning \code{network_data}. \code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed.\code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed and must contain four columns with column names "source" (driver) and "target" (target genes), "MI" (mutual information) and "spearman" (spearman correlation coefficient).
#' @param network2 network object, which is generated by \code{get.single.network} contaning \code{network_data}. \code{network_dat} is a data.frame, contains all the information of the network SJARACNe constructed and must contain four columns with column names "source" (driver) and "target" (target genes), "MI" (mutual information) and "spearman" (spearman correlation coefficient).
#' @examples
#' out_net_file1<-"path1_to_network_file/consensus_network_ncol_.txt"
#' out_net_file2<-"path2_to_network_file/consensus_network_ncol_.txt"
#' net1  <- get.single.network(network_file=out_net_file1)
#' net2  <- get.single.network(network_file=out_net_file2)
#' combined.network<-get.combined.network(network1=net1,network2=net2)
#' \dontrun{
#' }
#' @export
get.combined.network <- function(network1=NULL,network2=NULL){
  if(is.null(network1)){n_1<-NULL}else{
    n_1 <- network1$network_dat
    #n_1$type <-rep("net1",length(n_1$source))
  }
  if(is.null(network2)){n_2<-NULL}else{
    n_2 <- network2$network_dat
    #n_2$type <-rep("net2",length(n_2$source))
  }
  net_dat <- rbind(n_1,n_2)
  net_dat<-net_dat[order(-net_dat$MI),]
  net_dat<-net_dat[rownames(unique(net_dat[,c('source','target')])),] #unique pairs
  target_list  <- get_net2target_list(net_dat) #source get_net2target_list function from NetBID2
  igraph_obj   <- graph_from_data_frame(net_dat[,c('source','target')],directed=TRUE) ## add edge weight ???
  igraph_obj   <- set_edge_attr(igraph_obj,'weight',index=E(igraph_obj),value=net_dat[,'MI'])
  igraph_obj   <- set_edge_attr(igraph_obj,'sign',index=E(igraph_obj),value=sign(net_dat[,'spearman']))
  return(list(network_dat=net_dat,target_list=target_list,igraph_obj=igraph_obj))
}
#' @title combineNet2target_list
#' @description This function will create a combined target list of seed-partner driver pairs.
#' @param target_list list, the driver-to-target list object.\code{target_list} is generated by \code{get.combined.network} and contains all the driver gene regulons.
#' @param seed_name the driver genes will be used as seed for combination, must be included the the target_list
#' @param partner_name optional, the partner driver genes will be used as partner for combination. The default setting will use all the other drivers except the seed gene in the target_list.
#' @param return_type character,the type of return data. User can choose from "SINBA_target_list" and "SINBA_originalID". Default is "SINBA_originalID", which returns the combination originalID name for downstream activity calculation, "SINBA_target_list" returns the full combination target list.
#' @examples
#' out_net_file1<-"path1_to_network_file/consensus_network_ncol_.txt"
#' out_net_file2<-"path2_to_network_file/consensus_network_ncol_.txt"
#' net1  <- get.single.network(network_file=out_net_file1)
#' net2  <- get.single.network(network_file=out_net_file2)
#' combined.network<-get.combined.network(network1=net1,network2=net2)
#' sinba_target_list<-combineNet2target_list(target_list=combined.network$target_list,seed_name="LCK")
#' sinba_target_list<-combineNet2target_list(target_list=SINBA.par$combined.network$target_list,seed_name="LCK",partner_name=c("BCL2","FYN"),return_type="SINBA_target_list")
#' sinba_originalID<-combineNet2target_list(target_list=SINBA.par$combined.network$target_list,seed_name="LCK",partner_name=c("BCL2","FYN"),return_type="SINBA_originalID")
#' \dontrun{
#' }
#' @export
combineNet2target_list<-function(target_list=NULL,seed_name=NULL,partner_name=NULL,return_type="SINBA_originalID"){
  #check
  use_names<-names(target_list)
  use_seed_name<-intersect(use_names,seed_name); message(sprintf("Total %s genes will be used as seed,including: \n",length(use_seed_name)));message(cat(use_seed_name))

  if(!is.null(partner_name)){
    use_partner_name<-intersect(use_names,partner_name); message(sprintf("Total %s genes will be used as partner,including: \n",length(use_partner_name)));message(cat(use_partner_name))
  }else{
    use_partner_name<-use_names; message(sprintf("Total %s genes will be used as partner",(length(use_partner_name)-1)))
  }
  name_list<-lapply(use_seed_name, combine_names,target_list=target_list,partner_name=use_partner_name)
  all_source <- unique(unlist(name_list))
  if(return_type=="SINBA_originalID"){
    df<-data.frame(originalID=all_source)
    return(df)}else{
      all_target <- lapply(all_source, function(x,target_list) {
        s<-gsub("(.*)(:.*)","\\1",x)
        p<-gsub("(.*:)(.*)","\\2",x)
        n_s<-target_list[[s]]
        rownames(n_s)<-NULL
        n_p<-target_list[[p]]
        rownames(n_p)<-NULL
        l<-rbind(n_s,n_p)
        l<-l[order(-l$MI),]
        l<-l[!duplicated(l$target),]
        rownames(l)<-l$target
        return(l)
      },target_list=target_list)
      names(all_target) <- all_source
      return(all_target)
    }
}
#'
#' @title cal.SINBA.Activity
#' @description This function will Calculate Activity Value for Each Driver or seed-partner pairs.
#' @param target_list list, the driver-to-target list object.\code{target_list} is generated by \code{get.combined.network} and contains all the driver gene regulons.
#' @param SINBA_originalID character vector or the return data frame generated by \code{combineNet2target_list} with return_type="SINBA_originalID". \code{SINBA_originalID} contains the seed-partner driver pairs; seed gene name and partner gene name are seperated by colon. \code{SINBA_originalID} is optional.
#' @param cal_mat numeric matrix, the expression matrix of genes/transcripts.
#' @param es.method character, method applied to calculate the activity value. User can choose from c ("mean", "weightedmean", "maxmean", "absmean", "max" and "min"). Default is "weightedmean".
#' @param std logical, if TRUE, the expression matrix will be normalized by column. Default is TRUE.
#' @examples
#' sinba_originalID<-combineNet2target_list(target_list=SINBA.par$combined.network$target_list,seed_name="LCK",partner_name=c("BCL2","FYN","STAT3","PTCRA"),return_type="SINBA_originalID")
#' cal_mat<-Biobase::exprs(SINBA.par$cal.eset)
#' sinba_ac_mat<-cal.SINBA.Activity(SINBA_originalID=sinba_originalID,target_list=SINBA.par$combined.network$target_list,cal_mat=cal_mat)
#' single_ac_mat<-cal.SINBA.Activity(target_list=SINBA.par$combined.network$target_list,cal_mat=cal_mat)
#' \dontrun{
#' }
#' @export
cal.SINBA.Activity <- function(SINBA_originalID=NULL, target_list=NULL, cal_mat=NULL, es.method = 'weightedmean',std=TRUE) {
  ## mean, absmean, maxmean, weightedmean,max,min
  use_genes <- row.names(cal_mat)
  if(class(SINBA_originalID)=="data.frame"){SINBA_originalID<-as.character(SINBA_originalID$originalID)}
  if(base::length(use_genes)==0){
    message('No genes in the cal_mat, please check and re-try!');return(FALSE);
  }
  all_target <- target_list
  #all_target <- all_target[base::intersect(use_genes, names(all_target))] ## if the driver is not included in cal_mat but its target genes are included, will also calculate activity

  #z-normalize each sample
  if(std==TRUE) cal_mat <- apply(cal_mat, 2, do.std)

  if(is.null(SINBA_originalID)){
    ac.mat <-
      matrix(NA, ncol = ncol(cal_mat), nrow = base::length(all_target)) ## generate activity matrix, each col for sample, each row for source target
    for (i in 1:base::length(all_target)) {
      x <- names(all_target)[i]
      x1 <- all_target[[x]]
      x2 <- base::unique(base::intersect(rownames(x1), use_genes)) ## filter target by cal genes
      x1 <- x1[x2, ] ## target info
      target_num <- base::length(x2)
      if (target_num == 0)
        next
      if (target_num == 1){
        if (es.method == 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]
        if (es.method != 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]*x1$MI * sign(x1$spearman)
        next
      }
      if (es.method != 'weightedmean')
        ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE], 2, es, es.method)
      if (es.method == 'weightedmean') {
        weight <- x1$MI * sign(x1$spearman)
        ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE] * weight, 2, es, 'mean')
      }
    }
    rownames(ac.mat) <- names(all_target)
    colnames(ac.mat) <- colnames(cal_mat)
  }else{
    ac.mat <-
      matrix(NA, ncol = ncol(cal_mat), nrow = base::length(SINBA_originalID)) ## generate activity matrix, each col for sample, each row for source target
    for (i in 1:base::length(SINBA_originalID)) {
      x<-SINBA_originalID[i]
      s<-gsub("(.*)(:.*)","\\1",x)
      p<-gsub("(.*:)(.*)","\\2",x)
      n_s<-target_list[[s]]
      rownames(n_s)<-NULL
      n_p<-target_list[[p]]
      rownames(n_p)<-NULL
      x1<-rbind(n_s,n_p)
      x1<-x1[order(x1$MI,decreasing = T),]
      x1<-x1[!duplicated(x1$target),]
      rownames(x1)<-x1$target
      #x <- names(all_target)[i]
      #x1 <- all_target[[x]]
      x2 <- base::unique(base::intersect(rownames(x1), use_genes)) ## filter target by cal genes
      x1 <- x1[x2, ] ## target info
      target_num <- base::length(x2)
      if (target_num == 0)
        next
      if (target_num == 1){
        if (es.method == 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]
        if (es.method != 'weightedmean') ac.mat[i, ] <- cal_mat[x2,]*x1$MI * sign(x1$spearman)
        next
      }
      if (es.method != 'weightedmean')
        ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE], 2, es, es.method)
      if (es.method == 'weightedmean') {
        weight <- x1$MI * sign(x1$spearman)
        ac.mat[i, ] <- apply(cal_mat[x2,,drop=FALSE] * weight, 2, es, 'mean')
      }
    }
    rownames(ac.mat) <- SINBA_originalID
    colnames(ac.mat) <- colnames(cal_mat)
  }
  w1 <- which(is.na(ac.mat[,1])==FALSE)
  if(base::length(w1)==0){
    message('Fail in calculating activity, please check the ID type in cal_mat and target_list and try again !')
  }
  ac.mat <- ac.mat[w1,]
  return(ac.mat)
}
#'
#' @title seed2partners_ftest
#' @description This function will perform fisher exact test with seed and partner driver regulons.
#' SINBA_originalID,target_list,N,padjust_method=
#' @param target_list list, the driver-to-target list object. \code{target_list} is generated by \code{get.combined.network} and contains all the driver gene regulons.
#' @param SINBA_originalID character vector, the seed-partner driver pairs, seed gene name and partner gene name are seperated by colon. \code{SINBA_originalID} can be manually created or generated by \code{combineNet2target_list} with return_type="SINBA_originalID".
#' @param N numeric, total unique genes in the whole network or the target_list.
#' @param padjust_method character,the type of p value correction methods. User can choose from c ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default is "BH".
#' @return return a table of fisher exact test.
#' @examples
#' sinba_target_list<-combineNet2target_list(target_list=SINBA.par$combined.network$target_list,seed_name="LCK",partner_name=c("BCL2","FYN"),return_type="SINBA_target_list")
#' sinba_originalID<-combineNet2target_list(target_list=SINBA.par$combined.network$target_list,seed_name="LCK",partner_name=c("BCL2","FYN"),return_type="SINBA_originalID")
#' sinba_originalID<-as.character(sinba_originalID$originalID)
#' N<-length(unique(c(SINBA.par$combined.network$network_dat$source,SINBA.par$combined.network$network_dat$target)))
#'
#' f_df<-seed2partners_ftest(SINBA_originalID=sinba_originalID,target_list=SINBA.par$combined.network$target_list,N=N)
#' \dontrun{
#' }
#' @export
seed2partners_ftest<-function(SINBA_originalID=NULL,target_list=NULL,N,padjust_method="BH"){
  all_source <- unique(SINBA_originalID)
  ftest <- sapply(all_source, function(x,target_list,N) {
    s<-gsub("(.*)(:.*)","\\1",x)
    p<-gsub("(.*:)(.*)","\\2",x)
    n_s<-target_list[[s]][,"target"]
    n_p<-target_list[[p]][,"target"]
    genes<-unique(intersect(n_s,n_p))
    o.genes<-paste(genes,collapse = ",")
    n.overlap<-length(genes)
    n1<-length(n_s)
    n2<-length(n_p)
    pval.FET<-phyper(n.overlap,n1,N-n1,n2,lower.tail=FALSE)
    df<-data.frame(signature=x,n.overlap=n.overlap,Seed_targets=n1,Partner_targets=n2,N=N,pval.FET=pval.FET,o.genes=o.genes)
    df$OR.FET<-df[,"n.overlap"]*N/df[,"Seed_targets"]/df[,"Partner_targets"] #double check here change Seed_targets to Partner_targets
    return(df)},target_list=target_list,N=N,simplify = FALSE)
  df<-plyr::ldply(ftest,data.frame)
  df$FDR.BH.FET<-p.adjust(df$pval.FET,method=padjust_method)
  #df<-df[,-1]
  return(df)
}
#'
#' @title get.SINBA.DE.2G
#' @description  Differential Expression Analysis and Differential Activity Analysis Between 2 Sample Groups. \code{get.SINBA.DE.2G} is a function performs differential gene expression analysis and differential driver activity analysis between control group (parameter G0) and experimental group (parameter G1), using different methods.
#'
#' @param eset ExpressionSet class object, contains gene expression data or driver activity data.
#' @param G1 a vector of characters, the sample names of experimental group.
#' @param G0 a vecotr of characters, the sample names of control group.
#' @param G1_name character, the name of experimental group (e.g. "Male"). Default is "G1".
#' @param G0_name character, the name of control group (e.g. "Female"). Default is "G0".
#' @param verbose logical, if TRUE, sample names of both groups will be printed. Default is TRUE.
#' @param test.use character, users can choose one of the methods c("limma","t.test","wilcox"). Default is "limma"
#' @param random_effect a vector of characters, vector or factor specifying a blocking variable.Default is NULL, no random effect will be considered.
#'
#' @return
#' Return a data frame. Rows are genes/drivers, columns are "ID","Z-statistics","P.Value","adj.P.Val","logFC"
#'
#' @examples
#' phe_info <- Biobase::pData(SINBA.par$cal.eset)
#' phe_info$subgroup<-phe_info$das_sensitivity_category1
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`=="resistant")] # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`=="sensitive")] # get sample list for G1
#' ac_eset<-generate.SINBA.eset(exp_mat=ac_mat, phenotype_info=phe_info, feature_info=NULL, annotation_info="")
#' DE_gene_limma <- get.SINBA.DE.2G(eset=SINBA.par$cal.eset,G1=G1,G0=G0,G1_name="sensitive",G0_name='resistant',test.use="limma")

#' DA_gene_limma <- get.SINBA.DE.2G(eset=ac_eset,G1=G1,G0=G0,G1_name="sensitive",G0_name='resistant',test.use="limma")
#' \dontrun{
#' }
#' @export
#'
get.SINBA.DE.2G<-function(eset=NULL,G1=NULL,G0=NULL,G1_name=NULL,G0_name=NULL,verbose=TRUE,do.log2=F,test.use="limma",random_effect=NULL){
  all_input_para <- c('eset','G1','G0','verbose')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  check_res <- c(check_option('verbose',c(TRUE,FALSE),envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}

  cells.1<-G1
  cells.2<-G0
  features<-rownames(Biobase::fData(eset));
  object<-Biobase::exprs(eset)[features,c(cells.1,cells.2)];
  phe_info<-Biobase::pData(eset)[c(cells.1,cells.2),];

  comp_name<-sprintf("%s.Vs.%s",G1_name,G0_name)
  cat(comp_name,"\n")

  if(do.log2==T){
    data.1 <-log2(Matrix::rowMeans(object[features, cells.1, drop = FALSE])+1) #pseudocount 1
    data.2 <- log2(Matrix::rowMeans(object[features, cells.2, drop = FALSE])+1)
  }else{
    data.1 <-Matrix::rowMeans(object[features, cells.1, drop = FALSE])
    data.2 <- Matrix::rowMeans(object[features, cells.2, drop = FALSE])}
  total.diff <- (data.1 - data.2)

  de.results <- PerformDE(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    verbose = verbose,
    random_effect=random_effect
    #min.cells.feature = min.cells.feature,
    #latent.vars = latent.vars
  )
  diff.col<-"logFC"
  de.results[, diff.col] <- total.diff[rownames(x = de.results)]
  if(test.use!="limma"){
    logfc_sign<-ifelse(de.results$logFC==0,1,sign(de.results$logFC))#JJ modify
    de.results[,"P.Value"]<-de.results[rownames(x = de.results),"p_val"]
    z_val <- sapply(de.results$P.Value*logfc_sign,function(x)combinePvalVector(x,twosided = TRUE)[1])#JJ modify
    de.results[,"Z-statistics"]<-z_val
  }
  #de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  de.results$adj.P.Val = p.adjust(
    p = de.results$P.Value,
    method = "bonferroni",
    n = nrow(x = object)
  )
  de.results$ID<-rownames(de.results)
  #de.list<-list()
  #de.list[[comp_name]]<-de.results
  return(de.results[,c("ID","Z-statistics","P.Value","adj.P.Val","logFC")])
}
#'
#' ROC-based driver test
#' @title ROCTest
#' @description An AUC value of 1 means that expression values for this gene alone can perfectly classify the two groupings (i.e. Each of the samples in G1 exhibit a higher level than each of the samples in G2). An AUC value of 0 also means there is perfect classification, but in the other direction. A value of 0.5 implies that the gene has no predictive power to classify the two groups.
#'
#' @return Returns a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of
#' putative differentially expressed genes/driver.
#' @examples
#' pair.AUC<-list()
#' single.AUC<-list()
#' for(i in 1:dim(comps)[1]){
#'  comp_name <- sprintf('%s.Vs.%s',comps[i,1],comps[i,2]) ## each comparison must give a name !!!
#'  G0  <- rownames(phe_info)[which(phe_info$group==comps[i,2])] # get sample list for G0
#'  G1  <- rownames(phe_info)[which(phe_info$group==comps[i,1])] # get sample list for G1
#'  pair.AUC[[comp_name]]<-ROCTest(eset=SINBA.par$SINBA_AC.eset,G0=G0,G1=G1)
#'  single.AUC[[comp_name]]<-ROCTest(eset=SINBA.par$single_AC.eset,G0=G0,G1=G1)}
#' \dontrun{
#' }
#' @concept DE
#' @export
ROCTest <- function(
  eset,
  G0,
  G1,
  verbose = TRUE
) {
  data.use=Biobase::exprs(eset)
  cells.1=G1
  cells.2=G0
  to.return <- AUCMarkerTest(
    data1 = data.use[, cells.1, drop = FALSE],
    data2 = data.use[, cells.2, drop = FALSE],
    mygenes = rownames(x = data.use),
    print.bar = verbose
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  return(to.return)
}
#'
#' @title generate.AUCTable
#' @description The master table gathers all the results of \code{ROCTest}from multiple comparisons.
#' @param pair.AUC a list of results from \code{ROCTest} with seed-partner combined activities.
#' @param single.AUC a list of results from \code{ROCTest} with singe driver activities.
#' @param use_comp a vector of characters, the name of multiple comparisons. It will be used to name the columns of AUC table.
#' @examples
#' pair.AUC<-list()
#' single.AUC<-list()
#' for(i in 1:dim(comps)[1]){
#'  comp_name <- sprintf('%s.Vs.%s',comps[i,1],comps[i,2]) ## each comparison must give a name !!!
#'  G0  <- rownames(phe_info)[which(phe_info$group==comps[i,2])] # get sample list for G0
#'  G1  <- rownames(phe_info)[which(phe_info$group==comps[i,1])] # get sample list for G1
#'  pair.AUC[[comp_name]]<-ROCTest(eset=SINBA.par$SINBA_AC.eset,G0=G0,G1=G1)
#'  single.AUC[[comp_name]]<-ROCTest(eset=SINBA.par$single_AC.eset,G0=G0,G1=G1)}
#'  use_comp<-names(pair.AUC)
#' auc_tab<-generate.AUCTable(pair.AUC = pair.AUC,single.AUC = single.AUC,use_comp = use_comp)
#' \dontrun{
#' }
#' @export
#'
generate.AUCTable<-function(pair.AUC,single.AUC,use_comp){
  pair_IDs=rownames(pair.AUC[[1]])
  seed_IDs=gsub("(.*)(:.*)","\\1",pair_IDs)
  partner_IDs=gsub("(.*:)(.*)","\\2",pair_IDs)
  label_info <- data.frame('originalID'=pair_IDs,SINBA.seed=seed_IDs,SINBA.partner=partner_IDs,stringsAsFactors = F)

  combine_info <- lapply(use_comp,function(x){
    pair_IDs=rownames(pair.AUC[[x]])
    seed_IDs=gsub("(.*)(:.*)","\\1",pair_IDs)
    partner_IDs=gsub("(.*:)(.*)","\\2",pair_IDs)

    auc_pair<-pair.AUC[[x]][pair_IDs,c("myAUC","power")]
    names(auc_pair)<-paste0(c("AUC","AUCpower"),".",x,".SINBA")

    auc_seed<-single.AUC[[x]][seed_IDs,c("myAUC","power")]
    names(auc_seed)<-paste0(c("AUC","AUCpower"),".",x,".seed")

    auc_partner<-single.AUC[[x]][partner_IDs,c("myAUC","power")]
    names(auc_partner)<-paste0(c("AUC","AUCpower"),".",x,".partner")

    tmp<-cbind(auc_seed[,paste0(c("AUCpower"),".",x,".seed")],auc_partner[,paste0(c("AUCpower"),".",x,".partner")])
    single_power<-base::apply(tmp,1,function(x){return(base::max(x))})
    out<-cbind(auc_pair,auc_seed,auc_partner)
    out[,sprintf("delta_AUCpower.%s",x)]<-auc_pair[,paste0(c("AUCpower"),".",x,".SINBA")]-single_power
    rownames(out)<-pair_IDs
    return(out)
  })
  combine_info_out<- do.call(cbind,combine_info)
  delta_col<-grep("^delta_",colnames(combine_info_out),value = T)
  sinba_col<-grep(".SINBA$",colnames(combine_info_out),value = T)
  seed_col<-grep(".seed",colnames(combine_info_out),value = T)
  partner_col<-grep(".partner",colnames(combine_info_out),value = T)
  out<-cbind(label_info,combine_info_out[pair_IDs,c(delta_col,sinba_col,partner_col,seed_col)])
  return(out)
}
#'
#' @title generate.masterTable.SINBA
#' @description The master table gathers all the DE (differential expression analysis), DA (differential activity analysis of of single driver) and DA.SINBA (differential activity analysis of driver pairs) from multiple comparisons.
#'
#' @param use_comp a vector of characters, the name of multiple comparisons. It will be used to name the columns of master table.
#' @param DA.SINBA list, a list of DA comparisons of driver pairs activity, each comparison is a data.frame. The element name in the list must contain the name in \code{use_comp}.
#' @param DA list, a list of DA comparisons of single driver activity, each comparison is a data.frame. The element name in the list must contain the name in \code{use_comp}.
#' @param DE list, a list of DE comparisons, each comparison is a data.frame. The element name in the list must contain the name in \code{use_comp}.
#' @param ftest data.frame, the output of \code{seed2partners_ftest}, summarized the fisher exact test resutls of seed and partern drivers.
#' @param z_col character, name of the column in \code{DE} and \code{DA} contains the Z statistics. Default is "Z-statistics".
#' @param display_col character, name of the column in \code{DE} and \code{DA} need to be kept in the master table (e.g."P.Value","adj.P.Val" and "logFC"). Default is c("logFC","P.Value").
#' @param deltaZ_method character, name of the method to be used for calculate deltaZ value. Users can choose from "max" and "mean". Default value is "max".
#' @examples
#' use_comp<-names(SINBA.par$DA)
#' out<-generate.masterTable.SINBA(use_comp=use_comp,DA.SINBA=SINBA.par$SINBA.DA,DA=SINBA.par$DA,DE=SINBA.par$DE,ftest=SINBA.par$SINBA.ftest,z_col='Z-statistics',display_col=c('logFC','P.Value'),deltaZ_method="max")
#' \dontrun{
#' }
#' @export
#'
generate.masterTable.SINBA <- function(use_comp=NULL,DA.SINBA=NULL,DA=NULL,DE=NULL,ftest=NULL,z_col='Z-statistics',display_col=c('logFC','P.Value'),deltaZ_method="max"){
  if(is.null(use_comp)){message('No input for use_comp, please check and re-try!');return(FALSE)}
  if(is.null(DA.SINBA)){message('No input for DA.SINBA, please check and re-try!');return(FALSE)}
  if(is.null(DA)){message('No input for DA, please check and re-try!');return(FALSE)}
  if(is.null(DE)){message('No input for DE, please check and re-try!');return(FALSE)}
  if(is.null(ftest)){message('No input for ftest, please check and re-try!');return(FALSE)}
  if(is.null(z_col)){message('No input for z_col, please check and re-try!');return(FALSE)}
  if(length(setdiff(use_comp,names(DA.SINBA)))>0){message('%s not calculated, please check and re-try!',setdiff(use_comp,names(DA.SINBA)));return(FALSE)}

  funcType <- rep('SINBAnet',nrow(DA.SINBA[[1]]))
  rn <- rownames(DA.SINBA[[1]])
  #rn_label <- paste(rn,funcType,sep='_')
  #update ftest
  rownames(ftest)<-ftest$signature
  ftest<-ftest[rn,]
  SINBA_size <- ftest$Seed_targets+ftest$Partner_targets-ftest$n.overlap #modify from ftest
  label_info <- data.frame('originalID'=rn,'funcType'=funcType,'Size'=SINBA_size,stringsAsFactors=FALSE)
  seed<-gsub("(.*)(:.*)","\\1",label_info$originalID)
  partner<-gsub("(.*:)(.*)","\\2",label_info$originalID)
  label_info$SINBA.seed<-seed
  label_info$SINBA.partner<-partner
  rownames(label_info)<-label_info$originalID
  combine_info <- lapply(use_comp,function(x,deltaZ.method){
    DA_info<-DA[[x]][unique(c(seed,partner)),c(z_col,setdiff(display_col,z_col))]
    DE_info<-DE[[x]][unique(c(seed,partner)),c(z_col,setdiff(display_col,z_col))]
    DA_SINBA.info <- DA.SINBA[[x]][,c(z_col,setdiff(display_col,z_col))]

    colnames(DA_info) <- paste0(colnames(DA_info),'.',x,'_DA')
    colnames(DE_info) <- paste0(colnames(DE_info),'.',x,'_DE')
    colnames(DA_SINBA.info)<-paste0(colnames(DA_SINBA.info),'.',x,'_SINBA.DA')
    colnames(DA_info)[1] <- paste0('Z.',x,'_DA')
    colnames(DE_info)[1] <- paste0('Z.',x,'_DE')
    colnames(DA_SINBA.info)[1]<-paste0('Z.',x,'_SINBA.DA')
    DA_DE_combined<-cbind(DA_info,DE_info,stringsAsFactors=FALSE)
    DA_DE_combined$ID<-rownames(DA_DE_combined)
    DA_DE_info<-data.frame(seed_gene=seed,partner_gene=partner)
    DA_DE_info2<-dplyr::left_join(DA_DE_info,DA_DE_combined,by=c("partner_gene"="ID"))
    colnames(DA_DE_info2)[-c(1,2)]<-paste(colnames(DA_DE_info2)[-c(1,2)],"partner",sep=".")
    DA_DE_info1<-dplyr::left_join(DA_DE_info,DA_DE_combined,by=c("seed_gene"="ID"))
    colnames(DA_DE_info1)[-c(1,2)]<-paste(colnames(DA_DE_info1)[-c(1,2)],"seed",sep=".")
    out <- cbind(DA_SINBA.info,DA_DE_info1,DA_DE_info2,stringsAsFactors=FALSE)
    out<-out[,setdiff(colnames(out),c("seed_gene","partner_gene"))]
    if(deltaZ.method=="mean"){
      #DA_delta<-data.frame(da=(out[,paste0("Z.",x,"_SINBA.DA")]-(out[,paste0("Z.",x,"_DA.seed")]+out[,paste0("Z.",x,"_DA.partner")])/2))
      DA_delta<-data.frame(da=(out[,paste0("Z.",x,"_SINBA.DA")]/sqrt(2)-(out[,paste0("Z.",x,"_DA.seed")]+out[,paste0("Z.",x,"_DA.partner")])/2)) #JY suggested
    }
    if(deltaZ.method=="max"){
      sp_z<-out[,paste0("Z.",x,"_SINBA.DA")]
      s_z<-out[,paste0("Z.",x,"_DA.seed")]
      p_z<-out[,paste0("Z.",x,"_DA.partner")]
      tmp1<-data.frame(SP.Z=sp_z,S.Z=s_z,P.Z=p_z)
      tmp2<-base::sapply(tmp1[,1:3],function(x){X<-ifelse(abs(x)<1.96,0,x);return(as.numeric(X))},simplify = T)
      deltaZ<-apply(tmp2,1,function(x){tmp.val<-ifelse(abs(x[2])>=abs(x[3]),x[2],x[3])
      out.val=x[1]-tmp.val
      return(out.val)})
      DA_delta=data.frame(da=deltaZ)
    }
    out<-cbind(DA_delta,out)
    colnames(out)[1]<-paste0("delta_Z.",x,"_SINBA")
    rownames(out) <- rn
    out
  },deltaZ.method=deltaZ_method)
  combine_info_DA <- do.call(cbind,lapply(combine_info,function(x)x[rn,grep('_SINBA.DA$',colnames(x))]))
  combine_info_DE.p <- do.call(cbind,lapply(combine_info,function(x)x[rn,grep('_DE.partner$',colnames(x))]))
  combine_info_DA.p<-do.call(cbind,lapply(combine_info,function(x)x[rn,grep('_DA.partner$',colnames(x))]))
  combine_info_DE.s <- do.call(cbind,lapply(combine_info,function(x)x[rn,grep('_DE.seed$',colnames(x))]))
  combine_info_DA.s<-do.call(cbind,lapply(combine_info,function(x)x[rn,grep('_DA.seed$',colnames(x))]))
  deltaZ<-do.call(cbind,lapply(combine_info,function(x)x[rn,grep('^delta_Z.',colnames(x)),drop=F])) #one clumn data.fram
  ftest<-ftest[,-1]
  # put them together
  ms_tab <- cbind(label_info,ftest,deltaZ,combine_info_DA,combine_info_DA.p,combine_info_DA.s,combine_info_DE.p,combine_info_DE.s)
  #ms_tab<-ms_tab[-which(is.na(ms_tab$delta_Z.Sensitive.Vs.Resistant_SINBA)),] # some hub genes has no single driver activity which is removed from master table
  return(ms_tab)
}
#'
#' @title fitness_function
#' @description Fitness function to calculate the SUM of the selected submatrix
#'
#' @param solution a vector of n+m numbers.
#' @param matrix the original matrix to optimize.
#' @param n row numbers.
#' @param m column numbers.
#'
fitness_function <- function(solution, matrix, n, m) {
  # Extract rows and columns from the solution vector
  selected_rows <- solution[1:n]
  selected_cols <- solution[(n + 1):(n + m)]

  #selected_rows <- unique(selected_rows)
  #selected_cols <- unique(selected_cols)
  # Extract the submatrix
  submatrix <- matrix[selected_rows, selected_cols, drop = FALSE]
  # Return the sum of the submatrix
  return(sum(submatrix))
}
#'
#' @title find_max_sum_submatrix
#' @description Main function to find the maximum sum deltaZ score submatrix using GA to optimize the screening layout of seed and partner drivers.
#'
#' @param matrix a deltaZ score matrix for seed and driver layout optimization. It will be created by run \code{convert_deltaZ2_matrix}.
#' @param n row numbers.
#' @param m colum number
#' @param use_popSize Population size, GA::ga function parameter,default=50.
#' @param use_maxiter Maximum iterations, GA::ga function parameter, default=100.
#' @param use_run  Early stopping criterion,GA::ga function parameter, default=50.
#' @param use_suggestions Optional: initial solutions, GA::ga function parameter,default=NULL.
#' @examples
#'
#' generate_matrix2<-function(n, m) {
#'     out_mat<-matrix(sample(1:10, 10 * n * 10 * m, replace = TRUE), nrow = 10 * n, ncol = 10 * m)
#'     rownames(out_mat)<-paste0("R_",1:dim(out_mat)[1])
#'     colnames(out_mat)<-paste0("C_",1:dim(out_mat)[2])
#'     return(out_mat)
#' }
#' base_matrix<-generate_matrix2(n=5,m=10)
#' out<-find_max_sum_submatrix(matrix=base_matrix, n=4, m=3)
#' opt_mat<-out$submatrix
#' \dontrun{
#' }
#' @export
#'
find_max_sum_submatrix <- function(matrix, n, m,use_popSize = 50, use_maxiter = 100, use_run = 50, use_suggestions=NULL) {
  # Dimensions of the original matrix
  total_rows <- nrow(matrix)
  total_cols <- ncol(matrix)
  switched=FALSE
  if(total_rows>total_cols){
    switched=TRUE
    matrix<-t(matrix)
    total_rows <- nrow(matrix)
    total_cols <- ncol(matrix)
    tmp<-n
    n<-m
    m<-tmp
  }
   # Define lower and upper bounds
  lower <- c(rep(1, n + m))  # Row and column indices start from 1
  upper <- c(rep(total_rows, n), rep(total_cols, m))  # Separate bounds for rows and columns
  # Run the genetic algorithm
  result <- GA::ga(
    type = "permutation",        # Permutation-based optimization
    fitness =function(solution) fitness_function(solution, matrix=matrix, n=n, m=m),
    lower = c(1,1),                   # Lower bound for row and column indices
    upper = c(total_rows,total_cols),  # Upper bound for row and column indices
    popSize = use_popSize,                # Population size
    maxiter = use_maxiter,               # Maximum iterations
    run = use_run,                    # Early stopping criterion
    suggestions = use_suggestions,          # Optional: initial solutions
    seed = 123                   # Set seed for reproducibility
  )
  # Extract the best solution
  best_solution <- result@solution[1, ]
  selected_rows <- best_solution[1:n]
  selected_cols <- best_solution[(n + 1):(n + m)]
  # Extract the submatrix
  submatrix <- matrix[selected_rows, selected_cols, drop = FALSE]
  if(switched){
    submatrix<-t(submatrix)
  }
  # Return the submatrix and its sum
 return(list(submatrix = submatrix, max_sum = sum(submatrix)))
}
#' @title convert_deltaZ2_matrix
#' @description Convert a dataframe with selected columns(including the deltaZ score column) to a matrix.
#'
#' @param df a data frame with deltaZ scores,  it could be created by run \code{generate.masterTable.SINBA }.
#' @param use_row.col column name in the df to be used as rows in the output matrix.
#' @param use_column.col colum name in the df to be used as columns in the output matrix.
#' @param deltaZ.col column name with the deltaZ scores, the values will be used in the output matrix.
#' @param use_method the method used to aggregate the deltaZ values from same seed/partner drivers. Options: "max", "mean", "min" and "median". default="max"
#' @examples
#'
#' ms_tab<-SINBA.par$ms_tab
#' df<-convert_deltaZ2_matrix(df=ms_tab,use_row.col="SINBA.seed_geneSymbol",use_column.col="SINBA.partner_geneSymbol",deltaZ.col="delta_Z.G3.Vs.others_SINBA",use_method="max")
#' \dontrun{
#' }
#' @export
#'
convert_deltaZ2_matrix<-function(df,use_row.col,use_column.col,deltaZ.col,use_method="max"){
  if (!all(c(use_row.col,use_column.col,deltaZ.col)%in%colnames(df))) {
    stop("The input matrix must have the selected column names! ")
  }
  df<-df[,c(use_row.col,use_column.col,deltaZ.col)]
  names(df)<-c("x","y","z")
  mat_base <- with(df, tapply( z, list(x,y), FUN = use_method))
  return(as.matrix(mat_base))
}
###########################################################
#######################add-on from NETBID2#############################
#' @title generate.SINBA.eset

#' @description \code{generate.eset} generates ExpressionSet class object to contain and describe the high-throughput assays.
#' Users need to define its slots, which are expression matrix (required),
#' phenotype information and feature information (optional).
#' It is very useful when only expression matrix is available.
#'
#' @param exp_mat matrix, the expression data matrix. Each row represents a gene/transcript/probe, each column represents a sample.
#' @param phenotype_info data.frame, the phenotype information for all the samples in \code{exp_mat}.
#' In the phenotype data frame, each row represents a sample, each column represents a phenotype feature.
#' The row names must match the column names of \code{exp_mat}. If NULL, it will generate a single-column data frame.
#' Default is NULL.
#' @param feature_info data.frame, the feature information for all the genes/transcripts/probes in \code{exp_mat}.
#' In the feature data frame, each row represents a gene/transcript/probe and each column represents an annotation of the feature.
#' The row names must match the row names of \code{exp_mat}. If NULL, it will generate a single-column data frame.
#' Default is NULL.
#' @param annotation_info character, the annotation set by users for easier reference. Default is "".
#'
#' @return Return an ExressionSet object.
#'
#' @examples
#' mat1 <- matrix(rnorm(10000),nrow=1000,ncol=10)
#' colnames(mat1) <- paste0('Sample',1:ncol(mat1))
#' rownames(mat1) <- paste0('Gene',1:nrow(mat1))
#' eset <- generate.eset(exp_mat=mat1)
#'
#' \dontrun{
#' }
#' @export
generate.SINBA.eset <- function(exp_mat=NULL, phenotype_info=NULL, feature_info=NULL, annotation_info="") {
  #
  all_input_para <- c('exp_mat')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(base::min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  if(is.null(dim(exp_mat))==TRUE){
    exp_mat <- t(as.matrix(exp_mat));rownames(exp_mat) <- 'g1'
  }
  if (is.null(phenotype_info)) {
    phenotype_info <- data.frame(group = colnames(exp_mat), stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if (is.null(feature_info)) {
    feature_info <- data.frame(gene = rownames(exp_mat), stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  if((class(phenotype_info)=='character' | is.null(dim(phenotype_info))==TRUE) & is.null(names(phenotype_info))==TRUE){
    phenotype_info <- data.frame(group = phenotype_info, stringsAsFactors = FALSE)
    rownames(phenotype_info) <- colnames(exp_mat)
  }
  if((class(feature_info)=='character' | is.null(dim(feature_info))==TRUE) & is.null(names(feature_info))==TRUE){
    feature_info <- data.frame(gene = feature_info, stringsAsFactors = FALSE)
    rownames(feature_info) <- rownames(exp_mat)
  }
  eset <-
    new(
      "ExpressionSet",
      phenoData = new("AnnotatedDataFrame", phenotype_info),
      featureData = new("AnnotatedDataFrame", feature_info),
      annotation = annotation_info,
      exprs = as.matrix(exp_mat)
    )
  return(eset)
}
#'
#' Test for Intersection of Target Genes between Two Drivers
#'
#' \code{test.targetNet.overlap} performs Fisher's exact test to see whether the target genes from two drivers are significantly intersected.
#'
#'
#' @param source1_label character, the label of the first selected driver.
#' @param source2_label character, the label of the second selected driver.
#' @param target1 a vector of characters, the list of target genes from the first driver.
#' @param target2 a vector of characters, the list of target genes from the second driver.
#' @param total_possible_target numeric or a vector of characters. If input is numeric, it is the total number of possible target genes.
#' If input is a vector of characters, it is the background list of all possible target genes.
#'
#' @return Return statistics of the testing, including the \code{P.Value}, \code{Odds_Ratio} and \code{Intersected_Number}.
#'
#' @examples
#' source1_label <- 'test1'
#' target1 <- sample(paste0('G',1:1000),size=80)
#' source2_label <- 'test2'
#' target2 <- sample(paste0('G',1:1000),size=120)
#' test.targetNet.overlap(source1_label=source1_label,source2_label=source2_label,
#'                target1=target1,target2=target2,
#'                total_possible_target=paste0('G',1:1000))
#' \dontrun{
#' }
#' @importFrom NetBID2 test.targetNet.overlap
#' @noRd
#'
test.targetNet.overlap <- function(source1_label=NULL,source2_label=NULL,
                                   target1=NULL,target2=NULL,
                                   total_possible_target=NULL){
  #
  all_input_para <- c('source1_label','source2_label','target1','target2')
  check_res <- sapply(all_input_para,function(x)check_para(x,envir=environment()))
  if(min(check_res)==0){message('Please check and re-try!');return(FALSE)}
  #
  t1  <- base::unique(target1)
  t2  <- base::unique(target2)
  print(sprintf('%s has %d unique targets !',source1_label,base::length(t1)))
  print(sprintf('%s has %d unique targets !',source2_label,base::length(t2)))
  n11 <- base::length(base::intersect(t1,t2))
  n12 <- base::length(base::setdiff(t1,t2))
  n21 <- base::length(base::setdiff(t2,t1))
  if(class(total_possible_target) %in% c('integer','numeric')){
    n22 <- total_possible_target-base::length(union(t1,t2))
  }else{
    n22 <- base::length(base::setdiff(total_possible_target,c(t1,t2)))
  }
  mm  <- base::cbind(c(n11,n21),c(n12,n22))
  ft  <- fisher.test(mm)$p.value
  or  <- n11/n12/(n21/n22)
  rownames(mm) <- c(sprintf('In %s target',source1_label),sprintf('Not in %s target',source1_label))
  colnames(mm) <- c(sprintf('In %s target',source2_label),sprintf('Not in %s target',source2_label))
  print(mm)
  res <- c('P.Value'=ft,'Odds_Ratio'=or,'Intersected_Number'=mm[1,1])
  return(res)
}
###########################################################
#######################PLOT function#######################
#7.
#' @title draw.GSEA_SINBA
#' @description GSEA (Gene Set Enrichment Analysis) Plot for seed, partner and combined Seed-partner regulons.
#' @param rank_profile a named vector of numerics, the differential values (DE or DA) calculated from a sample comparison. Names of the vector must be gene names.The differential values can be "logFC", "t-statistics" or "Z-statistics".
#' @param use_genes_SINBA a vector of characters, a vector of genes to display. The genes the targe genes from seed-partner pair. Seed-partner regulon can be generated by \code{combineNet2target_list} with (return_type="SINBA_target_list")
#' @param use_genes_seed a vector of characters, a vector of genes to display. The genes the targe genes from the seed driver.
#' @param use_genes_partner a vector of characters, a vector of genes to display. The genes the targe genes from the partner driver.
#' @param use_direction a logical value, if TRUE, The target direction values need provided by user. The default value is TRUE.
#' @param use_direction_SINBA a vector of numeric 1s and -1s, 1 is positive regulation from the seed-partner pair, -1 is negative regulation from seed-partner pair. Users can get this vector by converting the signs of "spearman". If NULL, no regulation direction will be displayed. Default is NULL.
#' @param use_direction_seed a vector of numeric 1s and -1s, 1 is positive regulation from the seed, -1 is negative regulation from seed. Users can get this vector by converting the signs of "spearman". If NULL, no regulation direction will be displayed. Default is NULL.
#' @param use_direction_partner a vector of numeric 1s and -1s, 1 is positive regulation from the partner, -1 is negative regulation from partner. Users can get this vector by converting the signs of "spearman". If NULL, no regulation direction will be displayed. Default is NULL.
#' @param main character, an overall title for the plot. Default is "".
#' @param pdf_file character, the file path to save as PDF file. If NULL, no PDF file will be saved. Default is NULL.
#' @param annotation character, the annotation set by users for easier reference.
#' Normally the annotation is the P-value or other statistics to show the significance of the interested seed-partner driver pairs.If NULL, will perform a Kolmogorov-Smirnov test to get the significance value. Default is NULL.
#' @param annotation_cex numeric, giving the amount by which the text of annotation should be magnified relative to the default. Default is 1.2.
#' @param left_annotation character, annotation displayed on the left of the figure, representing left condition of the \code{rank_profile}. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure, representing right condition of the \code{rank_profile}. Default is "".
#' @return Return a logical value. If TRUE, the plot has been created successfully.
#' @examples
#' lck_fyn<-combineNet2target_list(target_list=SINBA.par$sinba.network$target_list,seed_name="LCK",partner_name="FYN",return_type="SINBA_target_list")
#' genes_SINBA<-lck_fyn[[1]]$target
#' genes_seed<-SINBA.par$sinba.network$target_list$LCK$target
#' genes_partner<-SINBA.par$sinba.network$target_list$FYN$target
#' dir_SINBA<-sign(lck_fyn[[1]]$spearman)
#' dir_seed<-sign(SINBA.par$sinba.network$target_list$LCK$spearman)
#' dir_partner<-sign(SINBA.par$sinba.network$target_list$FYN$spearman)
#' rank_profile <- SINBA.par$DE[[comp_name]]$`Z-statistics`
#' names(rank_profile)<-SINBA.par$DE[[comp_name]]$ID
#' rank_profile<-rank_profile[order(-rank_profile)]
#' draw.GSEA_SINBA(rank_profile=rank_profile,use_genes_SINBA=genes_SINBA,use_genes_seed=genes_seed,use_genes_partner=genes_partner,use_direction=TRUE, use_direction_SINBA=dir_SINBA,use_direction_seed=dir_seed,use_direction_partner=dir_partner,main="LCK_FYN",pdf_file=sprintf("%s/GSEA_LCK_FYN.pdf",plot.dir),annotation=NULL,annotation_cex=1.2,left_annotation="Sensitive",right_annotation="Resistant")
#'
#' \dontrun{
#' }
#'
#' @export
#'
draw.GSEA_SINBA <- function(rank_profile=NULL,use_genes_SINBA=NULL,use_genes_seed=NULL,use_genes_partner=NULL,use_direction=TRUE, use_direction_SINBA=NULL,use_direction_seed=NULL,use_direction_partner=NULL,main="",pdf_file=NULL,annotation=NULL,annotation_cex=1.2,left_annotation=NULL,right_annotation=NULL){
  #### start plot
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=10,height=10)
  }
  if(is.null(use_direction_SINBA)==FALSE){
    w1 <- which(use_genes_SINBA %in% names(rank_profile))
    w2 <- unique(use_direction[w1])
    if(length(w2)==1){
      if(w2==1) use_direction <- NULL;
    }
  }
  if(is.null(use_direction_SINBA)==FALSE){
    new_rank_profile <- c(rank_profile,-rank_profile)
    names(new_rank_profile) <- c(paste0('POS_',names(rank_profile)),paste0('NEG_',names(rank_profile)))
    rank_profile <- new_rank_profile
    #SINGBA
    use_genes_SINBA[which(use_direction_SINBA==1)] <- paste0('POS_',use_genes_SINBA[which(use_direction_SINBA==1)])
    use_genes_SINBA[which(use_direction_SINBA==-1)] <- paste0('NEG_',use_genes_SINBA[which(use_direction_SINBA==-1)])
    # #Seed
    use_genes_seed[which(use_direction_seed==1)] <- paste0('POS_',use_genes_seed[which(use_direction_seed==1)])
    use_genes_seed[which(use_direction_seed==-1)] <- paste0('NEG_',use_genes_seed[which(use_direction_seed==-1)])
    #partner
    use_genes_partner[which(use_direction_partner==1)] <- paste0('POS_',use_genes_partner[which(use_direction_partner==1)])
    use_genes_partner[which(use_direction_partner==-1)] <- paste0('NEG_',use_genes_partner[which(use_direction_partner==-1)])
  }
  rank_profile <- sort(rank_profile,decreasing = TRUE)
  r_len <- length(rank_profile)
  use_pos <- which(names(rank_profile) %in% use_genes_SINBA)
  use_pos_seed <- which(names(rank_profile) %in% use_genes_seed)
  use_pos_partner <- which(names(rank_profile) %in% use_genes_partner)
  layout(matrix(c(rep(6,10),5,4,3,2,rep(1,10)),ncol=1))
  # plot
  ## rank for all
  par(mar=c(10,6,0,2))
  mm <- max(abs(rank_profile))
  y1 <- seq(-mm,mm,length.out=7); y1 <- round(y1,1)
  unit <- r_len/10; unit <- round(unit/100)*100
  x1 <- seq(0,r_len,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,max(x1),-mm,mm))
  plot(rank_profile,col='grey',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(-mm,mm))
  polygon(x=c(0,1:r_len,r_len),y=c(0,rank_profile,0),col='grey',border=NA)
  if(is.null(left_annotation)==FALSE) text(0+r_len/100,mm,adj=0,left_annotation,col='red',xpd=TRUE)
  if(is.null(right_annotation)==FALSE) text(r_len-r_len/100,-mm,adj=1,right_annotation,col='blue',xpd=TRUE)
  #axis(side=2,at=y1,labels=y1,las=2);
  pp <- par()$usr
  segments(0,0,r_len,0,lwd=0.5)
  segments(0,min(y1),0,max(y1),lwd=1.5)
  text(-(pp[2]-pp[1])/50,y1,y1,adj=1,xpd=TRUE)
  segments(-(pp[2]-pp[1])/100,y1,0,y1,lwd=1.5)
  mtext(side=2,line = 2,'Ranked list metric (PreRanked)',cex=1.2)
  mtext(side=1,line = 3,'Rank in Ordered Dataset',cex=1.2)
  if(is.null(use_direction)==FALSE){
    axis(side=1,at=x1,labels=get_label_manual(x1/2))
  }else{
    axis(side=1,at=x1,labels=get_label_manual(x1))
  }
  # get zero cross
  w1 <- which.min(abs(rank_profile))
  abline(v=w1,lty=2,col='grey')
  if(is.null(use_direction)==FALSE){
    text(w1,-mm/4,sprintf('Zero cross at %s',round(w1/2)),adj=0.5)
  }else{
    text(w1,-mm/4,sprintf('Zero cross at %d',w1),adj=0.5)
  }
  if(is.null(use_direction)==FALSE){
    legend(w1,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment_SINBA','Enrichment_Seed','Enrichment_Partner','Hits_positive_direction','Hits_negative_direction','Ranking metric scores'),
           col=c('green','red','blue',pos_col,neg_col,'grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=1.0)
  }else{
    legend(w1,pp[3]-(pp[4]-pp[3])/4,lty=1,lwd=2,c('Enrichment_SINBA','Enrichment_Seed','Enrichment_Partner','Hits','Ranking metric scores'),
           col=c('green','red','blue','black','grey'),xpd=TRUE,horiz=TRUE,xjust=0.5,cex=1.0)
  }
  pm <- par()$usr
  ## get image bar
  par(mar=c(0,6,0,2))
  use_col <- z2col(rank_profile,sig_thre = 0,n_len = 30,blue_col='blue',red_col='red')
  image(x=as.matrix(1:r_len),col=use_col,bty='n',xaxt='n',yaxt='n',xlim=c(pm[1],pm[2])/r_len)
  abline(v=use_pos/r_len,col='grey')
  ## mark gene position;
  #SINBA_combined
  par(mar=c(0.5,6,0.5,2))
  plot(1,col='white',xlab="",ylab="",bty='n',xlim=c(1,r_len),xaxt='n',yaxt='n')
  if(is.null(use_direction_SINBA)==FALSE){
    use_pos_P <- which(names(rank_profile) %in% use_genes_SINBA[grep('POS',use_genes_SINBA)])
    use_pos_N <- which(names(rank_profile) %in% use_genes_SINBA[grep('NEG',use_genes_SINBA)])
    abline(v=use_pos_P,col=pos_col)
    abline(v=use_pos_N,col=neg_col)
  }else{
    abline(v=use_pos,col='green')
  }
  #Seed
  par(mar=c(0,6,0.5,2))
  plot(1,col='white',xlab="",ylab="",bty='n',xlim=c(1,r_len),xaxt='n',yaxt='n')
  if(is.null(use_direction_seed)==FALSE){
    use_pos_P <- which(names(rank_profile) %in% use_genes_seed[grep('POS',use_genes_seed)])
    use_pos_N <- which(names(rank_profile) %in% use_genes_seed[grep('NEG',use_genes_seed)])
    abline(v=use_pos_P,col=pos_col)
    abline(v=use_pos_N,col=neg_col)
  }else{
    abline(v=use_pos_seed,col='red')
  }
  #partner
  par(mar=c(0,6,0.5,2))
  plot(1,col='white',xlab="",ylab="",bty='n',xlim=c(1,r_len),xaxt='n',yaxt='n')
  if(is.null(use_direction)==FALSE){
    use_pos_P <- which(names(rank_profile) %in% use_genes_partner[grep('POS',use_genes_partner)])
    use_pos_N <- which(names(rank_profile) %in% use_genes_partner[grep('NEG',use_genes_partner)])
    abline(v=use_pos_P,col=pos_col)
    abline(v=use_pos_N,col=neg_col)
  }else{
    abline(v=use_pos_partner,col='blue')
  }
  ## GSEA ES
  par(mar=c(0.5,6,2,2))
  # get SINBA_ES score
  es_res_sinba <- get_ES(rank_profile,use_genes_SINBA)
  y2_sinba <- seq(min(es_res_sinba$RES),max(es_res_sinba$RES),length.out=7); y2_sinba <- round(y2_sinba,1)
  es_res_seed <- get_ES(rank_profile,use_genes_seed)
  y2_seed <- seq(min(es_res_seed$RES),max(es_res_seed$RES),length.out=7); y2_seed <- round(y2_seed,1)
  es_res_partner <- get_ES(rank_profile,use_genes_partner)
  y2_partner <- seq(min(es_res_partner$RES),max(es_res_partner$RES),length.out=7); y2_partner <- round(y2_partner,1)
  res<-c(es_res_sinba$RES,es_res_seed$RES,es_res_partner$RES)
  min_res<-min(res)
  max_res<-max(res)
  if(is.null(use_direction)==FALSE){
    plot(es_res_sinba$RES,col='green',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
         xlim=c(1,r_len),type='l',lwd=3,ylim=c(min_res,max_res),main=main,xpd=TRUE)
    lines(es_res_seed$RES,col='red',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
          xlim=c(1,r_len),lty=3,lwd=1,ylim=c(min(es_res_seed$RES),max(y2_seed)),main="",xpd=TRUE)
    lines(es_res_partner$RES,col='blue',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',xlim=c(1,r_len),lty=3,lwd=1,ylim=c(min(es_res_partner$RES),max(y2_partner)),main="",xpd=TRUE)
  }else{
    plot(es_res_sinba$RES,col='green',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
         xlim=c(1,r_len),type='l',lwd=3,ylim=c(min_res,max_res),main=main,xpd=TRUE)
    lines(es_res_seed$RES,col='red',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',
          xlim=c(1,r_len),lty=3,lwd=1,ylim=c(min(es_res_seed$RES),max(y2_seed)),main=main,xpd=TRUE)
    lines(es_res_partner$RES,col='blue',xaxt='n',yaxt='n',xlab="",ylab="",bty='n',xlim=c(1,r_len),lty=3,lwd=1,ylim=c(min(es_res_partner$RES),max(y2_partner)),main=main,xpd=TRUE)
  }
  pp <- par()$usr
  #abline(h=0)
  #axis(side=2,at=y2,label=y2,las=2)
  y2<-seq(min(c(y2_sinba,y2_seed,y2_partner)),max(c(y2_sinba,y2_seed,y2_partner)),length.out=7); y2<- round(y2,1)
  segments(0,0,r_len,0,lwd=1.5)
  segments(0,min(y2),0,max(y2),lwd=1.5)
  text(-(pp[2]-pp[1])/50,c(y2),c(y2),adj=1,xpd=TRUE)
  segments(-(pp[2]-pp[1])/100,c(y2),0,c(y2),lwd=1.5)
  mtext(side=2,line = 2,'Enrichment score (ES)',cex=1.2)
  # add annotation
  if(is.null(annotation)==TRUE){
    annotation <- sprintf("KS test p-value:%s",format(ks.test(rank_profile,rank_profile[use_genes_SINBA])$p.value,digits = 3,scientific = TRUE))
  }
  if(res[which.max(abs(res))]>0)
    text(r_len-r_len/50,max(y2),annotation,adj=1,cex=annotation_cex,xpd=TRUE)
  else
    text(0+r_len/50,min(y2)+(max(y2)-min(y2))/10,annotation,adj=0,cex=annotation_cex,xpd=TRUE)
  if(is.null(pdf_file)==FALSE){
    dev.off()
  }
  layout(1);
  return(TRUE)
}
#'
#' @title draw.NetBID_SINBA
#' @description GSEA (gene set enrichment analysis) plot for the Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA) analysis results.
#'
#' \code{draw.GSEA.NetBID.SINBA} will generate a GSEA plot for Synergy Inference by data-driven Network-based Bayesian Analysis (SINBA)
#' analysis results. SINBA calculates the synergistic effect between a seed driver and a partner driver.
#' The plot includes the GSEA plot for the seed driver and the partner driver independently and
#' the GSEA plot for the combination for the seed driver to each partner driver.
#' The statistics on the plot include the differentiated expression (DE), differentiated activity (DA) for each driver,
#' and the different Z (deltaZ) showing the difference between the combination of the seed and the partner driver to the sum of the original Z statistics.
#'
#' This is a plot function to draw GSEA for synergistic effect prediction between the seed driver and a list of partner drivers.
#' User need to input the differentiated expression information, and choose to display the target genes in one row or two rows, by selecting black color or red to blue color bar.
#'
#' @param DE data.frame,the differentiated expression results.
#' This data.frame could be generated by using \code{getDE.limma.2G} or \code{getDE.BID.2G}.
#' If user want to generate this data.frame by other strategies, the rownames must be the gene names or need one column to be the gene name
#' (set in \code{name_col}) and must contain the columns indicating the differentiated expression profile.
#' @param name_col character, the name of the column in \code{DE}, which contains the gene name. If NULL, will use the rownames of \code{DE}.
#' Default is NULL.
#' @param profile_col character, the name of the column in \code{DE}, which will be used as the differentiated expression profile.
#' If DE is created by \code{getDE.limma.2G} or \code{getDE.BID.2G}, this parameter could be 'logFC' or 't'.
#' @param profile_trend character, the choice of how to display the profile, from high/positive to low/negative ('pos2neg')
#' or low/negative to high/positive ('neg2pos').Default is 'pos2neg'.
#' @param seed_driver character, name for the seed driver.
#' @param partner_driver_list a vector of characters, name for the partner driver list.
#' @param seed_driver_label character, label for the seed driver displayed on the plot. Default is seed_driver.
#' @param partner_driver_label a vector of characters, label for the partner driver list displayed on the plot. Default is partner_driver_list
#' @param driver_DA_Z a vector of numeric values, the Z statistics of differentiated activity (DA) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param driver_DE_Z a vector of numeric values, the Z statistics of differentiated expression (DE) for the driver list.
#' Better to give name to the vector, otherwise will automatically use driver list (seed + partner) as the name.
#' @param target_list a list for the target gene information for the drivers. The names for the list must contain the driver in driver list (seed + partner)
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list} could be automatically generated by \code{get_net2target_list} by
#' running \code{get.SJAracne.network}.
#' @param DA_Z_merge a vector of numeric values, the Z statistics of differentiated activity (DA) for the combination of the seed driver to partner drivers saperately.
#' Better to give name to the vector, otherwise will automatically use partner driver list as the name.
#' @param target_list_merge a list for the target gene information for thecombination of the seed driver to partner drivers saperately.
#' The names for the list must contain the driver in partner_driver_list
#' Each object in the list must be a data.frame and should contain one column ("target") to save the target genes.
#' Strongly suggest to follow the NetBID2 pipeline, and the \code{target_list_merge} could be automatically generated by \code{merge_target_list}.
#' @param top_driver_number numeric, number for the top significant partner drivers to be displayed on the plot. Default is 10.
#' @param top_order character, choice of order pattern used to display the partner drivers. Two options,'merge' or 'diff'.
#' 'merge' means the partner drivers will be sorted by the combined Z statistics.
#' 'diff' means the partner drivers will be sorted by the delta Z statistics.
#' Default is 'merge'.
#' @param target_nrow numeric, number of rows for each driver display on the plot. Two options, 1 or 2.
#' If set to 1, the target genes' position on the profile will be displayed in one row.
#' If set to 2, the target genes' position on the profile will be displayed in two rows,
#' with positive regulated genes displayed on the first row and negative regulated genes displayed on the second row.
#' Default is 2.
#' @param target_col character, choice of color pattern used to display the targets. Two options,'black' or 'RdBu'.
#' If set to 'black', the lines will be colored in black.
#' If set to 'RdBu', the lines will be colored into Red to Blue color bar.
#' If \code{target_col_type} is set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If \code{target_col_type} is set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' with significant high set for red and low for blue. The significant threshold is set by \code{profile_sig_thre}.
#' Default is 'RdBu'.
#' @param target_col_type character, choice of the pattern used to display the color for target genes, only work when \code{target_col} is set as 'RdBu'.
#' Two options,'PN' or 'DE'.
#' If set as 'PN', the positive regulated genes will be colored in red and negative regulated genes in blue.
#' If set as 'DE', the color for the target genes is set according to its value in the differentiated expression profile,
#' Default is 'PN'.
#' @param left_annotation character, annotation displayed on the left of the figure representing left condition of the rank_profile. Default is "".
#' @param right_annotation character, annotation displayed on the right of the figure representing right condition of the rank_profile. Default is "".
#' @param main character, title for the plot. Default is "".
#' @param profile_sig_thre numeric, threshold for the absolute values in profile to be treated as significance.
#' Target genes without signifcant values in the profile will be colored in grey. Only work when \code{target_col_type} is set as "DE" and \code{target_col} is set as "RdBu".
#' Default is 0.
#' @param Z_sig_thre numeric, threshold for the Z statistics in \code{driver_DA_Z} and \code{driver_DE_Z} to be treated as signifcance.
#' Only signifcant values will have background color. Default is 1.64.
#' @param pdf_file character, file path for the pdf file to save the figure into pdf format.If NULL, will not generate pdf file. Default is NULL.
#'
#' @return logical value indicating whether the plot has been successfully generated
#'
#' @examples
#' \dontrun{
#' analysis.par <- list()
#' analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
#' NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
#' ms_tab <- analysis.par$final_ms_tab
#' sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
#'                                logFC_col='logFC.G4.Vs.others_DA',
#'                                Pv_col='P.Value.G4.Vs.others_DA',
#'                                logFC_thre=0.4,
#'                                Pv_thre=1e-7,
#'                                main='Volcano Plot for G4.Vs.others_DA',
#'                                show_label=FALSE,
#'                                label_type = 'origin',
#'                                label_cex = 0.5)
#' driver_list <- rownames(sig_driver)
#' ## choose seed driver and partner driver list
#' seed_driver <- driver_list[1]
#' part_driver <- ms_tab$originalID_label
#' ## get merge target
#' merge_target <- lapply(part_driver,function(x){
#'   m1 <- merge_target_list(driver1=seed_driver,driver2=x,
#'                           target_list=analysis.par$merge.network$target_list)
#' })
#' names(merge_target) <- part_driver
#' ## get activity matrix for the merge target network
#' ac_combine_mat <- cal.Activity(all_target=merge_target,
#'                                cal_mat=exprs(analysis.par$cal.eset),
#'                                es.method='weightedmean')
#' ## get DA for the combined drivers
#' comp_name <- 'G4.Vs.others'
#' G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')]
#' # get sample list for G0
#' G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')]
#' # get sample list for G1
#' DA_driver_combine <- getDE.limma.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## or use: DA_driver_combine <- getDE.BID.2G(eset=generate.eset(ac_combine_mat),
#'                                     G1=G1,G0=G0,
#'                                     G1_name='G4',G0_name='others')
#' ## prepare for SINBA input
#' ori_part_Z <- analysis.par$DA[[comp_name]][part_driver,'Z-statistics']
#' ori_seed_Z <- analysis.par$DA[[comp_name]][seed_driver,'Z-statistics']
#' DE <- analysis.par$DE[[comp_name]]
#' driver_DA_Z <- analysis.par$DA[[comp_name]][,'Z-statistics']
#' names(driver_DA_Z) <- rownames(analysis.par$DA[[comp_name]])
#' driver_DE_Z <- analysis.par$DE[[comp_name]][,'Z-statistics']
#' names(driver_DE_Z) <- rownames(analysis.par$DE[[comp_name]])
#' DA_Z_merge <- DA_driver_combine[,'Z-statistics']
#' names(DA_Z_merge) <- rownames(DA_driver_combine)
#' target_list_merge <- merge_target
#' seed_driver_label <- ms_tab[seed_driver,'gene_label']
#' partner_driver_list <- part_driver
#' profile_col <- 't'
#' partner_driver_label <- ms_tab[partner_driver_list,'gene_label']
#' target_list <- analysis.par$merge.network$target_list
##
#' draw.NetBID_SINBA(DE=DE,profile_col = profile_col,
#'                        seed_driver=seed_driver,
#'                        partner_driver_list=partner_driver_list,
#'                        seed_driver_label=seed_driver_label,
#'                        partner_driver_label=partner_driver_label,
#'                        driver_DA_Z=driver_DA_Z,driver_DE_Z=driver_DE_Z,
#'                        target_list=target_list,
#'                        DA_Z_merge=DA_Z_merge,
#'                        target_list_merge=target_list_merge,
#'                        top_driver_number=20,profile_trend='pos2neg',
#'                        top_order='merge',Z_sig_thre = 1.64,
#'                        target_nrow=1,target_col='RdBu',target_col_type='PN',
#'                        pdf_file=sprintf('%s/NetBID_GSEA_SINBA_demo1.pdf',
#'                        analysis.par$out.dir.PLOT))
#'}
#' @export
draw.NetBID_SINBA <- function(DE=NULL,name_col=NULL,profile_col=NULL,profile_trend='pos2neg',
                              seed_driver=NULL,partner_driver_list=NULL,
                              seed_driver_label=NULL,partner_driver_label=NULL,
                              driver_DA_Z=NULL,driver_DE_Z=NULL,target_list=NULL,
                              target_list_merge=NULL,
                              DA_Z_merge=NULL,
                              diff_Z=NULL,
                              top_driver_number=10,top_order='merge',target_nrow=2,target_col='RdBu',target_col_type='PN',
                              left_annotation="",right_annotation="",main="",
                              profile_sig_thre=0,Z_sig_thre=1.64,pdf_file=NULL){
  pos_col <- brewer.pal(12,'Paired')[8]
  neg_col <- brewer.pal(12,'Paired')[4]
  if(!profile_col %in% colnames(DE)){
    message(sprintf('%s not in colnames of DE, please check and re-try!',profile_col))
    return(FALSE)
  }

  driver_list <- c(seed_driver,partner_driver_list)
  driver_list_gene <- gsub('(.*)_.*','\\1',driver_list)
  show_label <- c(seed_driver_label,partner_driver_label)

  names(target_list_merge)<-gsub("(.*:)","",names(target_list_merge))
  #sinba_id<-names(target_list_merge)
  # get names
  if(is.null(names(driver_DA_Z))) names(driver_DA_Z) <- driver_list
  if(is.null(names(driver_DE_Z))) names(driver_DE_Z) <- driver_list
  if(is.null(names(show_label))) names(show_label) <- driver_list
  #if(is.null(names(DA_Z_merge))) names(DA_Z_merge) <- driver_list
  #
  ori_part_Z <- driver_DA_Z[partner_driver_list]
  ori_seed_Z <- driver_DA_Z[seed_driver]
  #diff_Z <- 2*DA_Z_merge[partner_driver_list]-(ori_part_Z+ori_seed_Z)
  names(diff_Z) <- gsub("(.*:)","",names(diff_Z)) #change combination name to partner name

  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list_gene]; names(driver_DE_Z) <- driver_list
  #DA_Z_merge  <- DA_Z_merge[sinba_id]
  names(DA_Z_merge)<-gsub("(.*:)","",names(DA_Z_merge)) #change combination name to partner name
  if(top_order=='merge'){
    if(length(partner_driver_list)>top_driver_number){
      #partner_driver_list <- partner_driver_list[order(abs(DA_Z_merge[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive part
      #partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)][1:top_driver_number]
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(DA_Z_merge[partner_driver_list],decreasing = TRUE)]
  }else{
    if(length(partner_driver_list)>top_driver_number){
      #partner_driver_list <- partner_driver_list[order(abs(diff_Z[partner_driver_list]),decreasing = TRUE)][1:top_driver_number]
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)][1:top_driver_number] ## only consider positive increase
    }
    if(profile_trend=='pos2neg')
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = FALSE)]
    else
      partner_driver_list <- partner_driver_list[order(diff_Z[partner_driver_list],decreasing = TRUE)]
  }
  driver_list <- c(seed_driver,partner_driver_list)
  show_label <- show_label[driver_list]
  driver_DA_Z <- driver_DA_Z[driver_list]
  driver_DE_Z <- driver_DE_Z[driver_list]
  diff_Z <- diff_Z[partner_driver_list]
  DA_Z_merge <- DA_Z_merge[partner_driver_list]
  if(is.null(name_col)==TRUE){
    DE <- cbind(DE[,setdiff(colnames(DE),'name')],name=rownames(DE),stringsAsFactors=FALSE)
    name_col <- 'name'
  }
  w1 <- which(is.na(DE[,profile_col])==FALSE)
  DE <- DE[w1,]
  if(profile_trend=='pos2neg') DE <- DE[order(DE[,profile_col],decreasing = TRUE),] else DE <- DE[order(DE[,profile_col],decreasing = FALSE),]
  DE_profile <- DE[,profile_col]
  DE_profile_name <- DE[,name_col]
  ##############################################
  ## calculate layout
  target_list <- lapply(target_list,function(x){return(x[which(x$target %in% DE_profile_name),])})
  target_list_merge <- lapply(target_list_merge,function(x){return(x[which(x$target %in% DE_profile_name),])})
  n_gene <- length(DE_profile)
  if(target_nrow==2){
    n_driver <- length(partner_driver_list)*4+2
    ratio1 <- ceiling(n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  } else {
    n_driver <- length(partner_driver_list)*2+1
    ratio1 <- ceiling(1.5*n_driver/15) ## profile height to rows
    ratio2 <- 4 ## width of profile to DA/DE
    rr <- 1
  }
  #
  if(is.null(pdf_file)==FALSE){
    pdf(pdf_file,width=(rr*2+ratio2)*2,height=(ratio1+rr)*2)
  }
  # get layout
  layout(matrix(c(rep(0,length.out=rr),rep(1,length.out=ratio2),rep(0,length.out=rr*1),
                  rep(c(rep(4,length.out=rr),rep(2,length.out=ratio2),rep(3,length.out=rr*1)),
                      length.out=ratio1*(ratio2+rr*2))),
                ncol=c(ratio2+rr*2),byrow=TRUE))
  ## plot 1
  par(mar=c(1.5,1.5,1.5,0))
  mm <- quantile(DE_profile,probs=c(0.0001,0.9999));
  mm <- max(abs(mm)); mm <- c(-mm,mm)
  y1 <- seq(mm[1],mm[2],length.out=5); y1 <- round(y1,1)
  unit <- n_gene/10; unit <- round(unit/100)*100
  x1 <- seq(0,n_gene,by=unit);x1 <- unique(x1); x1 <- c(x1,max(x1)+unit)
  par(usr=c(0,length(DE_profile),mm[1],mm[2]))
  plot(DE_profile,col='white',pch=16,xaxt='n',yaxt='n',xlab="",ylab="",bty='n',type='n',ylim=c(mm[1],mm[2]),main=main,cex.main=1.8)
  pp <- par()$usr; rr <- (pp[2]-pp[1])/n_gene
  polygon(x=c(pp[1],c(1:n_gene)*rr+pp[1],pp[2]),y=c(0,DE_profile,0),col='grey',border='grey',xpd=TRUE,lwd=0.3)
  if(profile_trend=='pos2neg'){
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[2]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[1]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }else{
    if(is.null(left_annotation)==FALSE) text(pp[1]+(pp[2]-pp[1])/100,mm[1]*0.8,adj=0,left_annotation,col=brewer.pal(9,'Reds')[6],xpd=TRUE,cex=1.2)
    if(is.null(right_annotation)==FALSE) text(pp[2]-(pp[2]-pp[1])/100,mm[2]*0.8,adj=1,right_annotation,col=brewer.pal(9,'Blues')[6],xpd=TRUE,cex=1.2)
  }
  axis(side=2,at=y1,labels=y1)
  mtext(side=2,line = 2.5,profile_col,cex=1)
  segments(pp[1],mm[1],pp[2],mm[1],xpd=TRUE)
  segments(x1*rr,mm[1]-(mm[2]-mm[1])/50,x1*rr,mm[1],xpd=TRUE)
  text(x1*rr,mm[1]-(mm[2]-mm[1])/25,get_label_manual(x1),adj=0.5,xpd=TRUE)
  ## plot 2
  par(mar=c(2,1.5,2,0))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,n_gene),xaxt='n',yaxt='n')
  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  rect(xleft = pp[1],xright=pp[2],ybottom = yy1[length(yy1)-target_nrow],ytop=yy1[length(yy1)],border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.3)) ## for seed rows
  rect(xleft = pp[1],xright=pp[2],ybottom = yy2,ytop=yy4,border=NA,col=get_transparent(brewer.pal(11,'Set3')[2],0.2)) ## for combine rows
  segments(x0=pp[1],x1=pp[2],y0=yy1,y1=yy1,lwd=0.2,col='light grey')
  segments(x0=pp[1],x1=pp[2],y0=yy3,y1=yy3,lwd=1,col='dark grey')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white')
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1,col='black')
  segments(x0=pp[1],x1=pp[2],y0=yy1[length(yy1)-target_nrow],y1=yy1[length(yy1)-target_nrow],lwd=1.5,col='black')
  # shorten yy1
  dyy <- yy1[2]-yy1[1]
  yy11 <- yy1-dyy*0.3
  # add columns
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]
  if(target_col_type=='DE'){
    cc <- z2col(DE_profile,sig_thre=profile_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[5:9],blue_col=brewer.pal(9,'Blues')[5:9],
                col_max_thre=max(abs(DE_profile)))
    #names(cc) <- names(DE_profile)
    cc[which(cc=='white')] <- 'light grey'
  }
  if(target_nrow==1){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    w0 <- which(DE_profile_name %in% t1$target)
    w1 <- w0*rr+pp[1]
    if(target_col=='black'){
      segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
    }else{
      if(target_col_type=='DE'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=cc[w0])
      }else{
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],lwd=1.5,
                 col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i],y1=yy11[2*i+1],lwd=1.5,
                   col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
        }
      }
    }
    # for combine
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      w0 <- which(DE_profile_name %in% t1$target)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over)
      w1_over <- w0_over*rr+pp[1]
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=cc[w0])
        }else{
          segments(x0=w1,x1=w1,y0=yy1[2*i-1],y1=yy11[2*i],lwd=1.5,col=c(neg_col,'white',pos_col)[sign(t1$spearman)+2])
        }
      }
      points(w1_over,rep((yy11[2*i]+yy1[2*i])/2,length.out=length(w1_over)),pch='*',col='black')
    }
  }
  ###################
  if(target_nrow==2){
    # for seed driver
    t1 <- use_target_list[[seed_driver]]
    t11 <- t1[which(t1$spearman>=0),]$target
    t12 <- t1[which(t1$spearman<0),]$target
    w0 <- which(DE_profile_name %in% t11)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-1],y1=yy1[length(yy1)],col=pos_col,lwd=1.5)
        }
      }
    }
    w0 <- which(DE_profile_name %in% t12)
    w1 <- w0*rr+pp[1]
    if(length(w1)>0){
      if(target_col=='black'){
        segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col='black',lwd=1)
      }else{
        if(target_col_type=='DE'){
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=cc[w0],lwd=1.5)
        }else{
          segments(x0=w1,x1=w1,y0=yy1[length(yy1)-2],y1=yy1[length(yy1)-1],col=neg_col,lwd=1.5)
        }
      }
    }
    # for each partner driver
    for(i in 1:length(partner_driver_list)){
      t1 <- use_target_list[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i],y1=yy11[4*i+1],col=pos_col,lwd=1.5)
          }
        }
      }
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-1],y1=yy11[4*i],col=neg_col,lwd=1.5)
          }
        }
      }
    }
    # for each partner driver + seed combination
    for(i in 1:length(partner_driver_list)){
      t1 <- target_list_merge[[partner_driver_list[[i]]]]
      t11 <- t1[which(t1$spearman>=0),]$target
      t12 <- t1[which(t1$spearman<0),]$target
      w0 <- which(DE_profile_name %in% t11)
      w1 <- w0*rr+pp[1]
      t_over <- intersect(target_list[[partner_driver_list[[i]]]]$target,target_list[[seed_driver]]$target)
      w0_over <- which(DE_profile_name %in% t_over) ## setdiff(t_over,names(DE_profile)) !!!
      w1_over <- w0_over*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-2],y1=yy11[4*i-1],col=pos_col,lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-1]+yy1[4*i-1])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
      w0 <- which(DE_profile_name %in% t12)
      w1 <- w0*rr+pp[1]
      if(length(w1)>0){
        if(target_col=='black'){
          segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col='black',lwd=1)
        }else{
          if(target_col_type=='DE'){
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=cc[w0],lwd=1.5)
          }else{
            segments(x0=w1,x1=w1,y0=yy1[4*i-3],y1=yy11[4*i-2],col=neg_col,lwd=1.5)
          }
        }
      }
      points(intersect(w1_over,w1),rep((yy11[4*i-2]+yy1[4*i-2])/2,length.out=length(intersect(w1_over,w1))),pch='*',col='black')
    }
    ####
  }
  ###################
  ## plot 3
  par(mar=c(2,0.5,2,2))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,3),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr
  rect(xleft=pp[1],xright=pp[2],ybottom=pp[3],ytop=pp[4])

  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=2,col='white',xpd=TRUE)
  segments(x0=pp[1],x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='dark grey',xpd=TRUE)
  abline(v=xx1)
  ## add text
  mm_min <- min(min(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*0.9,
                min(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*0.9,min(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*0.9)
  mm_min <- max(mm_min,Z_sig_thre)
  mm_max <- max(max(abs(driver_DA_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(driver_DE_Z[partner_driver_list]),na.rm=TRUE)*1.1,
                max(abs(diff_Z[partner_driver_list]),na.rm=TRUE)*1.1,max(abs(DA_Z_merge[partner_driver_list]),na.rm=TRUE)*1.1)
  c1 <- z2col(driver_DA_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c2 <- z2col(driver_DE_Z[driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c3 <- z2col(diff_Z[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  c4 <- z2col(DA_Z_merge[partner_driver_list],sig_thre=Z_sig_thre,n_len=100,red_col = brewer.pal(9,'Reds')[7],blue_col=brewer.pal(9,'Blues')[7],
              col_min_thre=mm_min,col_max_thre=mm_max)
  # for seed driver
  yy1<-yy3
  z1 <- driver_DA_Z[seed_driver]
  z2 <- driver_DE_Z[seed_driver]
  p1 <- get_z2p(z1)
  p2 <- get_z2p(z2)
  rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c1[1],border='dark grey',xpd=TRUE)
  rect(xright=xx1[3],xleft=xx1[2],ybottom=yy1[length(yy1)-1],ytop=yy1[length(yy1)],col=c2[1],border='dark grey',xpd=TRUE)
  text(x=sum(xx1[1:2])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p1,adj=0.5)
  text(x=sum(xx1[3:4])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,p2,adj=0.5)
  text(x=sum(xx1[2:3])/2,y=(yy1[length(yy1)-1]+yy1[length(yy1)])/2,'-',adj=0.5)
  # for partner driver
  for(i in 1:length(partner_driver_list)){
    z1 <- driver_DA_Z[partner_driver_list[i]]
    z2 <- driver_DE_Z[partner_driver_list[i]]
    z3 <- diff_Z[partner_driver_list[i]]
    z4 <- DA_Z_merge[partner_driver_list[i]]
    p1 <- get_z2p(z1)
    p2 <- get_z2p(z2)
    p3 <- format(z3,digits=3)
    p4 <- get_z2p(z4)
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i],ytop=yy1[2*i+1],col=c1[i+1],border='dark grey',xpd=TRUE) ## DA_Z
    rect(xleft=xx1[1],xright=xx1[2],ybottom=yy1[2*i-1],ytop=yy1[2*i],col=c4[i],border='dark grey',xpd=TRUE) ## merge_Z
    rect(xleft=xx1[3],xright=xx1[4],ybottom=yy2[i],ytop=yy2[i+1],col=c2[i+1],border='dark grey',xpd=TRUE) ## DE
    rect(xleft=xx1[2],xright=xx1[3],ybottom=yy2[i],ytop=yy2[i+1],col=c3[i],border='dark grey',xpd=TRUE) ## delta Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i]+yy1[2*i+1])/2,p1,adj=0.5) ## DA_Z
    text(x=sum(xx1[1:2])/2,y=(yy1[2*i-1]+yy1[2*i])/2,p4,adj=0.5) ## merge_Z
    text(x=sum(xx1[3:4])/2,y=(yy2[i]+yy2[i+1])/2,p2,adj=0.5) ## DE
    text(x=sum(xx1[2:3])/2,y=(yy2[i]+yy2[i+1])/2,p3,adj=0.5) ## delta z
  }
  textheight <- strheight('DA',units='user',cex=1.5)
  text(sum(xx1[1:2])/2,pp[4]+textheight,'DA',xpd=TRUE,cex=1.5)
  textheight <- strheight('DE',units='user',cex=1.5)
  text(sum(xx1[3:4])/2,pp[4]+textheight,'DE',xpd=TRUE,cex=1.5)
  textheight <- strheight('deltaZ',units='user',cex=1.5)
  text(sum(xx1[2:3])/2,pp[4]+textheight,'deltaZ',xpd=TRUE,cex=1.5)
  ## plot 4
  par(mar=c(2,6,2,0.2))
  plot(1,col='white',xlab="",ylab="",xlim=c(0,2),xaxt='n',yaxt='n',bty='n')
  pp <- par()$usr;rr <- (pp[2]-pp[1])/n_gene
  yy1 <- seq(from=pp[3],to=pp[4],length.out=n_driver+1) # separate each detail
  yy3 <- seq(from=pp[3],to=pp[4],length.out=length(partner_driver_list)*2+1+1) # separate each driver(pos/neg)
  yy2 <- yy1[seq(from=1,to=length(yy1),by=target_nrow*2)] # separate each driver combine
  yy4 <- yy2+(yy1[2]-yy1[1])*target_nrow
  xx1 <- seq(pp[1],pp[2],length.out=4)
  yy22 <- (yy2[1:(length(yy2)-1)]+yy2[2:length(yy2)])/2
  dyy22 <- yy22[2]-yy22[1]
  yy33 <- (yy3[1:(length(yy3)-1)]+yy3[2:length(yy3)])/2
  dyy33 <- yy33[2]-yy33[1]
  xleft <- pp[1]+(pp[2]-pp[1])*0.55
  tt <- pp[2]-xleft
  text(show_label[1],x=pp[1]+(pp[1]+pp[2])*0.53,y=(yy3[length(yy3)-1]+yy3[length(yy3)])/2,xpd=TRUE,adj=1,cex=1.2) ## label for seed
  text(show_label[2:length(show_label)],x=pp[1]+(pp[1]+pp[2])*0.52,y=yy22,xpd=TRUE,adj=1) ## label for partner
  # add target size
  use_target_list <- target_list[driver_list]
  use_merge_target_list <- target_list_merge[partner_driver_list]
  target_size <- do.call(rbind,lapply(use_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  merge_target_size <- do.call(rbind,lapply(use_merge_target_list,function(x){
    x1 <- length(which(x$spearman>=0))
    x2 <- length(which(x$spearman<0))
    c(x1,x2)
  }))
  # for seed driver
  if(target_nrow==2){
    mm <- max(merge_target_size)
    i <- length(yy1)-1
    rect(xleft=xleft,xright=xleft+target_size[1,1]/mm*tt,
         ybottom=yy1[i],ytop=yy1[i]+dyy22/2*0.35,col=pos_col,border=NA)
    rect(xleft=xleft,xright=xleft+target_size[1,2]/mm*tt,
         ytop=yy1[i],ybottom=yy1[i]-dyy22/2*0.35,col=neg_col,border=NA)
  }else{
    target_size <- rowSums(target_size)
    merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    i <- length(yy1)
    rect(xleft=xleft,xright=xleft+target_size[1]/mm*tt,
         ybottom=(yy1[i]+yy1[i-1])/2-dyy22*0.2/2,ytop=(yy1[i]+yy1[i-1])/2+dyy22*0.2/2,col='dark grey',border=NA)
  }
  # for partner
  if(target_nrow==2){
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1,1]/mm*tt,
           ybottom=yy33[2*i],ytop=yy33[2*i]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+target_size[i+1,2]/mm*tt,
           ytop=yy33[2*i],ybottom=yy33[2*i]-dyy33*0.35,col=neg_col,border=NA)
      # merge
      rect(xleft=xleft,xright=xleft+merge_target_size[i,1]/mm*tt,
           ybottom=yy33[2*i-1],ytop=yy33[2*i-1]+dyy33*0.35,col=pos_col,border=NA)
      rect(xleft=xleft,xright=xleft+merge_target_size[i,2]/mm*tt,
           ytop=yy33[2*i-1],ybottom=yy33[2*i-1]-dyy33*0.35,col=neg_col,border=NA)
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  #
  if(target_nrow==1){
    #target_size <- rowSums(target_size)
    #merge_target_size <- rowSums(merge_target_size)
    mm <- max(merge_target_size)
    for(i in 1:length(partner_driver_list)){
      rect(xleft=xleft,xright=xleft+target_size[i+1]/mm*tt,ybottom=yy33[i*2]-dyy33*0.2,ytop=yy33[i*2]+dyy33*0.2,col='dark grey',border=NA) #each
      rect(xleft=xleft,xright=xleft+merge_target_size[i]/mm*tt,ybottom=yy33[i*2-1]-dyy33*0.2,ytop=yy33[i*2-1]+dyy33*0.2,col='dark grey',border=NA) #merge
    }
    segments(x0=xleft,x1=pp[2],y0=pp[4],y1=pp[4],xpd=TRUE)
    sst <- round(seq(0,mm,length.out=3))
    ss <- sst*tt/mm+xleft
    segments(x0=ss,x1=ss,y0=pp[4],y1=pp[4]+(pp[4]-pp[3])/150,xpd=TRUE)
    text(x=ss,y=pp[4]+(pp[4]-pp[3])/100,srt=90,sst,xpd=TRUE,adj=0,cex=0.8)
    text('Target Size',x=(pp[1]+pp[2])*0.45,y=pp[4]+(pp[4]-pp[3])/50,adj=1,xpd=TRUE,cex=0.8)
  }
  ## add lines
  segments(x0=xleft,x1=pp[2],y0=yy2,y1=yy2,lwd=1.2,col='grey',xpd=TRUE)
  abline(v=xleft,col='grey')
  #abline(v=pp[2],col='grey')
  ## test for significant overlap
  total_possible_target <- unique(unlist(lapply(target_list,function(x)x$target)))
  for(i in 1:length(partner_driver_list)){
    res1 <- test.targetNet.overlap(seed_driver,partner_driver_list[i],
                                   target1=intersect(use_target_list[[seed_driver]]$target,DE_profile_name),
                                   target2=intersect(use_target_list[[partner_driver_list[i]]]$target,DE_profile_name),
                                   total_possible_target = total_possible_target)
    pv <- format(res1[1],digits=2,scientific=TRUE)
    ov <- round(res1[3])
    if(res1[1]<0.05){
      text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.8,col='dark red')
    }else{
      if(ov==0){
        text('No overlap',x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }else{
        text(sprintf('Overlap:%d, P value:%s',ov,pv),x=xleft-(pp[2]-pp[1])/10,y=yy33[i*2-1],adj=1,xpd=TRUE,cex=0.7)
      }
    }
  }
  ##
  if(is.null(pdf_file)==FALSE) dev.off()
  #layout(1);
  return(TRUE)
}

#'
#' @title draw.ridgeplot_SINBA
#' @description show target gene distribution profile of seed-partner, seed and parter driver.
#' @param Seed_gene_ID character, the seed gene name.
#' @param Partner_gene_ID character, the partner gene name.
#' @param target_list list, the driver-to-target list object. \code{target_list} is generated by \code{get.combined.network} and contains all the driver gene regulons.
#' @param DE, data frame, the DE or DA table.
#' @param ID_col, character, the ID_col must include gene name used in the target list.
#' @param Seed_show_label, character, the name of seed gene to show in the plot.
#' @param Partner_show_label, character, the name of partner gene to show in the plot.
#' @return ggplot object
#' @importFrom ggplot2 ggplot
#' @noRd
#'
#' @examples
#' Seed_gene_ID<-"LCK"
#' Partner_gene_ID<-"STAT3"
#' profile_col<-"t"
#' ID_col<-"ID"
#' Seed_show_label<-"LCK"
#' Partner_show_label<-"STAT3"
#' DE<-SINBA.par$DE[[comp_name]]
#'
draw.ridgeplot_SINBA<-function(Seed_gene_ID=NULL,Partner_gene_ID=NULL,target_list=NULL,DE=NULL,profile_col=NULL,ID_col=NULL,Seed_show_label=NULL,Partner_show_label=NULL){
  Seed_gene_targets<-DE[target_list[[which(names(target_list)==Seed_gene_ID)]]$target,c(ID_col,profile_col)]
  Seed_gene_targets<-Seed_gene_targets[complete.cases(Seed_gene_targets),]
  Seed_gene_targets$group<-Seed_show_label
  Seed_gene_targets[,profile_col]<-round(abs(Seed_gene_targets[,profile_col]),digits = 1)
  #Seed_gene_targets[,profile_col]<-scale(Seed_gene_targets[,profile_col],center = T,scale = T)[,1]
  #Seed_center<-attributes(scale(Seed_gene_targets[,profile_col],center = T,scale = T))[[2]]

  Partner_gene_targets<-DE[target_list[[which(names(target_list)==Partner_gene_ID)]]$target,c(ID_col,profile_col)]
  Partner_gene_targets<-Partner_gene_targets[complete.cases(Partner_gene_targets),]
  Partner_gene_targets$group<-Partner_show_label
  Partner_gene_targets[,profile_col]<-round(abs(Partner_gene_targets[,profile_col]),digits = 1)
  #Partner_gene_targets[,profile_col]<-scale(Partner_gene_targets[,profile_col],center = T,scale = T)[,1]
  #Partner_center<-attributes(scale(Partner_gene_targets[,profile_col],center = T,scale = T))[[2]]

  merge_SP_targets<-DE[unique(c(Seed_gene_targets[,ID_col],Partner_gene_targets[,ID_col])),c(ID_col,profile_col)]
  merge_SP_targets$group<-paste(Seed_show_label,Partner_show_label,sep = ":")
  merge_SP_targets[,profile_col]<-round(abs(merge_SP_targets[,profile_col]),digits = 1)
  #merge_SP_targets[,profile_col]<-scale(merge_SP_targets[,profile_col],center = T,scale = T)[,1]
  #merge_SP_center<-attributes(scale(merge_SP_targets[,profile_col],center = T,scale = T))[[2]]

  dat<-rbind.data.frame(Seed_gene_targets,Partner_gene_targets,merge_SP_targets)
  col<-brewer.pal(3,"Set3")
  names(dat)[2]<-"profile"
  g<-ggplot(dat, aes(x = profile, y = group, fill = group)) +
    geom_density_ridges()+
    #theme_ridges() +
    scale_fill_manual(values=col)+
    xlab(profile_col)+ylab("")+theme_classic()+
    theme(legend.position = "none")
  return(g)
}
#'
#' @title gg_color_hue
#' @description get default colors
#' @concept ggplot color
#' @importFrom grDevices hcl
#' @noRd
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @title draw.box_SINBA
#' @description draw.box_SINBA draws a scatter box plot to visualize one selected seed partner driver single and combined activity values across different phenotype subgroups of samples. Three side-by-side scatter box plots will be created. The left plot shows combined driver paire activity values in different phenotype subgroups, each point is a sample. The middle and right plot shows seed and partner driver activity values in different phenotype subgroups, each point is a sample.
#' @param ac_mat_sinba, numeric matrix, the activity values of the selected seed-partner driver across all samples.
#' @param ac_mat_single, numeric matrix, the activity values of single driver across all samples
#' @param SINBA_originalID, seed and parterne driver combination, is generated by \code{combineNet2target_list} with return_type="SINBA_originalID"
#' @param phenotype_info, a data frame with sample group information
#' @param group_col, character, the column name to be used for group
#' @param group_levels, a vector of characters, the order of the group
#' @param group_colors, a vector of characters, colors to be used for each group
#' @param main, character, the main title of the plot
#'
#' @examples
#' ac_mat_sinba<-exprs(SINBA.par$SINBA.ac.eset)
#' ac_mat_single<-exprs(SINBA.par$SINBA_single.ac.eset)
#' SINBA_originalID<-"LCK:STAT3"
#' phenotype_info<-pData(SINBA.par$SINBA.ac.eset)
#' phenotype_info$condition<-ifelse(phenotype_info$das_lc50>80,"Resistant","Sensitive")
#' g<-draw.box_SINBA(ac_mat_sinba=ac_mat_sinba,ac_mat_single=ac_mat_single,SINBA_originalID=SINBA_originalID,phenotype_info=phenotype_info,group_col="condition",group_levels=c("Sensitive","Resistant"),group_colors=c("red","blue"))
#' @export
draw.box_SINBA<-function(ac_mat_sinba=NULL,
                         ac_mat_single=NULL,
                         SINBA_originalID=NULL,
                         phenotype_info=NULL,
                         group_col=NULL,
                         group_levels=NULL,
                         group_colors=NULL,
                         main=""){
  if(missing(ac_mat_sinba) || missing(ac_mat_single) || missing (phenotype_info) || missing(group_col) || missing(SINBA_originalID)){
    stop("ac_mat_sinba,ac_mat_single,SINBA_originalID,phenotype_info and group_col are required input for the function!")
  }

  use_samples<-intersect(colnames(ac_mat_sinba),rownames(phenotype_info))
  use_samples<-intersect(colnames(ac_mat_single),use_samples)

  seed<-gsub(":.*","",SINBA_originalID);partner<-gsub(".*:","",SINBA_originalID)

  o1<-base::isTRUE(SINBA_originalID%in%rownames(ac_mat_sinba))
  o2<-base::isTRUE(all(c(seed,partner)%in%rownames(ac_mat_single)))
  if(!o1||!o2) {
    stop("SINBA_originalID=seed:partner; SINBA_originalID should be included in the rownames of ac_mat_sinba and seed partner should be included in the ac_mat_single.")
  }

  if(length(use_samples)==0) {
    stop("Column names of the exp_mat should have overlap with row names of phenotype_info.")
  }

  use_mat_sinba<-t(ac_mat_sinba[SINBA_originalID,use_samples,drop=F])
  use_mat_single_s<-t(ac_mat_single[seed,use_samples,drop=F])
  use_mat_single_p<-t(ac_mat_single[partner,use_samples,drop=F])
  df<-cbind(use_mat_sinba,use_mat_single_s,use_mat_single_p)
  df<-data.frame(df)
  colnames(df)<-c("SINBA","Seed","Partner")
  df$sampleID<-rownames(df)
  df<-reshape::melt(df,id="sampleID")

  use_phe<-data.frame(phenotype_info[use_samples, group_col,drop=F])
  names(use_phe)[1]<-"group"
  use_phe$sampleID<-rownames(use_phe)

  df<-dplyr::left_join(df,use_phe,by="sampleID")

  df$variable<-factor(df$variable,levels = c("SINBA","Seed","Partner"))

  if(!is.null(group_levels)) df$group<-factor(df$group,levels = group_levels) else df$group<-factor(df$group)

  if(!is.null(group_colors)) use_colors<-group_colors else use_colors<-gg_color_hue(length(group_levels))

  main_title<-ifelse(main=="",SINBA_originalID,main)

  p <- ggpubr::ggboxplot(df, x = "group", y = "value",
                         color = "group",
                         add = "jitter",
                         title = main_title,
                         outlier.shape=NA, add.params = list(size=1))+
    scale_colour_manual(name = "group",
                        values = use_colors)+
    guides(color=guide_legend(title="group"),fill=F)+xlab("")+ylab("Activity")+
    theme(legend.position = "none",
          legend.title = element_text(size=8,face = "bold"),
          legend.text = element_text(size=8,face = "plain"),
          title = element_text(size=12,face = "bold"),
          axis.text.x=element_text(size=12,face = "bold")
    )+facet_wrap(~variable,ncol = 3,scales = "free")

  return(p)
}
#' @title draw.ROC_SINBA
#' @description draw.ROC_SINBA draws ROC curves, including both seed and partner drive curves and the combination driver curve.
#' @param ac_mat_sinba, numeric matrix, the activity values of the selected seed-partner driver across all samples.
#' @param ac_mat_single, numeric matrix, the activity values of single driver across all samples
#' @param SINBA_originalID, seed and parterne driver combination, is generated by \code{combineNet2target_list} with return_type="SINBA_originalID"
#' @param phenotype_info, a data frame with sample group information
#' @param group_col, character, the column name to be used for group
#' @param G1_name, a vector of characters, the sample names of experimental group.
#' @param G0_name, a vector of characters, the sample names of NON-experimental group.
#' @param main, character, the main title of the plot
#' @param pdf_file, character, the file path to save as PDF file. If NULL, no PDF file will be save. Default is NULL.
#'
#' @examples
#' ac_mat_sinba<-exprs(SINBA.par$SINBA.ac.eset)
#' ac_mat_single<-exprs(SINBA.par$SINBA_single.ac.eset)
#' SINBA_originalID<-"LCK:NOTCH1"
#' phenotype_info<-pData(SINBA.par$SINBA.ac.eset)
#' phenotype_info$condition<-ifelse(phenotype_info$das_lc50>80,"Resistant","Sensitive")
#' draw.ROC_SINBA(ac_mat_sinba=ac_mat_sinba,ac_mat_single=ac_mat_single,SINBA_originalID=SINBA_originalID,phenotype_info=phenotype_info,group_col="condition",G1_name="Sensitive",G0_name="Resistant",pdf_file=sprintf("./ROC_%s.pdf",SINBA_originalID))
#' @importFrom ROCR prediction performance
#' @export
#'
draw.ROC_SINBA<-function(ac_mat_sinba=NULL,
                         ac_mat_single=NULL,
                         SINBA_originalID=NULL,
                         phenotype_info=NULL,
                         group_col=NULL,
                         G1_name=NULL,
                         G0_name=NULL,
                         main="",
                         pdf_file=NULL){
  if(missing(ac_mat_sinba) || missing(ac_mat_single) || missing (phenotype_info) || missing(group_col) || missing(SINBA_originalID) || missing(G1_name) || missing(G0_name)){
    stop("ac_mat_sinba,ac_mat_single,SINBA_originalID,phenotype_info, group_col, G1_name and G0_name are required input for the function!")
  }

  use_samples<-intersect(colnames(ac_mat_sinba),rownames(phenotype_info))
  use_samples<-intersect(colnames(ac_mat_single),use_samples)

  seed<-gsub(":.*","",SINBA_originalID);partner<-gsub(".*:","",SINBA_originalID)

  o1<-base::isTRUE(SINBA_originalID%in%rownames(ac_mat_sinba))
  o2<-base::isTRUE(all(c(seed,partner)%in%rownames(ac_mat_single)))
  if(!o1||!o2) {
    stop("SINBA_originalID=seed:partner; SINBA_originalID should be included in the rownames of ac_mat_sinba and seed partner should be included in the ac_mat_single.")
  }

  if(length(use_samples)==0) {
    stop("Column names of the exp_mat should have overlap with row names of phenotype_info.")
  }

  use_mat_sinba<-t(ac_mat_sinba[SINBA_originalID,use_samples,drop=F])
  use_mat_single_s<-t(ac_mat_single[seed,use_samples,drop=F])
  use_mat_single_p<-t(ac_mat_single[partner,use_samples,drop=F])
  df<-cbind(use_mat_sinba,use_mat_single_s,use_mat_single_p)
  df<-data.frame(df)
  colnames(df)<-c("SINBA","Seed","Partner")
  df$sampleID<-rownames(df)
  #df<-reshape::melt(df,id="sampleID")

  use_phe<-data.frame(phenotype_info[use_samples, group_col,drop=F])
  names(use_phe)[1]<-"group"
  use_phe$sampleID<-rownames(use_phe)

  df<-dplyr::left_join(df,use_phe,by="sampleID")
  df<-df[df$group%in%c(G1_name,G0_name),]
  df$condition<-ifelse(df$group==G1_name,1,0)

  pred <- prediction(df$SINBA,df$condition )
  pred2 <- prediction(df$Seed,df$condition)
  pred3<-prediction(df$Partner,df$condition)

  perf <- performance( pred, "tpr", "fpr" )
  perf2 <- performance(pred2, "tpr", "fpr")
  perf3 <- performance(pred3, "tpr", "fpr")

  main_title<-ifelse(main=="",SINBA_originalID,main)

  pdf(file = pdf_file,height = 5,width = 5)

  plot(perf, colorize = F,col="blue")
  plot(perf2, add = TRUE, colorize = F,col="red")
  plot(perf3, add = TRUE, colorize = F,col="green")
  segments(x0 = 0, y0 = 0, x1 = 1, y1=1,col="grey")
  legend(0.85, 0.15, legend=c("SINBA", "Seed","Partner"),
         col=c("blue","red","green"), lty=rep(1,3), cex=0.5)
  title(main = main_title)
  dev.off()
}
###########################################################
#######################DRUG annotation DGIdb###############
#' @title sourceDatabases
#' @description the available source database used in DGIdb Users can select partial of the source database
#' @examples
#' sourceDatabases()
#' @export
#'
sourceDatabases <- function() {
  sourceDatabases <- c("CGI","CIViC","COSMIC","CancerCommons","ChemblInteractions","ClearityFoundationBiomarkers","ClearityFoundationClinicalTrial","DTC","DoCM","DrugBank","FDA","GuideToPharmacology","JAX-CKB","MyCancerGenome","MyCancerGenomeClinicalTrial","NCI","OncoKB","PharmGKB","TALC","TEND","TTD","TdgClinicalTrial")
  return(sourceDatabases)
}

#' @title interactionTypes
#' @description the available interactionTypes used in DGIdb. Users can select specific interactionTypes
#' @examples
#' interactionTypes()
#' @export
#'
interactionTypes <- function() {
  interactionTypes <- c("activator","adduct","agonist","allosteric modulator","antagonist","antibody","antisense oligonucleotide","binder","blocker","chaperone","cleavage","cofactor","inducer","inhibitor","inhibitory allosteric modulator","inverse agonist","ligand","modulator","multitarget","n/a","negative modulator","other/unknown","partial agonist","partial antagonist","positive modulator","potentiator","product of","stimulator","substrate","suppressor","vaccine")
  return(interactionTypes)
}

#' @title geneCategories
#' @description the available geneCategories used in DGIdb. Users can select specific geneCategories
#' @examples
#' geneCategories()
#' @export
#'
geneCategories <- function() {
  geneCategories <- c("ABC TRANSPORTER","B30_2 SPRY DOMAIN","CELL SURFACE","CLINICALLY ACTIONABLE","CYTOCHROME P450","DNA DIRECTED RNA POLYMERASE","DNA REPAIR","DRUG METABOLISM","DRUG RESISTANCE","DRUGGABLE GENOME","ENZYME","EXCHANGER","EXTERNAL SIDE OF PLASMA MEMBRANE","FIBRINOGEN","G PROTEIN COUPLED RECEPTOR","GROWTH FACTOR","HISTONE MODIFICATION","HORMONE ACTIVITY","ION CHANNEL","KINASE","LIPASE","LIPID KINASE","METHYL TRANSFERASE","MYOTUBULARIN RELATED PROTEIN PHOSPHATASE","NEUTRAL ZINC METALLOPEPTIDASE","NUCLEAR HORMONE RECEPTOR","PHOSPHATIDYLINOSITOL 3 KINASE","PHOSPHOLIPASE","PROTEASE","PROTEASE INHIBITOR","PROTEIN PHOSPHATASE","PTEN FAMILY","RNA DIRECTED DNA POLYMERASE","SERINE THREONINE KINASE","SHORT CHAIN DEHYDROGENASE REDUCTASE","THIOREDOXIN","TRANSCRIPTION FACTOR","TRANSCRIPTION FACTOR BINDING","TRANSCRIPTION FACTOR COMPLEX","TRANSPORTER","TUMOR SUPPRESSOR","TYROSINE KINASE","UNKNOWN")
  return(geneCategories)
}

#' @title getResultSummary
#' @description organized the result from queryDgidbPost. Deprecated Use \code{getResultSummary_gene_V5()} instead.
#' @noRd
#'
getResultSummary <- function(gene, output) {
  .Deprecated("getResultSummary", "getResultSummary_gene_V5")
  # Row index with gene interaction information
  idx <- which(output$geneName == gene)
  # Prepare table of interactions (drug vs. interaction DB)

  tmp <- data.frame(drugName=output$interactions[[idx]]$drugName,
                    Sources=unlist(lapply(output$interactions[[idx]]$sources,function(x){return(paste(x,collapse = ","))})),
                    PMIDS=unlist(lapply(output$interactions[[idx]]$pmids,function(x){return(paste(x,collapse = ","))})),
                    InteractionTypes=unlist(lapply(output$interactions[[idx]]$interactionTypes,function(x){return(paste(x,collapse = ","))})),
                    Score=output$interactions[[idx]]$score,
                    Source_count=unlist(lapply(output$interactions[[idx]]$sources,function(x){return(length(x))})),
                    drugConceptId=output$interactions[[idx]]$drugConceptId)

  #result <- table(output[idx,]$interactions[[1]]$drugName,output[idx,]$interactions[[1]]$sources)
  # Expand matrix to all possible DBs, set multiple occurances of
  # interactions to one, and add gene and drug names
  #tmp <- data.frame(matrix(0, nrow = nrow(result), ncol = 3 + length(sources),dimnames = list(NULL,c('Gene', 'Drug', sources, 'Score'))),stringsAsFactors = FALSE)
  #tmp[,colnames(result)] <- result
  #tmp[tmp > 1] <- 1 # Remove double counts
  #tmp$Score <- rowSums(tmp)
  tmp$Gene <- rep(gene, nrow(tmp))
  #tmp$Drug <- rownames(result)
  # Determine type of interaction
  #resultType <- table(output[idx,]$interactions[[1]]$drugName,
  #                    output[idx,]$interactions[[1]]$interactionType)
  #if (nrow(tmp) == 1) {
  #    tmp$Type <- paste(colnames(resultType), collapse = ",")
  #} else {
  #    listResult <- lapply(split(resultType, seq(nrow(resultType))),
  #                         function(x, names) { names[x>0] },
  #                         colnames(resultType))
  #    tmp$Type <- sapply(listResult, paste, collapse=',')
  #}
  return(as.matrix(tmp))
}
#' @title getResultSummary_gene_V5
#' @description organized the result from queryDgidb_V5_gene
#' @noRd
getResultSummary_gene_V5 <- function(gene, output) {
    if(length(output)>1){
    tmp<- do.call("rbind",output)
  }else{
      tmp <- output[[1]]
    }
  tmp$queryGene <- rep(gene, nrow(tmp))
  return(tmp)
}
#' @title getResultSummary_drug
#' @description organized the result from queryDgidbPost_drug. Deprecated Use \code{getResultSummary_drug_V5()} instead.
#' @noRd
#'
getResultSummary_drug <- function(drug, output) {
   .Deprecated("getResultSummary_drug", "getResultSummary_drug_V5")
  # Row index with gene interaction information
  idx <- which(output$drugName == drug)
  # Prepare table of interactions (drug vs. interaction DB)

  tmp <- data.frame(geneName=output$interactions[[idx]]$geneName,
                    Sources=unlist(lapply(output$interactions[[idx]]$sources,function(x){return(paste(x,collapse = ","))})),
                    PMIDS=unlist(lapply(output$interactions[[idx]]$pmids,function(x){return(paste(x,collapse = ","))})),
                    InteractionTypes=unlist(lapply(output$interactions[[idx]]$interactionTypes,function(x){return(paste(x,collapse = ","))})),
                    Score=output$interactions[[idx]]$score,
                    Source_count=unlist(lapply(output$interactions[[idx]]$sources,function(x){return(length(x))})))

  tmp$Drug <- rep(drug, nrow(tmp))
  return(as.matrix(tmp))
}
#' @title getResultSummary_drug_V5
#' @description organized the result from queryDgidb_V5_drug
#' @param output output from query res$data$drugs$nodes$interactions
#' @param drug the used drug name
#' @noRd
getResultSummary_drug_V5 <- function(drug, output) {
  if(length(output)>1){
    tmp<- do.call("rbind",output)
  }else{
      tmp <- output[[1]]
    }
  tmp$Drug <- rep(drug, nrow(tmp))
  return(tmp)
}
#' @title queryDGIdb.SINBA
#' @description a wrapper to query DGIdb API with gene name. Deprecated Use \code{queryDgidb_V5_gene()} instead.
#' @return a list contains 1.a data frame("summary_result") of summarized gene,drug and interaction information; 2.raw query data("raw_result") from DGIdb (optional).
#' @param genes a character vector of genes for which drug interactions are queried.
#' @param source_Databases a character vector of source databses to be used. Users can call helper function sourceDatabases() to find all the available database. To query all available databases, skip argument (highly recommend).
#' @param gene_Categories a character vector of gene categories to be used. Users can call helper function geneCategories() to find all the available categories. To query all available categories, skip argument.
#' @param interaction_Types a character vector of gene drug interaction types to be used. Users can call helper function interactionTypes() to find all the available interactionTypes. To query all available interaction types, skip argument.
#' @param return_raw a logical value, default is FALSE, only return summarized information. If users want to check raw query data, please change to TRUE.
#' @noRd
#' @examples
#' genes<-c("ATM","BCL2L1","LCK","CDK4")
#' res<-queryDGIdb.SINBA(genes=genes,interaction_Types=c("inhibitor","antagonist"))
#' res<-queryDGIdb.SINBA(genes=genes,source_Databases =c("ChemblInteractions","TTD","OncoKB","MyCancerGenome","DrugBank","FDA"),interaction_Types=c("inhibitor","antagonist"))
#'
queryDGIdb.SINBA <- function(genes=NULL,
                             source_Databases = c("CGI","CIViC","COSMIC","CancerCommons","ChemblInteractions","ClearityFoundationBiomarkers","ClearityFoundationClinicalTrial","DTC","DoCM","DrugBank","FDA","GuideToPharmacology","JAX-CKB","MyCancerGenome","MyCancerGenomeClinicalTrial","NCI","OncoKB","PharmGKB","TALC","TEND","TTD","TdgClinicalTrial"),
                             gene_Categories = c("ABC TRANSPORTER","B30_2 SPRY DOMAIN","CELL SURFACE","CLINICALLY ACTIONABLE","CYTOCHROME P450","DNA DIRECTED RNA POLYMERASE","DNA REPAIR","DRUG METABOLISM","DRUG RESISTANCE","DRUGGABLE GENOME","ENZYME","EXCHANGER","EXTERNAL SIDE OF PLASMA MEMBRANE","FIBRINOGEN","G PROTEIN COUPLED RECEPTOR","GROWTH FACTOR","HISTONE MODIFICATION","HORMONE ACTIVITY","ION CHANNEL","KINASE","LIPASE","LIPID KINASE","METHYL TRANSFERASE","MYOTUBULARIN RELATED PROTEIN PHOSPHATASE","NEUTRAL ZINC METALLOPEPTIDASE","NUCLEAR HORMONE RECEPTOR","PHOSPHATIDYLINOSITOL 3 KINASE","PHOSPHOLIPASE","PROTEASE","PROTEASE INHIBITOR","PROTEIN PHOSPHATASE","PTEN FAMILY","RNA DIRECTED DNA POLYMERASE","SERINE THREONINE KINASE","SHORT CHAIN DEHYDROGENASE REDUCTASE","THIOREDOXIN","TRANSCRIPTION FACTOR","TRANSCRIPTION FACTOR BINDING","TRANSCRIPTION FACTOR COMPLEX","TRANSPORTER","TUMOR SUPPRESSOR","TYROSINE KINASE","UNKNOWN"),
                             interaction_Types = c("activator","adduct","agonist","allosteric modulator","antagonist","antibody","antisense oligonucleotide","binder","blocker","chaperone","cleavage","cofactor","inducer","inhibitor","inhibitory allosteric modulator","inverse agonist","ligand","modulator","multitarget","n/a","negative modulator","other/unknown","partial agonist","partial antagonist","positive modulator","potentiator","product of","stimulator","substrate","suppressor","vaccine"),
                             return_raw=F) {
  .Deprecated("queryDGIdb.SINBA", "queryDgidb_V5_gene")
  #,curatedOnly = c(FALSE, TRUE)) {
  if (missing(genes)) stop("Need to specify a vector of genes to query.")

  if (is.null(genes) || length(genes) == 0 ||
      !is.character(genes) || genes == "") {
    stop("Need to specify a non-empty vector of genes names.")
  }

  if (missing(source_Databases) |
      all(sourceDatabases() %in% source_Databases)) {
    databases <- NULL
  } else {
    databases <- match.arg(arg = source_Databases,
                           choices = sourceDatabases(),
                           several.ok = TRUE)
    databases <- paste(databases, collapse = ",")
  }
  if (missing(gene_Categories) |
      all(geneCategories() %in% gene_Categories)) {
    categories <- NULL
  } else {
    categories <- match.arg(arg = gene_Categories,
                            choices = geneCategories(),
                            several.ok = TRUE)
    categories <- paste(categories, collapse=",")
  }
  if (missing(interaction_Types) |
      all(interactionTypes() %in% interaction_Types)) {
    interactions <- NULL
  } else {
    interactions <- match.arg(arg = interaction_Types,
                              choices = interactionTypes(),
                              several.ok = TRUE)
    interactions <- paste(interactions, collapse = ",")
  }

  # Check internet connection
  tryCatch({
    msg <- ""
    r <- GET("https://dgidb.org/api/v2/interaction_types.json")
    if (status_code(r) != 200) {
      msg <- "DGIdb service not available."
    }
  }, error = function(err) {
    msg <- "Check internet connection"
  })
  if (msg != "")
    stop(msg)

  # Query DGIdb
  cat("Querying DGIDB...")
  use_genes<-paste(genes, collapse = ",")
  url <- "https://dgidb.org/api/v2/interactions.json"

  body <- list(genes = use_genes,
               interaction_sources = databases,
               gene_categories = categories,
               interaction_types = interactions)
  #source_trust_levels = trustLevel)
  body <- body[!sapply(body, is.null)]
  postRequest <- httr::POST(url = url, body = body, encode = 'multipart')
  text <- httr::content(postRequest, as = "text", encoding = "ISO-8859-1")
  if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
  if (identical(text, "")) stop("Query response was emtpy.")
  queryResult <- jsonlite::fromJSON(text, simplifyVector = TRUE)

  # queryResult <- queryDgidbPost(genes=genes,
  #                               interactionSources = databases,
  #                               geneCategories = categories,
  #                               interactionTypes = interactions)

  #,trustLevel = trustLevel)
  cat("done!\n")

  # Init result class: rDGIdbResult
  interactionList <- base::lapply(queryResult$matchedTerms$geneName,getResultSummary, queryResult$matchedTerms)
  mis_id<-which(unlist(base::lapply(interactionList,function(x){return(dim(x)[1])}))==0)
  if(length(mis_id)>0){
    msg_genes<-genes[mis_id]
    cat("No interaction drug was found for the following genes: ", msg_genes)
    tmp <- data.frame(base::do.call(rbind, interactionList[-mis_id]),
                      stringsAsFactors = FALSE)
  }else{
    tmp <- data.frame(base::do.call(rbind, interactionList),
                      stringsAsFactors = FALSE)
  }
  if(dim(tmp)[1]==0) {
    stop("No interaction gene drug pairs were found!")}else{
      tmp<-tmp[order(tmp$Source_count,decreasing=T),]
      rownames(tmp) <- 1:nrow(tmp)
      result <- list()
      result[["summary_result"]]<-tmp
      if (return_raw) result[["raw_result"]]<-queryResult
      return(result)
      # End of function queryDGIdb()
    }
}
#' @title queryDgidb_V5_gene
#' @description a wrapper to query DGIdb API with gene name. Use ghql GraphqlClient function to query DGIdb version 5 to find the interaction between drug and gene, information including: "interactionScore","interactionTypes","interactionAttributes" ,"publications","sources","drug.name","drug.conceptId","queryGene"
#' @importFrom  ghql GraphqlClient
#' @importFrom jsonlite fromJSON
#' @param genes character vector, the query gene names.
#' @export
#'
#' @examples
#' use_genes<-c("ATM","BCL2L1","LCK","CDK4")
#' res<-queryDgidb_V5_gene (genes=use_genes)
#'
queryDgidb_V5_gene <- function(genes) {
  dgidbV5_api<-"https://dgidb.org/api/graphql"
  out_lst<-lapply(genes,function(gene){
      ## Create a query class first
    conn <- ghql::GraphqlClient$new(url = dgidbV5_api)
    query_gene <- ghql::Query$new()
    #"BRAF","MAP2K2"
    query_gene$query('x', sprintf('{
      genes(names: ["%s"]) {
        nodes {
          interactions {
            drug {
              name
              conceptId
            }
            interactionScore
            interactionTypes {
              type
              directionality
            }
            interactionAttributes {
              name
              value
            }
            publications {
              pmid
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }',gene))
    ## Execute the query
    res <- conn$exec(query_gene$queries$x)
    # Convert the the output from raw to json format
    res <- jsonlite::fromJSON(res,flatten = TRUE,simplifyVector = TRUE)
    if(is.null(dim(res$data$genes$nodes$interactions[[1]]))){
      message(sprintf("Gene:%s not found in the database!",gene))}else{
        res_up<-getResultSummary_gene_V5(gene=gene,output=res$data$genes$nodes$interactions)
        return(res_up)}
        })
result<-do.call("rbind",out_lst[!is.null(out_lst)])
return(result)
}
#'
#' @title queryDGIdb_by_drug
#' @description a wrapper to query DGIdb API with drug name. Deprecated Use \code{queryDgidb_V5_drug()} instead.
#' @return a list contains 1.a data frame("summary_result") of summarized gene,drug and interaction information; 2.raw query data("raw_result") from DGIdb (optional).
#' @param drugs a character vector of drugs to be used to find targeting genes.
#' @param return_raw a logical value, default is FALSE, only return summarized information. If users want to check raw query data, please change to TRUE.
#' @noRd
#' @examples
#' res<-queryDGIdb_by_drug(drugs=c("dasatinib","venetoclax"))
#'
queryDGIdb_by_drug<-function(drugs=NULL, return_raw=F){
  .Deprecated("queryDGIdb_by_drug", "queryDgidb_V5_drug")
  url <- "https://dgidb.org/api/v2/interactions.json"
  body <- list(drugs = paste(drugs, collapse = ","))
  #source_trust_levels = trustLevel)
  body <- body[!sapply(body, is.null)]
  postRequest <- httr::POST(url = url, body = body, encode = 'multipart')
  text <- httr::content(postRequest, as = "text", encoding = "ISO-8859-1")
  if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
  if (identical(text, "")) stop("Query response was emtpy.")
  queryResult <- jsonlite::fromJSON(text, simplifyVector = TRUE)
  #queryResult<-queryDgidbPost_drug(drugs = drugs)
  interactionList <- lapply(queryResult$matchedTerms$drugName,getResultSummary_drug, queryResult$matchedTerms)
  tmp <- data.frame(do.call(rbind, interactionList),
                    stringsAsFactors = FALSE)
  tmp<-tmp[order(tmp$Source_count,decreasing=T),]
  rownames(tmp) <- 1:nrow(tmp)
  result <- list()
  result[["summary_result"]]<-tmp
  if (return_raw) result[["raw_result"]]<-queryResult
  return(result)
}
#' @title queryDgidb_V5_drug
#' @description a wrapper to query DGIdb API with drug name. Use ghql GraphqlClient function to query DGIdb version 5 to find the interaction between drug and gene. Informations including: "interactionScore","interactionTypes","interactionAttributes","publications", "sources","gene.name","gene.conceptId","gene.longName" ,"Drug".
#' @importFrom  ghql GraphqlClient
#' @importFrom jsonlite fromJSON
#' @param drugs character vector, the query drug names.
#' @export
#'
#' @examples
#' use_drugs<-c("Dasatinib","DOVITINIB")
#' res<-queryDgidb_V5_drug(drugs=use_drugs)
#'
queryDgidb_V5_drug <- function(drugs) {
  dgidbV5_api<-"https://dgidb.org/api/graphql"

  out_lst<-lapply(drugs,function(drug){
      ## Create a query class first
    conn <- GraphqlClient$new(url = dgidbV5_api)
    query_drug <- Query$new()
    #"Dasatinib","DOVITINIB"
    query_drug$query('x', sprintf('{
      drugs(names: ["%s"]) {
        nodes {
          interactions {
            gene {
              name
              conceptId
              longName
            }
            interactionScore
            interactionTypes {
              type
              directionality
            }
            interactionAttributes {
              name
              value
            }
            publications {
              pmid
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }',drug))
    ## Execute the query
    res <- conn$exec(query_drug$queries$x)
    # Convert the the output from raw to json format
    res <- jsonlite::fromJSON(res,flatten = TRUE,simplifyVector = TRUE)
    if(is.null(dim(res$data$drugs$nodes$interactions[[1]]))){
      message(sprintf("Drug:%s not found in the database!",drug))}else{
         res_up<-getResultSummary_drug_V5(drug=drug,output=res$data$drugs$nodes$interactions)
        return(res_up)}
        })
result<-do.call("rbind",out_lst[!is.null(out_lst)])
return(result)
}
#' @title queryDgidb_V5_drugAnnotation
#' @description Using ghql GraphqlClient function to query DGIdb version 5 to find the drug meta information, including: "name","conceptId","approved" ,"immunotherapy", "antiNeoplastic","drugAttributes" ,"drugApprovalRatings" ,"drugApplications".
#' @importFrom  ghql GraphqlClient
#' @importFrom jsonlite fromJSON
#' @param drugs character vector, the query drug names.
#' @export
#'
#' @examples
#' use_drugs<-c("Dasatinib","DOVITINIB","nondrug")
#' res<-queryDgidb_V5_drugAnnotation(drugs=use_drugs)
#'
queryDgidb_V5_drugAnnotation <- function(drugs) {
  dgidbV5_api<-"https://dgidb.org/api/graphql"
   out_lst<-lapply(drugs,function(drug){
      ## Create a query class first
    conn <- GraphqlClient$new(url = dgidbV5_api)
    query_drug <- Query$new()
    #"Dasatinib","DOVITINIB"
    query_drug$query('x', sprintf('{
    drugs(names: ["%s"]) {
      nodes {
        name
        conceptId
        approved
        immunotherapy
        antiNeoplastic
        drugAttributes {
          name
          value
        }
        drugApprovalRatings {
          rating
          source {
            sourceDbName
            sourceTrustLevel {
              level
            }
          }
        }
        drugApplications {
          appNo
        }
      }
    }
}',drug))
    ## Execute the query
    res <- conn$exec(query_drug$queries$x)
    # Convert the the output from raw to json format
    res <- jsonlite::fromJSON(res,flatten = TRUE,simplifyVector = TRUE)
    if(is.null(dim(res$data$drugs$nodes))){
      message(sprintf("Drug:%s not found in the database!",drug))}else{
         res_up<-res$data$drugs$nodes
        return(res_up)}
        })
result<-do.call("rbind",out_lst[!is.null(out_lst)])
return(result)
}
#######seed partner selection with drug database###########
#' @title Drugbank_drugphase_levels
#' @description the drug phases used in drugbank. Users can select sepecific drug phase
#' @examples
#' Drugbank_drugphase_levels()
#' @export
Drugbank_drugphase_levels <- function() {
  Drugbank_drugphase_levels <- c("approved", "experimental", "illicit", "investigational", "nutraceutical", "vet_approved", "withdrawn" )
  return(Drugbank_drugphase_levels)
}

#' @title Repurposing_HUB_drugphase_levels
#' @description the drug phases used in Repurposing HUB. Users can select sepecific drug phase
#' @examples
#' Repurposing_HUB_drugphase_levels()
#' @export
Repurposing_HUB_drugphase_levels <- function() {
  Repurposing_HUB_drugphase_levels <- c("Launched","Phase 1","Phase 1/Phase 2","Phase 2" ,"Phase 2/Phase 3","Phase 3","Preclinical","Withdrawn")
  return(Repurposing_HUB_drugphase_levels)
}

#' @title ChEMBL_drugphase_levels
#' @description the drug phases used in ChEMBL. Users can select sepecific drug phase
#' @examples
#' ChEMBL_drugphase_levels()
#' @export
ChEMBL_drugphase_levels <- function() {
  ChEMBL_drugphase_levels <- c("0","1","2","3","4")
  return(ChEMBL_drugphase_levels)
}

#' @title drugdatabase_names
#' @description the curated drug databases. Users can select sepecific drug database
#' @examples
#' drugdatabase_names()
#' @export
drugdatabase_names <- function() {
  drugdatabase_names <- c("BPdb","ChEMBL","cMAP","Drugbank_V5.1", "rDGIdb_1.8.0", "Repurposing_HUB")
  return(drugdatabase_names)
}

#' @title selectSP_driver_byDB
#' @description selectSP_driver_byDB is used to select driver genes based on their interaction drug information, including predicted blood-brain barrier penetration and drug developing phase.
#' @param genes a character vector of genes for which the targeting drug phase or blood brain barrier penetration information are queried.
#' @param select_type character, which drug information is required. User can choose from "BP_prediction" and "drug_phase".
#' @param database character, which database will be used. "BPdb" contains drug blood brain barrier penetrantion information, others contain drug developing phase information. User can choose from "BPdb", "ChEMBL","Repurposing_HUB" and "Drugbank_V5.1".
#' @param select_phase a character vector of drug phases will be used. For checking detailed phase information in each database, try to call helper functions: \code{ChEMBL_drugphase_levels()}; \code{Repurposing_HUB_drugphase_levels()}; \code{Drugbank_drugphase_levels()}.
#' @examples
#' drugdb.preload()
#' genes<-c("MDM2","CDK4","CDK6","ATM","ATR","BCL2","BCL2L1","HDAC6","STAT3","JAK2")
#' sp_df<-selectSP_driver_byDB(genes=genes,select_type="BP_prediction",database="BPdb")
#' sp_df<-selectSP_driver_byDB(genes=genes,select_type="drug_phase",database="Drugbank_V5.1")
#' sp_df<-selectSP_driver_byDB(genes=genes,select_type="drug_phase",database="Drugbank_V5.1",select_phase="approved")
#' sp_df<-selectSP_driver_byDB(genes=genes,select_type="drug_phase",database="ChEMBL")
#' sp_df<-selectSP_driver_byDB(genes=genes,select_type="drug_phase",database="Repurposing_HUB")
#' \dontrun{
#' }
#' @export
#'
selectSP_driver_byDB<-function(genes=NULL,
                               select_type=NULL,
                               database=NULL,
                               select_phase=NULL){
  if (missing(genes)) stop("Need to specify a vector of genes!")

  if (is.null(genes) || length(genes) == 0 ||
      !is.character(genes) || genes == "") {
    stop("Need to specify a non-empty vector of genes names.")
  }

  if (missing(select_type) || length(intersect(select_type,c("BP_prediction","drug_phase")))==0) {stop("Need to specify a select_type.")}

  check_type <- match.arg(arg = select_type,
                          choices = c("BP_prediction","drug_phase"),
                          several.ok = FALSE)

  if (missing(database) || length(intersect(database,c("ChEMBL","Repurposing_HUB","Drugbank_V5.1","BPdb")))==0) {stop("Need to choose one of the following databases: ChEMBL,Repurposing_HUB, Drugbank_V5.1, BPdb.") }

  use_database <- match.arg(arg = database,
                            choices = c("ChEMBL","Repurposing_HUB","Drugbank_V5.1","BPdb"),
                            several.ok = FALSE)


  if (use_database=="BPdb"&&check_type=="drug_phase" ) {
    stop("BPdb has NO drug phase information. Please change database or select_type!")}

  if(use_database=="ChEMBL"){
    if(is.null(select_phase)) use_phase<-ChEMBL_drugphase_levels() else
      use_phase <- match.arg(arg = select_phase,
                             choices = ChEMBL_drugphase_levels(),
                             several.ok = TRUE)}

  if(use_database=="Drugbank_V5.1"){
    if(is.null(select_phase)) use_phase<-Drugbank_drugphase_levels() else
      use_phase <- match.arg(arg = select_phase,
                             choices = Drugbank_drugphase_levels(),
                             several.ok = TRUE)}

  if(use_database=="Repurposing_HUB"){
    if(is.null(select_phase)) use_phase<-Repurposing_HUB_drugphase_levels() else
      use_phase <- match.arg(arg = select_phase,
                             choices = Repurposing_HUB_drugphase_levels(),
                             several.ok = TRUE)
  }

  if(check_type=="BP_prediction"){
    df<-SINBA_drug_obj$SINBA_drug_gene$BPdb
    df<-df[df$geneSymbol%in%genes,]
    out_df<-df[,c("geneSymbol","BBB","Drug_Name","ChEMBL","Drug_type")]
  }else{
    df<-SINBA_drug_obj$gene_with_drug_phase
    i<-intersect(which(df$geneSymbol%in%genes),which(df$Source==use_database))
    i<-intersect(i,which(df$Drug_group%in%use_phase))
    out_df<-df[i,c("geneSymbol","Drug_group","Drug_ID")]
  }
  rownames(out_df)<-NULL
  return(out_df)
}
##################preload drug database#####################
##
#
#' Preload drug database files into R workspace for SINBA
#'
#' \code{drugdb.preload} is a pre-processing function for SINBA. It preloads needed drug databases into R workspace,and saves it locally under data/ directory.
#'
#' Users need to choose the drug databse sources(e.g. "BPdb","ChEMBL","cMAP","Drugbank_V5.1", "rDGIdb_1.8.0","Repurposing_HUB")
#'
#' @param db  A vector of characters, users can choose from "BPdb","ChEMBL","cMAP","Drugbank_V5.1", "rDGIdb_1.8.0" and "Repurposing_HUB". If NULL, all the drug databases will be used. Default is NULL.
#' @param main.dir character, the main directory for SINBA package.
#' If NULL, will be \code{system.file(package = "SINBA")}. Default is NULL.
#' @param out.dir charater, the directory to save drug database information. If NULL, will be \code{system.file(package = "SINBA")}. Default is NULL.
#' @examples
#' drugdb.preload()
#'
#' \dontrun{
#' drugdb.preload()
#' }
#' @export
drugdb.preload <- function(db=NULL,main.dir=NULL,out.dir=NULL){
  if(is.null(db)==TRUE){
    use_db <- c("BPdb","ChEMBL","cMAP","Drugbank_V5.1", "rDGIdb_1.8.0","Repurposing_HUB")
    message(sprintf('db not set, will use all the available databases: %s!',base::paste(use_db, collapse = ",")))
  }else{
    use_db <- base::match.arg(arg = db,
                              choices = drugdatabase_names(),
                              several.ok = TRUE)
    message(sprintf('The following drug databases will be used: %s!',base::paste(use_db, collapse = ",")))
  }
  if(is.null(main.dir)==TRUE){
    main.dir <- sprintf("%s/data",system.file(package = "SINBA"))
    message(sprintf('main.dir not set, will use package directory: %s',main.dir))
  }
  if(is.null(out.dir)==TRUE){
    out.dir <- sprintf("%s/data",system.file(package = "SINBA"))
    message(sprintf('out.dir not set, will use package directory: %s',out.dir))
  }
  #get the whole drug databases
  #obj<-readRDS(sprintf("%s/drug_gene_database_2021-08-17.rds",main.dir))
  load(sprintf("%s/Drugdb_SINBA.rda",main.dir))
  SINBA_drug_obj<-list()
  tmp<-Drugdb_SINBA[[1]]
  SINBA_drug_obj[[1]]<-tmp[which(tmp$Source%in%use_db),]
  SINBA_drug_obj[[2]]<-list()
  tmp<-Drugdb_SINBA[[2]]
  SINBA_drug_obj[[2]]<-tmp[use_db]
  tmp<-Drugdb_SINBA[[3]]
  SINBA_drug_obj[[3]]<-tmp[which(tmp$Source%in%use_db),]
  names(SINBA_drug_obj)<-names(Drugdb_SINBA);rm(Drugdb_SINBA,tmp)

  RData.file <- sprintf('%s/Drugdb_SINBA.RData', out.dir)
  save(SINBA_drug_obj,file = RData.file)

  load(RData.file,.GlobalEnv)
  return(TRUE)
}
###########################################################
