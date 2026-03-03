# =====================================================================
# Eigengenes Calculation in R
#  Goal: Calculating the relationship between co-expression network modules and phenotypes.
#
#  Input:
#    - cultivar_gene_matrix.txt : table with columns [gene,S001-4DPA,S002-4DPA,S003-4DPA, ...,S400-20DPA]
#                where the prefix indicates the sample name and the suffix indicates the sampling stage.
#    - geneid_label.txt : gene-module membership in the co-expression network
#    - new_phenotype.txt : table with columns [smple_DPA,UHML,FS,Enlongation,UI]
#
#  Output:
#    - cultivar_correlation_matrix.csv # Correlations between genes
#
#  Dependencies:
#    - MEGENA impute cowplot patchwork tidyr
# =====================================================================

## --- Step 1: function initialization
ModuleHubProfileSingle = function(datexpr, no_pcs=10) {
  ngenes = dim(datexpr)[1]
  print("PC by ModuleHubProfileSingle .................. ")
  corrlmatrix = cor( t(datexpr), use = "pairwise.complete.obs")
  corrlmatrix = abs(corrlmatrix)
  diag(corrlmatrix)<- 0
  kin <- apply(corrlmatrix,2,sum, na.rm=TRUE) 
  orderK  = order(-kin)
  step = as.integer(ngenes/no_pcs)
  if(step<1){step=1; }
  selIdx= seq(from=1,to=ngenes, by=step)
  v = t(datexpr[selIdx[1:no_pcs], ])
  d = rep(0, no_pcs)
  ret   = NULL
  ret$v = v
  ret$d = d
  return(ret)
}

ModulePrinComps = function(datexpr,modules, min_modulesize=10) {
  couleur = modules
  datExpr = t(datexpr)
  
  no.pcs=10
  
  allmean = mean(datexpr, na.rm=T)
  modlevels <- names(couleur)
  
  listPCs = as.list(rep(NULL, no.pcs+1) )
  for (i in c(1:no.pcs) ){
    listPCs[[i]] = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }
  listPCs[[no.pcs+1]]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels)))
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels
  
  for(i in c(1:length(modlevels)) ){
    modulename = modlevels[i]
    restrict1 = colnames(datexpr) %in% couleur[[i]]
    restrict1= ifelse(is.na(restrict1), FALSE, restrict1)
    if(modlevels[i]=="grey"){next}
	
    datModule=t(datexpr[, restrict1])
    
    # check whether some samples have missing rate (missing value in column) >75%
    naDat = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff  = 0.75*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
      print("patch samples with missing rate >=0.75")
      naSampleIdx = c(1:length(nasBySample))[selSamples]
      for(x in naSampleIdx){
        #print(paste("Sample Idx=", x) )
        datModule[,x]=ifelse(is.na( datModule[,x]), allmean, datModule[,x])
      }
    }
  
    if(sum(restrict1)<min_modulesize){
      listPCs[[1]][,i] = datModule[1,]
      next
    }
    
    imputed=impute.knn(as.matrix(datModule))
    datModule =imputed$data
    
    datModule=t(scale(t(datModule)))
    
    #calculating PC
    svd1=tryCatch(svd(datModule), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)
    
    mtitle=paste("PCs of ", modulename," module", sep="")
    
    no.samples = dim(datModule)[2]
    
    actualpcs=min(dim(svd1$v)[2], no.pcs)

    listPCs[[ no.pcs+1] ][,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)
    
    # this is the first principal component, the ith column of the j-th element in the list
    for (j in c(1:actualpcs) ){       
      pcj=svd1$v[,j]
      
      # detect NAs
      jnas = is.na(pcj)
      if (sum(jnas)>0){
        break;
      }
      
      signhj=sign(sum(cor(pcj,  t(datModule))))
      
      if( !(signhj == 0 | is.na(signhj)) )  pcj=signhj* pcj
      listPCs[[j]][,i] = pcj
    }
    
  }
  
  PC.out = listPCs
  
  PC.res <- data.frame()
  l <- length(PC.out)
  for (i in 1:(l-1))
  {
    df <- data.frame(type = rep(paste("PC",i,sep = ""),nrow(PC.out[[i]])),id = colnames(datExpr),
                     as.data.frame(PC.out[[i]]))
    PC.res <- rbind.data.frame(PC.res,df)
    rownames(PC.out[[i]]) <- colnames(datExpr)
    rm(df)
  }
  
  PC.res <- rbind.data.frame(PC.res,
                             data.frame(type = rep("variance.explained",nrow(PC.out[[l]])),id = paste("PC",1:nrow(PC.out[[l]]),sep = ""),as.data.frame(PC.out[[l]])))
  rownames(PC.out[[l]]) <- paste("PC",1:nrow(PC.out[[l]]),sep = "")
  
  names(PC.out) <- c(paste("PC",1:(l-1),sep = ""),"variance.explained")
  
  return(PC.out)
  
}

## --- Step 2: eigengenes calculation
cultivar_all=read.table("cultivar_gene_matrix.txt",header=T,row.names=1)
datexpr=t(cultivar_all)

selected_4DPA <- datexpr[grep("4DPA", rownames(datexpr)),]
selected_8DPA <- datexpr[grep("8DPA", rownames(datexpr)),]
selected_12DPA <- datexpr[grep("12DPA", rownames(datexpr)),]
selected_16DPA <- datexpr[grep("16DPA", rownames(datexpr)),]
selected_20DPA <- datexpr[grep("20DPA", rownames(datexpr)),]

module_data <- read.table("geneid_label.txt")
colnames(module_data)=c("gene","label")
modules <- split(module_data$gene, module_data$label)

result_4DPA <- ModulePrinComps(datexpr = selected_4DPA, modules=modules)
result_8DPA <- ModulePrinComps(datexpr = selected_8DPA, modules=modules)
result_12DPA <- ModulePrinComps(datexpr = selected_12DPA, modules=modules)
result_16DPA <- ModulePrinComps(datexpr = selected_16DPA, modules=modules)
result_20DPA <- ModulePrinComps(datexpr = selected_20DPA, modules=modules)

## --- Step 3: correlation calculation
pheno=read.table("new_phenotype.txt",header=T,row.names=1)
row_names_pheno=rownames(pheno)
colnames(pheno)=c("FL","FS","FE","FU")

for (i in c("4DPA","8DPA","12DPA","16DPA","20DPA")) {
  eigengenes=get(paste0("result_",i))[[1]]
  # order
  row_names_eigengenes=rownames(eigengenes)
  both_sample <- intersect(row_names_eigengenes,row_names_pheno)
  
  filtered_pheno <- pheno[row.names(pheno) %in% both_sample, ]
  filtered_eigengene <- eigengenes[row.names(eigengenes) %in% both_sample, ]
  
  aligned_eigengene <- filtered_eigengene[order(rownames(filtered_eigengene)), ]
  aligned_phenotype <- filtered_pheno[order(rownames(filtered_pheno)), ]
  
  n_phenotypes=dim(aligned_phenotype)[2]
  n_eigengenes=dim(aligned_eigengene)[2]
  
  # pearson
  assign(paste0("correlation_matrix_",i),cor(aligned_eigengene, aligned_phenotype, method = "pearson"))
  p_value_matrix=matrix(ncol=n_phenotypes,nrow=n_eigengenes)
  
  #BH
  for (j in 1:n_phenotypes) {
    p_values <- numeric(n_eigengenes)
    for (x in 1:n_eigengenes) {
      cor_result <- cor.test(aligned_eigengene[, x], aligned_phenotype[, j], method = "pearson")
      p_values[x] <- cor_result$p.value
    }
    adjusted_p_values <- p.adjust(p_values, method = "BH")
    p_value_matrix[,j]=adjusted_p_values
  }
  
  colnames(p_value_matrix)=colnames(aligned_phenotype)
  rownames(p_value_matrix)=colnames(aligned_eigengene)
  result <- -log10(p_value_matrix)
  result=as.data.frame(result)
  result$moduleID=rownames(result)
  long_data_p <- pivot_longer(result, cols=-"moduleID",
                              names_to = "pheno", 
                              values_to = "p.adjust")
							  
  df=as.data.frame(get(paste0("correlation_matrix_",i)))
  df$moduleID=rownames(df)
  long_data_cor <- pivot_longer(df, cols=-"moduleID",
                                names_to = "pheno", 
                                values_to = "correlation")
								
  merged_df <- merge(long_data_cor, long_data_p, by = c("moduleID","pheno"))
  merged_df$DPA=i
  new_order <- c("moduleID", "pheno", "correlation", "p.adjust","DPA")
  ordered_df <- merged_df[, new_order]
  assign(paste0("merged_df_",i),ordered_df)
}

combined_df <- rbind(merged_df_4DPA, merged_df_8DPA, merged_df_12DPA, merged_df_16DPA, merged_df_20DPA)  
write.csv(combined_df,"cultivar_correlation_matrix.csv",quote = FALSE, row.names = FALSE)
}