# =====================================================================
#  MEGENA Construction in R
#  Goal: Constructing a co-expression network based on gene expression matrix.
#
#  Input:
#    - cultivar_gene_matrix.txt : table with columns [gene,S001-4DPA,S002-4DPA,S003-4DPA, ...,S400-20DPA]
#                where the prefix indicates the sample name and the suffix indicates the sampling stage.
#
#  Output:
#    - MEGENA_ijw_cultivar.txt # Correlations between genes
#    - MEGENA_el_cultivar.txt # Gene correlations included in the PFN plane
#    - MEGENA_module_summary_detail_cultivar.txt # Features of each gene in the co-expression network
#    - summary_overall_cultivar.txt # Features of each module in the co-expression network
#    - cultivar_*.txt # Details of each module (edges and nodes)
#
#  Dependencies:
#    - MEGENA
# =====================================================================

## --- Step 1: Parameter initialization
n.cores <- 2 
doPar <-TRUE 
method = "pearson" 
FDR.cutoff = 0.05 
module.pval = 0.05 
hub.pval = 0.05 
cor.perm = 10 
hub.perm = 100 
annot.table=NULL
id.col = 1
symbol.col= 2

## --- Step 2: Correlation calculation
data=read.table("cultivar_gene_matrix.txt",header=T,row.name=1,stringsAsFactors = F,check.names = F)
rho.out = calculate.rho.signed(data,n.perm = cor.perm,FDR.cutoff = FDR.cutoff,estimator = method,
                               use.obs = "na.or.complete",direction = "absolute",
                               rho.thresh = NULL,sort.el = TRUE)
#head(rho.out$signif.ijw)
write.table(rho.out$signif.ijw,"MEGENA_ijw_cultivar.txt",sep = '\t', row.names = FALSE, quote = FALSE)

## --- Step 3: PFN construction
signif.ijw=read.table("MEGENA_ijw_cultivar.txt",header=T,stringsAsFactors = F,check.names = F)
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}
el <- calculate.PFN(signif.ijw,doPar = doPar,num.cores = n.cores,keep.track = FALSE)
g <- graph.data.frame(el,directed = FALSE)
write.table(el,"MEGENA_el_cultivar.txt",sep = '\t', row.names = FALSE, quote = FALSE)

## --- Step 4: Gene clustering
MEGENA.output <- do.MEGENA(g,mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                           save.output = FALSE)
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## --- Step 5: Data aggregation
summary.output <- MEGENA.ModuleSummary(MEGENA.output,mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                       min.size = 10,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)
if (!is.null(annot.table))
{
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}
module.output <- module_convert_to_table(MEGENA.output,mod.pval = 0.05,
                                         hub.pval = 0.05,min.size = 10,max.size=vcount(g)/2)
write.table(module.output,"MEGENA_module_summary_detail_cultivar.txt",sep = '\t', row.names = FALSE, quote = FALSE)
#print(summary.output$modules)
#print(summary.output$module.table)
write.table(summary.output$module.table,"summary_overall_cultivar.txt",sep = '\t', row.names = FALSE, quote = FALSE)
for (i in 1:length(summary.output$module.table$module.id)) {
    sub_cluster <- summary.output$modules[[i]]
    df1 = el[el[,1] %in% sub_cluster,]
    df2 = df1[df1[,2] %in% sub_cluster,]
    write.table(df2,paste0("cultivar_",summary.output$module.table$module.id[i],".txt"),sep = '\t', row.names = FALSE, quote = FALSE)
}