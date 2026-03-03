# =====================================================================
#  aFETp Calculation in R
#  Goal: Calculating similarity based on FET matrix.
#
#  Input:
#    - ./fisher_matrix/*.txt : FET matrix
#
#  Output:
#    - summary_aFETp_table.txt # aFETp between modules of semi-wild cotton and cultivated cotton.
# =====================================================================

semi_module=vector()
cultivar_module=vector()
pval=vector()
file_names <- list.files(path = "fisher_matrix")
i=1
for (file_name in file_names) {
    contingency_table <- read.table(file_name))
    contingency_matrix <- matrix(contingency_table$V1,nrow=2,ncol=2,
                                 dimnames=list(c("in_semi","no_semi"),c("in_cultivar","no_cultivar")))
    result <- fisher.test(contingency_matrix)
    semi_module[i]=paste(strsplit(file_name,"_")[[1]][1:2],collapse="_")
    cultivar_module[i]=sub(".txt","",paste(strsplit(file_name,"_")[[1]][3:4],collapse="_"))
	##Remove module pairs with no intersection
    if (contingency_table[1,1]==0) {
      pval[i]=1
    }else {
        pval[i]=result$p.value
	}
    i=i+1
}

output=data.frame(semi_module=semi_module,cultivar_module=cultivar_module,P=pval)
bh_corrected <- p.adjust(output$P, method = "BH")
data1=data.frame(output,bh_corrected)
colnames(data1)=c("cultivar_module","semi_module","P","p.adjust")
write.csv(data1,"summary_aFETp_table.txt")
