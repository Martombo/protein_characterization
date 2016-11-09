library(topGO)

meta_analysis = TRUE
fisher = TRUE
avg_genes = FALSE
ontology = "MF"
geneID2GO = readMappings("../annotation_uniprot")
all = read.table("../dati")

Fisher.test <- function(p) {
	    Xsq <- -2*sum(log(p))
		    p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
			    return(p.val)
}

psumunif = function(x,n) {
	fun = function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)
	return(1/factorial(n) * sum(sapply(0:n, fun)))
}

run_topGO = function(gene_set, gene_sel, n_filter){
	GO = new("topGOdata", ontology=ontology, allGenes=gene_set, geneSel=gene_sel, annot=annFUN.gene2GO, gene2GO=geneID2GO)
	result_W01_fish = runTest(GO, algorithm = "weight01", statistic = "fisher")
	to_return = GenTable(GO, pval=result_W01_fish, topNodes=length(result_W01_fish@score))
	colnames(to_return) = paste0(colnames(to_return), as.character(n_filter))
	colnames(to_return)[1] = "GO.ID"
	return(to_return)
}

# analysis based on FC, all + 2 thresholds
nonnull_fc = all$V2[all$V2>0]
nonnull_fc = nonnull_fc[order(nonnull_fc, decreasing=T)]
n_nonnull = length(nonnull_fc)
gene_sel0 = function(x){return(x>0)}
threshold1 = nonnull_fc[round((n_nonnull*2)/3)]
gene_sel1 = function(x){return(x>threshold1)}
threshold2 = nonnull_fc[round((n_nonnull*1)/3)]
gene_sel2 = function(x){return(x>threshold2)}
gene_set = all$V2
names(gene_set) = all$V1
top_go0 = run_topGO(gene_set, gene_sel0, 0)
if(!meta_analysis){
	write.table(top_go0, file="topGO_results0", sep="\t", quote=F, row.names=F)
	stop()
}
top_go1 = run_topGO(gene_set, gene_sel1, 1)
top_go2 = run_topGO(gene_set, gene_sel2, 2)

# analysis based on p-value, 2 thresholds
nonnull_pval = all$V3[all$V3>0]
nonnull_pval = nonnull_pval[order(nonnull_pval, decreasing=T)]
n_nonnull = length(nonnull_pval)
gene_sel3 = function(x){return(x>nonnull_pval[round((n_nonnull*2)/3)])}
gene_sel4 = function(x){return(x>nonnull_pval[round((n_nonnull*1)/3)])}
gene_set = all$V3
names(gene_set) = all$V1
top_go3 = run_topGO(gene_set, gene_sel3, 3)
top_go4 = run_topGO(gene_set, gene_sel4, 4)

merged = merge(top_go0, top_go1, by="GO.ID")
merged = merge(merged, top_go2, by="GO.ID")
merged = merge(merged, top_go3, by="GO.ID")
merged = merge(merged, top_go4, by="GO.ID")
combined = data.frame(GO.ID=merged$GO.ID, Term=merged$Term0, Annotated=merged$Annotated0)
if(avg_genes){
	combined$Significant = apply(merged[,paste0("Significant",c(0,0,1,2,3,4))], 1, function(x) mean(as.double(x)))
	combined$Expected = apply(merged[,paste0("Expected",c(0,0,1,2,3,4))], 1, function(x) mean(as.double(x)))
}else{
	combined$Significant = merged$Significant0
	combined$Expected = merged$Expected0
}
if(fisher){
	combined$pval = apply(merged[,paste0("pval",c(0,0,1,2,3,4))], 1, function(x) Fisher.test(as.double(x)))
}else{
	combined$pval = apply(merged[,paste0("pval",c(0,0,1,2,3,4))], 1, function(x) psumunif(as.double(x), 6))
}
combined$pval[is.na(combined$pval)] = 1e-100
combined = combined[order(combined$pval, decreasing=F),]

write.table(combined, file="topGO_results", sep="\t", quote=F, row.names=F)

save.image()
