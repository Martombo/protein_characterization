library(topGO)

geneID2GO <- readMappings('annotation_uniprot')
all = read.table("dati")
gene_set = all$V2
names(gene_set) = all$V1
geneSel = function(x){return(x>0)}
GO = new("topGOdata", ontology="BP", allGenes=gene_set, geneSel=geneSel, annot=annFUN.gene2GO, gene2GO=geneID2GO)
result_W01_fish = runTest(GO, algorithm = "weight01", statistic = "fisher")
top_go = GenTable(GO, pval=result_W01_fish, topNodes=10000)

write.table(top_go, file="top_w01F", sep="\t", quote=F, row.names=F)

getSigGenesGO=function(go){
	SigGenes = sigGenes(GO)
	genes = genesInTerm(GO,whichGO=go)
	return(genes[[1]][genes[[1]] %in% sigGenes])
}

save.image()

#goGenes = getSigGenesGO("GO:0007155")
#write.table(goGenes,"adhesion",quote=F,col.names=F,row.names=F)
#goGenes = getSigGenesGO("GO:0007158")
#write.table(goGenes,"neuron",quote=F,col.names=F,row.names=F)
#goGenes = getSigGenesGO("GO:0030198")
#write.table(goGenes,"extracell",quote=F,col.names=F,row.names=F)
#goGenes = getSigGenesGO("GO:0007595")
#write.table(goGenes,"lactation",quote=F,col.names=F,row.names=F)
