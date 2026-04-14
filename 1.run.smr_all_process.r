library(TwoSampleMR)
#Rscript eqtl_tob37.r
#gwas category outcome
#需要修改对应的eqtl和DEGs的sheet名
#setwd("/Users/nullyann/Desktop/anding_smr/data/gwas_category_outcome/Microglia")
#setwd("/Users/nullyann/Desktop/anding_smr/data/gwas_category_outcome/Ast")
setwd("/Users/nullyann/Desktop/anding_smr/data/gwas_category_outcome/oligo")
library(openxlsx)
deg<-read.xlsx("/Users/nullyann/Desktop/anding_smr/关键的DEGs.xlsx",sheet=4,colNames=F) #Oligodendrocytes sheet4
#deg<-read.xlsx("/Users/nullyann/Desktop/anding_smr/关键的DEGs.xlsx",sheet=1,colNames=F) #Inhibitory sheet1
#deg<-read.xlsx("/Users/nullyann/Desktop/anding_smr/关键的DEGs.xlsx",sheet=3,colNames=F) #Ast sheet3
#deg<-read.xlsx("/Users/nullyann/Desktop/anding_smr/关键的DEGs.xlsx",sheet=2,colNames=F)  #sheet=2 microglia

gene_list<-toupper(deg$X1)

#读取exposure
path="/Users/nullyann/Desktop/anding_smr/data/eqtl_b37/Oligodendrocytes/"
fhlist<-list.files(path=path,pattern="controlsonly_associations_b37.csv")
#fhlist<-list.files(path="/Users/nullyann/Desktop/anding_smr/data/eqtl_b37/Microglia/",pattern="controlsonly_associations_b37.csv")
for( i in fhlist){
tmp<-read.csv(paste0(path,i))
a<-subset(tmp,p.value < 1e-05) # 5e-08,5e-06,1e-05
if( i == fhlist[1]){b<-a}
else{b<-rbind(b,a)}
}
b$se<-abs(b$beta/b$t.stat)
write.csv(b,file="exposure.csv")
table(gene_list %in% b$gene)
gin<-gene_list[gene_list %in% b$gene]
b<-subset(b,gene %in% gin)
gin

for( gs in gin){ #Ast LRP1B，RIMS1

write.csv(subset(b,gene==gs),file="exposure.csv",row.names=F)

eq<-"exposure.csv"
eq_clumped<-read_exposure_data(filename=eq,sep=",",snp_col="SNP", phenotype_col = "gene", 
beta_col="beta",se_col="se",effect_allele_col="effect_allele",other_allele="other_allele",clump=FALSE)
write.csv(eq_clumped,file="exposure_clumped.csv",row.names=F)
exposure_dat<-eq_clumped

c<-read.table("/Users/nullyann/Desktop/anding_smr/data/gwas/GCST90726344.tsv",header=T)
d<-merge(exposure_dat,c,by.x="SNP",by.y="rs_id")
d<-d[,c("SNP","effect_allele","other_allele","beta","standard_error","p_value","n","effect_allele_frequency")]
colnames(d)<-c("SNP","effect_allele","other_allele","beta","se","p","n","eaf")
write.csv(d,file="outcome.csv")

outcome_dat<-read_outcome_data(
    snps=eq_clumped$SNP,
    filename="outcome.csv",
    sep=",",eaf_col="eaf",
    snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",other_allele_col="other_allele",pval_col="p"
    )

dat<-harmonise_data(exposure_dat=exposure_dat,outcome_dat=outcome_dat)
write.csv(dat,file="harmonized_data.csv")

mr(dat)
write.xlsx(mr(dat),file=paste0(gs,"_mr.xlsx"),rowNames=T)
generate_odds_ratios(mr_res=mr(dat)) 

pdf(paste0(gs,"_mr_scatter_plot.pdf"),width=6,height=6)
plot(mr_scatter_plot(mr_results=mr(dat),dat))
dev.off()

pdf(paste0(gs,"_mr_forest_plot.pdf"),width=6,height=6)
plot(mr_leaveoneout_plot(leaveoneout_results=mr_leaveoneout(dat)))
dev.off()

write.xlsx(mr_pleiotropy_test(dat),file=paste0(gs,"_pleiotropy_test.xlsx"),rowNames=T)
saveRDS(dat,file=paste0(gs,"_harmonized_data.rds"))
}


library(biomaRt)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(GenomicRanges)

gene_symbol <- "GRM7"
rsids <- c("rs11131052", "rs4686092", "rs62244614", "rs6781770")

# ---- (A) gene symbol -> Entrez ID ----
eg <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol,
                            keytype = "SYMBOL", columns = "ENTREZID")
entrez <- unique(na.omit(eg$ENTREZID))[1]
stopifnot(!is.na(entrez))

# ---- (B) gene range (hg19) ----
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr_gene <- genes(txdb)[entrez]  # GRanges length 1
chr <- as.character(seqnames(gr_gene))
g_start <- start(gr_gene)
g_end   <- end(gr_gene)

pad <- 2e4
from <- max(1, g_start - pad)
to   <- g_end + pad

# ---- (C) rsID -> genomic position (GRCh37) ----
# 返回 GRanges（可能有些 rsID 查不到；或有多条映射）
snp_gr <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsids)
snp_gr$RSID <- names(snp_gr)

# 过滤到该基因所在 chr 且位于展示区间
snp_gr <- snp_gr[as.character(seqnames(snp_gr)) == chr]
snp_gr <- snp_gr[start(snp_gr) >= from & start(snp_gr) <= to]

# ---- (D) Tracks ----
itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
gtrack <- GenomeAxisTrack()
txtrack <- GeneRegionTrack(txdb, genome = "hg19", chromosome = chr,
                           start = from, end = to,
                           name = paste0(gene_symbol, " model"),
                           transcriptAnnotation = "symbol")
dtrack <- DataTrack(range = snp_gr, genome = "hg19", chromosome = chr,
                    name = "rsID", type = "p",  # points
                    groups = snp_gr$RSID,
                    legend = TRUE)

plotTracks(list(itrack, gtrack, txtrack, dtrack),
           from = from, to = to,
           main = paste0(gene_symbol, " (hg19/GRCh37) rsID positions"))


