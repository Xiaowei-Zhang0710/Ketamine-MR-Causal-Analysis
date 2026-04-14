suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/Users/nullyann/Desktop/anding_smr/data/gwas_category_outcome"  # 含 Ast/Microglia/Oligodendrocytes/Inhibitory 子目录
GWAS_FILE <- "/Users/nullyann/Desktop/anding_smr/data/gwas/GCST90726344.tsv"
BIM_FILE  <- "/Users/nullyann/Desktop/anding_smr/data/ldref/1000G_EUR/EUR.bim"

OUT_FULL <- file.path(ROOT, "gene_ct_full.tsv")
OUT_MIN  <- file.path(ROOT, "gene_ct.tsv")

cell_order <- c("Ast", "Microglia", "oligo", "Inhibitory")

msg <- function(...) cat(sprintf(...), "\n")

# --- rsID->chr,pos 映射：优先GWAS，缺失再用BIM ---
gwas <- fread(GWAS_FILE)
setnames(gwas, c("chromosome","base_pair_location","rs_id","p_value"),
         c("chr","pos","SNP","pval.outcome"), skip_absent = TRUE)
gwas_map <- unique(gwas[!is.na(SNP) & SNP!="", .(SNP=as.character(SNP), chr=as.integer(chr), pos=as.integer(pos), pval.outcome=as.numeric(pval.outcome))], by="SNP")

bim_map <- NULL
if (file.exists(BIM_FILE)) {
  bim <- fread(BIM_FILE, header=FALSE)
  bim_map <- unique(bim[!is.na(V2) & V2!="", .(SNP=as.character(V2), chr=as.integer(V1), pos=as.integer(V4))], by="SNP")
  bim_map <- bim_map[chr %in% 1:22]
}

coord_map <- copy(gwas_map)
setkey(coord_map, SNP)
if (!is.null(bim_map)) {
  missing <- setdiff(bim_map$SNP, coord_map$SNP)
  if (length(missing)>0) coord_map <- rbind(coord_map, bim_map[SNP %in% missing], fill=TRUE)
}
setkey(coord_map, SNP)

# --- 扫描每个celltype目录的 *_harmonized_data.rds ---
res <- list()

for (ct in cell_order) {
  d <- file.path(ROOT, ct)
  if (!dir.exists(d)) next
  rds_files <- list.files(d, pattern="_harmonized_data\\.rds$", full.names=TRUE)
  if (length(rds_files)==0) next
  msg("[Celltype] %s rds=%d", ct, length(rds_files))

  for (f in rds_files) {
    dat <- as.data.table(readRDS(f))

    if (!("SNP" %in% names(dat)) || !("exposure" %in% names(dat))) next
    gene <- as.character(dat$exposure[1])

    # 你这类rds通常有 pval.outcome / pval.exposure（没有就NA）
    if (!("pval.outcome" %in% names(dat))) dat[, pval.outcome := NA_real_]
    if (!("pval.exposure" %in% names(dat))) dat[, pval.exposure := NA_real_]

    # 选 lead SNP：优先用 outcome p 最小；如果全NA，再用 exposure p 最小；再不行就随便取第一条
    tmp <- dat[!is.na(SNP)]
    if (nrow(tmp)==0) next

    if (any(!is.na(tmp$pval.outcome))) {
      lead <- tmp[which.min(pval.outcome)]
      lead_type <- "GWAS_lead"
    } else if (any(!is.na(tmp$pval.exposure))) {
      lead <- tmp[which.min(pval.exposure)]
      lead_type <- "eQTL_lead"
    } else {
      lead <- tmp[1]
      lead_type <- "first"
    }

    lead_snp <- as.character(lead$SNP[1])
    lead_p_out <- suppressWarnings(as.numeric(lead$pval.outcome[1]))
    lead_p_exp <- suppressWarnings(as.numeric(lead$pval.exposure[1]))

    # 映射到 chr/pos
    mp <- coord_map[J(lead_snp), nomatch=0]
    if (nrow(mp)==0) next
    chr <- mp$chr[1]; pos <- mp$pos[1]

    res[[length(res)+1]] <- data.table(
      gene=gene,
      celltype=ct,
      SNP=lead_snp,
      chr=chr,
      pos=pos,
      pval_outcome=lead_p_out,
      pval_exposure=lead_p_exp,
      lead_source=lead_type
    )
  }
}

dt <- rbindlist(res, fill=TRUE)
if (nrow(dt)==0) stop("No entries produced; check paths and RDS content.")

# 去重：同gene+celltype保留（若重复，保留outcome p最小）
setorder(dt, gene, celltype, pval_outcome)
dt <- dt[!duplicated(dt, by=c("gene","celltype"))]

fwrite(dt, OUT_FULL, sep="\t")
fwrite(dt[, .(gene,celltype,chr,pos)], OUT_MIN, sep="\t")

msg("[DONE] gene_ct_full.tsv -> %s (rows=%d)", OUT_FULL, nrow(dt))
msg("[DONE] gene_ct.tsv      -> %s", OUT_MIN)