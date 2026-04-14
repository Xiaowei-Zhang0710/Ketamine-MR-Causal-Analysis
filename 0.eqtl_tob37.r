suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ============== CONFIG ==============
#CELLTYPE <- "Astrocytes"
#CELLTYPE<-"Microglia"
#CELLTYPE<-"Inhibitory_Neurons"
CELLTYPE<-"Oligodendrocytes"

BASE <- "/Users/nullyann/Desktop/anding_smr/data"
EQTL_DIR <- file.path(BASE, "eqtl", CELLTYPE)
OUT_DIR  <- file.path(BASE, "eqtl_b37", CELLTYPE)

BCFTOOLS <- "/opt/homebrew/bin/bcftools"
DBSNP_VCF <- file.path(BASE, "ref", "dbsnp", "dbsnp.v153.b37.vcf.gz")

PATTERN <- paste0("^", CELLTYPE, "_chr[0-9XYM]+_controlsonly_associations\\.csv$")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============== Helpers ==============
msg <- function(...) cat(sprintf(...), "\n")

assert_file <- function(x, what="file") {
  if (!file.exists(x)) stop(what, " not found: ", x)
}

is_acgt <- function(x) x %in% c("A","C","G","T")

# ============== Check inputs ==============
assert_file(BCFTOOLS, "bcftools")
assert_file(DBSNP_VCF, "dbSNP b37 vcf.gz")
assert_file(paste0(DBSNP_VCF, ".tbi"), "dbSNP b37 vcf.gz.tbi")

files <- list.files(EQTL_DIR, pattern = PATTERN, full.names = TRUE)
if (length(files) == 0) stop("No eQTL files matched: ", EQTL_DIR, " / ", PATTERN)
msg("[Input] eQTL files: %d", length(files))

# ============== 1) Collect all SNP IDs ==============
msg("==> Collecting all rsIDs ...")
all_snps <- unique(unlist(lapply(files, function(f) {
  x <- fread(f, select = "SNP")
  x$SNP
}), use.names = FALSE))

snps_txt <- file.path(OUT_DIR, paste0(CELLTYPE, "_all_snps.txt"))
writeLines(all_snps, snps_txt)
msg("[SNP] unique rsIDs: %d  -> %s", length(all_snps), snps_txt)

# ============== 2) Query dbSNP b37 ==============
msg("==> Query dbSNP(b37) for chr/pos/ref/alt ...")
tmp_vcf <- file.path(OUT_DIR, paste0(CELLTYPE, "_snps_dbsnp_b37.vcf.gz"))
tmp_tsv <- file.path(OUT_DIR, paste0(CELLTYPE, "_rsid_to_b37.tsv"))

cmd_view <- sprintf(
  "%s view -i 'ID=@%s' %s -Oz -o %s",
  shQuote(BCFTOOLS), shQuote(snps_txt), shQuote(DBSNP_VCF), shQuote(tmp_vcf)
)
status <- system(cmd_view)
if (status != 0) stop("bcftools view failed, status=", status)

status <- system(sprintf("%s index -t %s", shQuote(BCFTOOLS), shQuote(tmp_vcf)))
if (status != 0) stop("bcftools index failed, status=", status)

cmd_q <- sprintf(
  "%s query -f '%%ID\\t%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' %s > %s",
  shQuote(BCFTOOLS), shQuote(tmp_vcf), shQuote(tmp_tsv)
)
status <- system(cmd_q)
if (status != 0) stop("bcftools query failed, status=", status)

map <- fread(tmp_tsv, header = FALSE)
setnames(map, c("SNP","chr_b37","pos_b37","ref_b37","alt_b37"))
map[, `:=`(
  chr_b37 = gsub("^chr","", chr_b37),
  ref_b37 = toupper(ref_b37),
  alt_b37 = toupper(alt_b37)
)]

# 处理多等位 ALT
map[, alt_b37 := tstrsplit(alt_b37, ",", fixed = TRUE)[[1]]]

msg("[Map] rsIDs found in dbSNP b37: %d / %d", nrow(map), length(all_snps))

# ============== 3) Convert each chromosome file ==============
msg("==> Converting per-chromosome eQTL files -> b37 ...")

stats <- list()

for (f in files) {
  fn <- basename(f)
  msg("  - %s", fn)

  dt <- fread(f)
  # standardize
  dt[, `:=`(
    SNP = as.character(SNP),
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele)
  )]

  # drop non-ACGT alleles early (indels in this eQTL table are rare; if you want to keep them later, we can enhance)
  dt <- dt[is_acgt(effect_allele) & is_acgt(other_allele)]

  # merge b37 map
  m <- merge(dt, map, by="SNP", all.x=TRUE)

  n0 <- nrow(m)
  m <- m[!is.na(pos_b37) & !is.na(ref_b37) & !is.na(alt_b37)]
  n_map <- nrow(m)

  # allele match to b37 REF/ALT
  same <- (m$effect_allele == m$ref_b37 & m$other_allele == m$alt_b37)
  swap <- (m$effect_allele == m$alt_b37 & m$other_allele == m$ref_b37)

  # flip beta if swapped to align effect allele = b37 REF
  m$beta[swap] <- -as.numeric(m$beta[swap])

  m_keep <- m[same | swap]
  n_keep <- nrow(m_keep)
  n_swap <- sum(swap, na.rm=TRUE)

  # output columns: keep original + add b37 coord + overwrite chrom/position/effect_allele/other_allele to b37
  out <- m_keep
  out[, `:=`(
    chrom_hg38 = chrom,
    position_hg38 = position,
    chrom = paste0("chr", chr_b37),
    position = pos_b37,
    effect_allele = ref_b37,
    other_allele  = alt_b37
  )]

  out_file <- file.path(OUT_DIR, str_replace(fn, "\\.csv$", "_b37.csv"))
  fwrite(out, out_file)

  stats[[fn]] <- data.table(
    file=fn,
    rows_input=nrow(dt),
    rows_after_map=n_map,
    rows_keep=n_keep,
    beta_flipped=n_swap,
    keep_rate=ifelse(n0>0, n_keep/n0, NA_real_)
  )
}

stats_dt <- rbindlist(stats, fill=TRUE)
stats_file <- file.path(OUT_DIR, paste0(CELLTYPE, "_b37_convert_stats.tsv"))
fwrite(stats_dt, stats_file, sep="\t")

msg("[DONE] b37 eQTL written to: %s", OUT_DIR)
msg("[DONE] stats: %s", stats_file)
msg("Preview stats:")
print(head(stats_dt, 10))
