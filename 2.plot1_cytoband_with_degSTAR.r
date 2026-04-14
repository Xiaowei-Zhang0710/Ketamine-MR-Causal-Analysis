suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(grid)
})

# ===================== CONFIG =====================
INFILE   <- "/Users/nullyann/Desktop/anding_smr/data/gwas_category_outcome/gene_ct_full.tsv"
BIM_FILE <- "/Users/nullyann/Desktop/anding_smr/data/ldref/1000G_EUR/EUR.bim"
CYTOBAND_FILE <- "/Users/nullyann/Desktop/anding_smr/data/hg19_cytoBand.txt"

OUT_PDF <- "gene_on_chr_vertical_22inrow_SNPtick_DEGstars_FINAL.pdf"
OUT_PNG <- "gene_on_chr_vertical_22inrow_SNPtick_DEGstars_FINAL.png"

# 固定 celltype 顺序（用来排列星星/点）
preferred_order <- c("Ast", "Microglia", "oligo", "Inhibitory")

# ---------------- Layout ----------------
# dots（原来那列圆点）
Y_DOT_DY_MB  <- 0.55
X_DOT_OFFSET <- 0.30

# DEG 星星列（新增的一列）
X_DEG_OFFSET  <- 0.60    # 星星列离染色体中心更远（越大越靠左）
Y_STAR_DY_MB  <- 3  # ✅ 星星纵向间距（建议 0.75~1.0，越大越不挤）
STAR_SIZE     <- 0.7     # ✅ 星星大小（2.0~3.2）
STAR_STROKE   <- 0.8
STAR_SHAPE    <- 8       # 8 = asterisk/star

# gene label 在右侧
X_LABEL_OFFSET <- 0.48

# SNP tick（横线）
TICK_HALF_WIDTH <- 0.11

# 染色体宽度（cytoband画条的半宽）
CHR_HALF_WIDTH <- 0.12

# 输出质量
DPI <- 600

# ---------------- SNP tick style ----------------
SNP_TICK_COLOR <- "#007BFF"  # 鲜艳蓝
SNP_TICK_LWD   <- 1.45
SNP_TICK_ALPHA <- 0.98

# ---------------- Cytoband lighten ----------------
LIGHTEN_NO_SNP <- 0.60

# ---------------- spacing fix ----------------
# 如果 chr_i 上有基因/lead SNP，则 i 与 i+1 间距增加 GAP_DELTA
BASE_STEP <- 1.00
GAP_DELTA <- 0.70  # 你之前已经调大过；还不够可到 0.9

# ---------------- celltype colors ----------------
ct_colors <- c(
  "Ast"        = "#0072B2",
  "Microglia"  = "#D55E00",
  "oligo"      = "#009E73",
  "Inhibitory" = "#CC79A7"
)

# ===================== DEG mapping (from your handwritten table) =====================
# FRAS1: int
# ADGRB3: oligo, Mic
# ERC2: oligo, Mic, Ast
# NTM: oligo, Mic
# GRIN2A: oligo, Ast
# DAB1: oligo
# PACRG: oligo
# LRP1B: Ast
deg_map <- data.table(
  gene = c("FRAS1","ADGRB3","ERC2","NTM","GRIN2A","DAB1","PACRG","LRP1B"),
  deg_celltypes = c("int",
                    "oligo Mic",
                    "oligo Mic Ast",
                    "oligo Mic",
                    "oligo Ast",
                    "oligo",
                    "oligo",
                    "Ast")
)
deg_long <- deg_map[, .(celltype = unlist(strsplit(deg_celltypes, "\\s+"))), by = gene]
deg_long[, celltype := fifelse(celltype %in% c("Mic","MIC","micro","microglia"), "Microglia", celltype)]
deg_long[, celltype := fifelse(celltype %in% c("int","Int","INT","Inh","Inhib"), "Inhibitory", celltype)]
deg_long[, celltype := fifelse(tolower(celltype) %in% c("oligo","olig"), "oligo", celltype)]
deg_long[, celltype := fifelse(celltype %in% c("Ast","AST","astro"), "Ast", celltype)]
deg_long <- unique(deg_long[celltype %in% preferred_order])

# ===================== Helpers =====================
get_chr_sizes_from_bim <- function(bim_file) {
  bim <- fread(bim_file, header = FALSE)
  chr_sizes <- bim[, .(chr = as.integer(V1),
                       chr_len = max(as.integer(V4), na.rm = TRUE)), by = V1]
  chr_sizes <- chr_sizes[chr %in% 1:22]
  setorder(chr_sizes, chr)
  chr_sizes[, chr_len_mb := chr_len / 1e6]
  chr_sizes
}

stain_to_grey <- function(stain) {
  stain <- tolower(stain)
  out <- rep("grey90", length(stain))
  out[grepl("^gneg", stain)]    <- "grey96"
  out[grepl("^gpos100", stain)] <- "grey15"
  out[grepl("^gpos75", stain)]  <- "grey30"
  out[grepl("^gpos66", stain)]  <- "grey35"
  out[grepl("^gpos50", stain)]  <- "grey50"
  out[grepl("^gpos33", stain)]  <- "grey65"
  out[grepl("^gpos25", stain)]  <- "grey75"
  out[grepl("^gvar", stain)]    <- "grey80"
  out[grepl("^stalk", stain)]   <- "grey88"
  out[grepl("^acen", stain)]    <- "#D55E00"
  out
}

lighten_hex <- function(hex, frac = 0.6) {
  rgb <- grDevices::col2rgb(hex) / 255
  rgb2 <- rgb + (1 - rgb) * frac
  grDevices::rgb(rgb2[1,], rgb2[2,], rgb2[3,])
}

# ===================== Load data =====================
dt0 <- fread(INFILE)
stopifnot(all(c("gene","celltype","chr","pos") %in% names(dt0)))

dt0[, chr := as.integer(gsub("^chr","", as.character(chr)))]
dt0[, pos := as.integer(pos)]
dt0 <- dt0[chr %in% 1:22 & !is.na(pos) & gene != ""]
dt0[, pos_mb := pos / 1e6]

# normalize celltype strings in dt0
dt0[, celltype := as.character(celltype)]
dt0[celltype %in% c("Oligodendrocytes","Oligodendrocyte","Oligo","Oligo ","Oligodendrocytes "), celltype := "oligo"]
dt0[celltype %in% c("Astrocytes","Astrocyte","Astro"), celltype := "Ast"]
dt0[celltype %in% c("Microglia "), celltype := "Microglia"]
dt0[celltype %in% c("int","Int","Inhib"), celltype := "Inhibitory"]

# ===================== gene anchor (1 row per gene for positioning) =====================
# 若有 pval_outcome 则取最显著那条作为锚点，否则取最小pos
gene_anchor <- copy(dt0)
if ("pval_outcome" %in% names(gene_anchor)) {
  gene_anchor[, ptmp := suppressWarnings(as.numeric(pval_outcome))]
  setorder(gene_anchor, gene, ptmp, pos)
} else {
  setorder(gene_anchor, gene, pos)
}
gene_anchor <- gene_anchor[!duplicated(gene), .(gene, chr, pos, pos_mb)]

# ===================== chr sizes + has_snp =====================
chr_sizes <- if (file.exists(BIM_FILE)) get_chr_sizes_from_bim(BIM_FILE) else {
  tmp <- gene_anchor[, .(chr_len_mb = max(pos_mb, na.rm = TRUE)), by = chr][order(chr)]
  tmp[, chr_len := as.integer(chr_len_mb * 1e6)]
  tmp
}

# “有SNP/有基因锚点”的chr
chr_has <- gene_anchor[, .(has_snp = TRUE), by = chr]
chr_sizes <- merge(chr_sizes, chr_has, by="chr", all.x=TRUE)
chr_sizes[is.na(has_snp), has_snp := FALSE]
chr_sizes[, outline_col := ifelse(has_snp, "grey25", "grey80")]

# ===================== cumulative x spacing (correct) =====================
setorder(chr_sizes, chr)
x_center <- numeric(nrow(chr_sizes))
x_center[1] <- 1
if (nrow(chr_sizes) > 1) {
  for (i in 2:nrow(chr_sizes)) {
    extra <- ifelse(chr_sizes$has_snp[i-1], GAP_DELTA, 0)
    x_center[i] <- x_center[i-1] + BASE_STEP + extra
  }
}
chr_sizes[, x := x_center]

gene_anchor <- merge(gene_anchor, chr_sizes[,.(chr, x)], by="chr", all.x=TRUE)
gene_anchor[, x_label := x + X_LABEL_OFFSET]

# ===================== SNP tick data =====================
ticks <- gene_anchor[, .(
  x,
  y = pos_mb,
  x1 = x - TICK_HALF_WIDTH,
  x2 = x + TICK_HALF_WIDTH
)]

# ✅ Marker legend data (SNP dot)
marker_snp <- ticks[, .(x = x, y = y, mark = "SNP")]

# ===================== original dots column (optional, from gene_ct_full.tsv) =====================
dt_dots <- unique(dt0[!is.na(celltype) & celltype != "", .(gene, celltype)])
dt_dots <- dt_dots[celltype %in% preferred_order]
dt_dots <- merge(dt_dots, gene_anchor[,.(gene, x, pos_mb)], by="gene", all.x=TRUE)

dt_dots[, ct_rank := match(celltype, preferred_order)]
ct_n <- length(preferred_order)
dt_dots[, y_dot := pos_mb + ((ct_n - 1)/2 - (ct_rank - 1)) * Y_DOT_DY_MB]
dt_dots[, x_dot := x - X_DOT_OFFSET]

# ===================== DEG stars column =====================
deg_plot <- merge(deg_long, gene_anchor[,.(gene, x, pos_mb)], by="gene", all.x=TRUE)
deg_plot <- deg_plot[!is.na(x) & !is.na(pos_mb)]
deg_plot[, ct_rank := match(celltype, preferred_order)]

# ✅ 用更大的 Y_STAR_DY_MB，保证多celltype星星不挤
deg_plot[, y_star := pos_mb + ((ct_n - 1)/2 - (ct_rank - 1)) * Y_STAR_DY_MB]
deg_plot[, x_star := x - X_DEG_OFFSET]

# 用于 shape legend 的常量映射
deg_plot[, mark := "DEG"]

# ===================== Cytoband =====================
stopifnot(file.exists(CYTOBAND_FILE))
cb <- fread(CYTOBAND_FILE, header = FALSE)
setnames(cb, c("chrom","start","end","band","stain"))
cb[, chrom := gsub("^chr","", chrom)]
cb[, chr := suppressWarnings(as.integer(chrom))]
cb <- cb[chr %in% 1:22]
cb[, ymin := start / 1e6]
cb[, ymax := end / 1e6]
cb[, fill := stain_to_grey(stain)]

cb <- merge(cb, chr_sizes[,.(chr, has_snp, x)], by="chr", all.x=TRUE)
cb[has_snp == FALSE, fill := lighten_hex(fill, frac = LIGHTEN_NO_SNP)]

# ===================== Plot =====================
base_family <- "Helvetica"

p <- ggplot() +
  # Cytoband blocks
  geom_rect(
    data = cb,
    aes(xmin = x - CHR_HALF_WIDTH, xmax = x + CHR_HALF_WIDTH, ymin = ymin, ymax = ymax),
    fill = cb$fill, color = NA
  ) +
  # Chromosome outline
  geom_rect(
    data = chr_sizes,
    aes(xmin = x - CHR_HALF_WIDTH, xmax = x + CHR_HALF_WIDTH, ymin = 0, ymax = chr_len_mb),
    fill = NA, color = chr_sizes$outline_col, linewidth = 0.4
  ) +
  # SNP ticks (white halo + blue)
  geom_segment(
    data = ticks,
    aes(x = x1, xend = x2, y = y, yend = y),
    linewidth = SNP_TICK_LWD + 1.2, color = "white", alpha = 1, lineend = "butt"
  ) +
  geom_segment(
    data = ticks,
    aes(x = x1, xend = x2, y = y, yend = y),
    linewidth = SNP_TICK_LWD, color = SNP_TICK_COLOR, alpha = SNP_TICK_ALPHA, lineend = "butt"
  ) +
  # ✅ SNP dot (for marker legend + slightly clearer position)
  geom_point(
    data = marker_snp,
    aes(x = x, y = y, shape = mark),
    size = 1.7,
    color = "grey10",
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  # original dots column (celltype colored)
  geom_point(
    data = dt_dots,
    aes(x = x_dot, y = y_dot, color = celltype),
    size = 2.0, alpha = 0.85
  ) +
  # ✅ DEG stars column (celltype colored + marker legend)
  geom_point(
    data = deg_plot,
    aes(x = x_star, y = y_star, color = celltype, shape = mark),
    size = STAR_SIZE,
    stroke = STAR_STROKE,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  # Gene labels (stable white background)
  ggrepel::geom_label_repel(
    data = gene_anchor,
    aes(x = x_label, y = pos_mb, label = gene),
    direction = "y",
    hjust = 0,
    size = 3.4,
    fontface = "bold",
    family = base_family,
    color = "grey10",
    fill = alpha("white", 0.85),
    label.size = 0,
    label.r = unit(0.12, "lines"),
    min.segment.length = 0,
    segment.alpha = 0.25,
    segment.color = "grey40",
    box.padding = 0.20,
    point.padding = 0.12,
    max.overlaps = Inf
  ) +
  # Axes
  scale_x_continuous(
    breaks = chr_sizes$x,
    labels = paste0("chr", chr_sizes$chr),
    expand = expansion(mult = c(0.03, 0.20))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, " Mb"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  # Legends: color for celltype; shape for marker types
  scale_color_manual(values = ct_colors, name = "Cell type") +
  scale_shape_manual(
    name = "Marker",
    values = c("SNP" = 16, "DEG" = STAR_SHAPE)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 16, size = 3)),
    shape = guide_legend(order = 2, override.aes = list(color = "grey10"))
  ) +
  labs(x = NULL, y = "Genomic position (hg19/b37)") +
  theme_classic(base_size = 12, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line.x = element_blank(),
    plot.margin = margin(12, 60, 12, 90),  # 左侧留空间：星星列 + 点列
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(OUT_PDF, p, width = 26, height = 10, device = cairo_pdf, limitsize = FALSE)
ggsave(OUT_PNG, p, width = 26, height = 10, dpi = DPI, limitsize = FALSE)

print(p)