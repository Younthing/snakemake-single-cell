library(survival)
library(survminer)
library(reshape2)
library(SingleCellExperiment)
library(zellkonverter)
library(ggExtra)
library(ggpmisc)
library(MuSiC)
# TODO 添加多变量精确
dir.create("27-MUSIC反卷积")
load("/home/fanxi/projects/work/sct-project/sct-project-1/bulk_rna/1-TCGA数据下载/crc_rnarse.rda")
rse <- crc_rnarse
counts_exp = assay(rse, "unstranded")

sce <- readH5AD("results/adata_raw_CCL_LCL/nmf_subtype/anndata_nmf_subtype_harmony_celltypist_Macrophages.h5ad")

Est.prop <- music_prop(
  bulk.mtx = counts_exp, sc.sce = sce, clusters = "cNMF_cluster",
  samples = "patients_organ", verbose = F
)

df <- melt(Est.prop$Est.prop.weighted)
colnames(df) <- c("Sample", "Celltype", "Prop")
df$Celltype <- factor(df$Celltype, levels = levels(sce$cNMF_cluster))
p4 <- ggplot(df, aes(x = Sample, y = Prop, fill = Celltype)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.position = "bottom") +
  labs(x = "Samples", y = "Precent")
p4
ggsave(p4, file = "./27-MUSIC反卷积/music_est_prop.pdf", width = 8, height = 6)

surv <- read.csv("/home/fanxi/projects/work/sct-project/sct-project-1/bulk_rna/5-TCGA预后信息/临床变量rse_clinical.csv")



cells <- unique(sce$cNMF_cluster)
for (cell in cells) {
  try({
    df3 <- merge(data.frame(
      sample = rownames(Est.prop$Est.prop.weighted),
      EC = Est.prop$Est.prop.weighted[, cell]
    ), surv, by.x = "sample", by.y = "id")

    # 前后15%分组比较
    df3_15 <- df3[order(df3$EC), ]
    top_15 <- df3_15[1:floor(0.15 * nrow(df3_15)), ]
    bottom_15 <- df3_15[(nrow(df3_15) - floor(0.15 * nrow(df3_15)) + 1):nrow(df3_15), ]
    top_15$ec.group <- "Top 15%"
    bottom_15$ec.group <- "Bottom 15%"
    df3_subset_15 <- rbind(top_15, bottom_15)

    fit_15 <- survfit(Surv(os_time, os) ~ ec.group, data = df3_subset_15)
    ggsurvplot(fit_15,
      pval = TRUE,
      conf.int = FALSE,
      xlab = "Time in days",
      ggtheme = theme_light(),
      palette = c("#E7B800", "#00AFBB")
    )
    ggsave(paste0("./27-MUSIC反卷积/", cell, "_KM_OS_TopBottom15.pdf"), width = 8, height = 6)

    # 前后20%分组比较
    df3_20 <- df3[order(df3$EC), ]
    top_20 <- df3_20[1:floor(0.2 * nrow(df3_20)), ]
    bottom_20 <- df3_20[(nrow(df3_20) - floor(0.2 * nrow(df3_20)) + 1):nrow(df3_20), ]
    top_20$ec.group <- "Top 20%"
    bottom_20$ec.group <- "Bottom 20%"
    df3_subset_20 <- rbind(top_20, bottom_20)

    fit_20 <- survfit(Surv(os_time, os) ~ ec.group, data = df3_subset_20)
    ggsurvplot(fit_20,
      pval = TRUE,
      conf.int = FALSE,
      xlab = "Time in days",
      ggtheme = theme_light(),
      palette = c("#E7B800", "#00AFBB")
    )
    ggsave(paste0("./27-MUSIC反卷积/", cell, "_KM_OS_TopBottom20.pdf"), width = 8, height = 6)

    # 前后30%分组比较
    df3_30 <- df3[order(df3$EC), ]
    top_30 <- df3_30[1:floor(0.3 * nrow(df3_30)), ]
    bottom_30 <- df3_30[(nrow(df3_30) - floor(0.3 * nrow(df3_30)) + 1):nrow(df3_30), ]
    top_30$ec.group <- "Top 30%"
    bottom_30$ec.group <- "Bottom 30%"
    df3_subset_30 <- rbind(top_30, bottom_30)

    fit_30 <- survfit(Surv(os_time, os) ~ ec.group, data = df3_subset_30)
    ggsurvplot(fit_30,
      pval = TRUE,
      conf.int = FALSE,
      xlab = "Time in days",
      ggtheme = theme_light(),
      palette = c("#E7B800", "#00AFBB")
    )
    ggsave(paste0("./27-MUSIC反卷积/", cell, "_KM_OS_TopBottom30.pdf"), width = 8, height = 6)
  })
}
