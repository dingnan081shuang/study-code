# Final Version of Metabolomics R Analysis
rm(list=ls())
library(tidyverse)
library(ropls)
library(openxlsx)

# Set parameters
options(scipen = 9999)

# Load and process data
data <- read.csv("RPLVSCON_final_metabolomics.csv", header = TRUE, row.names = 1)
data <- as.data.frame(t(data[,-1]))
exprSet <- data
ex <- as.matrix(exprSet)

# Check if log transformation is needed
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  # Log2 transform
  exprSet <- log2(ex)
  print("log2 transformation finished")
} else {
  print("log2 transformation not needed")  
}

exprSet <- as.data.frame(exprSet)
exprSet <- as.data.frame(t(exprSet))
# Remove outlier samples
exprSet <- exprSet[-c(187,188,189,190,191), ]
Group <- c(rep("Control", 122), rep("PL", 66))
exprSet <- data.frame(Group, exprSet)
data <- exprSet
save(data, file = "outputfile/data.Rdata")

# Separate features and labels
x <- data[, -1]
y <- data$Group

# OPLS-DA modeling
plsda <- opls(x, y,
              predI = 1, orthoI = 1,
              log10L = FALSE,
              crossvalI = nrow(x),
              scaleC = "pareto",
              fig.pdfC = "outputfile/plsda.pdf",
              permI = 200)

# Extract VIP values
vip <- plsda@vipVn %>% as.data.frame() %>%
  rename("VIP" = ".") %>%
  rownames_to_column(var = 'variables')

# Function to calculate FC, logFC, and p-values
calc_fc <- function(df, group1, group2) {
  fc <- rowMeans(df[, group2]) / rowMeans(df[, group1])
  log_fc <- log2(fc)
  t_test_p <- apply(df, 1, function(row) {
    t.test(row[group1], row[group2])$p.value
  })
  return(data.frame(FC = fc, logFC = log_fc, p.value = t_test_p))
}

group1_indices <- which(data$Group == "Control")
group2_indices <- which(data$Group == "PL")

fc_results <- calc_fc(t(x), group1_indices, group2_indices)
fc_results <- data.frame(rownames(fc_results), fc_results)
colnames(fc_results)[1] <- 'variables'

# Merge results
dat_vip <- vip %>%
  left_join(fc_results, by = "variables")
dat_vip %>% write_csv("outputfile/opls-daResult1.csv")
dat_vip %>% write_csv("outputfile/opls-daResult.csv")
save(dat_vip, file = "outputfile/allDiff.Rdata")

# Calculate score data for OPLS-DA plot
plsdaScore <- data.frame(
  t1 = plsda@scoreMN,
  to1 = plsda@orthoScoreMN
) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  rename(
    "t1" = "p1",
    "to1" = "o1"
  ) %>%
  rownames_to_column(var = 'sample') %>%
  mutate(group = data$Group)

t1Weight <- sprintf("%.1f%%", plsda@modelDF[1, 1] * 100)
to1Weight <- sprintf("%.1f%%", plsda@modelDF[2, 1] * 100)

R2X <- plsda@modelDF[1, 1] + plsda@modelDF[2, 1]
R2Y <- plsda@modelDF[1, 3] + plsda@modelDF[2, 3]
Q2Y <- plsda@modelDF[1, 6] + plsda@modelDF[2, 6]

subTitle <- paste0("R2X=", R2X, "  R2Y=", R2Y, "  Q2Y=", Q2Y)

# Plot OPLS-DA
oplsdaFig <- ggplot(plsdaScore, aes(x = t1, y = to1, color = group)) +
  geom_point(aes(shape = group), size = 0.8) +
  stat_ellipse(aes(fill = group), alpha = 0.2, geom = "polygon") +
  theme_classic() +
  scale_fill_manual(values = c("#17becf", "#d62728"),
                    labels = c("Control", "PL")) +
  scale_color_manual(values = c("#17becf", "#d62728"),
                     labels = c("Control", "PL")) +
  labs(title = "OPLS-DA",
       subtitle = subTitle,
       x = paste0("t1 (", t1Weight, ")"),
       y = paste0("to1 (", to1Weight, ")")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    plot.subtitle = element_text(hjust = 0.5, size = rel(0.6))
  )

# Save plot
ggsave("outputfile/oplsda.png", oplsdaFig, width = 4, height = 4)
ggsave("outputfile/oplsda.pdf", oplsdaFig, width = 4, height = 4, onefile = FALSE)

# Plot OPLS-DA with sample labels
oplsdaFigWithLabels <- ggplot(plsdaScore, aes(x = t1, y = to1, color = group)) +
  geom_point(aes(shape = group), size = 2) +
  geom_text(aes(label = sample), hjust = 0.5, vjust = -0.5, size = 3, show.legend = FALSE) +
  stat_ellipse(aes(fill = group), alpha = 0.2, geom = "polygon") +
  theme_classic() +
  scale_fill_manual(values = c("#17becf", "#d62728"),
                    labels = c("Control", "PL")) +
  scale_color_manual(values = c("#17becf", "#d62728"),
                     labels = c("Control", "PL")) +
  labs(title = "OPLS-DA with Sample Labels",
       x = paste0("t1 (", t1Weight, ")"),
       y = paste0("to1 (", to1Weight, ")")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)),
    plot.subtitle = element_text(hjust = 0.5, size = rel(0.6))
  )

# Save labeled plot
ggsave("outputfile/oplsda_with_labels.png", oplsdaFigWithLabels, width = 8, height = 6)

# Calculate Z-score to identify outliers
plsdaScore$z_t1 <- scale(plsdaScore$t1)
plsdaScore$z_to1 <- scale(plsdaScore$to1)

# Identify outliers with Z-score threshold > 3
threshold <- 3
outliers <- plsdaScore %>%
  filter(abs(z_t1) > threshold | abs(z_to1) > threshold)

# Print outliers
print(outliers)

# Save final list of differential metabolites with VIP > 1.0, FC > 1.2, and p-value < 0.05
diffMetabolites <- dat_vip %>% 
  filter(p.value < 0.05) %>%
  filter(VIP > 1) %>%
  filter(abs(logFC) > 0.585) %>%
  rownames_to_column(var = 'variables')

save(diffMetabolites, file = "outputfile/diffMetabolites.Rdata")
write_csv(diffMetabolites, "outputfile/diffMetabolites.csv")
