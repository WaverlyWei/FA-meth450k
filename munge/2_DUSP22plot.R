library(here)
library(ggplot2)
library(reshape)

dusp22 <- read.csv(here("results","Desktop/SuperFund/ordered_dusp22_pheno.csv"))

# row 35th: FA
# col 2: 73: target ids 
betas <- t(dusp22[,c(2:73)])
# name the column by probe IDs
# the last column is FA
colnames(betas) <- dusp22[,1]
                
# target ID vs. FA values
FA <- melt(betas[,35], id.var = "row.names")
# Var1: target ID
FA$Var1 <- rownames(FA)

beta <- melt(betas[,-35], id.var = "pos")
colnames(beta)[1] <- "Var1"

# merge by target ID
merged <- merge(beta, FA, by = "Var1", all = TRUE)
colnames(merged) <- c("TargetID", "ProbeID", "Value", "FA")

# 74: probe ID 91: pos
site_pos <- dusp22[-35,c(74,91)]
colnames(site_pos) <- c("ProbeID", "Position")
merged2 <- merge(merged, site_pos, by = "ProbeID", all = TRUE)
merged2 <- arrange(merged2, Position)
merged2$Pos_ID <- paste(merged2$ProbeID, merged2$Position)

pdf(file = here::here("Desktop/SuperFund/results", "DUSP22_plot.pdf"))
merged2 %>%
  group_by(FA, Pos_ID) %>%
  ggplot(., aes(x = Pos_ID, y = Value, fill = as.factor(FA))) +
  geom_boxplot(outlier.size = 1) +
  scale_x_discrete("ProbeID_Position") +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  scale_fill_manual(values=c("#009E73", "#D55E00"),
                    name="",
                    labels=c("FA Controls", "FA Exposed")) +
  labs(x = "CpG Probe ID and Genomic Location",
       y = "Beta Value (% DNA Methylation)",
       title = "Probes within DUSP22 Gene")  +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
