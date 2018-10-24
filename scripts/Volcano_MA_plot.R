library(scales)
library(ggplot2)
library(ggrepel)

# install.packages("scales")
# install.packages("ggplot2")
# install.packages("ggrepel")


#Load the dspi file (output from SUPPA diffSplice)
dpsi <- read.table(file="~/TRA2_diffSplice.dpsi",sep="\t")
colnames(dpsi) <- c("dPSI","p_value")
#Load the psi file (output from SUPPA generateEvents)
psi_events <- read.table(file="~/TRA2_events.psi",sep="\t")
colnames(psi_events) <- c("CTRL1","CTRL2","CTRL3","KD1","KD2","KD3")


#Load the tpm of the events (output from SUPPA diffSplice, activating flag --save_tpm_events)
event_TPM <- read.table(file="~/TRA2_diffSplice_avglogtpm.tab",sep="\t",header = TRUE)
colnames(event_TPM) <- c("event","mean_TPM")


#Merge dpsi and psi
merge1 <- merge(dpsi,psi_events,by="row.names")
merge2 <- merge(merge1,event_TPM,by.x="Row.names",by.y="event")
rownames(merge2) <- merge2$Row.names
final_table <- merge2
final_table <- final_table[,-1]
final_table <- final_table[!is.nan(final_table$dPSI),]
final_table$cpval <- p.adjust(final_table$p_value, method = "bonferroni")
final_table$log10pval <- -log10(final_table$p_value)
final_table$sig <- "not sig"
final_table[final_table$p_value < 0.05,]$sig <- "sig"
# final_table[, c(2:13)] <- sapply(final_table[, c(2:13)], function(x) as.numeric(gsub(",", ".", x)))
final_table$logRNAc <-final_table$mean_TPM

final_table <- cbind(final_table,rownames(final_table))
colnames(final_table)[14] <- "Name"
label <- "ENSG00000149554"

#Volcano plot

jpeg(file = "~/volcano_TRA2.jpeg")

p <- ggplot(final_table, aes(x=dPSI, y=log10pval, color=sig))
p + geom_point() + geom_label_repel(aes(label=ifelse(Name=="ENSG00000149554;SE:chr11:125496728-125497502:125497725-125499127:+","CHEK1",''))) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="solid", size=1) +
  geom_hline(yintercept=1.3, size=1) +
  xlab(expression(~Delta~PSI)) + ylab("-log10(p-value)") + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="none",
    legend.text=element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)) +
  scale_color_manual(values=c("sig" = "darkblue", "not sig" = "gray40", "nan" = "gray80")) + 
  scale_x_continuous(breaks=pretty_breaks(n=5), limits = c(-1, 1)) +
  scale_y_continuous(breaks=pretty_breaks(n=5), limits = c(0, 4.5)) +
  guides(fill = guide_legend(reverse = FALSE))

dev.off()





#MA plot (filtering out genes with average logRNAc value below 0)
jpeg(file = "~/MA_TRA2.jpeg")

p <- ggplot(final_table[final_table$logRNAc > 0,], aes(x=logRNAc, y=dPSI, color=sig))
p + geom_point() + geom_label_repel(aes(label=ifelse(Name=="ENSG00000149554;SE:chr11:125496728-125497502:125497725-125499127:+","CHEK1",''))) +
  scale_color_manual(values=c("sig" = "darkblue", "not sig" = "gray40", "nan" = "gray80")) + 
  labs(title="...", x="Average transcript abundance", y=expression(~Delta~PSI)) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title=element_blank())

dev.off()

