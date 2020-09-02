# Load the mammalian data.
dat = read.table("results/hg19_to_mm9.out", stringsAsFactors=FALSE, header=FALSE, sep="\t", col.names=c("chrom", "start", "end", "name", "peak", "label", "mappedChrom", "mappedStart", "mappedEnd", "mappedPeak"))
# Draw the pie chart.
pdf("results/hg19-mm9.pie.pdf")
pie(c(length(which(dat$label == "ortholog")), length(which(dat$label == "gain_hg19")), length(which(dat$label == "loss_mm9"))), labels=c("Ortholog","Human Gain","Mouse Loss"), main="Human to Mouse")
dev.off()

# Load the primate data.
dat = read.table("results/hg19_to_panTro4.chain.out", stringsAsFactors=FALSE, header=FALSE, sep="\t", col.names=c("chrom", "start", "end", "name", "peak", "label", "mappedChrom", "mappedStart", "mappedEnd", "mappedPeak"))
# Draw the pie chart.
pdf("results/hg19-panTro4.pie.pdf")
pie(c(length(which(dat$label == "ortholog")), length(which(dat$label == "gain_hg19")), length(which(dat$label == "loss_panTro4"))), labels=c("Ortholog","Human Gain","Chimp Loss"), main="Human to Chimpanzeee")
dev.off()

# Load the invertebrate data.
dat = read.table("results/dm3_to_droSim1.out", stringsAsFactors=FALSE, header=FALSE, sep="\t", col.names=c("chrom", "start", "end", "name", "peak", "label", "mappedChrom", "mappedStart", "mappedEnd", "mappedPeak"))
# Draw the pie chart.
pdf("results/dm3-droSim1.pie.pdf")
pie(c(length(which(dat$label == "ortholog")), length(which(dat$label == "gain_dm3")), length(which(dat$label == "loss_droSim1"))), labels=c("Ortholog","D. melanogaster Gain","D. simulans Loss"), main="D. melanogaster to D. simulans")
dev.off()
