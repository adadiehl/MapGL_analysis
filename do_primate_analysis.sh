# This shell script will execute commands to reproduce the primate analysis presented in:
# MapGL: Inferring evolutionary gain and loss of short genomic sequence features by
# phylogenetic maximum parsimony. Diehl, AG, and Boyle, AP. 2020. BMC Bioinformatics.

# See README for software requirements and versions

# Retrieve CTCF ChIP-seq data from ENCODE and merge peaks. (URLs correct as of 2 Septemeber 2020)
cd data/ChIP-seq
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF085HTY.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF002CEL.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF738TKN.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF002DDJ.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF002DBD.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF710VEH.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF096AKZ.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF002DAJ.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF963PJY.bed.gz
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF123VEZ.bed.gz
bedtools merge -i <(zcat ENCFF* | bedtools sort) -c 2,3,5,7,8,9,10 -o collapse -delim "," | awk '{ N=split($4, A, ","); split($5, B, ","); split($6, C, ","); split($7, D, ","); split($8, E, ","); split($9, F, ","); split($10, G, ","); SCORE = 0; SVAL = 0; PVAL = 0; QVAL = 0; for (i=1; i<=N; i++) { SCORE += C[i]; SVAL += D[i]; PVAL += E[i]; QVAL += F[i] }; TOT = 0; for (i=1; i<=N; i++) { TOT+=(A[i]+G[i]) }; printf "%s\t%d\t%d\t.\t%d\t.\t%.6f\t%.6f\t%.6f\t%d\n", $1, $2, $3, SCORE/N, SVAL/N, PVAL/N, QVAL/N, (TOT/N)-$2 }' > CTCF.hg19.merged.narrowPeak

# Retrieve 100-way Newick tree from UCSC (URL correct as of 2 Septemeber 2020)
cd ../phylo
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/hg19.100way.nh

# Prune out chosen primate species and draw a rendering with draw_tree...
tree_doctor -P hg19,panTro4,rheMac3,macFas5,calJac3,saiBol1 hg19.100way.nh > hg19.primates.small.nh
draw_tree hg19.primates.small.nh > ../../results/hg19.primates.small.ps

# Download alignment chains and chrom.sizes. (URLs correct as of 2 Septemeber 2020)
cd ../chains
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToPanTro4.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToRheMac3.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMacFas5.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToCalJac3.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToSaiBol1.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/panTro4/bigZips/panTro4.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/rheMac3/bigZips/rheMac3.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/calJac3/bigZips/calJac3.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.chrom.sizes
wget http://hgdownload.soe.ucsc.edu/goldenPath/saiBol1/bigZips/saiBol1.chrom.sizes
for spp2 in panTro4 rheMac3 macFas5 calJac3 saiBol1 ; do ln -s $(printf "hg19To%s.over.chain.gz\n" $(echo $spp2 | sed 's/^./\U&\E/')) hg19.$spp2.over.chain.gz ; done
cd ../../

# Run MapGL on the human ChIP-seq data...
# This was run within an Anaconda environment using Python 3.7.6 and MapGL version 1.2.0
mapGL.py data/ChIP-seq/CTCF.hg19.merged.narrowPeak data/phylo/hg19.primates.small.nh hg19 panTro4 data/chains/hg19.panTro4.over.chain.gz data/chains/hg19.rheMac3.over.chain.gz data/chains/hg19.macFas5.over.chain.gz data/chains/hg19.calJac3.over.chain.gz data/chains/hg19.saiBol1.over.chain.gz > results/hg19_to_panTro4.out
