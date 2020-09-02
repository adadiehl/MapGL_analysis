# This shell script will execute commands to reproduce the invertebrate analysis presented in:
# MapGL: Inferring evolutionary gain and loss of short genomic sequence features by
# phylogenetic maximum parsimony. Diehl, AG, and Boyle, AP. 2020. BMC Bioinformatics.

# See README for software requirements and versions

# Download D. melanogaster CTCF ChIP-seq data (URL correct as of 2 Septemeber 2020)
cd data/ChIP-seq
wget https://www.encodeproject.org/files/ENCFF123VEZ/@@download/ENCFF123VEZ.bed.gz
gunzip ENCFF123VEZ.bed.gz

# Retrieve invertebrate 26-way Newick tree from UCSC (URL correct as of 2 Septemeber 2020)
cd ../phylo
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/dm6.27way.nh

# Prune out chosen invertebrate species and draw a rendering with draw_tree...
tree_doctor -P dm6,droSim1,droAna3,droPse3,anoGam1 dm6.27way.nh | sed 's/dm6/dm3/' | sed 's/droPse3/dp3/' > dm3.nh
draw_tree dm3.nh > ../../results/dm3.ps

# Download alignment chains and chrom.sizes. (URLs correct as of 2 Septemeber 2020)
cd ../chains
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDroSim1.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDroAna3.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDp3.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToAnoGam1.over.chain.gz
for spp2 in droSim1 droAna3 dp3 anoGam1; do ln -s $(printf "dm3To%s.over.chain.gz\n" $(echo $spp2 | sed 's/^./\U&\E/')) dm3.$spp2.over.chain.gz ; done
cd ../../

# Run MapGL on the drosophila ChIP-seq data...
# This was run within an Anaconda environment using Python 3.7.6 and MapGL version 1.2.0
mapGL.py data/ChIP-seq/ENCFF123VEZ.bed data/phylo/dm3.nh dm3 droSim1 data/chains/dm3.droSim1.over.chain.gz data/chains/dm3.droAna3.over.chain.gz data/chains/dm3.dp3.over.chain.gz data/chains/dm3.anoGam1.over.chain.gz > results/dm3_to_droSim1.out 
