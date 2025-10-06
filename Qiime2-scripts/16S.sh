############### Qiime Scripts used for 16S rRNA metabarcoding dataset ###############
### All scripts were ran on the Georgia Advanced Computing Resource Center (GACRC) at UGA
### These scripts should be used on a QIIME2 artifact that contained the all demultiplexed files 
### All steps were ran individually in the format Job Submission Resource Parameters + Step #

##### Example Job Submission Resource Parameters
#!/bin/bash
#SBATCH --job-name="Importing_16S_Data"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=24:00:00
#SBATCH -e Importing_16S_Data.err-%N
#SBATCH -o Importing_16S_Data.out-%N

# Step 1 Importing Data
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/16s/16s_manifest.txt
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca.qza

# Script
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $INPUT \
  --output-path $OUTPUT \
  --input-format PairedEndFastqManifestPhred33V2



# Step 2 Use cutadapt to remove primers
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca.qza
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca-trimmed.qza

# Sript
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $INPUT \
  --p-adapter-r TTACCGCGGCKGCTGRCACACAATTACCATA \
  --p-adapter-f ATTAGAWACCCBNGTAGTCCGGCTGGCTGACT \
  --p-error-rate 0.2 \
  --o-trimmed-sequences $OUTPUT



# Step 3 Summarize Demultiplexed Sequences
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca.qza
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca.qzv

# Script
qiime demux summarize \
  --i-data $INPUT \
  --o-visualization $OUTPUT



# Step 4 Denoise with DADA2
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca-trimmed.qza
OUTPUT=/home/gas51591/boba/16s/analysis/data-sequences-ca-trimmed-dada2.qza

# Script
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $INPUT \
  --p-trim-left-f 0 \
  --p-trunc-len-f 250 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 200 \
  --p-n-threads 12 \
  --o-representative-sequences $OUTPUT \
  --o-table /home/gas51591/boba/16s/ca_analysis/16s-dada2-table-ca.qza \
  --o-denoising-stats /home/gas51591/boba/16s/ca_analysis/16s-dada2-stats-ca.qza



# Step 5 Exporting DADA2 Stats
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime metadata tabulate \
  --m-input-file /home/gas51591/boba/16s/ca_analysis/16s-dada2-stats-ca.qza \
  --o-visualization /home/gas51591/boba/16s/ca_analysis/16s-dada2-stats-ca.qzv

qiime feature-table summarize \
  --i-table /home/gas51591/boba/16s/ca_analysis/16s-dada2-table-ca.qza \
  --o-visualization /home/gas51591/boba/16s/ca_analysis/16s-dada2-table-ca.qzv \
  --m-sample-metadata-file /home/gas51591/boba/16s/metadata/16s_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data /home/gas51591/boba/16s/ca_analysis/data-sequences-ca-trimmed-dada2.qza \
  --o-visualization /home/gas51591/boba/16s/ca_analysis/data-sequences-ca-tabulated.qzv



# Step 6 Import Taxonomy Files
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/tax
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime tools import \
  --type FeatureData[Sequence] \
  --input-path /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.fasta \
  --output-path /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.qza

qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /home/gas51591/boba/16s/tax/SILVA138-nr99_taxonomy_16S.txt \
  --output-path /home/gas51591/boba/16s/tax/SILVA138-nr99_ref_taxonomy_16S.qza



# Step 7 BLAST Taxonomy 
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/tax
OUTPUT=/home/gas51591/boba/16s/analysis

# Script
qiime feature-classifier classify-consensus-blast \
  --i-query /home/gas51591/boba/16s/ca_analysis/data-sequences-ca-trimmed-dada2.qza \
  --i-reference-reads /home/gas51591/boba/16s/tax/SILVA138-nr99_sequences_16S.qza \
  --i-reference-taxonomy /home/gas51591/boba/16s/tax/SILVA138-nr99_ref_taxonomy_16S.qza \
  --p-maxaccepts 10 \
  --p-min-consensus .51 \
  --p-perc-identity .80 \
  --o-classification /home/gas51591/boba/16s/ca_analysis/ca-16s-rep-sequences-taxonomy.qza \
  --o-search-results /home/gas51591/boba/16s/ca_analysis/ca-16s-rep-sequences-tax-search-results.qza



# Step 8 Taxonomic Visualization 
# Activate QIIME2/2024.2-amplicon and add Path Variables 
module load QIIME2/2024.2-amplicon
INPUT=/home/gas51591/boba/16s/analysis/ca-16s-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/16s/analysis/ca-16s-rep-sequences-taxonomy.qzv

# Script
qiime metadata tabulate \
  --m-input-file $INPUT \
  --o-visualization $OUTPUT



# Step 9 Create Tree for Phylogentic Diversity Analysis 
# Activate QIIME2-2024.2 in Miniconda3 Environment 
module load Miniconda3
source activate /home/gas51591/conda-envs/envs/QIIME2-2024.2
INPUT=/home/gas51591/boba/16s/analysis/18S-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/16s/analysis/

# Script
qiime phylogeny  align-to-tree-mafft-fasttree \
  --i-sequences /home/gas51591/boba/16s/ca_analysis/data-sequences-ca-trimmed-dada2.qza \
  --o-alignment  /home/gas51591/boba/16s/ca_analysis/ca-aligned-16S-rep-seqs.qza \
  --o-masked-alignment /home/gas51591/boba/16s/ca_analysis/ca-16s-masked-aligned-rep-seqs.qza \
  --o-tree /home/gas51591/boba/16s/ca_analysis/ca-unrooted-16s-tree.qza \
  --o-rooted-tree /home/gas51591/boba/16s/ca_analysis/ca-rooted-16s-tree.qza