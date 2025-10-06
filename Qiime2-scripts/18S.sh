############### Qiime Scripts used for 18S rRNA metabarcoding dataset ###############
### All scripts were ran on the Georgia Advanced Computing Resource Center (GACRC) at UGA
### These scripts should be used on a QIIME2 artifact that contained the all demultiplexed files 
### All steps were ran individually in the format Job Submission Resource Parameters + Step #

##### Example Job Submission Resource Parameters
#!/bin/bash
#SBATCH --job-name="Importing_18S_Data"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=24:00:00
#SBATCH -e Importing_18S_Data.err-%N
#SBATCH -o Importing_18S_Data.out-%N

# Step 1 Importing Data
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/data-clean
OUTPUT=/home/gas51591/boba/updated-database/analysis/data-sequences.qza

# Script
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $INPUT \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path $OUTPUT



# Step 2 Summarize Demultiplexed Sequences
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/updated-database/analysis/data-sequences.qza
OUTPUT=/home/gas51591/boba/updated-database/analysis/data-sequences.qzv

# Script
qiime demux summarize \
--i-data $INPUT \
--o-visualization $OUTPUT



# Step 3 Denoise with DADA2
# Activate QIIME2/2023.2_picrust2 and add Path Variables 
module load QIIME2/2023.2_picrust2
INPUT=/home/gas51591/boba/updated-database/analysis/data-sequences.qza
OUTPUT=/home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qza

# Script
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $INPUT \
--p-trim-left-f 0 \
--p-trunc-len-f 209 \
--p-trim-left-r 0 \
--p-trunc-len-r 230 \
--o-representative-sequences $OUTPUT \
--o-table /home/gas51591/boba/updated-database/analysis/18S-dada2-table.qza \
--o-denoising-stats /home/gas51591/boba/updated-database/analysis/18S-dada2-stats.qza



# Step 4 Exporting DADA2 Stats
# Activate QIIME2/2023.2_picrust2 and add Path Variables 
module load QIIME2/2023.2_picrust2
INPUT=/home/gas51591/boba/updated-database/analysis
OUTPUT=/home/gas51591/boba/updated-database/analysis

# Script
qiime metadata tabulate \
--m-input-file /home/gas51591/boba/updated-database/analysis/18S-dada2-stats.qza \
--o-visualization /home/gas51591/boba/updated-database/analysis/18S-dada2-stats.qzv

qiime feature-table summarize \
--i-table /home/gas51591/boba/updated-database/analysis/18S-dada2-table.qza \
--o-visualization /home/gas51591/boba/updated-database/analysis/18S-dada2-table.qzv \
--m-sample-metadata-file /home/gas51591/boba/analysis/seagrasses_metadata_habitat_location.tsv

qiime feature-table tabulate-seqs \
--i-data /home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qza \
--o-visualization /home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qzv



# Step 5 BLAST Taxonomy 
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qza
OUTPUT=/home/gas51591/boba/updated-database/analysis/18S-rep-sequences-taxonomy.qza

# Script
qiime feature-classifier classify-consensus-blast \
--i-query /home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qza \
--i-reference-reads /home/gas51591/boba/updated-database/tax/SILVA138-nr99_fixed_strings_custom_seq_sequences-May.qza \
--i-reference-taxonomy /home/gas51591/boba/updated-database/tax/SILVA138-nr99_fixed_strings_custom_seq_taxonomy-May.qza \
--p-maxaccepts 1 \
--p-perc-identity .90 \
--o-classification /home/gas51591/boba/updated-database/analysis/18S-rep-sequences-taxonomy.qza \
--o-search-results /home/gas51591/boba/updated-database/analysis/18S-rep-sequences-tax-search-results.qza



# Step 6 Taxonomic Visualization 
# Activate QIIME2/2024.2 and add Path Variables 
module load QIIME2/2024.2
INPUT=/home/gas51591/boba/updated-database/analysis/18S-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/updated-database/analysis/18S-rep-sequences-taxonomy.qzv

# Script
qiime metadata tabulate \
--m-input-file $INPUT \
--o-visualization $OUTPUT


# Step 7 Create Tree for Phylogentic Diversity Analysis 
# Activate QIIME2-2024.2 in Miniconda3 Environment 
module load Miniconda3
source activate /home/gas51591/conda-envs/envs/QIIME2-2024.2
INPUT=/home/gas51591/boba/updated-database/analysis/18S-rep-sequences-taxonomy.qza
OUTPUT=/home/gas51591/boba/updated-database/analysis/

# Script
qiime phylogeny  align-to-tree-mafft-fasttree \
--i-sequences /home/gas51591/boba/updated-database/analysis/denoised-rep-sequences.qza \
--o-alignment  /home/gas51591/boba/updated-database/analysis/aligned-18S-rep-seqs.qza \
--o-masked-alignment /home/gas51591/boba/updated-database/analysis/masked-aligned-rep-seqs.qza \
--o-tree /home/gas51591/boba/updated-database/analysis/unrooted-18S-tree.qza \
--o-rooted-tree /home/gas51591/boba/updated-database/analysis/rooted-18S-tree.qza