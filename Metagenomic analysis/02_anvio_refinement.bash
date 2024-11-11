# anvio script to refine EX MAGs found in metagenomic studies

# test with one sample: PRJNA713414

## Workflow for anvio 7.1
## This workflow is based on the assumption that you already did your bin reassembly and want to refine your target bins.

# Prepare files
# Create a concatenated assembly by combining all your reassembled bins into one file
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG" #your WDIR should be WDIR="/storage/hdd5/lwu/Metagenomes/Sulfurimonas/MG/SulfMG", with this i think it should run as it is here now

# Change fasta header for anvio while creating a concatenated assembly file 
mkdir -p $WDIR/intermediate_results/Bin_final/anvio_refine/input $WDIR/intermediate_results/Bin_final/anvio_refine/output
ls -1 $WDIR/intermediate_results/Bin_final/bins/*.fa | xargs -n 1 basename | while read SID
do
   sed -e '/^>/s/\./_/' -e "/^>/s/[^[:alnum:]>-]/_/g"  $WDIR/intermediate_results/Bin_final/bins/${SID}
done >> $WDIR/intermediate_results/Bin_final/anvio_refine/input/concat_assembly.fasta

# create a list with all contigs and corresponding bins
ls -1 $WDIR/intermediate_results/Bin_final/bins/  | xargs -n 1 basename | sed 's/\.fa//'  | while read SID
do
   grep '^>' $WDIR/intermediate_results/Bin_final/bins/${SID}".fa" | sed  -e "/^>/s/[^[:alnum:]>-]/_/g"  -e 's/^>//' -e "s/$/\t${SID}/" -e 's/\./_/g'
done >> $WDIR/intermediate_results/Bin_final/anvio_refine/input/binning_results.txt


PROJECT="MAC_metaG"

# Run mapping of all clean reads on the concatenated bin assembly
conda activate bowtie2-2.5.3
bowtie2-build $WDIR/intermediate_results/Bin_final/anvio_refine/input/concat_assembly.fasta $WDIR/intermediate_results/Bin_final/anvio_refine/input/contigs
bowtie2 --threads 100 -x $WDIR/intermediate_results/Bin_final/anvio_refine/input/contigs  -1 $WDIR/intermediate_results/clean_reads/all_R1.fastq -2 $WDIR/intermediate_results/clean_reads/all_R2.fastq -S $WDIR/intermediate_results/Bin_final/anvio_refine/input/${PROJECT}".sam "
samtools sort -m 3G -o $WDIR/intermediate_results/Bin_final/anvio_refine/input/${PROJECT}".bam" $WDIR/intermediate_results/Bin_final/anvio_refine/input/${PROJECT}".sam"
samtools index $WDIR/intermediate_results/Bin_final/anvio_refine/input/${PROJECT}".bam"


## Create a contigs DB
# While creating a contigs database k-mer frequencies are caclulated (default 4, change by using --kmer-size),
# contigs longer than 20,000 bp are split into smaller contigs (change size by using --split-length), 
# ORFs are called by using Prodigal (skip by using --skip-gene-calling)and contig splitting will cut mindfully when
# using Prodigal, so that no geens are cut in the middle.

conda activate anvio-7.1
anvi-gen-contigs-database -f $WDIR/intermediate_results/Bin_final/anvio_refine/input/concat_assembly.fasta -o $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -n 'Contigs database ${PROJECT}' -T 80

## Run anvio hmm profiles
# With anvi'o there are already some hmm models supplied, constituting of multiple published bacterial single-copy gene collections
# Further you can run your own HMM profile (--hmm-profile-dir) or only a specific installed HMM profile (--installed-hmm-profile)

anvi-run-hmms -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db --also-scan-trnas -T 80

# Run taxonomy
anvi-run-ncbi-cogs -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db --sensitive -T 80 

# Get taxonomy
anvi-run-scg-taxonomy -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -T 80
-------------
## Profiling BAM files
# The anvi'o profile db contains specific information about contigs.
# Profiling a bam file creates a single profile reporting properties for each contig in a sample based on mapping.
# Each profile.db links to a contig.db

anvi-profile -i $WDIR/intermediate_results/Bin_final/anvio_refine/input/${PROJECT}".bam" -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db --output-dir $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db --sample-name ${PROJECT} --min-contig-length 2500 --profile-SCVs -T 80

# Import own binning
# binning_results.txt should contain a TAB-delimited file containing infos about which contig belongs to what bin.
anvi-import-collection $WDIR/intermediate_results/Bin_final/anvio_refine/input/binning_results.txt -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} --contigs-mode

# Estimate completeness of bins
anvi-estimate-genome-completeness -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} 

# Estimate taxonomy of bins
anvi-estimate-scg-taxonomy -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} --compute-scg-coverages

# Estimate metabolism (2 h / 1,5 h)
anvi-run-kegg-kofams -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -T 80  
anvi-estimate-metabolism -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} -O ${PROJECT} --include-metadata --metagenome-mode


# This part you can run outside of a screen. It is important, that you ran your gtdbtk classification before, so that you know the bins you want to refine
# in your case it probably makes sense to refine all of them, just so it is done.
# Here I am only refining all archaeal and Clostridia bins, as these are of interest for me
# Refine MAGs
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
conda activate anvio-7.1
cd $WDIR
# use this to see a list of all bins. Remember, that you changed your bin names slightly at the beginning, as anvio required '_' instead of '.' in the bin names
anvi-show-collections-and-bins -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db


# Now let's get into it. You refine the bin manually in the web browser. What you need to look for is differences in the GC content of single contigs.
# If you see some contigs having a much higher or much lower GC these very likely don't belong to your bin.
# Anvio offers you a real-time completeness and contamination calculation, so you can track if removing contigs changed anything.
anvi-refine -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} --server-only -P 8080 -b M4_d458_spades_bin_17_orig

http://localhost:8080

anvi-summarize -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -o $WDIR/intermediate_results/Bin_final/anvio_refine/summary -C ${PROJECT}

# some high quality MAGs still had quite high contamination, therefore refien again

# This part you can run outside of a screen. It is important, that you ran your gtdbtk classification before, so that you know the bins you want to refine
# in your case it probably makes sense to refine all of them, just so it is done.
# Here I am only refining all archaeal and Clostridia bins, as these are of interest for me
# Refine MAGs
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
conda activate anvio-7.1
cd $WDIR
# use this to see a list of all bins. Remember, that you changed your bin names slightly at the beginning, as anvio required '_' instead of '.' in the bin names
anvi-show-collections-and-bins -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db

# Co_bin_88_orig, G3_d458_spades_bin_29_indper, Co_bin_49_orig

# Now let's get into it. You refine the bin manually in the web browser. What you need to look for is differences in the GC content of single contigs.
# If you see some contigs having a much higher or much lower GC these very likely don't belong to your bin.
# Anvio offers you a real-time completeness and contamination calculation, so you can track if removing contigs changed anything.
anvi-refine -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -C ${PROJECT} --server-only -P 8080 -b G3_d458_spades_bin_29_indper

http://localhost:8080

anvi-summarize -p $WDIR/intermediate_results/Bin_final/anvio_refine/profile_db/PROFILE.db -c $WDIR/intermediate_results/Bin_final/anvio_refine/output/CONTIGS.db -o $WDIR/intermediate_results/Bin_final/anvio_refine/summary_v2 -C ${PROJECT}


cp summary_v2/bin_by_bin/G3_d458_spades_bin_29_indper_refined/G3_d458_spades_bin_29_indper_refined-contigs.fa ../../../Clostridia/additional_bins/
 
conda activate checkm2-1.0.2
checkm2 predict --threads 80 -x fa --input $WDIR/${TAXA}/additional_bins --output-directory $WDIR/${TAXA}/additional_bins/checkm2_out --allmodels

