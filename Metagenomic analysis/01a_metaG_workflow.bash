# New workflow for metagenomics

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

tar -xvf $WDIR/X204SC23102800-Z01-F001.tar

# copy directories with raw reads into WDIR

# create files with sample names

nano sample_names.txt

E4_d458
G3_d458
L3_d458
M4_d458
Van_BES

mkdir logfiles

while read line
do
  mkdir -p $WDIR/${line}/raw_reads $WDIR/${line}/raw_reads_clean $WDIR/${line}/intermediate_results/clean_reads $WDIR/${line}/final_results
done < $WDIR/sample_names.txt

# save read files in raw_reads directory and unzip

while read line
do
  mv $WDIR/${line}/*.fq.gz $WDIR/${line}/raw_reads
  gunzip $WDIR/${line}/raw_reads/*".fq.gz"
done < $WDIR/sample_names.txt

# rename read files to Read_1.fq and Read_2.fq

while read line
do
  cat $WDIR/${line}/raw_reads/*1.fq > $WDIR/${line}/raw_reads/Read_1.fq
  cat $WDIR/${line}/raw_reads/*2.fq > $WDIR/${line}/raw_reads/Read_2.fq
done < $WDIR/sample_names.txt


# quality check using fastqc
module load fastqc/0.11.9

# Parameter
THREADS=80

while read line
do
  fastqc -o $WDIR/${line}/raw_reads/ -f fastq -t $THREADS $WDIR/${line}/raw_reads/"Read_1.fq" $WDIR/${line}/raw_reads/"Read_2.fq" >> $WDIR/logfiles/${line}_fastqc_raw_reads.log 2>&1
done < $WDIR/sample_names.txt


# multiqc

REPORT="MAC_MG"

conda activate multiqc
cd $WDIR
multiqc -d -dd -1 */raw_reads/ -n ${REPORT}"_raw" -o MultiQC/raw_reads/

# higher duplication levels in E4, G3 and L3
# per sequence GC content failed for E4 and L3
# per base sequence content shows wobbly bases in first 10 bp


# quality trimming
module load bbmap/38.86
module load fastqc/0.11.9
     
# Parameter
k=23
MINK=11
THREADS=80
TRIMQ=20  					
MIN=100   					
ftl=10
maq=10

while read line
do
  /storage/hdd2/mmaeke/Metagenomics/BBMap_trimming_olorin.sh ${WDIR} ${line} $k $MINK $THREADS $TRIMQ $MIN $ftl $maq
done < $WDIR/sample_names.txt


# repeat quality check
conda activate multiqc

REPORT="MAC_MG_clean"

cd $WDIR
multiqc -d -dd -1 */intermediate_results/clean_reads/ -n ${REPORT}"_clean" -o ./MultiQC/clean_reads/


# run dereplication
module load bbmap/38.86

# Parameter
DEDUPE=t
OPTICAL=f
DUPESUBS=0
DUPEDIST=40
deletetemp=t
THREADS=60

while read line
do
  clumpify.sh in=$WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_1.fq in2=$WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_2.fq out=$WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq out2=$WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq dedupe=$DEDUPE optical=$OPTICAL dupesubs=$DUPESUBS dupedist=$DUPEDIST deletetemp=t threads=$THREADS >> $WDIR/logfiles/${line}_clumpify.log 2>&1
done < $WDIR/sample_names.txt



# new quality check
module load fastqc/0.11.9

THREADS=100

while read line
do
  fastqc -o $WDIR/${line}/intermediate_results/clean_reads/ -f fastq -t $THREADS $WDIR/${line}/intermediate_results/clean_reads/"final_pure_reads_dedupe_1.fq" $WDIR/${line}/intermediate_results/clean_reads/"final_pure_reads_dedupe_2.fq" >> $WDIR/logfiles/${line}_fastqc_dedupe_reads.log 2>&1
done < $WDIR/sample_names.txt

# run multiqc for report
conda activate multiqc

REPORT="MAC_MG"

cd $WDIR
multiqc -d -dd -1 */intermediate_results/clean_reads/final_pure_reads_dedupe* -n ${REPORT}"_dedupe" -o ./MultiQC/dedupe_reads/


# nonpareil
conda activate nonpereil

cd $WDIR

# Parameter
THREADS=60 

while read line
do
  mkdir -p $WDIR/${line}/intermediate_results/clean_reads/nonpareil
  nonpareil -s $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -T kmer -f fastq -b $WDIR/${line}/intermediate_results/clean_reads/nonpareil/${line}_1 -t $THREADS >> $WDIR/logfiles/${line}"_nonpareil_R1.log" 2>&1
  nonpareil -s $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -T kmer -f fastq -b $WDIR/${line}/intermediate_results/clean_reads/nonpareil/${line}_2 -t $THREADS >> $WDIR/logfiles/${line}"_nonpareil_R2.log" 2>&1
done < $WDIR/sample_names.txt

# to parse nonpareil output copy all nonpareil .npo files into new directory
mkdir -p $WDIR/Nonpareil_overview
cp $WDIR/*/intermediate_results/clean_reads/nonpareil/*.npo $WDIR/Nonpareil_overview

module load R/4.2.2

$WDIR/nonpareil_parse.R -d $WDIR/Nonpareil_overview -s $WDIR/Nonpareil_overview/nonpareil_summary.txt -p $WDIR/Nonpareil_overview/nonpareil_plot.pdf


# Assembly

#megahit
conda activate megahit-1.2.9
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

# Parameter
MEMORY=0.5
THREADS=60

cd $WDIR

while read line
do
  mkdir -p ${line}/intermediate_results/Assembly
  megahit -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -m $MEMORY -t $THREADS --presets meta-sensitive -o $WDIR/${line}/intermediate_results/Assembly/Megahit  >> $WDIR/logfiles/${line}"_megahit.log" 2>&1
  rm -r $WDIR/${line}/intermediate_results/Assembly/Megahit/intermediate_contigs
done < $WDIR/sample_names.txt

conda deactivate

#spades
conda activate spades

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

cd $WDIR

# Parameter
MEMORY=400
THREADS=80

while read line
do
  spades.py -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq  -o $WDIR/${line}/intermediate_results/Assembly/SPAdes --meta -t $THREADS -m $MEMORY >> $WDIR/logfiles/${line}"_SPAdes.log" 2>&1
  rm -r $WDIR/${line}/intermediate_results/Assembly/SPAdes/K*  
done < $WDIR/sample_names.txt

while read line
do
   rm -r $WDIR/${line}/intermediate_results/Assembly/SPAdes/K*  
done < $WDIR/sample_names.txt

# mapping
conda activate anvio-7.1

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

cd $WDIR

while read line 
do
  anvi-script-reformat-fasta $WDIR/${line}/intermediate_results/Assembly/Megahit/final.contigs.fa -o $WDIR/${line}/intermediate_results/Assembly/Megahit/contigs_fixed.fa -l 1000 -r $WDIR/${line}/intermediate_results/Assembly/Megahit/fixed_contigs_report.txt --simplify-names >> $WDIR/logfiles/${line}"_megahit_simplify_fasta.log" 2>&1
  anvi-script-reformat-fasta $WDIR/${line}/intermediate_results/Assembly/SPAdes/scaffolds.fasta -o $WDIR/${line}/intermediate_results/Assembly/SPAdes/contigs_fixed.fa -l 1000 -r $WDIR/${line}/intermediate_results/Assembly/SPAdes/fixed_contigs_report.txt --simplify-names >> $WDIR/logfiles/${line}"_SPAdes_simplify_fasta.log" 2>&1
done < $WDIR/sample_names.txt

conda deactivate

# bowtie2
conda activate bowtie2-2.5.3
conda install -c bioconda samtools # into bowtie environment

# Parameter
THREADS=80
cd $WDIR
while read line 
do
  mkdir -p $WDIR/${line}/intermediate_results/Mapping/bowtie_megahit $WDIR/${line}/intermediate_results/Mapping/bowtie_spades
    /storage/hdd2/mmaeke/Metagenomics/Index_bowtie_olorin.sh ${line} 80
done < $WDIR/sample_names.txt

#Run read mapping for each combination of index and sample
cd $WDIR

while read line 
do
  while read SID 
  do
    /storage/hdd2/mmaeke/Metagenomics/Mapping_bowtie_olorin.sh ${line} ${SID} $THREADS
  done < $WDIR/sample_names.txt
done < $WDIR/sample_names.txt

cd $WDIR
while read IID 
do
  while read SID 
  do
    samtools sort -o ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".bam" ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".sam"
    samtools index ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".bam"
    
	samtools sort -o ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".bam" ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".sam"
    samtools index ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".bam"
  done < $WDIR/sample_names.txt
done < $WDIR/sample_names.txt


# Remove intermediate files
while read line 
do
  rm $WDIR/${line}/intermediate_results/Mapping/bowtie_megahit/contigs* 
  rm $WDIR/${line}/intermediate_results/Mapping/bowtie_spades/contigs*
  rm $WDIR/${line}/intermediate_results/Mapping/bowtie_megahit/*sam 
  rm $WDIR/${line}/intermediate_results/Mapping/bowtie_spades/*sam
done < $WDIR/sample_names.txt


# run metabat2 and CONCOCT binning
conda activate metawrap-env

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

while read line
do
  mkdir -p $WDIR/${line}/intermediate_results/Binning
done < $WDIR/sample_names.txt

mkdir -p $WDIR/CoAssembly/Binning

# Parameter
MINCONTIG=2000 
CUTOFF=2500    
CHUNK=10000
OVERLAP=0
THREADS=100

while read line
do
  cd $WDIR
  /storage/hdd2/mmaeke/Metagenomics/Binning_CONCOCT_olorin.bash  $WDIR/${line}/intermediate_results/Assembly $WDIR/${line}/intermediate_results/Binning $WDIR/${line}/intermediate_results/Mapping ${line} $MINCONTIG $CHUNK $OVERLAP $THREADS >> $WDIR/logfiles/${line}"_CONCOCT.log" 2>&1

  /storage/hdd2/mmaeke/Metagenomics/Binning_Metabat2_olorin.bash  $WDIR/${line}/intermediate_results/Assembly $WDIR/${line}/intermediate_results/Binning $WDIR/${line}/intermediate_results/Mapping ${line} $CUTOFF $THREADS >> $WDIR/logfiles/${line}"_metabat2.log" 2>&1
done < $WDIR/sample_names.txt

# run vamb

mkdir -p $WDIR/intermediate_results/vamb/SPAdes/Mapping $WDIR/intermediate_results/vamb/Megahit/Mapping $WDIR/intermediate_results/vamb/coassembly/Mapping

# create file lists as input for concatenate.py
conda activate vamb-4.1.3

ASSEMBLY_spades=$(ls -1 $WDIR/*/intermediate_results/Assembly/SPAdes/contigs_fixed.fa )
ASSEMBLY_megahit=$(ls -1 $WDIR/*/intermediate_results/Assembly/Megahit/contigs_fixed.fa )


python /opt/moep/vamb/4.1.3/vamb/src/concatenate.py $WDIR/intermediate_results/vamb/SPAdes/Mapping/catalogue.fasta $ASSEMBLY_spades --nozip
python /opt/moep/vamb/4.1.3/vamb/src/concatenate.py $WDIR/intermediate_results/vamb/Megahit/Mapping/catalogue.fasta $ASSEMBLY_megahit --nozip

conda deactivate

# create index
conda activate bowtie2-2.5.3

bowtie2-build $WDIR/intermediate_results/vamb/SPAdes/Mapping/catalogue.fasta $WDIR/intermediate_results/vamb/SPAdes/Mapping/contigs
bowtie2-build $WDIR/intermediate_results/vamb/Megahit/Mapping/catalogue.fasta $WDIR/intermediate_results/vamb/Megahit/Mapping/contigs

# run mapping

NCPU=80
while read line 
do
   bowtie2 --threads $NCPU -x $WDIR/intermediate_results/vamb/SPAdes/Mapping/contigs -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -S $WDIR/intermediate_results/vamb/SPAdes/Mapping/${line}".sam" 
   
   samtools view $WDIR/intermediate_results/vamb/SPAdes/Mapping/${line}".sam" -F 3584 -b --threads 40 -o $WDIR/intermediate_results/vamb/SPAdes/Mapping/${line}".bam"
   samtools sort $WDIR/intermediate_results/vamb/SPAdes/Mapping/${line}".bam" -o $WDIR/intermediate_results/vamb/SPAdes/Mapping/${line}"_sorted.bam"
   
   bowtie2 --threads $NCPU -x $WDIR/intermediate_results/vamb/Megahit/Mapping/contigs -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -S $WDIR/intermediate_results/vamb/Megahit/Mapping/${line}".sam" 
   samtools view $WDIR/intermediate_results/vamb/Megahit/Mapping/${line}".sam" -F 3584 -b --threads 40 -o $WDIR/intermediate_results/vamb/Megahit/Mapping/${line}".bam"
   samtools sort $WDIR/intermediate_results/vamb/Megahit/Mapping/${line}".bam" -o $WDIR/intermediate_results/vamb/Megahit/Mapping/${line}"_sorted.bam"
done < $WDIR/sample_names.txt


# run binning with vamb
conda activate vamb-4.1.3

vamb bin default -l 24 -n 512 512 --outdir $WDIR/intermediate_results/vamb/Megahit/Binning --fasta $WDIR/intermediate_results/vamb/Megahit/Mapping/catalogue.fasta --bamfiles $WDIR/intermediate_results/vamb/Megahit/Mapping/*sorted.bam -o  --minfasta 500000

vamb bin default -l 24 -n 512 512 --outdir $WDIR/intermediate_results/vamb/SPAdes/Binning --fasta $WDIR/intermediate_results/vamb/SPAdes/Mapping/catalogue.fasta --bamfiles $WDIR/intermediate_results/vamb/SPAdes/Mapping/*sorted.bam -o C --minfasta 500000


# sort bins based on sample names into corresponding directories
# Sample order
# E4_d458
# G3_d458
# L3_d458
# M4_d458
# Van_BES

SAMPLE_1="E4_d458"
SAMPLE_2="G3_d458"
SAMPLE_3="L3_d458"
SAMPLE_4="M4_d458"
SAMPLE_5="Van_BES"

while read line
do
   mkdir -p $WDIR/${line}/intermediate_results/Binning/vamb_megahit $WDIR/${line}/intermediate_results/Binning/vamb_spades
done < $WDIR/sample_names.txt

cp $WDIR/intermediate_results/vamb/Megahit/Binning/bins/S1* $WDIR/$SAMPLE_1/intermediate_results/Binning/vamb_megahit
cp $WDIR/intermediate_results/vamb/SPAdes/Binning/bins/S1* $WDIR/$SAMPLE_1/intermediate_results/Binning/vamb_spades
cp $WDIR/intermediate_results/vamb/Megahit/Binning/bins/S2* $WDIR/$SAMPLE_2/intermediate_results/Binning/vamb_megahit
cp $WDIR/intermediate_results/vamb/SPAdes/Binning/bins/S2* $WDIR/$SAMPLE_2/intermediate_results/Binning/vamb_spades
cp $WDIR/intermediate_results/vamb/Megahit/Binning/bins/S3* $WDIR/$SAMPLE_3/intermediate_results/Binning/vamb_megahit
cp $WDIR/intermediate_results/vamb/SPAdes/Binning/bins/S3* $WDIR/$SAMPLE_3/intermediate_results/Binning/vamb_spades
cp $WDIR/intermediate_results/vamb/Megahit/Binning/bins/S4* $WDIR/$SAMPLE_4/intermediate_results/Binning/vamb_megahit
cp $WDIR/intermediate_results/vamb/SPAdes/Binning/bins/S4* $WDIR/$SAMPLE_4/intermediate_results/Binning/vamb_spades
cp $WDIR/intermediate_results/vamb/Megahit/Binning/bins/S5* $WDIR/$SAMPLE_5/intermediate_results/Binning/vamb_megahit
cp $WDIR/intermediate_results/vamb/SPAdes/Binning/bins/S5* $WDIR/$SAMPLE_5/intermediate_results/Binning/vamb_spades

# Bin refinement
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
cd $WDIR

while read line
do
  mkdir -p $WDIR/${line}/intermediate_results/Refined_Bins/megahit
  mkdir -p $WDIR/${line}/intermediate_results/Refined_Bins/spades
  mkdir -p $WDIR/CoAssembly/Refined_Bins
done < $WDIR/sample_names.txt

module load metawrap
conda activate metawrap-env

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
cd $WDIR

# Parameter
THREADS=80
COMPLETENESS=60
CONTAMINATION=10

# single samples
while read line
do
  metawrap bin_refinement -o $WDIR/${line}/intermediate_results/Refined_Bins/spades -t $THREADS -A $WDIR/${line}/intermediate_results/Binning/concoct_spades/concoct_output/fasta_bins -B $WDIR/${line}/intermediate_results/Binning/metabat_spades -C $WDIR/${line}/intermediate_results/Binning/vamb_spades -c $COMPLETENESS -x $CONTAMINATION >> $WDIR/logfiles/${line}"_bin_refinement_spades.log" 2>&1

  metawrap bin_refinement -o $WDIR/${line}/intermediate_results/Refined_Bins/megahit -t $THREADS -A $WDIR/${line}/intermediate_results/Binning/concoct_megahit/concoct_output/fasta_bins -B $WDIR/${line}/intermediate_results/Binning/metabat_megahit -C $WDIR/${line}/intermediate_results/Binning/vamb_megahit -c $COMPLETENESS -x $CONTAMINATION >> $WDIR/logfiles/${line}"_bin_refinement_megahit.log" 2>&1
done < $WDIR/sample_names.txt



#########
# processing of the coassembly
mkdir -p $WDIR/CoAssembly/Assembly_megahit

#list
READ_1=$(ls -m $WDIR/*/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq | sed 's/[[:space:]]//g')
READ_2=$(ls -m $WDIR/*/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq | sed 's/[[:space:]]//g')

conda activate megahit-1.2.9
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

cd $WDIR

# Parameter
MEMORY=0.75
THREADS=100

#meta-sensitive
megahit -1 /storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/E4_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/G3_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/L3_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/M4_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/Van_BES/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 /storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/E4_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/G3_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/L3_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/M4_d458/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq,/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG/Van_BES/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -m $MEMORY -t $THREADS --presets meta-sensitive -o $WDIR/CoAssembly/Assembly_megahit/metasensitive  >> $WDIR/logfiles/"CoAssembly_metasensitive.log" 2>&1


# readmapping with bowtie2
conda activate anvio-7.1

anvi-script-reformat-fasta $WDIR/CoAssembly/Assembly_megahit/metasensitive/final.contigs.fa -o $WDIR/CoAssembly/Assembly_megahit/metasensitive/contigs_fixed.fa -l 1000 -r $WDIR/CoAssembly/Assembly_megahit/metasensitive/fixed_contigs_report.txt --simplify-names >> $WDIR/logfiles/"CoAssembly_metasensitive_simplify_fasta.log" 2>&1
conda deactivate

conda activate bowtie2-2.5.3
mkdir -p $WDIR/CoAssembly/Mapping

while read line
do
  bowtie2-build $WDIR/CoAssembly/Assembly_megahit/metasensitive/contigs_fixed.fa $WDIR/CoAssembly/Mapping/contigs >> $WDIR/logfiles/"CoAssembly_metasensitive_indexing.log" 2>&1
  /storage/hdd2/mmaeke/Metagenomics/PhD/Mapping_bowtie_CoAssembly.bash ${line} $THREADS  >> $WDIR/logfiles/"CoAssembly_metasensitive_"${line}"_mapping.log" 2>&1
done < $WDIR/sample_names.txt

conda deactivate

# binning coassembly
conda activate metawrap-env
mkdir -p $WDIR/CoAssembly/Binning
# concoct
/storage/hdd2/mmaeke/Metagenomics/PhD/CONCOCT_olorin_CoAssembly.bash  $WDIR/CoAssembly/Assembly_megahit/metasensitive $WDIR/CoAssembly/Binning $WDIR/CoAssembly/Mapping  $MINCONTIG $CHUNK $OVERLAP $THREADS >> $WDIR/logfiles/CoAssembly_CONCOCT.log 2>&1  
# metabat
/storage/hdd2/mmaeke/Metagenomics/PhD/Metabat2_olorin_CoAssembly.bash  $WDIR/CoAssembly/Assembly_megahit/metasensitive $WDIR/CoAssembly/Binning $WDIR/CoAssembly/Mapping  $CUTOFF $THREADS >> $WDIR/logfiles/CoAssembly_metabat2.log 2>&1  

# vamb binning
mkdir -p $WDIR/CoAssembly/Mapping_vamb

conda activate vamb-4.1.3

python /opt/moep/vamb/4.1.3/vamb/src/concatenate.py $WDIR/CoAssembly/Mapping_vamb/catalogue.fasta $WDIR/CoAssembly/Assembly_megahit/metasensitive/contigs_fixed.fa --nozip

conda deactivate

# create index
conda activate bowtie2-2.5.3

bowtie2-build $WDIR/CoAssembly/Mapping_vamb/catalogue.fasta $WDIR/CoAssembly/Mapping_vamb/contigs

# run mapping
NCPU=80
while read line 
do
   bowtie2 --threads $NCPU -x $WDIR/CoAssembly/Mapping_vamb/contigs -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -S $WDIR/CoAssembly/Mapping_vamb/${line}".sam" 
   
   samtools view $WDIR/CoAssembly/Mapping_vamb/${line}".sam" -F 3584 -b --threads 60 -o $WDIR/CoAssembly/Mapping_vamb/${line}".bam"
   samtools sort $WDIR/CoAssembly/Mapping_vamb/${line}".bam" --threads 60 -o $WDIR/CoAssembly/Mapping_vamb/${line}"_sorted.bam"
done < $WDIR/sample_names.txt

# run binning with vamb
conda activate vamb-4.1.3
vamb bin default -l 24 -n 512 512 --outdir $WDIR/CoAssembly/Binning/vamb --fasta $WDIR/CoAssembly/Mapping_vamb/catalogue.fasta --bamfiles $WDIR/CoAssembly/Mapping_vamb/*sorted.bam -o C --minfasta 500000

# Bin refinement Co-assembly
metawrap bin_refinement -o $WDIR/CoAssembly/Refined_Bins -t $THREADS -A $WDIR/CoAssembly/Binning/concoct/concoct_output/fasta_bins -B $WDIR/CoAssembly/Binning/metabat -C $WDIR/CoAssembly/Binning/vamb/bins -c $COMPLETENESS -x $CONTAMINATION >> $WDIR/logfiles/"CoAssembly_bin_refinement.log" 2>&1
########


# Remove working files metawrap bin refinement
while read line
do
 rm -r $WDIR/${line}/intermediate_results/Refined_Bins/megahit/work_files $WDIR/${line}/intermediate_results/Refined_Bins/megahit/fasta_bins $WDIR/${line}/intermediate_results/Refined_Bins/megahit/metabat_megahit $WDIR/${line}/intermediate_results/Refined_Bins/megahit/vamb_megahit
 rm -r $WDIR/${line}/intermediate_results/Refined_Bins/spades/work_files $WDIR/${line}/intermediate_results/Refined_Bins/spades/fasta_bins $WDIR/${line}/intermediate_results/Refined_Bins/spades/metabat_spades $WDIR/${line}/intermediate_results/Refined_Bins/spades/vamb_spades
 rm -r $WDIR/CoAssembly/intermediate_results/Refined_Bins/work_files $WDIR/CoAssembly/intermediate_results/Refined_Bins/metabat $WDIR/CoAssembly/intermediate_results/Refined_Bins/fasta_bins $$WDIR/CoAssembly/intermediate_results/Refined_Bins/vamb
done < $WDIR/sample_names.txt

# Bin dereplication for Coassembly and single sample assemblies
#These following commands will copy and rename the bins into a new directory called `input_bins`

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
mkdir -p $WDIR/intermediate_results/Bin_dereplication/input_bins
cd $WDIR/intermediate_results/Bin_dereplication/input_bins

ls -1 $WDIR/CoAssembly/Refined_Bins/metawrap_60_10_bins/ | xargs -n1 basename | while read line
do
  sed "/^>/s/^>/>Co_/" $WDIR/CoAssembly/Refined_Bins/metawrap_60_10_bins/${line} > $WDIR/intermediate_results/Bin_dereplication/input_bins/Co_${line}
done

cat $WDIR/sample_names.txt | while read SID
do
 ls -1 $WDIR/${SID}/intermediate_results/Refined_Bins/spades/metawrap_60_10_bins/ | xargs -n1 basename | while read line
 do
	sed "/^>/s/^>/>${SID}_spades_/" $WDIR/${SID}/intermediate_results/Refined_Bins/spades/metawrap_60_10_bins/${line} > $WDIR/intermediate_results/Bin_dereplication/input_bins/${SID}_spades_${line}
 done
done

cat $WDIR/sample_names.txt | while read SID
do
 ls -1 $WDIR/${SID}/intermediate_results/Refined_Bins/megahit/metawrap_60_10_bins/ | xargs -n1 basename | while read line
 do
	sed "/^>/s/^>/>${SID}_megahit_/" $WDIR/${SID}/intermediate_results/Refined_Bins/megahit/metawrap_60_10_bins/${line} > $WDIR/intermediate_results/Bin_dereplication/input_bins/${SID}_megahit_${line}
 done
done

cd $WDIR/intermediate_results/Bin_dereplication/input_bins/
ls -d "$PWD"/* > ../bin_list.txt

# combine all checkm stats files as input for genome info
# create file with adding Co_assembly bins
sed -e '1s/bin/genome/' -e '2,$ s/^/Co_/' $WDIR/CoAssembly/Refined_Bins/metawrap_60_10_bins.stats | cut -f 1-3  >  $WDIR/intermediate_results/Bin_dereplication/checkm_stats_tmp.csv 
# now add stats for all other bin refinements
while read line
do 
   sed "2,$ s/^/${line}_spades_/" $WDIR/${line}/intermediate_results/Refined_Bins/spades/metawrap_60_10_bins.stats | sed '1d' | cut -f 1-3  >> $WDIR/intermediate_results/Bin_dereplication/checkm_stats_tmp.csv
   sed "2,$ s/^/${line}_megahit_/" $WDIR/${line}/intermediate_results/Refined_Bins/megahit/metawrap_60_10_bins.stats | sed '1d' | cut -f 1-3  >> $WDIR/intermediate_results/Bin_dereplication/checkm_stats_tmp.csv  
done < $WDIR/sample_names.txt

# run dereplication (drep v3.2.2)
cat $WDIR/intermediate_results/Bin_dereplication/checkm_stats_tmp.csv  | awk  'NR>=2 {$1=$1".fa"}1' | sed -e 's/\t/,/g' -e 's/ /,/g' > $WDIR/intermediate_results/Bin_dereplication/checkm_stats.csv
rm $WDIR/intermediate_results/Bin_dereplication/checkm_stats_tmp.csv

# parameters
THREADS=100
MINLENGTH=50000
COMPLETENESS=60
CONTAMINATION=10
PANI=0.90
SANI=0.95
MINOVERLAP=0.5

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

conda activate drep-3.0.0 # runs drep 3.2.2
dRep dereplicate -p $THREADS -g $WDIR/intermediate_results/Bin_dereplication/bin_list.txt --length $MINLENGTH -comp $COMPLETENESS -con $CONTAMINATION --S_algorithm ANImf -pa $PANI -sa $SANI -nc $MINOVERLAP  --genomeInfo $WDIR/intermediate_results/Bin_dereplication/checkm_stats.csv $WDIR/intermediate_results/Bin_dereplication/drep_out >> $WDIR/logfiles/"dRep_dereplication.log" 2>&1
conda deactivate


# Adjust dereplication of drep
# run coverm for single sample assemblies

conda activate coverm-0.6.1
while read line
do
   coverm genome -b $WDIR/${line}/intermediate_results/Mapping/bowtie_megahit/*.bam -m mean --min-read-percent-identity 95 --min-read-aligned-length 50 --exclude-supplementary -d $WDIR/${line}/intermediate_results/Refined_Bins/megahit/metawrap_60_10_bins -x fa -t 80 > $WDIR/${line}/intermediate_results/Refined_Bins/megahit/${line}"_bins_megahit_quant_tmp.txt"

   coverm genome -b $WDIR/${line}/intermediate_results/Mapping/bowtie_spades/*.bam -m mean --min-read-percent-identity 95 --min-read-aligned-length 50 --exclude-supplementary -d $WDIR/${line}/intermediate_results/Refined_Bins/spades/metawrap_60_10_bins -x fa -t 80 > $WDIR/${line}/intermediate_results/Refined_Bins/spades/${line}"_bins_spades_quant_tmp.txt"
done < $WDIR/sample_names.txt

# run coverm for coassembly
coverm genome -b $WDIR/CoAssembly/Mapping/*.bam -m mean --min-read-percent-identity 95 --min-read-aligned-length 50 --exclude-supplementary -d $WDIR/CoAssembly/Refined_Bins/metawrap_60_10_bins -x fa -t 80 > $WDIR/CoAssembly/Refined_Bins/Co_bins_quant_tmp.txt

# adjust files and add sample name to each line
while read line
do
   sed "2,$ s/^/${line}_megahit_/" $WDIR/${line}/intermediate_results/Refined_Bins/megahit/${line}"_bins_megahit_quant_tmp.txt" > $WDIR/${line}/intermediate_results/Refined_Bins/megahit/${line}"_bins_quant.txt"
   sed "2,$ s/^/${line}_spades_/" $WDIR/${line}/intermediate_results/Refined_Bins/spades/${line}"_bins_spades_quant_tmp.txt" > $WDIR/${line}/intermediate_results/Refined_Bins/spades/${line}"_bins_quant.txt"
done < $WDIR/sample_names.txt

sed '2,$ s/^/Co_/' $WDIR/CoAssembly/Refined_Bins/Co_bins_quant_tmp.txt  >  $WDIR/CoAssembly/Refined_Bins/Co_bins_quant.txt


# run reassembly with adjsuted dRep bins based on Christiane Hassenrück's workflow. This workflow here starts after the dereplication
# for this you can find a snakemake workflow in /storage/hdd2/mmaeke/Metagenomics/PhD/metaG_reassembly, containing all neccessary config files, rules and scripts.
# The config file and asset files of your sample paths you need to adjust yourself

# create your asset file
# Define the output file
output_file="sample_list.txt"
# create a sample_list.txt file holding sample name and paths to raw reads (these won't be used lateron, but are required as input)
while IFS= read -r line || [[ -n "$line" ]]; do
    # Construct paths for read 1 and read 2
    read1_path=$WDIR/${line}/raw_reads/Read_1.fq
    read2_path=$WDIR/${line}/raw_reads/Read_2.fq

    # Print the sample name and paths to the output file
    echo -e "$line\t$read1_path\t$read2_path" >> /storage/hdd2/mmaeke/Metagenomics/PhD/metaG_reassembly/assets/${output_file}
done < $WDIR/sample_names.txt

conda activate avamb
snakemake --configfile /storage/hdd2/mmaeke/Metagenomics/PhD/metaG_reassembly/config/config.yaml --snakefile /storage/hdd2/mmaeke/Metagenomics/PhD/metaG_reassembly/Snakefile -j 12 --use-conda 

# gtdb_tk for taxonomic classification of final bins
conda activate gtdbtk-2.1.0
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

gtdbtk classify_wf --genome_dir $WDIR/intermediate_results/Bin_final/bins --out_dir $WDIR/intermediate_results/Bin_final/gtdbtk_out --cpus 80 -x fa  >> $WDIR/logfiles/"final_bins_classify_gtdbtk.log" 2>&1

# coverm for relative abundance of bins
# mapping
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

mkdir -p $WDIR/intermediate_results/Bin_final/coverm/output $WDIR/intermediate_results/Bin_final/coverm/input
cat $WDIR/intermediate_results/Bin_final/bins/*.fa > $WDIR/intermediate_results/Bin_final/coverm/input/concat_bins.fa 

conda activate bowtie2-2.5.3
bowtie2-build $WDIR/intermediate_results/Bin_final/coverm/input/concat_bins.fa  $WDIR/intermediate_results/Bin_final/coverm/input/contigs

while read line
do
 bowtie2 --threads 80 -x $WDIR/intermediate_results/Bin_final/coverm/input/contigs -1 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq -2 $WDIR/${line}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq -S $WDIR/intermediate_results/Bin_final/coverm/input/${line}".sam"
 samtools sort -o $WDIR/intermediate_results/Bin_final/coverm/input/${line}".bam" $WDIR/intermediate_results/Bin_final/coverm/input/${line}".sam" >> $WDIR/logfiles/"final_bins_mapping.log" 2>&1
 samtools index $WDIR/intermediate_results/Bin_final/coverm/input/${line}".bam" >> $WDIR/logfiles/"final_bins_mapping.log" 2>&1
done < $WDIR/sample_names.txt


# coverm relative abundance
conda activate coverm-0.6.1

coverm genome -b $WDIR/intermediate_results/Bin_final/coverm/input/*.bam -m relative_abundance --min-read-percent-identity 95 --min-read-aligned-length 50 --exclude-supplementary -d $WDIR/intermediate_results/Bin_final/bins -x fa -t 80 > $WDIR/intermediate_results/Bin_final/coverm/output/bin_rel_abund_coverm.txt
coverm genome -b $WDIR/intermediate_results/Bin_final/coverm/input/*.bam -m mean --min-read-percent-identity 95 --min-read-aligned-length 50 --exclude-supplementary -d $WDIR/intermediate_results/Bin_final/bins -x fa -t 80 > $WDIR/intermediate_results/Bin_final/coverm/output/bin_mean_cov_coverm.txt