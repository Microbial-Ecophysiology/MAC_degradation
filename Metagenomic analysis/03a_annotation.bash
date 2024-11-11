# Clostridia analysis
# now that all Clostridia bins were refined, we are starting all Clostridia related analyses. For this I am copying all Clostridia bins to the new directory

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
TAXA="Clostridia"

mkdir -p $WDIR/${TAXA}/bins_all
# create bins.txt containing a list of all MAGs

while read SID
do
  cp $WDIR/intermediate_results/Bin_final/anvio_refine/summary/bin_by_bin/${SID}*/*contigs.fa $WDIR/${TAXA}/bins
done < $WDIR/${TAXA}/bins.txt

# rerun checkm2 for bins
conda activate checkm2-1.0.2
checkm2 predict --threads 60 -x fa --input $WDIR/${TAXA}/bins_all --output-directory $WDIR/${TAXA}/checkm2_out --allmodels

# select only MAGs with completeness >80%, contamination <5%
awk -v FS="\t" -v OFS="\t" '
    # Check if the Completeness_Model_Used (column 5) is "Gradient Boost"
    $5 == "Gradient Boost (General Model)" {
        # If the model used is "Gradient Boost", check Completeness_General (column 2)
        if ($2 >= 80 && $3 <= 5) {
            print $0
        }
    }
    # Check if the Completeness_Model_Used (column 5) is "Neural Network"
    $5 == "Neural Network (Specific Model)" {
        # If the model used is "Neural Network", check Completeness_Specific (column 4)
        if ($4 >= 80 && $3 <= 5) {
            print $0
        }
    }
' "$WDIR/${TAXA}/checkm2_out/quality_report.tsv" > "$WDIR/${TAXA}/quality_filtered_MAGs.txt"

# 37 MAGs (of 46) meet quality statistics

# copy bins to directory bins_filtered
mkdir -p $WDIR/${TAXA}/bins_filtered

cat $WDIR/${TAXA}/quality_filtered_MAGs.txt | cut -f1 | while read SID
do
   cp $WDIR/${TAXA}/bins_all/${SID}.fa  $WDIR/${TAXA}/bins_filtered/${SID}.fa
done


# rerun gtdb_tk after refinement
conda activate gtdbtk-2.1.0
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"

gtdbtk classify_wf --genome_dir $WDIR/${TAXA}/bins_filtered --out_dir $WDIR/${TAXA}/gtdbtk_out --cpus 40 -x fa  >> $WDIR/logfiles/"Clostridia_bins_classify_gtdbtk.log" 2>&1



# run annotation of all bins

# gene prediction
mkdir -p $WDIR/${TAXA}/annotation/prodigal_out

conda activate anvio-7.1
ls -1 $WDIR/${TAXA}/bins_filtered| sed 's/\.fa//' | while read SID
do
  prodigal -i $WDIR/${TAXA}/bins_filtered/${SID}".fa" -o $WDIR/${TAXA}/annotation/prodigal_out/${SID}_coords.gbk -a $WDIR/${TAXA}/annotation/prodigal_out/${SID}.faa
done


### Run diamond scan with KEGG against MAGs of interest ###
# In /storage/hdd2/mmaeke/databases/KEGG/release_20231231 already present are scan db and pepunit db (both of which are required)

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
TAXA="Clostridia"
DB="/storage/hdd2/mmaeke/databases/KEGG/release_20231231"

conda activate checkm2-1.0.2

mkdir -p $WDIR/${TAXA}/annotation/kegg/tmp
mkdir -p $WDIR/${TAXA}/annotation/kegg/kegg_annotation_dmnd/
mkdir -p $WDIR/${TAXA}/annotation/selfblast/
mkdir -p $WDIR/${TAXA}/annotation/kegg/kegg_annotation_parsed/

THREADS=40

# scan DB
# use input from prodigal

ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond/kegg_genes.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out $WDIR/${TAXA}/annotation/kegg/${line}".kegg.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/${TAXA}/annotation/kegg/tmp
  diamond makedb --in $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa -d $WDIR/${TAXA}/annotation/kegg/kegg_annotation_dmnd/cds_db_${line}
  diamond blastp -d $WDIR/${TAXA}/annotation/kegg/kegg_annotation_dmnd/cds_db_${line}.dmnd -q $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa -k 1 -o $WDIR/${TAXA}/annotation/selfblast/${line}".selfblast" -f 6
done
 
# parse annotation
module load R/4.2.2
ls -d $WDIR/${TAXA}/annotation/kegg/*.kegg.annotation | sed "s/\.kegg.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_kegg.R -b $WDIR/${TAXA}/annotation/kegg/${line}".kegg.annotation" -s $WDIR/${TAXA}/annotation/selfblast/${line}".selfblast" -m $DB/mapping_info -c 0.4 -t 120 -o $WDIR/${TAXA}/annotation/kegg/kegg_annotation_parsed_new/${line}".kegg.parsed"
done


# NR annotation ###
DB="/storage/hdd6/DB/Diamond/NR_07022024"
conda activate checkm2-1.0.2 # diamond 2.0.15

mkdir -p $WDIR/${TAXA}/annotation/NR/tmp
mkdir -p $WDIR/${TAXA}/annotation/NR/NR_annotation_dmnd/
mkdir -p $WDIR/${TAXA}/annotation/NR/NR_annotation_parsed/
mkdir -p $WDIR/${TAXA}/annotation/selfblast/
THREADS=30

ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond_nr.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out $WDIR/${TAXA}/annotation/NR/${line}".NR.annotation" --evalue 1e-10 --sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames salltitles -b 10 --threads $THREADS --tmpdir $WDIR/${TAXA}/annotation/NR/tmp
done

# parse annotation
module load R/4.2.2
ls -d $WDIR/${TAXA}/annotation/NR/*.NR.annotation | sed "s/\.NR.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_NR.R -b $WDIR/${TAXA}/annotation/NR/${line}".NR.annotation" -s $WDIR/annotation/selfblast/${line}".selfblast"  -c 0.4 -o $WDIR/annotation/NR/NR_annotation_parsed_new/${line}".NR.parsed"
done


# Annotation using the Ocean Genome Atlas
DB="/storage/hdd6/DB/Diamond/OM-RGC_v2"

conda activate checkm2-1.0.2

mkdir -p $WDIR/${TAXA}/annotation/OM-RGCv2/tmp
mkdir -p $WDIR/${TAXA}/annotation/OM-RGCv2/OM-RGCv2_annotation_dmnd/
mkdir -p $WDIR/${TAXA}/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/
THREADS=60

ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/OM-RGC_v2_ref.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out $WDIR/${TAXA}/annotation/OM-RGCv2/${line}".OM-RGCv2.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/${TAXA}/annotation/OM-RGCv2/tmp
done

# parse annotation
ls -d $WDIR/${TAXA}/annotation/OM-RGCv2/*.OM-RGCv2.annotation | sed "s/\.OM-RGCv2.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_OM-RGCv2.R -b $WDIR/${TAXA}/annotation/OM-RGCv2/${line}".OM-RGCv2.annotation" -s $WDIR/${TAXA}/annotation/selfblast/${line}".selfblast" -m $DB/OM-RGC_v2.tsv -c 0.4 -t 60 -o $WDIR/${TAXA}/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed"
done


# Prediction of peptidases ###
# new MEROPS was released, therefore run new MEROPS and SignalP test, diamond 2.0.15
conda activate checkm2-1.0.2

mkdir -p $WDIR/${TAXA}/annotation/MEROPS/scan $WDIR/${TAXA}/annotation/MEROPS/pepunit $WDIR/${TAXA}/annotation/MEROPS/tmp $WDIR/${TAXA}/annotation/MEROPS/peptide_blast $WDIR/${TAXA}/annotation/MEROPS/out_merops_parsed

DB="/storage/hdd6/DB/MEROPS/v12.4/diamond"
THREADS=60

# scan DB
# use input from prodigal

ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
 do
  diamond blastp --db $DB/merops_scan.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out $WDIR/${TAXA}/annotation/MEROPS/scan/${line}".scan.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
done


ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa| sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/merops_pepunit.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out  $WDIR/${TAXA}/annotation/MEROPS/pepunit/${line}".pepunit.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
done


# check for false positives by only keeping those hits from pepunit that were also present in scan.db
ls -1 $WDIR/${TAXA}/annotation/MEROPS/pepunit | sed 's/\.pepunit.blastout//' | while read line
do
 cat $WDIR/${TAXA}/annotation/MEROPS/pepunit/${line}".pepunit.blastout" | cut -f1 | grep -F -f - <(cut -f1 $WDIR/${TAXA}/annotation/MEROPS/scan/${line}".scan.blastout") | uniq  > $WDIR/${TAXA}/annotation/MEROPS/tmp/${line}"_hits"
done


ls -1 $WDIR/${TAXA}/annotation/MEROPS/tmp | sed 's/\_hits//' | while read line
do
 cat $WDIR/${TAXA}/annotation/MEROPS/tmp/${line}"_hits" | while read IID
 do
  grep ${IID} $WDIR/${TAXA}/annotation/MEROPS/pepunit/${line}.pepunit.blastout >> $WDIR/${TAXA}/annotation/MEROPS/peptide_blast/${line}".peptide.blastout"
 done
done


# parse files ###
module load R/4.2.2
ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_merops.R -f $WDIR/${TAXA}/annotation/MEROPS/scan/${line}".scan.blastout" -p $WDIR/${TAXA}/annotation/MEROPS/pepunit/${line}".pepunit.blastout" -m /storage/hdd6/DB/MEROPS/v12.4/metadata_pepunit.txt -s $WDIR/${TAXA}/annotation/selfblast/${line}".selfblast" -c 0.4 -o $WDIR/${TAXA}/annotation/MEROPS/out_merops_parsed/${line}".merops.parsed"
done


# Transporter using TCDB database 
DB="/storage/hdd6/DB/TCDB/v122023/diamond"
conda activate checkm2-1.0.2


mkdir -p $WDIR/${TAXA}/annotation/TCDB/tmp
mkdir -p $WDIR/${TAXA}/annotation/TCDB/TCDB_annotation_parsed/

THREADS=100

ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/tcdb.dmnd --query $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa --out $WDIR/${TAXA}/annotation/TCDB/${line}".TCDB.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/${TAXA}/annotation/TCDB/tmp
done

# parse
module load R/4.2.2
ls -d $WDIR/${TAXA}/annotation/TCDB/*.TCDB.annotation | sed "s/\.TCDB.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_TCDB.R -b $WDIR/${TAXA}/annotation/TCDB/${line}".TCDB.annotation" -s $WDIR/${TAXA}/annotation/selfblast/${line}".selfblast" -m /storage/hdd6/DB/TCDB/v122023 -c 0.4 -t 120 -o $WDIR/${TAXA}/annotation/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed"
done


# Prediction of extracellular peptidases
# run signalp on samples
conda activate signalp-6.0g_fast

mkdir -p $WDIR/${TAXA}/annotation/signalp
cd $WDIR/${TAXA}/annotation/signalp
ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//"  | xargs -n 1 basename  | while read line
do
  signalp6 -fasta $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa -format none -org other -od $WDIR/${TAXA}/annotation/signalp/${line}_signalp
done


# dbCAN annotation
DB="/storage/hdd6/DB/dbCAN/v3.0.7"

conda activate dbcan-3.0.7

mkdir -p $WDIR/${TAXA}/annotation/dbCAN
cd $WDIR/${TAXA}/annotation/dbCAN
ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//"  | xargs -n 1 basename |  while read line
do
   #sed '/^\\#\\#FASTA$/,$d'  $WDIR/Annotation/MetaErg/output/${line}/data/all.gff > $WDIR/dbCAN/tmp_${line}.gff
   run_dbcan $WDIR/${TAXA}/annotation/prodigal_out/${line}.faa protein --dbCANFile $DB/dbCAN.txt --out_dir $WDIR/${TAXA}/annotation/dbCAN/${line} --db_dir $DB >> $WDIR/logfiles/dbCAN_${line}.log 2>&1
done
rm $WDIR/${TAXA}/annotation/dbCAN/tmp*



# combine results from annotations for all MAGs

# kegg annotation
ls -1 $WDIR/${TAXA}/annotation/kegg/kegg_annotation_parsed_new/*best.txt | sed 's/\.kegg.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/kegg/kegg_annotation_parsed_new/${line}".kegg.parsed_best.txt" >> $WDIR/${TAXA}/annotation/concat_kegg_annotation_best.txt
done

# NR
ls -1 $WDIR/${TAXA}/annotation/NR/NR_annotation_parsed_new/*best.txt | sed 's/\.NR.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/NR/NR_annotation_parsed_new/${line}".NR.parsed_best.txt" >> $WDIR/${TAXA}/annotation/concat_NR_annotation_best.txt
done

# OM-RGC-v2
ls -1 $WDIR/${TAXA}/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/*best.txt | sed 's/\.OM-RGCv2.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed_best.txt" >> $WDIR/${TAXA}/annotation/concat_OM-RGCv2_annotation_best.txt
done

# merops
ls -1 $WDIR/${TAXA}/annotation/MEROPS/out_merops_parsed/*best.txt | sed 's/\.merops.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/MEROPS/out_merops_parsed/${line}".merops.parsed_best.txt" >> $WDIR/${TAXA}/annotation/concat_merops_annotation_best.txt
done

# signalp
ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/signalp/${line}"_signalp"/prediction_results.txt >> $WDIR/${TAXA}/annotation/concat_signalp_annotation.txt
done

# TCDB
ls -1 $WDIR/${TAXA}/annotation/TCDB/TCDB_annotation_parsed/*best.txt | sed 's/\.TCDB.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed_best.txt" >> $WDIR/${TAXA}/annotation/concat_TCDB_annotation_best.txt
done

# dbCAN
ls -d $WDIR/${TAXA}/annotation/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/${TAXA}/annotation/dbCAN/${line}/overview.txt >> $WDIR/${TAXA}/annotation/concat_dbCAN_annotation.txt
done



# 16S genes
## Run 16S tree with raxml

# start with extracting 16S sequences from new MAGs

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/MAC_metaG"
TAXA="Clostridia"

conda activate barrnap-0.9

THREADS=60

ls -1 $WDIR/${TAXA}/bins_filtered/*fa | sed 's/\.fa//' | while read line
do
 barrnap --threads $THREADS --kingdom bac --outseq ${line}.gff --quiet < ${line}".fa"
done

# this gives me all present rRNA genes including 5S and 23S
# now I need to get only 16S genes

mv $WDIR/${TAXA}/bins_filtered/*.gff $WDIR/${TAXA}/annotation/prodigal_out/
ls -1 $WDIR/${TAXA}/annotation/prodigal_out/*.gff | xargs -n1 basename | sed  's/\.gff//' | while read line
do
 grep -A1 '^>16S' $WDIR/${TAXA}/annotation/prodigal_out/${line}.gff | sed "/^>/s/^>/>${line}_/"   >> $WDIR/${TAXA}/annotation/concat_16S_MAGs.gff
done

# only 5 16S genes
