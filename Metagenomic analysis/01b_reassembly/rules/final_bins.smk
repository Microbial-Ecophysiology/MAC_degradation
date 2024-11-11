rule fastani_check:
	input:
		winning = config["adir"] + "/intermediate_results/Bin_dereplication/drep_out/data_tables/Wdb.csv",
		cluster = config["adir"] + "/intermediate_results/Bin_dereplication/drep_out/data_tables/Cdb.csv"
	output:
		check_clusters = config["adir"] + "/intermediate_results/Bin_dereplication/fastani_check/check_clusters.txt",
		fastani_check = config["adir"] + "/intermediate_results/Bin_dereplication/fastani_check/check_summary.txt"
	params:
		bindir = config["adir"] + "/intermediate_results/Bin_dereplication/input_bins",
		anidir = config["adir"] + "/intermediate_results/Bin_dereplication/fastani_check"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/drep.yaml"
	shell:
		"""
		# extract winning genomes from initial clustering from cluster DB
		cut -d',' -f1 {input.winning} | sed '1d' | grep -F -f - {input.cluster} | tr ',' '\\t' | sort -k1,1 | cut -f1,2,6 > {output.check_clusters}
		# only look at comparisons within the same primary cluster
		cut -f3 {output.check_clusters} | sort | uniq -d | while read line
		do
		  awk -v FS="\\t" -v OFS="\\t" -v cl="$line" '$3 == cl' {output.check_clusters} > {params.anidir}/check_cl_${{line}}.txt
		  cut -f1 {params.anidir}/check_cl_${{line}}.txt | sed "s!^!{params.bindir}/!" > {params.anidir}/tmp
		  fastANI --ql {params.anidir}/tmp --rl {params.anidir}/tmp -o {params.anidir}/fastani_cl_${{line}}.txt --fragLen 1000
		  rm {params.anidir}/tmp
		  sed "s/^/${{line}}\\t/" {params.anidir}/fastani_cl_${{line}}.txt
		done > {output.fastani_check}
		"""

checkpoint adjust_dereplication:
	input:
		co_quant = config["adir"] + "/CoAssembly/Refined_Bins/Co_bins_quant.txt",
		sample_quant_spades = expand(config["adir"] + "/{sample}/intermediate_results/Refined_Bins/spades/{sample}_bins_quant.txt", sample = SAMPLES),
		sample_quant_megahit = expand(config["adir"] + "/{sample}/intermediate_results/Refined_Bins/megahit/{sample}_bins_quant.txt", sample = SAMPLES),
		genomeInfo = config["adir"] + "/intermediate_results/Bin_dereplication/drep_out/data_tables/genomeInfo.csv",
		winning = config["adir"] + "/intermediate_results/Bin_dereplication/drep_out/data_tables/Wdb.csv",
		cluster = config["adir"] + "/intermediate_results/Bin_dereplication/drep_out/data_tables/Cdb.csv",
		fastani_check = config["adir"] + "/intermediate_results/Bin_dereplication/fastani_check/check_summary.txt"
	output:
		winning_adj = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/Wdb_adjusted.csv",
		cluster_adj = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/Cdb_adjusted.csv",
		workspace = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/drep_adjustment.Rdata",
		done = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/done",
		gendir = directory(config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes")
	params:
		script = config["wdir"] + "/scripts/adjust_drep.R",
		bindir = config["adir"] + "/intermediate_results/Bin_dereplication/input_bins"
	conda: config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -c {input.cluster} -w {input.winning} -g {input.genomeInfo} -f {input.fastani_check} -q "{input.co_quant} {input.sample_quant_spades} {input.sample_quant_megahit}" -a 95 -A 0.5 -d 0.8 -C {output.cluster_adj} -W {output.winning_adj} -r {output.workspace}
		mkdir -p {output.gendir}
		cut -d',' -f1 {output.winning_adj} | sed '1d' | while read line
		do
		  cp {params.bindir}/$line {output.gendir}/$line
		done
		touch {output.done}
		"""

rule cat_clean_reads:
	input:
		clean_R1 = expand(config["adir"] + "/{sample}/intermediate_results/clean_reads/final_pure_reads_dedupe_1.fq", sample = SAMPLES),
		clean_R2 = expand(config["adir"] + "/{sample}/intermediate_results/clean_reads/final_pure_reads_dedupe_2.fq", sample = SAMPLES)
	output:
		all_R1 = config["adir"] + "/intermediate_results/clean_reads/all_R1.fastq",
		all_R2 = config["adir"] + "/intermediate_results/clean_reads/all_R2.fastq"
	shell:
		"""
		cat {input.clean_R1} > {output.all_R1}
		cat {input.clean_R2} > {output.all_R2}
		"""

rule bin_reassembly_index:
	input:
		binlist = config["adir"] + "/intermediate_results/Bin_dereplication/bin_list.txt",
	output:
		allbins = config["adir"] + "/intermediate_results/Bin_reassembly/bin_index/all_bins.fa",
		binindex = config["adir"] + "/intermediate_results/Bin_reassembly/bin_index/all_bins.fa.sa"
	params:
		bindir = config["adir"] + "/intermediate_results/Bin_reassembly/original_bins"
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.bindir}
		cat {input.binlist} | while read line
		do
		  BIN=$(echo $line | xargs basename | sed 's/\\.fa//')
		  sed "/^>/s/^>/>${{BIN}}_/" $line > {params.bindir}/${{BIN}}.fa
		done 
		cat {params.bindir}/*.fa > {output.allbins}
		bwa index {output.allbins}
		"""

# remove supplementary and secondary alignments from bam file based on sam flag
rule bin_reassembly_mapping:
	input:
		all_R1 = config["adir"] + "/intermediate_results/clean_reads/all_R1.fastq",
		all_R2 = config["adir"] + "/intermediate_results/clean_reads/all_R2.fastq",
		allbins = config["adir"] + "/intermediate_results/Bin_reassembly/bin_index/all_bins.fa",
		binindex = config["adir"] + "/intermediate_results/Bin_reassembly/bin_index/all_bins.fa.sa"
	output:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads/done"
	params:
		script_dir = config["metawrap_scripts"],
		indir = config["adir"] + "/intermediate_results/Bin_reassembly/original_bins",
		outdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads"
	#threads: config["parallel_threads"]
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.outdir}
		ulimit -n 10000
		# future development:  maybe allow redundant mappings here... (at the moment, this is disabled with the samtools flag filter)
		bwa mem -t {threads} {input.allbins} {input.all_R1} {input.all_R2} | samtools view -F 256 - | {params.script_dir}/filter_reads_for_bin_reassembly.py {params.indir} {params.outdir} 2 5
		touch {output.done}
		"""

# include rmdup (seqkit) to remove duplicates in final output fastq
rule bin_reassembly_collect:
	input:
		winning = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/Wdb_adjusted.csv",
		cluster = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/Cdb_adjusted.csv",
		done = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads/done"
	output:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done"
	params:
		readdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads",
		collectdir = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly"
	conda: config["wdir"] + "/envs/kraken2.yaml"
	shell:
		"""
		mkdir -p {params.collectdir}
		sed '1d' {input.winning} | while read line
		do
		  REP=$(echo $line | cut -d',' -f1 | sed 's/\\.fa//')
		  CLU=$(echo $line | cut -d',' -f2)
		  rm -f {params.collectdir}/${{REP}}.strict_1.fastq {params.collectdir}/${{REP}}.strict_2.fastq {params.collectdir}/${{REP}}.permissive_1.fastq {params.collectdir}/${{REP}}.permissive_2.fastq
		  grep -w "$CLU" {input.cluster} | cut -d',' -f1 | sed 's/\.fa//' | while read fa
		  do
		    cat {params.readdir}/${{fa}}.strict_1.fastq | seqkit rmdup -n >> {params.collectdir}/${{REP}}.strict_1.fastq
		    cat {params.readdir}/${{fa}}.strict_2.fastq | seqkit rmdup -n >> {params.collectdir}/${{REP}}.strict_2.fastq
		    cat {params.readdir}/${{fa}}.permissive_1.fastq | seqkit rmdup -n >> {params.collectdir}/${{REP}}.permissive_1.fastq
		    cat {params.readdir}/${{fa}}.permissive_2.fastq | seqkit rmdup -n >> {params.collectdir}/${{REP}}.permissive_2.fastq
		  done
		done
		touch {output.done}
		"""

rule bin_reassembly_strict:
	input:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done",
		fa = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes/{cbin}.fa"
	output:
		fa = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.strict.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly",
		tmpdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_spades",
		spadesdir = config["adir"] + "/intermediate_results/Bin_reassembly/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.strict --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.strict_1.fastq -2 {params.readdir}/{params.bin}.strict_2.fastq -o {params.spadesdir}/{params.bin}.strict || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.strict/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.strict/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.strict/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.strict/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_strict_individual:
	input:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done",
		fa = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes/{cbin}.fa"
	output:
		fa = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.indstr.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads",
		tmpdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_spades",
		spadesdir = config["adir"] + "/intermediate_results/Bin_reassembly/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.indstr --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.strict_1.fastq -2 {params.readdir}/{params.bin}.strict_2.fastq -o {params.spadesdir}/{params.bin}.indstr || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.indstr/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.indstr/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.indstr/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.ind_str/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_permissive:
	input:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done",
		fa = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes/{cbin}.fa"
	output:
		fa = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.permissive.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly",
		tmpdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_spades",
		spadesdir = config["adir"] + "/intermediate_results/Bin_reassembly/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.permissive --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.permissive_1.fastq -2 {params.readdir}/{params.bin}.permissive_2.fastq -o {params.spadesdir}/{params.bin}.permissive || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.permissive/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.permissive/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_permissive_individual:
	input:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done",
		fa = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes/{cbin}.fa"
	output:
		fa = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.indper.fa"
	params:
		script_dir = config["metawrap_scripts"],
		readdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_reads",
		tmpdir = config["adir"] + "/intermediate_results/Bin_reassembly/tmp_spades",
		spadesdir = config["adir"] + "/intermediate_results/Bin_reassembly/reassemblies",
		bin = "{cbin}"
	threads: config["annotation_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.spadesdir}
		# if spades fails, just copy the original assembly (this should avoid breaking the snakemake workflow)
		spades.py -t {threads} -m 150 --tmp {params.tmpdir}/tmp_{params.bin}.indper --careful --untrusted-contigs {input.fa} -1 {params.readdir}/{params.bin}.permissive_1.fastq -2 {params.readdir}/{params.bin}.permissive_2.fastq -o {params.spadesdir}/{params.bin}.indper || sed '/^>/s/^[^_]*_/>/' {input.fa} > {params.spadesdir}/{params.bin}.indper/scaffolds.fasta
		if [[ -s {params.spadesdir}/{params.bin}.indper/scaffolds.fasta ]]
		then
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.indper/scaffolds.fasta > {output.fa}
		else
		  {params.script_dir}/rm_short_contigs.py 1000 {params.spadesdir}/{params.bin}.indper/contigs.fasta > {output.fa}
		fi
		"""

rule bin_reassembly_orig:
	input:
		done = config["adir"] + "/intermediate_results/Bin_reassembly/reads_for_reassembly/done",
		fa = config["adir"] + "/intermediate_results/Bin_dereplication/drep_adj/dereplicated_genomes/{cbin}.fa"
	output:
		fa = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.orig.fa"
	shell:
		"""
		cp {input.fa} {output.fa}
		"""

# I am using the wildcard name cbin here (candidate bin), to not cause any conflicts with the wildcard bin used later
def aggregate_strict(wildcards):
	checkpoint_output = checkpoints.adjust_dereplication.get().output['gendir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	cbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.strict.fa", cbin = cbin)

def aggregate_permissive(wildcards):
	checkpoint_output = checkpoints.adjust_dereplication.get().output['gendir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	cbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.permissive.fa", cbin = cbin)

def aggregate_indstr(wildcards):
	checkpoint_output = checkpoints.adjust_dereplication.get().output['gendir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	cbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.indstr.fa", cbin = cbin)

def aggregate_indper(wildcards):
	checkpoint_output = checkpoints.adjust_dereplication.get().output['gendir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	cbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.indper.fa", cbin = cbin)

def aggregate_orig(wildcards):
	checkpoint_output = checkpoints.adjust_dereplication.get().output['gendir']
	fasta_files = sorted(glob.glob(os.path.join(checkpoint_output, '*.fa')))
	cbin = [re.sub('\.fa$', '', os.path.basename(x)) for x in fasta_files]
	return expand(config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins/{cbin}.orig.fa", cbin = cbin)

rule bin_reassembly_checkm:
	input:
		strict = aggregate_strict,
		permissive = aggregate_permissive,
		indstr = aggregate_indstr,
		indper = aggregate_indper,
		orig = aggregate_orig
	output:
		stats = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins.stats"
	params:
		script_dir = config["metawrap_scripts"],
		tmpdir = config["adir"] + "/tmp",
		checkmdir = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins.checkm",
		indir = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/metawrap.yaml"
	shell:
		"""
		mkdir -p {params.tmpdir}
		checkm lineage_wf -x fa {params.indir} {params.checkmdir} -t {threads} --tmpdir {params.tmpdir} --pplacer_threads 10
		{params.script_dir}/summarize_checkm.py {params.checkmdir}/storage/bin_stats_ext.tsv | (read -r; printf "%s\\n" "$REPLY"; sort) > {output.stats}
		"""

# settings for parameters a and b for calculating score to select best bin have been empirically chosen
# they should not be changed independently from minimum completeness and maximum contamination settings
rule bin_reassembly_best:
	input:
		stats = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins.stats"
	output:
		stats = config["adir"] + "/intermediate_results/Bin_reassembly/best_bins.stats"
	params:
		script = config["wdir"] + "/scripts/select_best_bins.R",
		comp = config["completeness"],
		cont = config["contamination"]
	conda: config["wdir"] + "/envs/r.yaml"
	shell:
		"""
		{params.script} -s {input.stats} -o {output.stats} -a 5 -b 1.2 -c {params.comp} -r {params.cont}
		"""

checkpoint collect_final_bins:
	input:
		checkm = config["adir"] + "/intermediate_results/Bin_reassembly/best_bins.stats"
	output:
		done = config["adir"] + "/intermediate_results/Bin_final/done",
		outdir = directory(config["adir"] + "/intermediate_results/Bin_final/bins")
	params:
		indir = config["adir"] + "/intermediate_results/Bin_reassembly/reassembled_bins"
	shell:
		"""
		mkdir -p {output.outdir}
		cut -f1 {input.checkm} | sed '1d' | while read line
		do 
		  sed -e "/^>/s/^>/>${{line}}_/" -e '/^>/s/[-\\.]/_/g' {params.indir}/${{line}}.fa > {output.outdir}/${{line}}.fa
		done
		touch {output.done}
		"""


rule run_checkm2:
	input:
		done = config["adir"] + "/intermediate_results/Bin_final/done"
	output:
		checkm2 = config["adir"] + "/intermediate_results/Bin_final/checkm2_out/quality_report.tsv"
	params:
		bindir = config["adir"] + "/intermediate_results/Bin_final/bins",
		outdir = config["adir"] + "/intermediate_results/Bin_final/checkm2_out",
		tmpdir = config["adir"] + "/tmp"
	threads: config["parallel_threads"]
	conda: config["wdir"] + "/envs/checkm2.yaml"
	shell:
		"""
		mkdir -p {params.tmpdir}
		checkm2 predict --threads {threads} -x fa --allmodels --tmpdir {params.tmpdir} --input {params.bindir} --output-directory {params.outdir}
		"""
