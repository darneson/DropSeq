#to submit:
#  snakemake dirs
#  snakemake -j 100 --cluster "qsub {params.sge_opts}" --forceall
#  have to make log directories to submit --> otherwise wont work

#generate flowchart:
# snakemake --dag Output/{701,702,703,705,706}/gene_exon_tagged.dge.txt.gz | dot -Tsvg > flowchart.svg

#to get all files in a directory
# IDS, = glob_wildcards("thedir/{id}.fastq")

#get samples
configfile: "config.yaml"
ALL_SAMPLES = config["bam_files"]

#Make log directories
DIRS = ['logs/','benchmarks/','logs/tag_cells','logs/tag_molecules','logs/filter_bam','logs/trim_starting_sequence','logs/trim_poly_a','logs/sam_to_fastq',
'logs/star','logs/sort_aligned','logs/merge_bam','logs/gene_exon_tagged','logs/synth_err','logs/tag_hist','logs/kneeplot','logs/dig_expr','logs/report']

for dir in DIRS:
    if not os.path.exists(dir):
        os.makedirs(dir)

##### email/job title #####
email = "darneson@ucla.edu"
job_label = "TBI"

##### Important Parameters #####
Number_Barcodes = "2000"   #Double number of expected cells for detect bead synthesis errors
Number_Core_Barcodes = "1000"

##### reference files #####
STARREFDIR = '/u/home/d/darneson/nobackup-xyang123/DropSeq/Star_Index/mm10'
REF_SEQ = '/u/home/d/darneson/shared-project/datasets/Drop_Seq/mm10_reference/mm10.fasta'
REF_FLAT = '/u/home/d/darneson/shared-project/datasets/Drop_Seq/mm10_reference/mm10.refFlat'

##### TOOLS #####
DropSeqTools = "/u/home/d/darneson/shared-project/tools/Drop_Seq/Drop-seq_tools-1.12/"
TagBamWithReadSequenceExtended = DropSeqTools+"TagBamWithReadSequenceExtended"
FilterBAM = DropSeqTools+"FilterBAM"
TrimStartingSequence = DropSeqTools+"TrimStartingSequence"
PolyATrimmer = DropSeqTools+"PolyATrimmer"
TagReadWithGeneExon = DropSeqTools+"TagReadWithGeneExon"
DetectBeadSynthesisErrors = DropSeqTools+"DetectBeadSynthesisErrors"
BAMTagHistogram = DropSeqTools+"BAMTagHistogram"
DigitalExpression = DropSeqTools+"DigitalExpression"
picard_jar=DropSeqTools+'/3rdParty/picard/picard.jar'
STAR='/u/home/d/darneson/shared-project/tools/Drop_Seq/STAR-2.5.0c/bin/Linux_x86_64/STAR'

rule all:
    input:
        # expand("Output/{sample}/aligned.sorted.bam", sample = config["bam_files"])
        "report.html"

# Stage 1: pre-alignment tag and trim
rule tag_cells:
    input:
        bam="/u/home/d/darneson/nobackup-xyang123/DropSeq/TBI_Project/FastqToSam/TBI_{sample}_fastq_to_bam.bam",
    output:
        summary="Output/{sample}/unaligned_tagged_Cellular.bam_summary.txt",
        tagged_cell="Output/{sample}/unaligned_tagged_Cell.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/tag_cells -e logs/tag_cells"
    benchmark:
        "benchmarks/{sample}.tag_cells.benchmark.txt"
    shell:
        """
        {TagBamWithReadSequenceExtended} BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 \
        INPUT={input.bam} SUMMARY={output.summary} OUTPUT={output.tagged_cell}
        """
rule tag_molecules:
    input:
        tagged_cell="Output/{sample}/unaligned_tagged_Cell.bam"
    output:
        summary="Output/{sample}/unaligned_tagged_Molecular.bam_summary.txt",
        tagged_molecule="Output/{sample}/unaligned_tagged_CellMolecular.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:35:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/tag_molecules -e logs/tag_molecules"
    benchmark:
        "benchmarks/{sample}.tag_molecules.benchmark.txt"
    shell:
        """
        {TagBamWithReadSequenceExtended} BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 \
        INPUT={input.tagged_cell} SUMMARY={output.summary} OUTPUT={output.tagged_molecule} 
        """
rule filter_bam:
    input:
        tagged_molecule="Output/{sample}/unaligned_tagged_CellMolecular.bam"
    output:
        tagged_filtered="Output/{sample}/unaligned_tagged_filtered.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/filter_bam -e logs/filter_bam"
    benchmark:
        "benchmarks/{sample}.filter_bam.benchmark.txt"
    shell:
        """
        {FilterBAM} TAG_REJECT=XQ \
        INPUT={input.tagged_molecule} OUTPUT={output.tagged_filtered}
        """
rule trim_starting_sequence: 
    input:
        tagged_filtered="Output/{sample}/unaligned_tagged_filtered.bam"
    output:
        summary="Output/{sample}/adapter_trimming_report.txt",
        trimmed_smart="Output/{sample}/unaligned_tagged_trimmed_smart.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/trim_starting_sequence -e logs/trim_starting_sequence"
    benchmark:
        "benchmarks/{sample}.traim_starting_sequence.benchmark.txt"
    shell:
        """
        {TrimStartingSequence} SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 \
        INPUT={input.tagged_filtered} OUTPUT_SUMMARY={output.summary} OUTPUT={output.trimmed_smart}
        """
rule trim_poly_a:
    input:
        trimmed_smart="Output/{sample}/unaligned_tagged_trimmed_smart.bam"
    output:
        summary="Output/{sample}/polyA_trimming_report.txt",
        unmapped_bam="Output/{sample}/unaligned_mc_tagged_polyA_filtered.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/trim_poly_a -e logs/trim_poly_a"
    benchmark:
        "benchmarks/{sample}.trim_poly_a.benchmark.txt"
    shell:
        """
        {PolyATrimmer} MISMATCHES=0 NUM_BASES=6 \
        INPUT={input.trimmed_smart} OUTPUT_SUMMARY={output.summary} OUTPUT={output.unmapped_bam}
        """
# Stage 2: alignment
rule sam_to_fastq:
    input:
        unmapped_bam="Output/{sample}/unaligned_mc_tagged_polyA_filtered.bam"
    output:
        fastq="Output/{sample}/unaligned_mc_tagged_polyA_filtered.fastq"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/sam_to_fastq -e logs/sam_to_fastq"
    benchmark:
        "benchmarks/{sample}.sam_to_fastq.benchmark.txt"
    shell:
        """
        java -Xmx500m -jar {picard_jar} SamToFastq \
        INPUT={input.unmapped_bam} FASTQ={output.fastq}
        """
rule star_alignment:
    input:
        fastq="Output/{sample}/unaligned_mc_tagged_polyA_filtered.fastq",
        starref=STARREFDIR
    output:
        prefix="Output/{sample}/star.",
        aligned_sam="Output/{sample}/star.Aligned.out.sam"
    params: sge_opts="-cwd -j y -pe shared 10 -l h_data=8G,highmem,h_rt=1:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/star -e logs/star"
    benchmark:
        "benchmarks/{sample}.star_alignment.benchmark.txt"
    shell:
        """
        {STAR} --genomeDir {input.starref} --runThreadN 10 \
        --readFilesIn {input.fastq} --outFileNamePrefix {output.prefix}
        """
#Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
rule sort_aligned:
    input:
        aligned_sam="Output/{sample}/star.Aligned.out.sam"
    output:
        aligned_sorted_bam="Output/{sample}/aligned.sorted.bam",
        temp_dir="Output/{sample}/"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/sort_aligned -e logs/sort_aligned"
    benchmark:
        "benchmarks/{sample}.sort_aligned.benchmark.txt"
    shell:
        """
        java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar {picard_jar} SortSam \
        INPUT={input.aligned_sam} OUTPUT={output.aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR={output.temp_dir}
        """
#Stage 4: merge and tag aligned reads
rule merge_bam:
    input:
        aligned_sorted_bam="Output/{sample}/aligned.sorted.bam",
        unmapped_bam="Output/{sample}/unaligned_mc_tagged_polyA_filtered.bam",
        reference_seq=REF_SEQ
    output:
        merged="Output/{sample}/merged.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/merge_bam -e logs/merge_bam"
    benchmark:
        "benchmarks/{sample}.merge_bam.benchmark.txt"
    shell:
        """
        java -Xmx4000m -jar {picard_jar} MergeBamAlignment INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false \
        REFERENCE_SEQUENCE={input.reference_seq} UNMAPPED_BAM={input.unmapped_bam} ALIGNED_BAM={input.aligned_sorted_bam} OUTPUT={output.merged}
        """
rule tag_with_gene_exon:
    input:
        merged="Output/{sample}/merged.bam",
        refflat=REF_FLAT
    output:
        gene_exon_tagged="Output/{sample}/star_gene_exon_tagged.bam"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/gene_exon_tagged -e logs/gene_exon_tagged"
    benchmark:
        "benchmarks/{sample}.tag_exon_with_gene.benchmark.txt"
    shell:
        """
        {TagReadWithGeneExon} TAG=GE CREATE_INDEX=true \
        INPUT={input.merged} ANNOTATIONS_FILE={input.refflat} O={output.gene_exon_tagged}
        """
rule synth_err:
    input:
        gene_exon_tagged="Output/{sample}/star_gene_exon_tagged.bam"
    output:
        gene_exon_tagged_clean="Output/{sample}/gene_exon_tagged_clean.bam",
        stats="Output/{sample}/my.synthesis_stats.txt",
        summary="Output/{sample}/my.synthesis_stats.summary.txt"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/synth_err -e logs/synth_err"
    benchmark:
        "benchmarks/{sample}.synth_err.benchmark.txt"
    shell:
        """
        {DetectBeadSynthesisErrors} NUM_BARCODES={Number_Barcodes} PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
        I={input.gene_exon_tagged} O={output.gene_exon_tagged_clean} OUTPUT_STATS={output.stats} SUMMARY={output.summary}
        """
rule tag_hist:
    input:
        gene_exon_tagged_clean="Output/{sample}/gene_exon_tagged_clean.bam"
    output:
        cell_readcounts="Output/{sample}/cell_readcounts.txt.gz"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/tag_hist -e logs/tag_hist"
    benchmark:
        "benchmarks/{sample}.tag_hist.benchmark.txt"
    shell:
        """
        {BAMTagHistogram} TAG=XC \
        I={input.gene_exon_tagged_clean} O={output.cell_readcounts}
        """
rule kneeplot:
    input:
        cell_readcounts="Output/{sample}/cell_readcounts.txt.gz"
    output: 
        pdf="Output/{sample}/Kneeplot.pdf"
    params: sge_opts="-cwd -j y -l h_data=2G,h_rt=0:10:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/kneeplot -e logs/kneeplot"
    benchmark:
        "benchmarks/{sample}.knee_plot.benchmark.txt"
    shell:
        """
        module load R;
        Rscript scripts/KneePlot.R {input.cell_readcounts} {output.pdf}
        """
rule dig_expr:
    input:
        pdf="Output/{sample}/Kneeplot.pdf",
        gene_exon_tagged_clean="Output/{sample}/gene_exon_tagged_clean.bam"
    output:
        dge="Output/{sample}/gene_exon_tagged.dge.txt.gz",
        summary="Output/{sample}/gene_exon_tagged.dge.summary.txt"
    params: sge_opts="-cwd -j y -l h_data=8G,h_rt=0:45:00 -M "+email+" -m bea -N "+job_label+"{sample} -o logs/dig_expr -e logs/dig_expr"
    benchmark:
        "benchmarks/{sample}.dig_expr.benchmark.txt"
    shell:
        """
        {DigitalExpression} NUM_CORE_BARCODES={Number_Core_Barcodes} \
        INPUT={input.gene_exon_tagged_clean} OUTPUT={output.dge} SUMMARY={output.summary}
        """
#can use this to make more sophisticated reports
#http://pythonhosted.org/sphinxleash/
# https://bitbucket.org/snakemake/snakemake/issues/86/embed-images-in-report
rule report:
    input:
        Workflow="flowchart.svg",
        T1=expand("Output/{sample}/unaligned_tagged_Cellular.bam_summary.txt", sample=config["bam_files"]),
        T2=expand("Output/{sample}/unaligned_tagged_Molecular.bam_summary.txt", sample=config["bam_files"])
    output:
        "report.html"
    params: sge_opts="-cwd -j y -l h_data=2G,h_rt=0:10:00 -M "+email+" -m bea -N "+job_label+"_report -o logs/report -e logs/report"
    script:
        "scripts/report.py"