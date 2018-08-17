Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It starts from raw sequencing data and outputs digital gene expression matrices (DGEs).

[WE ARE CURRENTLY USING AN UPDATED VERSION HERE](https://github.com/darneson/dropSeqPipeDropEST)

Updates
------------------
We have since moved on from using this method and now use an adapted version of [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) which we have found to be much more sophisticated and flexible than our [snakemake](https://snakemake.readthedocs.io/en/stable/) implementation of [Drop-seq Tools](http://mccarrolllab.com/download/1276/). 

We have made a number of changes to the original workflow:

1). We have adapted the [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) workflow to work on a [SGE Cluster Environment](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html), specifically the [UCLA Hoffman2 Cluster](https://www.hoffman2.idre.ucla.edu/computing/sge/).

2). We have made the workflow compatible with the new pipeline [dropEST](https://github.com/hms-dbmi/dropEst) which is capable of handling reads aligning to intronic regions and results in a large boost in cell number and genes/UMIs per cell in our hands.

If you would like to use this updated [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow, please go [here](https://github.com/darneson/dropSeqPipeDropEST).

Requirements
------------------
[snakemake](https://snakemake.readthedocs.io/en/stable/) is the workflow management which wraps [Drop-seq Tools](http://mccarrolllab.com/dropseq/) which runs everything

[Drop-seq Tools](http://mccarrolllab.com/dropseq/) can be downloaded from [here](http://mccarrolllab.com/download/1276/) and a manual of the various functions is available [here](http://mccarrolllab.com/wp-content/uploads/2016/03/Drop-seqAlignmentCookbookv1.2Jan2016.pdf)

[STAR](https://github.com/alexdobin/STAR) is used to map the reads.

Changes to config file
------------------
In order to run the workflow, you must modify the [config.yaml](https://github.com/darneson/DropSeq/blob/master/config.yaml) file with your sample information. Simply provide the "SampleName" for each sample for read1 of your paired end fastq files.

Changes to Snakefile
------------------
The following changes must be made to the [Snakefile](https://github.com/darneson/DropSeq/blob/master/Snakefile)

1). If you are running in a [SGE Cluster Environment](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) and you would like email notifications of your job status, add your email to the **email** parameter.

2). If you would like to give your cluster jobs a specific prefix, you can specify that with **job_label**

3). Provide a path to your STAR Index for your reference genome with: **STARREFDIR**

4). Give the paths to your reference genome fasta file: **REF_SEQ** and gtf file: **GTF_FILE**

5). Provide the path to the Drop-seq Tools directory: **DropSeqTools**

6). Provide the path to your STAR executable: **STAR**

7). If necessary, provide the path to your particular java version executable: **JAVA_distro**

Running the Workflow
------------------
Run the snakemake workflow with the following command:
```
snakemake -j 100 --cluster "qsub {params.sge_opts}"
```
