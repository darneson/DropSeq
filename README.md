Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It starts from raw sequencing data and outputs digital gene expression matrices (DGEs).

[CURRENT VERSION HERE](https://github.com/darneson/dropSeqPipeDropEST)

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

