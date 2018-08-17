Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It allows to go from raw data of your Single Cell RNA seq experiment until the final count matrix with QC plots along the way.

This is the tool we use in our lab to improve our wetlab protocol as well as provide an easy framework to reproduce and compare different experiments with different parameters.

It uses STAR to map the reads. It is usable for any single cell protocol using two reads where the first one holds the Cell and UMI barcodes and the second read holds the RNA. Here is a non-exhausitve list of compatible protocols:

* Drop-Seq
* SCRB-Seq
* 10x Genomics
* DroNc-seq

This package is trying to be as user friendly as possible. One of the hopes is that non-bioinformatician can make use of it without too much hassle. It will still require some command line execution, this is not going to be fully interactive package.
