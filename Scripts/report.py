from snakemake.utils import report

# with open(snakemake.input.T1) as vcf:
#     n_calls = sum(1 for l in vcf if not l.startswith("#"))

report("""
TBI Drop-seq Workflow Documentation
===================================

Data was processed using Drop-seq Tools (version) Star (version), etc.

Workflow (see Workflow_).
Cellular Barcodes (see Table T1_).
Molecular Barcodes (see Table T2_).
""", snakemake.output[0],metadata="Author: Douglas Arneson (darneson@ucla.edu)", **snakemake.input)


