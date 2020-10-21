Introduction
============

MetaCortex is an assembler for metagenomic, or environmental sequence data. It is based onnthe consensus version of the Cortex assembler (cortex_con) developed by Mario Caccamo and Zamin Iqbal.

Key Concepts
------------

Though it is possible to go straight from FASTA or FASTQ files to a complete meta-assembly, MetaCortex (like Cortex) implements an intermediate binary file format which enables the parallelisation of the process of converting raw reads into a de Bruijn graph structure. As well as speeding up the overal assembly process by sharing out the reads across multiple processing cores, the approach also makes it possible to carru out big assemblies in lower memory environments.

Figure 1 illustrates the typical approach taken with MetaCortex, Reads are separated into files -0 either by hand, or beacuse the sequencing instrument has produced multiple files. MetaCortex is then executed in two stages: Firstly, in parallel, to import each reads file, create a de Bruijn graph and output a CTX file. Then secondly, to combine the separate CTX files and output contigs and sequence graphs.


