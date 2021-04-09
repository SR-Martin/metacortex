Usage
==========

If the dataset you are dealing with is small enough and/or you have enough memory and time, you don't need to create intermediate CTX files and can go from FASTQ to contigs in one step::

	echo "file1.fastq" > allreads.txt
	echo "file2.fastq" >> allreads.txt
	echo "file3.fastq" >> allreads.txt
	metacortex_31 -k 31 -n 23 -b 65 -i allreads.txt -t fastq -f contigs.fa -l log.txt

Each time you run MetaCortex, you need to specify a `file of files', which is simply a plain text file that provides a list of the input files. In the above example, the first three lines create this file, then the second line invokes MetaCortex.

To write out sequence graphs for the contigs::

	echo "file1.fastq" > allreads.txt
	echo "file2.fastq" >> allreads.txt
	echo "file3.fastq" >> allreads.txt
	metacortex_31 -k 31 -n 23 -b 65 -i allreads.txt -t fastq -f contigs.fa -l log.txt -A MCC -G

To use the subtractive walk algorithm::

	echo "file1.fastq" > allreads.txt
	echo "file2.fastq" >> allreads.txt
	echo "file3.fastq" >> allreads.txt
	metacortex_31 -k 31 -n 23 -b 65 -i allreads.txt -t fastq -f contigs.fa -l log.txt -A SW

To use th perfect path algorithm::

	echo "file1.fastq" > allreads.txt
	echo "file2.fastq" >> allreads.txt
	echo "file3.fastq" >> allreads.txt
	metacortex_31 -k 31 -n 23 -b 65 -i allreads.txt -t fastq -f contigs.fa -l log.txt -A PP

Creating intermediate CTX files
-------------------------------

If you have many input files and you wish to process them separately, individual FASTQ files can be converted into CTX files (binary representations of the de Bruijn graph) using the following command::

	echo "file1.fastq" > file1.txt
	metacortex_k31 -k 31 -n 23 -b 65 -i file1.txt -t fastq -o file1.ctx

As before, the first line just creates a file of files. The ``-o`` option specifies the name of the binary output file.

Merging CTX files and writing contigs
-------------------------------------

Once you have a set of CTX files, these can be merged and contigs output. A typical command will look like the following::

	echo "file1.ctx" > allfiles.txt
	echo "file2.ctx" >> allfiles.txt
	echo "file3.ctx" >> allfiles.txt
	metacortex_k31 -k 31 -n 23 -b 65 -i allfiles.txt -t binary -o all.ctx -f contigs.fa -g 100 -l log.txt

Again, we start by making a file of files - now containing all the CTX files. This time we specify 'binary' for the ``-t`` option to tell Cortex to expect CTX files. 


Options
-------

Below is a list of options, split into input options, output options, and algorithm parameters.

+----------------------------+--------------------------------------------------------+
| Input Option <argument>    | Description                                            |
+============================+========================================================+
| ``-k <int>``               | The kmer size to be used for the de Bruijn graph.      |
+----------------------------+--------------------------------------------------------+
| ``-n <int>``               | The hash table width.                                  |
+----------------------------+--------------------------------------------------------+
| ``-b <int>``               | The hast table height.                                 |
+----------------------------+--------------------------------------------------------+
| ``--max-db-mem <int>``     | Max memory to use for hash table. Can specify K/M/G/T  |
|                            | e.g. "8G"                                              |
+----------------------------+--------------------------------------------------------+
| ``-i <filename>``          | The name of an input file of files.                    |
+----------------------------+--------------------------------------------------------+
| ``-t <input type>``        | Type of input, either "binary", "fastq" or "fasta".    |
+----------------------------+--------------------------------------------------------+     


+----------------------------+--------------------------------------------------------+
| Output Option <argument>   | Description                                            |
+============================+========================================================+
| ``-o <ctx file>``          | Name of file to write hash table to.                   |
+----------------------------+--------------------------------------------------------+
| ``-f <fasta file>``        | Name of file to write contigs to.                      |
+----------------------------+--------------------------------------------------------+
| ``-l <log file>``          | Name of log file.                                      |
+----------------------------+--------------------------------------------------------+
| ``-G``                     | Write sequence graph in fastg and GFA format. Must be  |
|                            | used with MCC algorithm. Filename is taken from ``-f`` |
+----------------------------+--------------------------------------------------------+


+----------------------------+--------------------------------------------------------+
| Parameter <argument>       | Description                                            |
+============================+========================================================+
| ``-A <algorithm>``         | The graph traversal algorithm to use, must be one of   |
|                            | "MCC"(default), "SW", "PP" or "GS".                    |
+----------------------------+--------------------------------------------------------+
| ``-g <min contig length>`` | Minimum contig length to output. Default 1.            |
+----------------------------+--------------------------------------------------------+
| ``-C <min path coverage>`` | Minimum value for coverage along paths. Default 2.     |
+----------------------------+--------------------------------------------------------+
| ``-W <SW_delta>``          | Value to use as delta in algorithm SW. Must be between |
|                            | 0 and 1. Default 0.8.                                  |
+----------------------------+--------------------------------------------------------+ 