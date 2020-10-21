Usage
==========

If the dataset you are dealing with is small enough and/or you have enough memory and time, you don't need to create intermediate CTX files and can go from FASTQ to contigs in one step: ::
	echo "file1.fastq" > allfiles.txt
	echo "file2.fastq" >> allfiles.txt
	echo "file3.fastq" >> allfiles.txt
	cortex_con_31 -k 31 -n 23 -b 65 -i allfiles.txt -t fastq -o all.ctx -f contigs.fa -g 100 -l log.txt

Each time you run MetaCortex, you need to specify a `file of files', which is simply a plain text file that provides a list of the input files. In the above example, the first three lines create this file, then the second line invokes MetaCortex. The options have the following meanings:

+-------------------------+--------------------------------------------------------+
| Option <argument>       | Description                                            |
+=========================+========================================================+
| ``-k <int>``            | The kmer size to be used for the de Bruijn graph.      |
+-------------------------+--------------------------------------------------------+
| ``-n <int>``            |  The hash table width.                                 |
+-------------------------+--------------------------------------------------------+
| ``-b <int>``            |  The hast table height.                                |
+-------------------------+--------------------------------------------------------+
| ``-i <filename>``       |  The name of an input file of files.                   |
+-------------------------+--------------------------------------------------------+
| ``-t <input type>``     | Type of input, either binary, fastq or fasta.          |
+-------------------------+--------------------------------------------------------+

Creating intermediate CTX files
-------------------------------

If you have many input files and you wish to process them separately, individual FASTQ files can be converted into CTX files (binary representations of the de Bruijn graph) using the following command: ::
	echo "file1.fastq" > file1.txt
	metacortex_k31 -k 31 -n 23 -b 65 -i file1.txt -t fastq -o file1.ctx

As before, the first line just creates a file of files. The new ``-o`` option specifies the name of the binary output file.

Merging CTX files and writing contigs
-------------------------------------

Once you have a set of CTX files, these can be merged and contigs output. A typical command will look like the following: ::
	echo "file1.ctx" > allfiles.txt
	echo "file2.ctx" >> allfiles.txt
	echo "file3.ctx" >> allfiles.txt
	metacortex_k31 -k 31 -n 23 -b 65 -i allfiles.txt -t binary -o all.ctx -f contigs.fa -g 100

Again, we start by making a file of files - now containing all the CTX files. This time we specify `binary' for the ``-t`` option to tell Cortex to expect CTX files. 