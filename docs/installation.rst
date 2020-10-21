Installation
==============

A decision needs to be made at compile-time on the largest k-mer size which you wish a given executable to support and this value may be either 31, 63, 95, 127, 160 or 192 nucleotides. Selecting this at compile-time allows for much more efficient use of memory during program execution. Often users will compile a number of different versions of the code, selecting at run-time the most appropriate one to use. The maximum kmer size may be specified as a parameter to the make command: ::
	make MAXK=31 MAC=1 metacortex

This command should be executed in the directory containing the Makefile. When the build process completes, the executable may be found in the ``bin`` directory as ``metacortex_kX`` where ``X`` is the maximum kmer size chosen. If building on Mac OS, it is necessary to specify an additional option: ::
	make MAXK=31 MAC=1 metacortex


