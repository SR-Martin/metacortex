/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk)
 *
 * MetaCortex is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MetaCortex is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MetaCortex.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 *
 * This file is modified from source that was part of CORTEX. The
 * original license notice for that is given below.
 *
 ************************************************************************
 *
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 *
 ************************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <getopt.h>
#include <cmd_line.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <binary_kmer.h>
#include <element.h>
#include <path.h>
#include <err.h>
#include <logger.h>

void usage(void)
{
    printf("\nusage: metacortex [-h] [--input file_of_files] [--mem_height n] [--dump_binary bin_output] [--input_format fastq|fasta|binary] [--output_contigs contigs.fa] \n");
    printf("\n" \
    "   [-h|--help] = This help screen.\n" \
    "   [-a|--remove_bubbles] = Removes the bubbles in the graph.\n"\
    "   [-A|--algorithm ALG] = Traversal algorithm for contigs (MCC | SW | PP | GS).\n"
    "   [-b|--mem_width INT] = Size of hash table buckets (default 100).\n" \
    "   [-c|--tip_clip INT] = Clips the tips in the graph, the argument defines the max length for the tips.\n" \
    "   [-C|--high_confidence INT] = Only outputs contigs which have no nodes with coverage below THESHOLD.\n" \
    "   [-e|--output_coverages] = Print coverages for contigs/supernodes in a different file with _cov suffix.\n" \
    "   [-d|--ouput_contigs FILENAME] = Fasta file with all the contigs (after applying all specified actions on graph).\n" \
    "   [-f|--ouput_supernodes FILENAME] = Fasta file with all the supernodes (after applying all specified actions on graph).\n" \
    "   [-g|--min_contig_length INT] = minimum contig length produced.\n" \
    "   [-G|--gfa] = Output gfa2 and fastg sequence files.\n"
    "   [-i|--input FILENAME] = File of filenames to be processed (start and end read is optional, format <filename>  <start read index>  <end read index> ).\n" \
    "   [-k|--kmer_size INT] = Kmer size (default 21), it has to be an odd number.\n" \
    "   [-M|--multiple_subgraph_contigs] = Find multiple contigs for each disconnected subgraph during MCC.\n"
    "   [-n|--mem_height INT] = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n" \
    "   [-o|--dump_binary FILENAME] = Dump binary for graph in file (after applying all specified actions on graph).\n" \
    "   [-O|--hash_output_file ] = Dumps the whole graph into a file. Read with the input_format hash. The file stores the information required to restore the hash table, hence mem_height and mem_width don't have any effect."
    "   [-p|--quality_score_offset] = Fastq quality offset. Default 33. Use 63 for illumina.\n"\
    "   [-q|--quality_score_threshold INT] = Filter for quality scores in the input file, any k-mer wiht a base with quality in the threshold or smaller is not considered (default 0).\n"\
    "   [-q|--quality_score_threshold ] = The minimiun phred value for a base to be considered in an assembly.\n"\
    "   [-t|--input_format FORMAT] = File format for input (binary | fasta | fastq | hash ).\n" \
    "   [-u|--remove_seq_errors] = remove sequence of kmers induced by errors. Equivalent to --remove_low_coverage_kmers 1\n"\
    "   [-W|--SW_delta FLOAT] = minimum normalised coverage difference for SW.\n"
    "   [-z|--remove_low_coverage_kmers INT] = Filter for kmers with coverage in the threshold or smaller.\n" \
    "   [-Z|--max_read_len] = Maximum read length over all input files.\n"\
	"	[--max_db_mem] = Maximum memory to be allocated for db-graph, in bytes or with K|M|G|T appended.\n"
    "\n");
    printf("\n");
}

int default_opts(CmdLine * c)
{
    //core parameters
    c->verbose = false;
    c->kmer_size = 21;
    c->bucket_size = 100;
    c->number_of_buckets_bits = 10;

    //----------------
    //actions on/off
    //----------------
    //input
    c->input_file = false; //it is not present until we see it
    c->input_file_format_known = false;
    c->health_check_binary = false;
	c->input_reference_known = false;

    //cleaning
    c->tip_clip = false;
    c->remove_bubbles = false;
    c->low_coverage_path_clip = false;
    c->low_coverage_node_clip = false;
    c->remove_low_coverage_supernodes = false;
    c->remove_spurious_links = false;
    c->multiple_subgraph_contigs = false;

    //output
    c->high_confidence = false;
    c->dump_binary = false;
    c->output_fasta = false;
    c->detect_bubbles = false;
    c->output_coverages = false;
    c->dump_hash = false;
    c->graphviz = false;
    c->print_uncertain_as_n = true;
    c->output_log = false;
    c->min_contig_length = 0;
    c->path_coverage_minimum = 1;
    c->gfa_and_fastg = false;

    //-----------
    //parameters
    //-----------
    //threads
    c->threads = 1;
    c->max_double_y_complexity = PATH_MAX_DOUBLE_Y_COMPLEXITY;

    //input
    //c->input_file_format required
    c->quality_score_threshold = 0;
    //c->input_filename required
    //c->qual_filename required
    c->quality_score_offset = 33;
    c->max_read_len = 1000;
    //c->binary_version //TODO

    //cleaning
    //c->node_coverage_threshold required
    c->max_node_edges = 9; //TODO
    c->delta_coverage = 100; //TODO
    c->subtractive_walk_delta = 1.0;
    c->linked_list_max_size = 0; //TODO
    c->tip_length = 100; //TODO
    c->tip_clip_iterations = 100;
    //c->remove_low_coverage_supernodes_threshold required;
    //c->bubble_max_depth required
    //c->bubble_max_length required
    c->remove_spurious_links_min_coverage = 40;
    c->remove_spurious_links_max_difference = 1;

    //output
    //c->output_ctx_filename required
    //c->output_fasta_filename required
    //c->output_graphviz_filename required
    c->singleton_length = 100;
    c->algorithm = METACORTEX_CONSENSUS;
    c->max_length=200000;
    c->min_subgraph_size=0;
    c->hashtable_max_mem = 0;

    return 1;
}

CmdLine parse_cmdline(int argc, char *argv[], int unit_size)
{
    int i;
    printf("Command: ");
    for (i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    printf("\n");
    printf("Unit size: %i\n", unit_size);
    printf("Node size: %lu\n", sizeof(dBNode));

    CmdLine cmd_line;
    default_opts(&cmd_line);

    int opt;
    int longopt_index = 0;

    static struct option long_options[] = {
        {"remove_bubbles", no_argument, NULL, 'a'},
        {"mem_width", required_argument, NULL, 'b'},
        {"tip_clip", required_argument, NULL, 'c'},
        {"high_confidence", required_argument, NULL, 'C'},
        {"output_contigs", required_argument, NULL, 'd'},
        {"output_coverages", no_argument, NULL, 'e'},
        {"output_supernodes", required_argument, NULL, 'f'},
        {"min_contig_length", required_argument, NULL, 'g'},
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"health_check",no_argument,NULL,'j'},
        {"kmer_size", required_argument, NULL, 'k'},
        {"log_file", required_argument, NULL, 'l'},
        {"min_subgraph_size", required_argument, NULL, 'm'},
        {"mem_height", required_argument, NULL, 'n'},
        {"dump_binary", required_argument, NULL, 'o'},
        {"quality_score_offset", required_argument, NULL, 'p'},
        {"quality_score_threshold", required_argument, NULL, 'q'},
        {"max_node_edges", required_argument, NULL, 'r'},  // using still available letters
        {"remove_low_coverage_supernodes", required_argument, NULL,'s'},
        {"input_format", required_argument, NULL, 't'},
        {"remove_seq_errors", no_argument, NULL, 'u'},
        {"verbose", no_argument, NULL, 'v'},
        {"detect_bubbles", required_argument, NULL, 'w'},
        {"singleton_length", required_argument, NULL, 'x'},
        {"linked_list_max_size", required_argument, NULL, 'y'},  // using still available letters
        {"remove_low_coverage_kmers", required_argument, NULL, 'z'},
        {"algorithm",required_argument,NULL,'A'},
        {"gfa", no_argument, NULL, 'G'},
     	{"input_reference", required_argument, NULL, 'H'},
  	{"output_kmer_coverage", required_argument, NULL, 'J'},
        {"remove_spurious_links",required_argument,NULL,'L'},
        {"multiple_subgraph_contigs",no_argument,NULL,'M'},
        {"hash_output_file", required_argument, NULL, 'O'},
        {"tip_clip_iterations", required_argument, NULL, 'P'},
        {"delta_coverage", required_argument, NULL, 'R'},  // using still available letters
        {"threads", required_argument, NULL, 'T'},
        {"graphviz", required_argument, NULL, 'V'},
        {"SW_delta", required_argument, NULL, 'W'},
		{"max_db_mem", required_argument, NULL, 0},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv,
                              "ab:c:d:ef:g:hi:jk:l:m:n:o:p:q:r:s:t:uvw:x:y:z:A:B:C:D:E:FGH:I:J:K:L:MN:O:P:R:TV:W:Z:",
                              long_options, &longopt_index)) >= 0)
    {
        //Parse the default options
        switch (opt) {
        	case 0:
        		if(strcmp(long_options[longopt_index].name, "max_db_mem") == 0)
        		{
					if (optarg)
					{
						unsigned long long hashtable_max_mem = parse_mem_string(optarg);
						if(hashtable_max_mem == 0)
						{
							errx(1, "[--max_db_mem] Could not understand string %s \n Please use format <number><K/M/G/T>", optarg);
						}
						else
						{
							printf("Max DB size is %llu bytes\n", hashtable_max_mem);
							cmd_line.hashtable_max_mem = hashtable_max_mem;
						}
					}
        		}
        		break;
            case 'a':
                cmd_line.remove_bubbles = true;
                break;

            case 'b':	//bucket size
                if (optarg == NULL)
                    errx(1, "[-b | --mem_width] option requires int argument [hash table bucket size]");
                cmd_line.bucket_size = atoi(optarg);
                if (cmd_line.bucket_size == 0 || cmd_line.bucket_size > SHRT_MAX)	//check that -b is not bigger than max_short
                    errx(1, "[-b | --mem_width] option requires 'short' argument bigger than 0");
                break;

            case 'c':	//clip tip
                if (optarg == NULL)
                    errx(1,
                         "[-c | --clip_tip] option requires int argument [max length of tips]");
                cmd_line.tip_length = atoi(optarg);

                if (cmd_line.tip_length <= 0)
                    errx(1,
                         "[-c | --clip_tip] option requires int argument bigger than 0");
                cmd_line.tip_clip = true;
                break;

            case 'd': //dump paths from the tree decomposion of the graph
                cmd_line.output_fasta = true;
                if (optarg==NULL)
                    errx(1,"[-d | --output_contigs] option requires a filename");

                if (strlen(optarg)<LENGTH_FILENAME)
                {
                    strcpy(cmd_line.output_fasta_filename,optarg);
                }
                else
                {
                    errx(1,"[-d | --output_contigs] filename too long [%s]",optarg);
                }

                if (access(optarg,F_OK)==0){
                    errx(1,"[-d | --output_contigs] filename [%s] exist!",optarg);
                }

                cmd_line.algorithm=Y_WALK;
                break;

            case 'e':	//output coverages
                cmd_line.output_coverages = true;
                break;

            case 'f':	//output of supernodes (more conservative paths)
                cmd_line.output_fasta = true;
                if (optarg == NULL)
                    errx(1,
                         "[-f | --output_supernodes] option requires a filename");

                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_fasta_filename, optarg);
                } else {
                    errx(1,"[-f | --output_supernodes] filename too long [%s]",optarg);
                }

                if (access(optarg, F_OK) == 0) {
                    errx(1,"[-f | --output_supernodes] filename [%s] exist!",optarg);
                }

                break;
            case 'g':
                if (optarg == NULL) {
                    errx(1,"[-g | --min_contig_length] option requires int argument");
                }
                cmd_line.min_contig_length = atoi(optarg);
                if(cmd_line.min_contig_length < 0)
                {
                    errx(1,"[-g | --min_contig_length] option must be non-negative!");
                }
                break;

            case 'h':
                usage();
                exit(0);
                break;


            case 'i': //file of filenames
                if (optarg == NULL)
                    errx(1,"[-i | --input] option requires a filename [file of filenames]");

                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.input_filename, optarg);
                } else {
                    errx(1, "[-i | --input] filename too long [%s]", optarg);
                }

                if (access(optarg, R_OK) == -1) {errx(1,"[-i | --input] filename [%s] cannot be accessed", optarg);
                }
                cmd_line.input_file = true;
                break;

            case 'j'://health check of binary
                cmd_line.health_check_binary = true;
                break ;

            case 'k':	//kmer size
                if (optarg == NULL){
                    errx(1,"[-k | --kmer_size] option requires int argument [kmer size]");
                }
                cmd_line.kmer_size = atoi(optarg);

                if (cmd_line.kmer_size == 0){
                    errx(1,"[-k | --kmer_size] option requires int argument bigger than 0");
                }
                if(cmd_line.kmer_size % 2 == 0){
                    errx(1, "[-k | --kmer_size] option can't be even");
                }
                if (cmd_line.kmer_size > (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * 32)) {
                    errx(1, "[-k | --kmer_size] The kmer size is greater than the maximum kmer_size. ");
                }
                break;

            case 'l':
                if (optarg == NULL) {
                    errx(1, "[-l | --log_file] option requires a filename");
                }
                if (strlen(optarg) < LENGTH_FILENAME) {
                    cmd_line.output_log = true;
                    strcpy(cmd_line.log_filename, optarg);
                    log_start(optarg);
                    printf("Log file: %s\n", optarg);
                } else {
                    errx(1, "[-l | --log_file] filename too long");
                }
                break;

            case 'm':
                if (optarg == NULL)
                    errx(1,"[-m | --min_subgraph_size] option requires int argument");
                cmd_line.min_subgraph_size = atoi(optarg);
                break;

            case 'n':	//number of buckets
                if (optarg == NULL)
                    errx(1,"[-n | --mem_height] option requires int argument [hash table number of buckets in bits]");
                cmd_line.number_of_buckets_bits = atoi(optarg);
                break;

            case 'o':	//output of binary ctx
                cmd_line.dump_binary = true;
                if (optarg == NULL)
                    errx(1, "[-o | --output_ctx] option requires a filename");
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_ctx_filename, optarg);
                } else {
                    errx(1, "[-o | --output_ctx] filename too long [%s]",optarg);
                }
                if (access(optarg, F_OK) == 0) {
                    errx(1, "[-o | --output_ctx] filename [%s] exists!",optarg);
                }
                break;

            case 'p':	//quality_score_offset
                if (optarg == NULL){
                    errx(1, "[-p | --quality_score_offset] option requires int argument");
                }
                cmd_line.quality_score_offset = atoi(optarg);

                if (cmd_line.quality_score_offset == 0){
                    errx(1,"[-p | --quality_score_offset] option requires int argument bigger than 0");
                }
                break;

            case 'q':	//quality threshold
                if (optarg == NULL){
                    errx(1, "[-q | --quality_score_threshold] option requires int argument [quality score threshold]");
                }
                cmd_line.quality_score_threshold = atoi(optarg);

                if (cmd_line.quality_score_threshold == 0){
                    errx(1, "[-q | --quality_score_threshold] option requires int argument bigger than 0");
                }
                break;

            case 'r':
                if (optarg == NULL){
                    errx(1,"[-r ] max_node_edges option requires int argument");
                }
                cmd_line.max_node_edges = atoi(optarg);
                break;

            case 's': //remove_low_coverage_supernodes
                cmd_line.remove_low_coverage_supernodes=true;
                if (optarg==NULL)
                    errx(1, "[-s | --remove_low_coverage_supernodes] option requires int argument [coverage threshold]");
                cmd_line.remove_low_coverage_supernodes_threshold= atoi(optarg);

                if (cmd_line.remove_low_coverage_supernodes_threshold <= 0)
                    errx(1, "[-s | --remove_low_coverage_supernodes] option requires intint argument bigger than 0");
                break;

            case 't':
                if (optarg == NULL){
                    errx(1, "[-t  | --input_format binary | fasta | fastq | hash ] option requires a file type");
                }

                cmd_line.input_file_format_known=true;

                if(strcmp(optarg, "binary") == 0){
                    cmd_line.input_file_format = CTX;
                }else if(strcmp(optarg, "fasta") == 0){
                    cmd_line.input_file_format = FASTA;
                }else if(strcmp(optarg, "fastq") == 0){
                    cmd_line.input_file_format = FASTQ;
                }else if(strcmp(optarg, "hash") == 0){
                    cmd_line.input_file_format = HASH;
                }else  {
                    fprintf(stderr, "[-t  | --input_format  binary | fasta | fastq | hash ] invalid option %s ]", optarg);
                    exit(-1);
                }
                break;

            case 'u': //remove_seq_errors --> this is consistent with an option on cortex_var -- it removes low coverate supernodes fixing the length tp k+1 and the coverage at 1
                cmd_line.remove_low_coverage_supernodes=true;
                cmd_line.remove_low_coverage_supernodes_threshold=1;
                break;

            case 'v':
                cmd_line.verbose = true;
                break;
            case 'x':
                if (optarg == NULL) {
                    errx(1, "[-x | --singleton_length] option requires argument");
                }
                cmd_line.singleton_length = atoi(optarg);
                if(optarg == 0){
                    errx(1,  "[-x | --singleton_length] requires an integer greater than 0");
                }
                break;

            case 'y':
                if (optarg == NULL){
                    errx(1,"[-y ] linked_list_max_size option requires int argument");
                }
                cmd_line.linked_list_max_size = atoi(optarg);
                break;

            case 'z':	//node coverage threshold
                if (optarg == NULL){
                    errx(1, "[-z | --node_coverage_threshold] option requires int argument [node coverate cut off]");
                }
                cmd_line.node_coverage_threshold = atoi(optarg);
                cmd_line.low_coverage_node_clip = true;

                if (cmd_line.node_coverage_threshold == 0){
                    errx(1, "[-z | --remove_low_coverage_kmers] option requires int argument bigger than 0");
                }
                break;

            case 'A':             
                if(strcmp(optarg, "MCC") == 0)
                {
                    cmd_line.algorithm = METACORTEX_CONSENSUS;
                }
                else if(strcmp(optarg, "SW") == 0)
                {
                    cmd_line.algorithm = SUBTRACTIVE_WALK;
                }
                else if(strcmp(optarg, "PP") == 0)
                {
                    cmd_line.algorithm = PERFECT_PATH;
                }
                else if(strcmp(optarg, "GS") == 0)
                {
                    cmd_line.algorithm = GRAPH_STATS;
                }
                else
                {
                    errx(1, "[-A  | --algorithm  MCC | SW | PP | GS ] invalid option %s ]", optarg);
                }
                break;

            case 'C':
                if (optarg == NULL) {
                    errx(1, "[-C | --high_confidence] option called without value, using default threshold (2)");
                    cmd_line.path_coverage_minimum=1;
                }
                else{
                  cmd_line.path_coverage_minimum = atoi(optarg);
                }               
                cmd_line.high_confidence=true;
                break;

            case 'G':
               cmd_line.gfa_and_fastg = true;
               break;
            case 'J':
                cmd_line.output_kmer_coverage_know = true;
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_kmer_coverage, optarg);
                } else {
                    errx(1, "[-J | --output_kmer_coverage] filename too long [%s]", optarg);
                }
                break;

            case 'L':
                if (optarg == NULL) {
                    errx(1, "[-L | --remove_spurious_links] option requires argument [max difference,min coverage]");
                } else {
                    sscanf(optarg, "%d,%d", &cmd_line.remove_spurious_links_max_difference, &cmd_line.remove_spurious_links_min_coverage);
                    if (cmd_line.remove_spurious_links_max_difference < 1) {
                        errx(1, "[-L | --remove_spurious_links] max difference must be 1 or more.");
                    } else if (cmd_line.remove_spurious_links_min_coverage < 1) {
                        errx(1, "[-L | --remove_spurious_links] min coverage must be 1 or more.");
                    } else {
                        log_and_screen_printf("Warning: Spurious links algorithm incomplete.\n");
                        cmd_line.remove_spurious_links = true;
                    }
                }
                break;
            case 'M':
                cmd_line.multiple_subgraph_contigs = true;
                break;
            case 'N':
                cmd_line.print_uncertain_as_n = false;
                break;
            case 'O':
                cmd_line.dump_hash = true;
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_hash_filename, optarg);
                } else {
                    errx(1, "[-O | --hash_output_file] filename too long [%s]", optarg);
                }
                break;
            case  'P':
                if (optarg == NULL)
                    errx(1,
                         "[-P | --tip_clip_iterations] option requires int argument ");
                cmd_line.tip_clip_iterations = atoi(optarg);

                if (cmd_line.tip_length <= 0)
                    errx(1,
                         "[-P | --tip_clip_iterations] option requires int argument bigger than 0");

                break;
            case 'Q':
                errx(1, "454 quality files not implemented yet. Use the regular fasta format. \n");
                if (optarg == NULL){
                    errx(1, "[-Q  | --qual_file] option requires a filename [file of filenames]");
                    exit(-1);
                }

                printf("Qual file %s\n", optarg);
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.qual_filename, optarg);
                } else {
                    errx(1, "[-Q | --qual_file] filename too long [%s]", optarg);
                }

                if (access(optarg, R_OK) == -1) {
                    errx(1, "[-Q | --qual_file ] filename [%s] cannot be accessed", optarg);
                }
                break;

            case 'R':
                if (optarg == NULL){
                    errx(1,"[-R ] delta_coverage option requires int argument (as a percentage)");
                }
                cmd_line.delta_coverage = atoi(optarg);
                cmd_line.delta_coverage=cmd_line.delta_coverage/100;
                break;
                
            case 'T':	//Number of threads
                if (optarg == NULL) {
                    errx(1, "[-T | --threads INT] option requires int argumen");
                }
                cmd_line.threads = atoi(optarg);
                if (cmd_line.threads <= 0) {
                    errx(1, "[-T | --threads INT] option requires int argument bigger than 0");
                }
                int tmp = cmd_line.threads  - 1;
                int check = cmd_line.threads & tmp;
                if (check != 0 ) {
                    errx(1, "[-T | --threads INT] has to be a power of two. ");
                }
                break;
            case 'V':
                if (optarg == NULL) {
                    errx(1, "[-V | --graphviz] option requires a filename [file of filenames]");
                }
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.output_graphviz_filename, optarg);
                } else {
                    errx(1, "[-V | --graphviz] filename too long [%s]", optarg);
                }
                cmd_line.graphviz = true;
                printf("Graphviz file: %s\n", optarg);
                break;
            case 'H':
               	if (optarg == NULL) {
                   errx(1, "[-H | --input_reference ] option requires a filename [fasta file reference]");
                }
                if (strlen(optarg) < LENGTH_FILENAME) {
                    strcpy(cmd_line.input_reference, optarg);
                } else {
                    errx(1, "[-H | --input_reference] filename too long [%s]", optarg);
                }
                if (access(optarg, R_OK) == -1) {
                                         errx(1,"[-H | --input_reference] filename [%s] cannot be accessed", optarg);
                }
                printf("Reference file: %s\n", optarg);
                            cmd_line.input_reference_known = true;
                break;
            case 'W':	//SW_delta
                if (optarg == NULL) {
                    errx(1, "[-W | --SW_delta FLOAT] option requires float argument");
                }
                cmd_line.subtractive_walk_delta = atof(optarg);
                if (cmd_line.subtractive_walk_delta < 0 || cmd_line.subtractive_walk_delta > 1 ) {
                    errx(1, "[-W | --SW_delta FLOAT] must be between 0 and 1.");
                }
                break;
            case 'Z'://max read length --> use by the parser to read fasta/fastq files.
                if (optarg==NULL) {
                    errx(1,"[-Z | --max_read_len] option requires int argument [maximum read length in input]");
                }

                cmd_line.max_read_len = atoi(optarg);

                if (cmd_line.max_read_len <= 0)
                    errx(1,"[-Z | --max_read_len] option requires int argument bigger than 0");
                break;

            default:
                // Error already displayed by getopt_long
                exit(1);
                break;
        }
    }

    // We've already sent this to the screen, but the log wasn't open then...
    if (cmd_line.output_log) {
        log_and_screen_printf("\nMetaCortex ");
        log_and_screen_printf(METACORTEX_VERSION);
        log_and_screen_printf("\n\n");

        log_printf("Command: ");
        for (i = 0; i < argc; i++) {
            log_printf("%s ", argv[i]);
        }
        log_printf("\nUnit size: %i\n", unit_size);
    }

    //check if input file format is known
    if (!cmd_line.input_file_format_known) {
        errx(1, "file format not defined [ --input_format fasta | fastq | binary | hash]");
    }
    //check input file
    if (!cmd_line.input_file) {
        errx(1, "input file required [option -i | --input_file]");
    }

    if (cmd_line.number_of_buckets_bits < 0) {
        printf("args error -b %i -h %i \n", cmd_line.bucket_size,cmd_line.number_of_buckets_bits);
        errx(1, "memory configuration error - revise options -b -n or --max_db_mem -b");
    }

    if(cmd_line.input_file_format == ROCHE){
        if(strlen(cmd_line.qual_filename) == 0){
            printf("454 format requires a quality file [option -Q | --qual_file FILENAME] \n");
        }
    }

    //check kmer_size is odd
    if (cmd_line.kmer_size % 2 == 0) {
        errx(1, "[-k | --kmer_size] is even [%i]!", cmd_line.kmer_size);
    }

    return cmd_line;
}

unsigned long long parse_mem_string(char* mem_string)
{
	if(!isdigit(*mem_string))
	{
		return 0;
	}
	else
	{
		unsigned long long value;
		char* end;
		value = strtol(mem_string, &end, 10);
		if(strlen(end) > 1)
		{
			return 0;
		}
		else
		{
			unsigned long long multiplier = 0;
			if(strlen(end) == 0)
			{
				multiplier = 1;
			}
			else if(*end == 'K')
			{
				multiplier = 1000;
			}
			else if(*end == 'M')
			{
				multiplier = 1000 * 1000;
			}
			else if(*end == 'G')
			{
				multiplier = 1000 * 1000 * 1000;
			}
			else if(*end == 'T')
			{
				multiplier = 1000ul * 1000ul * 1000ul * 1000ul;
			}
			return value * multiplier;
		}
	}
}
