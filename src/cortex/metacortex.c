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
#include <time.h>
#include <string.h>
#include <locale.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <file_reader.h>
#include <perfect_path.h>
#include <branches.h>
#include <y_walk.h>
#include <logger.h>
#include <assert.h>
#include <cleaning.h>
#include <cmd_line.h>
#include "metagraphs.h"
#include "graph_stats.h"
#include "bubble_find.h"
#include "subtractive_walk.h"

typedef struct {
    char* filename;
    int colour;
    int pair;
} ReadFileDescriptor;

#define FILE_LIST_SIZE 1024
#define TOTAL_MAX_LENGTH 10000  // bubble_find.c
#define BUBBLE_MAX_DEPTH 10     // bubble_find.c

void write_graphviz_file(char *filename, dBGraph * db_graph);
void timestamp();

char *remove_file_extension(char *filename)
{
	int end = strlen(filename) - 1;
	int start = end - 6;
	int i;

	if (DEBUG)
		printf("[remove_file_extension] In:  %s\n", filename);

	if (start < 0)
		start = 0;

	for (i = end; i >= start; i--) {
		if (filename[i] == '.') {
			filename[i] = 0;
			break;
		}
	}

	if (DEBUG)
		printf("[remove_file_extension] Out: %s\n", filename);

	return filename;
}

boolean add_read_file(char* filename, int colour, int pair, int* n_files, ReadFileDescriptor* list[])
{
    boolean already_seen = false;
    FILE *fp;
    int i;

    // Check if file already in list
    for (i=0; i<*n_files; i++) {
        if (strcmp(filename, list[i]->filename) == 0) {
            already_seen = true;
            printf("Error: file %s is included more than once.\n", filename);
            exit(-1);
        }
    }

    // Check file exists
    fp = fopen(filename, "r");
    if (fp) {
        fclose(fp);
    } else {
        printf("Error: can't open file %s\n", filename);
        exit(-1);
    }

    // Add to list if not already seen
    if (!already_seen) {
        if (*n_files == FILE_LIST_SIZE) {
            printf("Error: maximum file list size reached. Continuing might be dangerous.\n");
            exit(-1);
        } else {
            list[*n_files] = malloc(sizeof(ReadFileDescriptor));
            if (!list[*n_files]) {
                printf("Error: Couldn't allocate memory for file descriptor\n");
                exit(-1);
            }

            list[*n_files]->filename = malloc(strlen(filename)+1);
            if (!list[*n_files]->filename) {
                printf("Error: Couldn't allocate memory for filename\n");
                exit(-1);
            }
            strcpy(list[*n_files]->filename, filename);
            list[*n_files]->colour = colour;
            list[*n_files]->pair = pair;

            *n_files = *n_files + 1;
        }
    }

    return already_seen;
}

void output_basic_info(CmdLine cmd_line)
{
    log_and_screen_printf("Max k: %i\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)-1);
    log_and_screen_printf("max_node_edges: %i\n", (cmd_line.max_node_edges));
    log_and_screen_printf("delta_coverage: %f\n", (cmd_line.delta_coverage));
    if (cmd_line.input_file_format == FASTQ) {
        log_and_screen_printf("Quality score offset: %i", cmd_line.quality_score_offset);
        if (cmd_line.quality_score_offset == 33) {
            log_and_screen_printf(" (Sanger format)");
        } else if (cmd_line.quality_score_offset == 64) {
            log_and_screen_printf(" (Solexa/Illumina format)");
        } else {
            log_and_screen_printf(" (Unknown format)");
        }
        log_and_screen_printf("\nQuality score threshold: %i\n", cmd_line.quality_score_threshold);
    }
}

int main(int argc, char **argv)
{
	setlocale (LC_ALL, "");
    FILE *fp_fnames;
    dBGraph *db_graph = NULL;
    short kmer_size;
    CmdLine cmd_line;
    boolean using_colours = false;
    ReadFileDescriptor* file_list[FILE_LIST_SIZE];
    long long seq_length = 0;
    int n_file_list = 0;
    int i;
    //int total_max_length = TOTAL_MAX_LENGTH;
    //int bubble_max_depth = BUBBLE_MAX_DEPTH;

    log_and_screen_printf("\nMetaCortex ");
    log_and_screen_printf(METACORTEX_VERSION);
    log_and_screen_printf("\n");
	log_and_screen_printf("Compiled on %s at %s \n\n", __DATE__, __TIME__);

    //command line arguments
    cmd_line = parse_cmdline(argc, argv, sizeof(Element));

    output_basic_info(cmd_line);

    if(cmd_line.input_file_format == HASH) {
        db_graph = hash_table_read_dumped_memory(cmd_line.input_filename);
    } else {
        fp_fnames = fopen(cmd_line.input_filename, "r");	//open file of file names
        kmer_size = cmd_line.kmer_size;
       // DEBUG = cmd_line.verbose;

        timestamp();
        log_and_screen_printf("\nInput file of filenames: %s\n", cmd_line.input_filename);
        log_and_screen_printf("Kmer size: %'d Hash table size (%'d bits): %'d Hash table bucket size: %'d Total size: %'qd\n",
                              cmd_line.kmer_size,
                              cmd_line.number_of_buckets_bits,
                              1 << cmd_line.number_of_buckets_bits,
                              cmd_line.bucket_size,
                              ((long long)1 << cmd_line.number_of_buckets_bits) * cmd_line.bucket_size);
        fflush(stdout);

        //Create the de Bruijn graph/hash table
        db_graph = hash_table_new(cmd_line.number_of_buckets_bits, cmd_line.bucket_size, 25, cmd_line.kmer_size);

        if (db_graph == NULL) {
            printf("Please free up memory and re-run.\n");
            exit(-1);
        }

        //TODO
        if (cmd_line.threads > 0) {
        	db_graph->number_of_threads = cmd_line.threads;
        }

        //Zone to initalise the buffers
        path_array_initialise_buffers(kmer_size);
        timestamp();
        log_and_screen_printf("\nTable created: %'d\n", 1 << cmd_line.number_of_buckets_bits);
        db_graph_print_status(db_graph);
        fflush(stdout);
        int count_file = 0;
        db_graph->loaded_kmers = 0;	//total sequence length

        //Go through all the files, loading data into the graph

        boolean all_entries_are_unique = false;

        if (cmd_line.input_file_format == CTX) {
            all_entries_are_unique = true;
        }

        long long bad_reads = 0;

        while (!feof(fp_fnames)) {
            short colour = 0;
            short pair = 0;
            char filename[1024];
            char line[1024];
            count_file++;

            line[0] = 0;
            if (fgets(line, 1024, fp_fnames)) {
                if (strlen(line) > 1) {
                    sscanf(line, "%s %hd %hd\n", filename, &colour, &pair);
                    add_read_file(filename, colour, pair, &n_file_list, file_list);
                }
            }
        }

        fclose(fp_fnames);

        for (i=0; i<n_file_list; i++) {
            short colour = file_list[i]->colour;
            char* filename = file_list[i]->filename;


            if (DEBUG) {
                log_and_screen_printf("\nNew file: %s colour %'d\n", filename, colour);
                fflush(stdout);
            }

            fflush(stdout);

            if ((colour < 0) || (colour >= NUMBER_OF_COLOURS)) {
                printf("Colour is out of allowable range (maximum number of colours: %'d).", NUMBER_OF_COLOURS);
                fflush(stdout);
                exit(-1);
            }

            if (colour > 0) {
                using_colours = true;
            }

            switch (cmd_line.input_file_format) {

                case CTX:
                    log_and_screen_printf("\nReading ctx file %'d: %s\n", i+1, filename);
				    fflush(stdout);
                    seq_length = load_binary_from_filename_into_graph(filename, db_graph, colour, all_entries_are_unique);

                    all_entries_are_unique = false;
                    break;

                case FASTQ:
                    if (colour != 0) {
                        printf("\n*** Warning: Colour should be 0 for FASTQ files.\n");
                    }
                    log_and_screen_printf("\nReading fastq file %'d: %s\n", i+1, filename);
				    fflush(stdout);
                    seq_length = load_fastq_from_filename_into_graph(filename, colour, &bad_reads, cmd_line.quality_score_threshold, 5000, cmd_line.quality_score_offset, db_graph);
                    break;

                case FASTA:
                    if (colour != 0) {
                        printf("\n*** Warning: Colour should be 0 for FASTA files.\n");
                    }
                    log_and_screen_printf("\nReading fasta file %'d: %s\n", i+1, filename);
				    fflush(stdout);
                    seq_length = load_fasta_from_filename_into_graph(filename, colour, &bad_reads, 5000, db_graph);
                    break;
                default:
                    printf("Unknown file format. ");
                    fflush(stdout);
                    break;
            }

            db_graph->loaded_kmers += seq_length;
            timestamp();
            log_and_screen_printf("\nRead of file %'d complete. Total kmers: %'lld Bad reads: %'qd Seq length: %'qd Total seq length: %'qd\n\n",
                                  i+1, hash_table_get_unique_kmers(db_graph), bad_reads, seq_length, db_graph->loaded_kmers);
            hash_table_print_stats(db_graph);

            fflush(stdout);


        }
    }
	path_array_initialise_buffers(db_graph->kmer_size);
    db_graph->max_double_y_complexity = cmd_line.max_double_y_complexity;
	fflush(stdout);

	if (cmd_line.low_coverage_node_clip) {
		timestamp();
		log_and_screen_printf("\nRemove low coverage nodes (<=%'d) \n", cmd_line.node_coverage_threshold);
		fflush(stdout);
		cleaning_remove_low_coverage(cmd_line.node_coverage_threshold, db_graph);
        hash_table_print_stats(db_graph);

	}

    if (cmd_line.remove_low_coverage_supernodes) {
		timestamp();
		log_and_screen_printf("\nRemoving paths with low coverage...\n");
		fflush(stdout);
		int p = cleaning_prune_low_coverage_path(cmd_line.remove_low_coverage_supernodes_threshold, cmd_line.tip_length, db_graph);
		log_and_screen_printf("%'d nodes removed\n", p);
        hash_table_print_stats(db_graph);
	}

	if (cmd_line.tip_clip) {
		timestamp();
		log_and_screen_printf("\nClip tips\n");
		fflush(stdout);
        log_and_screen_printf("%'d tips clipped\n", cleaning_remove_tips(cmd_line.tip_length, cmd_line.tip_clip_iterations ,db_graph));
        hash_table_print_stats(db_graph);
	}

	if (cmd_line.remove_bubbles) {
		timestamp();
		log_and_screen_printf("\nRemoving bubbles\n");
		fflush(stdout);
		cmd_line.bubble_max_length = db_graph->kmer_size * 10 + 1;
		//TODO: we should defined the bubble depth here too
		int kmers_removed_1 = cleaning_remove_bubbles(cmd_line.bubble_max_length, db_graph);
		log_and_screen_printf("%'d kmers removed\n", kmers_removed_1);
        hash_table_print_stats(db_graph);
        if (cmd_line.tip_clip) {
            timestamp();
            log_and_screen_printf("\nClip tips\n");
            fflush(stdout);
            log_and_screen_printf("%'d tips clipped\n", cleaning_remove_tips(cmd_line.tip_length, cmd_line.tip_clip_iterations ,db_graph));
            hash_table_print_stats(db_graph);
        }
    }

    if(cmd_line.dump_hash){
		timestamp();
		log_and_screen_printf("\nDumping hash table to file: %s\n", cmd_line.output_hash_filename);
		hash_table_dump_memory(cmd_line.output_hash_filename, db_graph);
		fflush(stdout);
	}

    hash_table_print_stats(db_graph);

    if (cmd_line.dump_binary) {
        timestamp();
        log_and_screen_printf("\n");
        // If using colours, we dump a separate file for each colour. Otherwise, it's just one.
        if (!using_colours) {
            log_and_screen_printf("Dumping graph: %s\n", cmd_line.output_ctx_filename);
            db_graph_dump_binary(cmd_line.output_ctx_filename, &db_node_check_flag_not_pruned, db_graph);
            fflush(stdout);
        } else {
            int c;
            for (c = 0; c < NUMBER_OF_COLOURS; c++) {
                char fname[1024];
                char suffix[64];
                char *extension;

                sprintf(suffix, "_c%d.ctx", c);
                strcpy(fname, cmd_line.output_ctx_filename);
                extension = strstr(fname, ".ctx");
                if (!extension) {
                    extension = strstr(fname, ".CTX");
                }

                if (extension) {
                    strcpy(extension, suffix);
                } else {
                    strcat(fname, suffix);
                }

                log_and_screen_printf("Dumping graph %s - ", fname);
                db_graph_dump_binary_by_colour(fname, &db_node_check_flag_not_pruned, c, db_graph);
            }
            fflush(stdout);
        }
    }

    // Write graphviz file?
    if (cmd_line.graphviz) {
      timestamp();
      write_graphviz_file(cmd_line.output_graphviz_filename, db_graph);
      fflush(stdout);
    }

    if (cmd_line.output_fasta) {
        db_graph->path_coverage_minimum = cmd_line.path_coverage_minimum;
        switch (cmd_line.algorithm) {
            case PERFECT_PATH:
                log_and_screen_printf("\nDumping supernodes: %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);

                perfect_path_print_paths(cmd_line.output_fasta_filename,
                                         cmd_line.max_length, cmd_line.singleton_length,
                                         cmd_line.output_coverages,
                                         db_graph);

                break;
            case Y_WALK:
                log_and_screen_printf("\nDumping supernodes (Y_WALK): %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);
                y_walk_print_paths(cmd_line.output_fasta_filename,
                                   cmd_line.max_length, cmd_line.singleton_length, cmd_line.gfa_and_fastg,
                                   cmd_line.output_coverages, cmd_line.graphviz, db_graph);
                log_and_screen_printf("Dump complete\n");
                break;
            case BRANCHES:
                log_and_screen_printf("\nDumping supernodes (branches): %s\n", cmd_line.output_fasta_filename);
                fflush(stdout);
                //branches_get_path(cmd_line.output_fasta_filename, cmd_line.max_length, cmd_line.output_coverages, db_graph);
                break;
            case METACORTEX_CONSENSUS:
                log_and_screen_printf("\nDumping subgraph consensus contigs: %s\n", cmd_line.output_fasta_filename);
                metacortex_find_subgraphs(db_graph, cmd_line.output_fasta_filename, cmd_line.min_subgraph_size, cmd_line.min_contig_length, 
                                                    cmd_line.multiple_subgraph_contigs, cmd_line.gfa_and_fastg);
                break;
            case GRAPH_STATS:
                // if --high_confidence option set by user, alter output name accordingly
                if(cmd_line.high_confidence) {
                    char temp[LENGTH_FILENAME];
                    strcpy(temp, cmd_line.output_fasta_filename);
                    sprintf(cmd_line.output_fasta_filename, "%s.high_conf", temp);
                }

                log_and_screen_printf("\nSearching graph for stats...\n");
                find_subgraph_stats(db_graph, cmd_line.output_fasta_filename,
                  cmd_line.min_subgraph_size, cmd_line.min_contig_length,
                  cmd_line.max_node_edges, cmd_line.delta_coverage,
                  cmd_line.linked_list_max_size,
                  cmd_line.multiple_subgraph_contigs,
                  cmd_line.gfa_and_fastg);


                /* Put all of this into a separate command line call (BUBBLEFIND)*/
                //db_graph_identify_branches(1000, db_graph);
                //log_and_screen_printf("\nSearching graph for branches and bubbles...\n");
                //db_graph_walk_branches(cmd_line.output_fasta_filename, total_max_length, db_graph->kmer_size * 10 + 1, bubble_max_depth, db_graph);

                //db_graph_walk_branches(char *filename, int total_max_length, int bubble_max_length, int bubble_max_depth, dBGraph * db_graph)
                //metacortex_find_subgraphs(db_graph, cmd_line.output_fasta_filename, cmd_line.min_subgraph_size, cmd_line.min_contig_length, cmd_line.multiple_subgraph_contigs);
                break;
            case SUBTRACTIVE_WALK:

                log_and_screen_printf("\nPerforming iterative subtractive walk..\n");
                log_and_screen_printf("Using value %f for determining covering variants.\n", cmd_line.subtractive_walk_delta); 
                subtractive_walk(db_graph, cmd_line.output_fasta_filename, cmd_line.min_contig_length, cmd_line.subtractive_walk_delta);
                break;
                
            default:
                log_and_screen_printf("Algorithm not implemented \n");
                break;
        }
    }

    // Stop program terminating, so XCode leaks tool can report!
    //printf("press char to continue");
    //char c = getchar();
    //printf("Character pressed %c\n", c);

    timestamp();
    log_and_screen_printf("\n\nDONE\n\n");

    return 0;
}

void print_graph(dBGraph * db_graph)
{
	void print_graphviz(dBNode * node) {
		if (node != NULL && element_get_coverage_all_colours(node) > 0) {
			BinaryKmer tmp, co;

			short kmer_size = db_graph->kmer_size;
			binary_kmer_assignment_operator(tmp, node->kmer);
			binary_kmer_reverse_complement(&tmp, kmer_size, &co);
			char seq[kmer_size + 1], seqNext[kmer_size + 1],
			    seq1[kmer_size + 1];
			binary_kmer_to_seq(&tmp, kmer_size, seq1);
			char *print = db_node_check_for_any_flag(node,
								 STARTING_FORWARD
								 |
								 BRANCH_NODE_FORWARD
								 |
								 BRANCH_NODE_REVERSE
								 |
								 END_NODE_FORWARD
								 |
								 END_NODE_REVERSE
								 | X_NODE) ?
			    "ellipse" : "point";

			char *color = db_node_check_for_any_flag(node,
								 STARTING_FORWARD
								 )
			    ? "turquoise4" : db_node_check_for_any_flag(node,
									X_NODE)
			    ? "red" : db_node_check_for_any_flag(node,
								 BRANCH_NODE_FORWARD
								 |
								 BRANCH_NODE_REVERSE)
			    ? "orange" : db_node_check_for_any_flag(node,
								    END_NODE_FORWARD
								    |
								    END_NODE_REVERSE)
			    ? "royalblue" : "black";
			printf("%s [label=%s, shape=%s, color=%s]\n", seq1,
			       seq1, print, color);
			binary_kmer_to_seq(&tmp, kmer_size, seq);
			binary_kmer_left_shift(&tmp, 2, kmer_size);
			Orientation o = forward;
			Edges e = db_node_get_edges_for_orientation(node, o);
			Edges e2 =
			    db_node_get_edges_for_orientation(node, reverse);
			/*
			 */
			void hasEdge1(Nucleotide n) {
				Edges et = e;
				BinaryKmer bk;
				Key k = &bk;
				if (((et >> n) & 1) == 1) {
					binary_kmer_modify_base(&tmp, n,kmer_size, 0);
					binary_kmer_to_seq(element_get_key(&tmp, kmer_size, k), kmer_size, seqNext);
					printf("%s -> %s [ label = \"%c\", color=\"%s\"];\n",
					       seq, seqNext,
					       binary_nucleotide_to_char(n),
					       "blue");
				}
			}
			nucleotide_iterator(&hasEdge1);
			binary_kmer_assignment_operator(tmp, co);
			binary_kmer_left_shift(&tmp, 2, kmer_size);
			o = reverse;
			e = e2;
			void hasEdge2(Nucleotide n) {
				Edges et = e;
				BinaryKmer bk;
				Key k = &bk;
				if (((et >> n) & 1) == 1) {
					binary_kmer_modify_base(&tmp, n,
								kmer_size, 0);
					binary_kmer_to_seq(element_get_key
							   (&tmp, kmer_size, k),
							   kmer_size, seqNext);
					printf
					    ("%s -> %s [ label = \"%c\", color=\"%s\"];\n",
					     seq, seqNext,
					     binary_nucleotide_to_char(n),
					     "red");
				}
			}
			nucleotide_iterator(&hasEdge2);

		}
	}
	hash_table_traverse(&print_graphviz, db_graph);
}

void write_graphviz_file(char *filename, dBGraph * db_graph)
{
	db_graph_write_graphviz_file(filename, db_graph);
}

void timestamp() {
	time_t ltime = time(NULL);
	log_and_screen_printf("\n-----\n%s",asctime(localtime(&ltime)));
	fflush(stdout);
}
