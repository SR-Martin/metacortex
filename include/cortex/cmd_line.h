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

/************************************************************************
 * cmd_line.h
 *
 * Manipulation of command line.
 ************************************************************************/

#ifndef CMD_LINE_H_
#define CMD_LINE_H_

#include <stdio.h>
#include <global.h>
#include <file_format.h>

typedef enum{
    PERFECT_PATH = 0,
    BRANCHES = 1,
    Y_WALK = 2,
    BUBBLES = 3,
    READ_PAIR = 5,
    METACORTEX_CONSENSUS = 6,
    GRAPH_STATS = 7,
    SUBTRACTIVE_WALK = 8,           
}TraverseAlgorithm;

typedef struct
{
    //core parameters
    boolean verbose;
    int kmer_size;
    int bucket_size;
    int number_of_buckets_bits;
    
    //----------------
    //actions on/off
    //----------------
    //input
    boolean input_file; //if it is present
    boolean input_file_format_known;
    boolean health_check_binary;
    boolean input_reference_known;
    
    //cleaning
    boolean tip_clip;
    boolean remove_bubbles;
    boolean low_coverage_path_clip;
    boolean low_coverage_node_clip;
    boolean remove_low_coverage_supernodes;
    boolean remove_spurious_links;
    
    //output
    boolean dump_binary;
    boolean high_confidence;
    boolean output_fasta;
    boolean dump_hash;
    boolean detect_bubbles;
    boolean output_coverages;
    boolean graphviz;
    boolean print_uncertain_as_n;
    boolean output_log;
    boolean output_kmer_coverage_know;
    boolean multiple_subgraph_contigs;
    boolean gfa_and_fastg;
    
    //-----------
    //parameters
    //-----------
    //core parameters
    int threads;
    short max_double_y_complexity;
    
    //input
    FileFormat input_file_format;
    int quality_score_threshold;
    char input_filename[LENGTH_FILENAME];
    char qual_filename[LENGTH_FILENAME];
    char input_reference[LENGTH_FILENAME];
    
    int quality_score_offset;
    int max_read_len;
    int binary_version;
    
    //cleaning
    int node_coverage_threshold;
    int max_node_edges;
    float delta_coverage;
    float subtractive_walk_delta;
    int linked_list_max_size;
    int tip_length;
    int remove_low_coverage_supernodes_threshold;
    int bubble_max_depth;
    int bubble_max_length;
    int remove_spurious_links_min_coverage;
    int remove_spurious_links_max_difference;
    int path_coverage_minimum;
    int tip_clip_iterations;
    
    //output
    char output_ctx_filename[LENGTH_FILENAME];
    char output_fasta_filename[LENGTH_FILENAME];
    char output_graphviz_filename[LENGTH_FILENAME];
    char output_hash_filename[LENGTH_FILENAME];
    char log_filename[LENGTH_FILENAME];
    char output_kmer_coverage[LENGTH_FILENAME];
    char output_reference_coverage_file[LENGTH_FILENAME];
    int max_length;
    int singleton_length;
    int min_subgraph_size;
    int min_contig_length;
    TraverseAlgorithm algorithm;

    unsigned long hashtable_max_mem;

} CmdLine;

CmdLine parse_cmdline( int argc, char* argv[],int unit_size);
int default_opts(CmdLine *);
unsigned long long parse_mem_string(char* mem_string);

#endif /* CMD_LINE_H_ */
