/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk) and
 *     Samuel Martin (samuel.martin@earlham.ac.uk)
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
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>
#include <libgen.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "cleaning.h"
#include "coverage_walk.h"
#include "dB_graph.h"
#include "element.h"
#include "graph_stats.h"
#include "logger.h"
#include "metacortex.h"
#include "metagraphs.h"
#include "subtractive_walk.h"



#define SUBTRACTIVE_WALK_QUEUE_SIZE 20000000 // 10000000


// ----------------------------------------------------------------------
// Work through graph, count coverage, X, Y nodes
// ----------------------------------------------------------------------
void subtractive_walk(dBGraph * graph, char* consensus_contigs_filename, int min_contig_size, float delta_coverage)
{
    FILE* fp_contigs_fasta;
    int counter = 0;
    int min_path_size = min_contig_size - graph->kmer_size;
    min_contig_size = min_contig_size > graph->kmer_size + 1 ? min_contig_size : graph->kmer_size + 1;
    
    Queue* graph_queue = queue_new(SUBTRACTIVE_WALK_QUEUE_SIZE, sizeof(dBNode*));
    if (!graph_queue) {
        log_and_screen_printf("Couldn't get memory for graph queue.\n");
        exit(-1);
    }

    Path *simple_path = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
   // Path *path_fwd = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    //Path *path_rev = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);


   /* Open contigs file */
    fp_contigs_fasta = fopen(consensus_contigs_filename, "w");
    if (!fp_contigs_fasta) {
        log_and_screen_printf("ERROR: Can't open contig file.\n%s\n", consensus_contigs_filename);
        exit(-1);
    }   

    db_graph_reset_flags(graph);

    // Hash table iterator to walk graphs, produce paths
    void traversal_for_contigs(dBNode * node) 
    {
        uint32_t coverage = element_get_coverage_all_colours(node);     
        if(coverage >= graph->path_coverage_minimum)
        {
            
            dBNode* seed_node = NULL;

            /* Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point */
            //log_printf("Growing graph from node\n");
            graph_queue->number_of_items = 0;
            int nodes_in_graph = grow_graph_from_node(node, &seed_node, graph, graph_queue, MAX_EXPLORE_BRANCHES);
            
            if(nodes_in_graph > min_path_size && seed_node)
            {
                node = seed_node;
                coverage_walk_get_path(node, forward, NULL, graph, simple_path, false);
               // coverage_walk_get_path(node, reverse, NULL, graph, path_rev, false);

               // path_reverse(path_fwd, simple_path);
               // path_append(simple_path, path_rev);

                simple_path->id = counter++;

                if(simple_path->length > min_path_size)
                {
                	log_printf("Write path of size %d\n", simple_path->length);

                    // don't print coverage statistics, as these will change as the algorithm proceeds.
                    path_to_fasta_with_statistics(simple_path, fp_contigs_fasta, 0, 0, 0, false);
                }

                    
				uint32_t min_cov = element_get_coverage_all_colours(seed_node);
				int min_index = -1;
				for(int i = 0; i < simple_path->length; i++)
				{
					dBNode* current_node = simple_path->nodes[i];
					uint32_t cov = element_get_coverage_all_colours(current_node);
					if(cov <= min_cov)
					{
						min_cov = cov;
						min_index = i;
					}
				}
				int levels [simple_path->length];
				//float deltas [simple_path->length];

				int current_level = 1;
				levels[min_index] = current_level;
				uint32_t last_cov = min_cov;
				for(int i = min_index + 1; i < simple_path->length; i++)
				{
					dBNode* current_node = simple_path->nodes[i];
					uint32_t cov = element_get_coverage_all_colours(current_node);
					int diff = cov - last_cov;
					int denominator = cov > last_cov ? cov : last_cov;
					float delta = (float)diff/denominator;
					//deltas[i] = delta;
					if(delta > delta_coverage)
					{
						current_level++;
					}
					else if(delta < -delta_coverage && current_level > 1)
					{
						current_level--;
					}
					levels[i] = current_level;
					last_cov = cov;
				}
				current_level = 1;
				last_cov = min_cov;
				for(int i = min_index - 1; i >= 0; i--)
				{
					dBNode* current_node = simple_path->nodes[i];
					uint32_t cov = element_get_coverage_all_colours(current_node);
					int diff = cov - last_cov;
					int denominator = cov > last_cov ? cov : last_cov;
					float delta = (float)diff/denominator;
					//deltas[i] = delta;
					if(delta > delta_coverage)
					{
						current_level++;
					}
					else if(delta < -delta_coverage && current_level > 1)
					{
						current_level--;
					}
					levels[i] = current_level;
					last_cov = cov;
				}


				//debug hist
				/*
				if(simple_path->length > 100000)
				{
					char* filename;
					asprintf(&filename, "node_%qd.hist", simple_path->id);
					FILE* hist_file = fopen(filename, "w");
					int last_coverage = 0;
					for(int n = 0; n < simple_path->length; n++)
					{
						dBNode* current_node = simple_path->nodes[n];
						Orientation current_orientation = simple_path->orientations[n];
						uint32_t coverage = element_get_coverage_all_colours(current_node);
						int num_edges = db_node_edges_count_all_colours(current_node, current_orientation);
						fprintf(hist_file, "%u\t%f\t%i\t%i\n", coverage, deltas[n], levels[n], num_edges);
						last_coverage = coverage;
					}
					fclose(hist_file);
				}
				*/

				//Now subtract the coverages:
				// level 1 -> 0
				// all others get reduced by last min_cov;
				last_cov = min_cov;
				int i = 0;
				while( i < simple_path->length)
				{
					dBNode* current_node = simple_path->nodes[i];
					if(levels[i] == 1)
					{
						last_cov = element_get_coverage_all_colours(current_node);
						current_node->coverage[0] = 0;
						i++;
					}
					else
					{
						// subtract an amount from the nodes with more then 1 covering
						// variant. The amount to subtract is a lerp between the
						// last level 1 and the next level 1.
						int j = i + 1;
						while(j < simple_path->length && levels[j] > 1)
						{
							j++;
						}
						int next_cov = last_cov;
						if( j < simple_path->length)
						{
							next_cov = element_get_coverage_all_colours(simple_path->nodes[j]);
						}
						if( i == 0 )
						{
							last_cov = next_cov;
						}

						int diff = j - i;
						int increment = (next_cov - last_cov) / diff;
						for(int k = i; k < j; k++)
						{
							current_node = simple_path->nodes[k];
							int subtract = last_cov + ((k-i) * increment);
							if(current_node->coverage[0]  > subtract)
							{
								current_node->coverage[0] -= subtract;
							}
							else
							{
								current_node->coverage[0] = 0;
							}
						}
						i = j;
					}
                }
             
                /* Reset paths */
                path_reset(simple_path);
               // path_reset(path_fwd);
                //path_reset(path_rev);
            }
            
            while (graph_queue->number_of_items > 0) 
            {
                dBNode* queue_node = (dBNode*)queue_pop(graph_queue);
                db_node_action_unset_flag(queue_node, VISITED);
            }
        }
    } // traversal_for_contigs

    db_graph_reset_flags(graph);
    log_and_screen_printf("Full traversal started...");
    hash_table_traverse(&traversal_for_contigs, graph);
    log_and_screen_printf("DONE\n");
    
    db_graph_reset_flags(graph);
    fclose(fp_contigs_fasta);

    queue_free(graph_queue);
}
