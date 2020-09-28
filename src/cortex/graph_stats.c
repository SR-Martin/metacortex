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
#include "global.h"
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "logger.h"
#include "graph_tools.h"
#include "graph_formats.h"
#include "node_queue.h"
#include "coverage_walk.h"
#include "perfect_path.h"
#include "graph_stats.h"
#include "cleaning.h"
#include "metacortex.h"
#include "metagraphs.h"
#include "report_output.h"

#define COVERAGE_BINS 10
#define COVERAGE_BIN_SIZE 1
#define MAX_BRANCHES 10
#define GRAPH_LOG10_LIMIT 10 // little hacky to do this here, because it needs to match size of subgraph_dist in graph_stats.h
#define NUM_BEST_NODES 5

void timestamp_gs() {
    time_t ltime = time(NULL);
    log_printf("\n-----\n%s",asctime(localtime(&ltime)));
    fflush(stdout);
}

///////////////////// literally pasting in Richard's linked-list code for the moment. Unsure of memory issues.
typedef struct _TopItem {
    dBNode * ptr;
    int value;
    struct _TopItem* next;
    struct _TopItem* prev;
} TopItem;

TopItem *start = 0;
// int linked_list_max_size = 1; // should be user defined - top N coverage nodes
int current_size = 0;
void add_item(dBNode* ptr, int value, int linked_list_max_size)
{
    //printf("Adding %s with %d\n", ptr, value);
    TopItem *ti = malloc(sizeof(TopItem));
    if (ti == 0) {
        printf("Error: can't get memory for TopItem!\n");
        exit(1);
    }

    ti->ptr = ptr;
    ti->value = value;
    ti->next = 0;
    ti->prev = 0;

    if ((start == 0) && (linked_list_max_size>0)) {
        start = ti;
        current_size++;
    } else {
        // If value is > than start of list, we need to find where to insert it
        if ((value > start->value) || (current_size < linked_list_max_size)) {
            // Find insertion point
            TopItem* current = start;
            while ((current->next != 0) && (value > current->value)) {
                current = current->next;
            }

            // If value still > current item, then we met the end of the list.
            if (value > current->value) {
                //printf("Inserting after %s\n", current->ptr);
                ti->next = 0;
                ti->prev = current;
                current->next = ti;
                current_size++;
            } else {
                //printf("Inserting before %s\n", current->ptr);

                if (current->prev != 0) {
                    current->prev->next = ti;
                }
                ti->next = current;
                ti->prev = current->prev;

                current->prev = ti;
                current_size++;

                if (current == start) {
                    start = ti;
                }
            }

            // Is list too big now?
            if (current_size > linked_list_max_size) {
                //printf("Removing %s\n", start->ptr);
                TopItem* temp = start;
                start = temp->next;
                free(temp);
                current_size--;
            }

        } else {
            free(ti);
        }
    }
}

void clear_list(dBGraph* graph)
{
    TopItem* current = start;
    TopItem* past = 0;
    char* seq = calloc(256, 1);

    log_printf("HIGH COVERAGE KMERS\n");
    while(current != 0) {
        binary_kmer_to_seq(&(current->ptr->kmer), graph->kmer_size, seq);
        log_printf("%s\n", seq);

        cleaning_prune_db_node(current->ptr, graph);
        past = current;
        current = current->next;
        free(past);
        graph->unique_kmers--;
    }
}

////////////////////// end of linked list functions


/*----------------------------------------------------------------------*
 * Function: grow_graph_from_node_stats                                 *
 * Purpose: takes start node, walks a complete graph from there         *
 *          does not produce paths, or output contigs, just stats       *
 * Params:                                                              *
 * Returns:                                                             *
 *----------------------------------------------------------------------*/
int grow_graph_from_node_stats(dBNode* start_node, dBNode** best_node, dBGraph* graph, Queue* graph_queue, GraphInfo* nodes_in_graph, float delta)
{
    Queue* nodes_to_walk;
    dBNode* node;
    int orientation;
    int depth;
    int best_edges[NUM_BEST_NODES];
    int i;
    *best_node = 0;
    for(i=0; i<NUM_BEST_NODES; i++){
        best_edges[i]=0;
    }
    char* seq = calloc(256, 1);
    float delta_coverage;

    // Nucleotide iterator, used to walk all possible paths from a node
    void walk_if_exists(Nucleotide n) {
        //if (debug) printf("Trying nucleotide %i\n", n);
        int end_orientation;

        // If there is an edge in any colour for this nucleotide...
        if (db_node_edge_exist_any_colour(node, n, orientation)) {

            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode * next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, graph);
            if (!next_node) {
                log_and_screen_printf("Error: Something went wrong with db_graph_get_next_node\n");
                exit(-1);
            }

            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                pathStep first_step;
                Path * new_path;
                dBNode* end_node;
                i = 0;

                // Get path
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                first_step.flags = 0;
                new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
                if (!new_path) {
                    log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(-1);
                }

                db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);

                // check for path coverage here
                uint32_t starting_coverage = element_get_coverage_all_colours(node);

                double path_coverage=0;
                uint32_t min_coverage=0; uint32_t max_coverage=0; // required for path_get_statistics()
               	path_get_statistics(&path_coverage, &min_coverage, &max_coverage, new_path);

                delta_coverage = delta * (float) starting_coverage;
                if (delta_coverage<1){
                    delta_coverage=1;
                }

                min_coverage = starting_coverage - (int) delta_coverage;
                if (min_coverage <1){
                    min_coverage=1;
                }

                max_coverage = starting_coverage + (int) delta_coverage;
                if (((path_coverage >= min_coverage) && (path_coverage <= max_coverage)) || best_node == NULL)  {
                    // Add end node to list of nodes to visit
                    end_node = new_path->nodes[new_path->length-1];
                    end_orientation = new_path->orientations[new_path->length - 1];
                    if (!db_node_check_flag_visited(end_node)) {
                        if (!db_node_is_blunt_end_all_colours(end_node, new_path->orientations[new_path->length-1])) {
                            if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                                log_and_screen_printf("Queue too large. Ending. (WALK)\n");
                                exit(1);
                            }
                        }
                    }


                    // check nodes in path now
                    // only really need to check final node as it's a perfect path
                    // is it blunt? has it been seen before?
                    // things I'm dropping for now - loop detection, branch + loop, catching too large a bubble
                    if (db_node_is_blunt_end_all_colours(end_node, end_orientation)) {
                        // DO NOTHING WITH THIS
                        //db_graph_check_and_add_path(merged_path, patharray);
                    }
                    if (db_node_check_flag_visited(end_node)) {
                        // need to count back from here to original branching point?
                        log_printf("\nBUBBLE FOUND, path length\t%i\n", new_path->length);
                        // length of path here? not perfect - if bubble structure is complex, it will only report on the most immediate perfect path size.
                        // end_node
                        binary_kmer_to_seq(&(end_node->kmer), graph->kmer_size, seq);
                        log_printf("BUBBLE FOUND at kmer %s\n", seq);
                    }


                    // Now go through all nodes, look for best and mark all as visited
                    for (i=0; i<new_path->length; i++) {
                        if (!db_node_check_flag_visited(new_path->nodes[i])) {
                            uint32_t this_coverage = element_get_coverage_all_colours(new_path->nodes[i]);
                            int this_FOR_edges = db_node_edges_count_all_colours(new_path->nodes[i], forward);
                            int this_REV_edges = db_node_edges_count_all_colours(new_path->nodes[i], reverse);


                            // add node degrees to 2D array of all degrees in subgraph
                            nodes_in_graph->node_degree[this_FOR_edges][this_REV_edges]++;

                            // if this is the new best node update the other bests
                            if ((this_coverage > nodes_in_graph->best_coverage[0]) ||
                                ((this_coverage == nodes_in_graph->best_coverage[0]) && ((this_FOR_edges + this_REV_edges) < best_edges[0])))
                            {
                                best_edges[0] = (this_FOR_edges + this_REV_edges);
                                *best_node = new_path->nodes[i];
                            }

                            db_node_action_set_flag_visited(new_path->nodes[i]);
                            queue_push(graph_queue, new_path->nodes[i]);
                            nodes_in_graph->total_size++;
                        }
                    } // path->length loop
                }
                else{
                    cleaning_prune_db_node(new_path->nodes[0], graph);
                    // NOTE best_node - needs to be checked here. Don't want to return NULL
                }

                // Clean up
                path_destroy(new_path);
            }
        }
    } // walk_if_exists

    // Start a queue of nodes to walk
    //log_and_screen_printf("Allocating %d Mb to store queue information (max %d nodes, when full each node could be %d)...\n", ((METACORTEX_QUEUE_SIZE * sizeof(QueueItem*)) / 1024) / 1024, METACORTEX_QUEUE_SIZE, sizeof(QueueItem));
    nodes_to_walk = node_queue_new(METACORTEX_QUEUE_SIZE);
    if (!nodes_to_walk) {
        log_and_screen_printf("Couldn't get memory for node queue.\n");
        exit(-1);
    }

    // Add start node to list of nodes to visit
    if (queue_push_node(nodes_to_walk, start_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(-1);
    }


    if (!db_node_check_flag_visited(start_node)) {
        db_node_action_set_flag_visited(start_node);
        nodes_in_graph->total_size++;
    }

    // Now keep visiting nodes and walking paths
    while (nodes_to_walk->number_of_items > 0) {
        // Take top node from list
        node = queue_pop_node(nodes_to_walk, &depth);

        // Look at all paths out from here
        orientation = forward;
        nucleotide_iterator(&walk_if_exists);
        orientation = reverse;
        nucleotide_iterator(&walk_if_exists);
    }

    queue_free(nodes_to_walk);

    return 0;
}


// ----------------------------------------------------------------------
// Work through graph, count coverage, X, Y nodes
// ----------------------------------------------------------------------
void find_subgraph_stats(dBGraph * graph, char* consensus_contigs_filename,
                         int min_subgraph_kmers, int min_contig_size, int max_node_edges, float delta_coverage,
                         int linked_list_max_size, int walk_paths, boolean gfa_fastg_output)
{
    FILE* fp_analysis;
    FILE* fp_report;
    FILE* fp_degrees;
    FILE* fp_contigs_fasta;
    FILE* fp_contigs_fastg;
    FILE* fp_contigs_gfa;
    long int Contig_Branches[MAX_BRANCHES];
    char* seq = calloc(256, 1);
    long int total_nodes = 0;
    int i;  int j;
    int counter= 0;
    int min_distance = 0; //10 * (graph->kmer_size);  // NOTE: needs to be a cmd_line option

    char cwd[1024];

    if (getcwd(cwd, sizeof(cwd)) != NULL){
        // do NOTHING
    }
    else{
        log_and_screen_printf("CWD command returned NULL\n");
    }

    char*  graph_wd = NULL; 

    Path *simple_path = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_fwd = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_rev = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);

    GraphInfo* nodes_in_graph = calloc(1,sizeof(GraphInfo));
    new_GraphInfo(nodes_in_graph);

    // array to bin coverage 0-5, 5-10, 10-15..95-100
    long int Coverage_Dist[COVERAGE_BINS*COVERAGE_BIN_SIZE]; // will this work?
    char fastg_filename[256];
    char gfa_filename[256];
    char analysis_filename[256];
    char degrees_filename[256];

    Queue* graph_queue;
    for(i=0;i<MAX_BRANCHES;i++){
        Contig_Branches[i]=0;
    }
    // Initialise Coverage_Dist  int i;
    for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE);i++){
        Coverage_Dist[i]=0;
    }


    /* Open contigs file */
    fp_contigs_fasta = fopen(consensus_contigs_filename, "w");
    if (!fp_contigs_fasta) {
        log_and_screen_printf("ERROR: Can't open contig file.\n%s\n", consensus_contigs_filename);
        exit(-1);
    }

    remove_file_extension(consensus_contigs_filename);

    /* Open the analysis file */
    sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
    fp_analysis = fopen(analysis_filename, "w");
    if (!fp_analysis) {
        log_and_screen_printf("ERROR: Can't open analysis file.\n");
        exit(-1);
    }


    /* Open fastg contigs file */
    if(gfa_fastg_output)
    {
        sprintf(fastg_filename, "%s.fastg", consensus_contigs_filename);
        fp_contigs_fastg = fopen(fastg_filename, "w");
        if (!fp_contigs_fastg) {
            log_and_screen_printf("ERROR: Can't open contig (fastg) file.\n%s\n", fastg_filename);
            exit(-1);
        }
        // write the header
        fprintf(fp_contigs_fastg, "#FASTG:begin;");
        fprintf(fp_contigs_fastg, "\n#FASTG:version=1.0:assembly_name=\"%s\";", consensus_contigs_filename);

        /* Open gfa contigs file */
        sprintf(gfa_filename, "%s.gfa", consensus_contigs_filename);
        fp_contigs_gfa = fopen(gfa_filename, "w");
        if (!fp_contigs_gfa) {
            log_and_screen_printf("ERROR: Can't open contig (gfa) file.\n%s\n", gfa_filename);
            exit(-1);
        }
        
        // write the header
        fprintf(fp_contigs_gfa, "H\n");

    }

    /* Open the sugraph degree file */
    sprintf(degrees_filename, "%s.degrees", consensus_contigs_filename);
    fp_degrees = fopen(degrees_filename, "w");
    if (!fp_degrees) {
        log_and_screen_printf("ERROR: Can't open degrees file.\n");
        exit(-1);
    }

    // check for graphs dir existance
    if (basename(consensus_contigs_filename)==consensus_contigs_filename){
        log_and_screen_printf("(Relative path for contig output given, prefixing CWD)\n");
        // returns '.' which breaks other paths later on
        int size = 8 + strlen(cwd);
        graph_wd = calloc(size, sizeof(char));
        sprintf(graph_wd, "%s/graphs/", cwd);
    }
    else{
        // dirname modifies 'consensus_contigs_filename' on some platforms, shifted in here to avoid that
        int size = 8 + strlen(dirname(consensus_contigs_filename));
        graph_wd = calloc(size, sizeof(char));
        
        sprintf(graph_wd, "%s/graphs/", dirname(consensus_contigs_filename));
        sprintf(analysis_filename, "%s%s.tex", graph_wd, basename(consensus_contigs_filename));
    }

    mkdir(graph_wd, 777);

    log_and_screen_printf("graphs dir\t%s\n", graph_wd);
    log_and_screen_printf("graphs\t%s\n", analysis_filename);
    
    free(graph_wd);

    /* Open the DIGEST file */
    fp_report = fopen(analysis_filename, "w");
    if (!fp_report) {
        log_and_screen_printf("ERROR: Can't open analysis (DIGEST) file.\n\t%s\n", analysis_filename);
        exit(-1);
    }

    // header line for degrees file
    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            fprintf(fp_degrees,"for[%d]rev[%d]\t", i, j);
        }
    }
    fprintf(fp_degrees,"total\n");

    db_graph_reset_flags(graph);

    // called by stats_traversal(), walks paths from a branch
    void find_path_length_with_first_edge_all_colours(Nucleotide n, dBNode * node, int * path_length, Orientation orientation) {
        if (db_node_edge_exist_any_colour(node, n, orientation)) {
            pathStep first_step;
            Path * new_path;
            new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
            first_step.node = node;
            first_step.orientation = orientation;
            first_step.label = n;
            first_step.flags = 0;

            db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);
            new_path->length=1;
            * path_length += new_path->length;
            path_destroy(new_path);
            
            // is this essentially:
            //  if (db_node_edge_exist_any_colour(node, n, orientation)) {
            //      * path_length += 1;
            //  }
            // ???
        }
    }

    // Hash table iterator to label nodes
    void stats_traversal(dBNode * node) {
        //if (!db_node_check_flag_visited(node)) {
        uint32_t this_coverage = element_get_coverage_all_colours(node) - 1;
        int edges_forward= db_node_edges_count_all_colours(node, forward);
        int edges_reverse = db_node_edges_count_all_colours(node, reverse);
        int all_edges = edges_forward + edges_reverse;
        int local_distance = 1;
        int orientation;
        if (this_coverage<0) {
            log_and_screen_printf("Error: Coverage is <1 in the graph?\n");
            exit(-1);
        }


        // GRAPH DENSITY ESTIMATES
        // hash_table_traverse if edges_forward+edges reverse>2
        //   if(db_node_edge_exist_any_colour)(node, n, orientation)
        //   perfect path each edge
        //   count length
        //   average length > kmer?

        Contig_Branches[all_edges-1]++;

        if ((all_edges>2) && (all_edges<=max_node_edges)){

            if (linked_list_max_size){
                add_item(node, this_coverage, linked_list_max_size);
            }

            //log_and_screen_printf("\nWalking branch node...\n");
            // Look at all paths out from here

            orientation = forward;
            int i;
            for (i = 0; i < 4; i++) {
                find_path_length_with_first_edge_all_colours(i, node, &local_distance, orientation);
            }

            orientation = reverse;
            // nucleotide_iterator(&walk_if_exists);
            for (i = 0; i < 4; i++) {
                find_path_length_with_first_edge_all_colours(i, node, &local_distance, orientation);
            }

            local_distance = local_distance / all_edges;
        }
        else{
            // could check
            local_distance = min_distance;
        }

        // PARTITIONING - REMOVE EXTREMELY BRANCHED NODES
        if ((all_edges<=max_node_edges) && (local_distance>=min_distance))
        {
            if(this_coverage>COVERAGE_BINS*COVERAGE_BIN_SIZE-1)
            {
                this_coverage = COVERAGE_BINS*COVERAGE_BIN_SIZE-1;
            }

            Coverage_Dist[this_coverage]++;
            total_nodes++;

            // Look for Y shape branch forward orientation
            // The nodes at the top of the Y should contain different colours
            if (edges_forward > 1
                && edges_reverse == 1) {
                db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
            }
            // Look for Y shape branch reverse orientation
            if (edges_reverse > 1
                && edges_forward == 1) {
                db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
            }
            // Look for X-shaped branch
            if (edges_reverse > 1
                && edges_forward > 1) {
                db_node_action_set_flag(node, X_NODE);
            }
        }
        else{
            log_and_screen_printf("\nPruning node:\tedges\t%i\tdistance\t%i\n", all_edges, local_distance);
            cleaning_prune_db_node(node, graph);
        }
    } // stats_traversal()

    graph_queue = node_queue_new(METACORTEX_QUEUE_SIZE);
    if (!graph_queue) {
        log_and_screen_printf("Couldn't get memory for graph queue.\n");
        exit(-1);
    }
    /* Initialise temporaray path array buffers */
    path_array_initialise_buffers(graph->kmer_size);

    // Hash table iterator to initially walk all
    void explore_graph_size(dBNode * node) {
        if(db_node_check_for_any_flag(node, PRUNED | VISITED) == false){

            initialise_GraphInfo(nodes_in_graph);
            // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
            log_printf("\nGrowing graph from node");
            graph_queue->number_of_items = 0;

            timestamp_gs();
            log_printf("\n");

            // now with a subgraph, walk the graph counting degrees by graph
            explore_subgraphs(node, graph, nodes_in_graph);

            if (nodes_in_graph->total_size ==1) {
                // ignore; pruned node
            }
            else if (nodes_in_graph->total_size) {
                // print out the size of the current subgraph
                log_printf("graph size\t%i\n",nodes_in_graph->total_size);
                fprintf(fp_analysis, "%i\t%i\t",nodes_in_graph->branch_nodes,nodes_in_graph->total_size);
                binary_kmer_to_seq(&nodes_in_graph->highest_cov_in_subgraph, graph->kmer_size, seq);
                fprintf(fp_analysis, "%s\n", seq);

                // update graph wide stats
                print_degree_stats(nodes_in_graph, fp_degrees);
                if (nodes_in_graph->total_size>nodes_in_graph->largest_subgraph) {
                    nodes_in_graph->largest_subgraph=nodes_in_graph->total_size;
                }
                nodes_in_graph->branch_nodes_total=nodes_in_graph->branch_nodes_total+nodes_in_graph->branch_nodes;
                nodes_in_graph->num_subgraphs++;
                i=log10(nodes_in_graph->total_size);
                if(i>=GRAPH_LOG10_LIMIT){
                    i=GRAPH_LOG10_LIMIT-1;
                }
                nodes_in_graph->subgraph_dist[i]++;
                if(nodes_in_graph->total_size>min_subgraph_kmers){
                    nodes_in_graph->num_subgraphs_2k++;
                }
            }
            if (nodes_in_graph->branch_nodes>(MAX_BRANCHES-1)){
                nodes_in_graph->branch_nodes=MAX_BRANCHES-1;
            }
        }
    } // end of &explore_graph_size

    // Hash table iterator to walk graphs, produce paths
    void traversal_for_contigs(dBNode * node) {
        if(db_node_check_for_any_flag(node, PRUNED | VISITED) == false){

            dBNode* seed_node;
            initialise_GraphInfo(nodes_in_graph);

            // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
            log_printf("\nGrowing graph from node");
            graph_queue->number_of_items = 0;

            timestamp_gs();
            log_printf("\n");

            // now with a subgraph, walk the graph counting degrees by graph
            // - this sets VISITED flag as true for many of the nodes in the graph.
            grow_graph_from_node_stats(node, &seed_node, graph, graph_queue, nodes_in_graph, delta_coverage);

            if (nodes_in_graph->total_size ==1) {
                // ignore; pruned node
                cleaning_prune_db_node(node, graph);
                db_node_action_set_flag(node, VISITED);
                log_printf("\t[Singleton pruned.]\n");
            }
            else if (seed_node == NULL) {
                printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph->total_size);
            } else if (nodes_in_graph->total_size) {
                int kmer_size = graph->kmer_size;
                char kmer_string[kmer_size + 1];
                binary_kmer_to_seq(&seed_node->kmer, kmer_size, kmer_string);
                log_printf("Seed node: %s\n", kmer_string);
                /* enough nodes to bother with? If so, get consensus contig */
                if (walk_paths && (nodes_in_graph->total_size >= min_subgraph_kmers)) {

                    // should be a perfect path? might be two paths though, if we started in the middle
                    // NOTE: unecessary coverage element but repeating the whole path finding without coverage
                    //  is more work than necessary I think. See what processing time it changes?
                    coverage_walk_get_path(seed_node, forward, NULL, graph, path_fwd, true);
                    coverage_walk_get_path(seed_node, reverse, NULL, graph, path_rev, true);

                    path_reverse(path_fwd, simple_path);
                    path_append(simple_path, path_rev);

                    simple_path->id = counter;
                    if (simple_path->length > (min_contig_size - graph->kmer_size)) {
                        log_printf("Write path of size %d\n", simple_path->length);
                        log_printf("graph size\t%i\n",nodes_in_graph->total_size);

                        // could save the path walking again here if needed, hold these figures in path structure
                        double average_coverage=0;
                        uint32_t min_coverage=0;
                        uint32_t max_coverage=0;
                        path_get_statistics(&average_coverage, &min_coverage, &max_coverage, simple_path);
                        // NOTE: decision - minimum cov or average cov dictates confidence threshold met?
                        // Output for alternative formats
                        path_to_fasta(simple_path, fp_contigs_fasta);
                                               
                        if(gfa_fastg_output)
                        {
                            path_to_gfa2_and_fastg(simple_path,graph,fp_contigs_gfa, fp_contigs_fastg);
                        }
                        counter++;
                    } else {
                        log_printf("Didn't write path of size %d\n", simple_path->length);
                    }

            
                    /*	HERE - INSTEAD OF 'VISITING' AN X-NODE, REMOVE EDGES FROM BEFORE
                     AND AFTER IT ON CURRENT PATH */
                    
                    //unset VISITED flags from grow_graph_from_node_stats
                    dBNode* queue_node;
                    while (graph_queue->number_of_items > 0) {
                        queue_node = (dBNode*)queue_pop(graph_queue);
                        db_node_action_unset_flag(queue_node, VISITED);
                    }
                    

                    // Now disconnect path from other nodes and mark path as visited, so it's not visited again //
                    for (i=0; i<simple_path->length; i++) {
                        if (simple_path->step_flags[i] & X_NODE)
                        {
                            log_printf("[x-NODE IN PATH RETAINED]\n");
                            // leave branching nodes in path, cleaning_prune_db_node
                            //   will remove reciprocal edges to adjacent nodes in path
                        }
                        else // for non-branching nodes (including Y-branches)
                        {
                            cleaning_prune_db_node(simple_path->nodes[i], graph);
                            db_node_action_set_flag(simple_path->nodes[i], VISITED);
                        }

                    }


                    /* Reset paths */
                    path_reset(simple_path);
                    path_reset(path_fwd);
                    path_reset(path_rev);

                } else  {
                    log_printf("  Number of nodes (%i) too small. Not outputting contig.\n", nodes_in_graph->total_size);
                } // end size graph size check


            } else {
                // catch graph size of zero? Not sure why this happens - grow-graph must be failing
                log_printf("graph size of zero?\n");
            }
            if (nodes_in_graph->branch_nodes>(MAX_BRANCHES-1)){
                nodes_in_graph->branch_nodes=MAX_BRANCHES-1;
            }
        }
    } // traversal_for_contigs

    // check each node in the graph, FLAG X&Y nodes (mark all nodes as visited)
    timestamp_gs();
    log_and_screen_printf("Stats traversal started...");
    hash_table_traverse(&stats_traversal, graph);
    log_and_screen_printf("DONE\n");

    timestamp_gs();
    // first line for stats output file
    fprintf(fp_analysis, "\n#Subgraph sizes\n");
    log_and_screen_printf("Graph size traversal started...");
    // first travesal - build subgraphs out, produce stats
    hash_table_traverse(&explore_graph_size, graph);
    log_and_screen_printf("DONE\n");

    // clear linked list
    log_and_screen_printf("Unique kmers before clearing:\t %lld\n", graph->unique_kmers);
    log_and_screen_printf("linked_list_max_size:\t %i\n", linked_list_max_size);
    if (linked_list_max_size){
        clear_list(graph);
    }
    log_and_screen_printf("Unique kmers after clearing:\t %lld\n", graph->unique_kmers);

    timestamp_gs();
    db_graph_reset_flags(graph);
    // second travesal - build subgraphs out, produce contigs
    log_and_screen_printf("Full traversal started...");
    hash_table_traverse(&traversal_for_contigs, graph);
    log_and_screen_printf("DONE\n");
    
    // write the fastg footer
    if(gfa_fastg_output)
    {
        fprintf(fp_contigs_fastg, "\n#FASTG:end;");
        fclose(fp_contigs_fastg);
        fclose(fp_contigs_gfa);
    }
    
    fclose(fp_analysis);
    fclose(fp_degrees);

    // Output graph wide stats (coverage)



    // run R script to produce figures for report
    // will this work? initialising 'cmd' like this?

    //char cmd = printf("Rscript %s %s", <path_to_src>/degree_plots.R, degrees_filename);
    //system(cmd);  // potential problems with this apparently? is permissions are an initialiseAlignmentSummaryFile

    char command[1024];
    //char r_script_path[]="/home/aylingm/grimoire/metacortex/";
    char * r_script_path=getenv("R_ENV_PATH");

    if (r_script_path==NULL){
        log_and_screen_printf("\nR_ENV_PATH not set, skipping graphs step...\n\n");
    }
    else{
        printf("\nPATH : %s\n",r_script_path);

        if (cwd != NULL){
            sprintf(command, "Rscript %sdegree_plots.R %s/%s", r_script_path, cwd, degrees_filename);
            log_and_screen_printf("\n%s\n", command);
            int systemRet = system(command);
            if(systemRet == -1){
                // The system method failed
                log_and_screen_printf("Failed call to system?\n");
            }
            log_and_screen_printf("\n");
        }
        else{
            log_and_screen_printf("CWD command returned NULL\n");
        }
    }

    writeLaTeXHeader(fp_report, consensus_contigs_filename);
    writeLaTeXreport(fp_report, (int) MAX_BRANCHES, (int) COVERAGE_BIN_SIZE, \
                     (int) COVERAGE_BINS, (int) GRAPH_LOG10_LIMIT, (int) NUM_BEST_NODES, &Contig_Branches[0], \
                     &Coverage_Dist[0], graph, nodes_in_graph);
    fclose(fp_report);
    writeLaTeXreport_to_log_and_screen((int) MAX_BRANCHES, (int) COVERAGE_BIN_SIZE, \
                                       (int) COVERAGE_BINS, (int) GRAPH_LOG10_LIMIT, (int) NUM_BEST_NODES, &Contig_Branches[0], \
                                       &Coverage_Dist[0], graph, nodes_in_graph);

    // memory issue - analysis_filename is being stomped on at some point
    log_and_screen_printf("\nanalysis filename\t%s\n", analysis_filename);
    sprintf(command, "pdflatex -interaction=nonstopmode %s", analysis_filename);
    log_and_screen_printf("\n%s\n", command);
    int systemRet = system(command);
    if(systemRet == -1){
        // The system method failed
        log_and_screen_printf("Failed call to system?\n");
    }

    db_graph_reset_flags(graph);
}

void print_degree_stats(GraphInfo * info, FILE* fp_degrees){
    int i;  int j;
    int total_nodes=info->total_size;

    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            fprintf(fp_degrees,"%f\t",  ((float) info->node_degree[i][j]) / (float) total_nodes);
        }
    }
    fprintf(fp_degrees,"%d\n", total_nodes);
}

void initialise_GraphInfo(GraphInfo * info){
    int i, j;
    info->total_size = 0;
    info->branch_nodes = 0;
    info->end_nodes = 0;
    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            info->node_degree[i][j]=0;
        }
    }
    info->highest_cov = 0;
    for (i=0; i<NUM_BEST_NODES; i++) {
        info->best_coverage[i]=0;
    }
}

void new_GraphInfo(GraphInfo * info){
    int i;
    info->largest_subgraph = 0;
    info->num_subgraphs = 0;
    info->num_subgraphs_2k = 0;
    info->simple_bubbles = 0;
    info->branch_nodes_total=0;
    for(i=0;i<GRAPH_LOG10_LIMIT;i++){
        info->subgraph_dist[i]=0;
    }
}

int explore_subgraphs(dBNode* start_node, dBGraph* graph, GraphInfo* nodes_in_graph){
    Queue* nodes_to_walk;
    dBNode* node;
    int orientation;
    int depth;
    int best_edges[NUM_BEST_NODES];
    int i;
    for(i=0; i<NUM_BEST_NODES; i++){
        best_edges[i]=0;
    }

    // Nucleotide iterator, used to walk all possible paths from a node
    void walk_if_exists(Nucleotide n) {
        // If there is an edge in any colour for this nucleotide...
        if (db_node_edge_exist_any_colour(node, n, orientation)) {

            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode * next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, graph);
            if (!next_node) {
                log_and_screen_printf("Error: Something went wrong with db_graph_get_next_node\n");
                exit(-1);
            }

            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited(next_node)) {
                pathStep first_step;
                Path * new_path;
                dBNode* end_node;
                i = 0;

                // Get path
                first_step.node = node;
                first_step.orientation = orientation;
                first_step.label = n;
                first_step.flags = 0;
                new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
                if (!new_path) {
                    log_and_screen_printf("ERROR: Not enough memory to allocate new path.\n");
                    exit(-1);
                }

                db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);

                double path_coverage=0;
                uint32_t min_coverage=0; uint32_t max_coverage=0; // required for path_get_statistics()
                path_get_statistics(&path_coverage, &min_coverage, &max_coverage, new_path);

                if (1)  {
                    // Add end node to list of nodes to visit
                    end_node = new_path->nodes[new_path->length-1];
                    if (!db_node_check_flag_visited(end_node)) {
                        if (!db_node_is_blunt_end_all_colours(end_node, new_path->orientations[new_path->length-1])) {
                            if (queue_push_node(nodes_to_walk, end_node, depth+1) == NULL) {
                                log_and_screen_printf("Queue too large. Ending. (WALK)\n");
                                exit(1);
                            }
                        }
                    }

                    // Now go through all nodes, look for best and mark all as visited
                    for (i=0; i<new_path->length; i++) {
                        if (!db_node_check_flag_visited(new_path->nodes[i])) {
                            uint32_t this_coverage = element_get_coverage_all_colours(new_path->nodes[i]);
                            int this_FOR_edges = db_node_edges_count_all_colours(new_path->nodes[i], forward);
                            int this_REV_edges = db_node_edges_count_all_colours(new_path->nodes[i], reverse);

                            // add node degrees to 2D array of all degrees in subgraph
                            nodes_in_graph->node_degree[this_FOR_edges][this_REV_edges]++;

                            if (this_coverage>nodes_in_graph->highest_cov){
                                nodes_in_graph->highest_cov=this_coverage;
                                binary_kmer_assignment_operator(nodes_in_graph->current_kmer,new_path->nodes[i]->kmer);
                                binary_kmer_assignment_operator(nodes_in_graph->highest_cov_in_subgraph,nodes_in_graph->current_kmer);
                            }

                            // if this is better than the lowest 'good' node (top five coverage)
                            if ((this_coverage > nodes_in_graph->best_coverage[NUM_BEST_NODES-1]) ||
                                ((this_coverage == nodes_in_graph->best_coverage[NUM_BEST_NODES-1]) && ((this_FOR_edges + this_REV_edges) > best_edges[NUM_BEST_NODES-1])))
                            {
                                // sort algorithm - because I sort as I build array, no need to make more than one pass
                                int temp_cov=0;
                                binary_kmer_initialise_to_zero(&(nodes_in_graph->temp_kmer));
                                // yes, this is the same as above.
                                binary_kmer_assignment_operator(nodes_in_graph->current_kmer,new_path->nodes[i]->kmer);
                                // seed_node->kmer

                                int j=0;
                                while(this_coverage){
                                    if(j>=NUM_BEST_NODES){
                                        this_coverage=0;
                                    }
                                    else if (this_coverage>nodes_in_graph->best_coverage[j]||
                                             ((this_coverage == nodes_in_graph->best_coverage[j]) && ((this_FOR_edges + this_REV_edges) > best_edges[j]))){
                                        temp_cov = nodes_in_graph->best_coverage[j];
                                        nodes_in_graph->best_coverage[j] = this_coverage;
                                        this_coverage=temp_cov;

                                        // recycle temp_cov for one line
                                        temp_cov=best_edges[j];
                                        best_edges[j] = (this_FOR_edges + this_REV_edges);
                                        // set the current edge count for the rest of the loop
                                        this_FOR_edges=temp_cov;
                                        this_REV_edges=0;
                                        temp_cov=0;

                                        binary_kmer_assignment_operator(nodes_in_graph->temp_kmer,nodes_in_graph->kmer[j]);
                                        binary_kmer_assignment_operator(nodes_in_graph->kmer[j],nodes_in_graph->current_kmer);
                                        binary_kmer_assignment_operator(nodes_in_graph->current_kmer,nodes_in_graph->temp_kmer);
                                        j++;
                                    }
                                    else{
                                        j++;
                                    }
                                }
                            }

                            if (db_node_check_for_any_flag(new_path->nodes[i], BRANCH_NODE_FORWARD)){
                                nodes_in_graph->branch_nodes++;
                            }
                            else if (db_node_check_for_any_flag(new_path->nodes[i], BRANCH_NODE_REVERSE)){
                                nodes_in_graph->branch_nodes++;
                            }
                            else if (db_node_check_for_any_flag(new_path->nodes[i], X_NODE)){
                                nodes_in_graph->branch_nodes++;
                            }

                            db_node_action_set_flag_visited(new_path->nodes[i]);

                            nodes_in_graph->total_size++;
                        }
                    } // path->length loop
                }
                else{
                    cleaning_prune_db_node(new_path->nodes[0], graph);
                    // NOTE best_node - needs to be checked here. Don't want to return NULL
                }

                // Clean up
                path_destroy(new_path);
            }
        }
    } // walk_if_exists

    // Start a queue of nodes to walk
    //log_and_screen_printf("Allocating %d Mb to store queue information (max %d nodes, when full each node could be %d)...\n", ((METACORTEX_QUEUE_SIZE * sizeof(QueueItem*)) / 1024) / 1024, METACORTEX_QUEUE_SIZE, sizeof(QueueItem));
    nodes_to_walk = node_queue_new(METACORTEX_QUEUE_SIZE);
    if (!nodes_to_walk) {
        log_and_screen_printf("Couldn't get memory for node queue.\n");
        exit(-1);
    }

    // Add start node to list of nodes to visit
    if (queue_push_node(nodes_to_walk, start_node, 0) == NULL) {
        log_and_screen_printf("Queue too large. Ending.\n");
        exit(-1);
    }


    if (!db_node_check_flag_visited(start_node)) {
        db_node_action_set_flag_visited(start_node);
        nodes_in_graph->total_size++;
    }

    // Now keep visiting nodes and walking paths
    while (nodes_to_walk->number_of_items > 0) {
        // Take top node from list
        node = queue_pop_node(nodes_to_walk, &depth);

        // Look at all paths out from here
        orientation = forward;
        nucleotide_iterator(&walk_if_exists);
        orientation = reverse;
        nucleotide_iterator(&walk_if_exists);
    }

    queue_free(nodes_to_walk);

    return 0;
}
