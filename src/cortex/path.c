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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <seq.h>
#include <binary_kmer.h>
#include <element.h>
#include <path.h>
#include <math.h>
#include <logger.h>
#include "dB_graph.h"
#include "gfa_segment.h"

int fastg_recursion_level;

Path *path_new(int max_length, short kmer_size)
{
    Path *path = calloc(1, sizeof(Path));
    if(path == NULL){
        fprintf(stderr, "[path_new]Unable to allocate path\n");
        exit(-1);
    }
    path->nodes = calloc(max_length, sizeof(dBNode *));
    path->orientations = calloc(max_length, sizeof(Orientation));
    path->labels = calloc(max_length + 1, sizeof(Nucleotide));
    path->seq = calloc(max_length + 1 + kmer_size, sizeof(char));
    path->in_nodes = calloc(PATH_MAX_IN_NODES, sizeof(int));
    path->step_flags = calloc(max_length, sizeof(Flags));
    path->max_virtual_length = max_length;
    if(path->nodes == NULL){
        fprintf(stderr, "[path_new]Unable to allocate nodes\n");
        exit(-1);
    }

    if(path->orientations == NULL){
        fprintf(stderr, "[path_new]Unable to allocate orientations\n");
        exit(-1);
    }

    if(path->labels == NULL){
        fprintf(stderr, "[path_new]Unable to allocate labels\n");
        exit(-1);
    }

    if(path->seq == NULL){
        fprintf(stderr, "[path_new]Unable to allocate seq\n");
        exit(-1);
    }

    if(path->in_nodes == NULL){
        fprintf(stderr, "[path_new]Unable to allocate in_nodes\n");
        exit(-1);
    }

    path->in_nodes_capacity = PATH_MAX_IN_NODES;
    path->in_nodes_count = 0;
    path->max_length = max_length;
    flags_action_unset_flag(IS_CYCLE, &(path->flags));
    flags_action_set_flag(PRINT_FIRST, &(path->flags));
    path->kmer_size = kmer_size;
    //	path->depth = -1;
    path->header = 0;
    path->used = false;
    path->subpath_id = 0;

    if ((!path->nodes) || (!path->orientations) || (!path->labels) || (!path->seq)) {
        free(path);
        path = 0;
    }

    return path;
}

void path_destroy(Path * path)
{
   //if (path == NULL) {
    //    fprintf(stderr, "[path_destroy] The path is null");
    //    exit(-1);
    //}
    free(path->nodes);
    free(path->orientations);
    free(path->labels);
    free(path->seq);
    free(path->in_nodes);
    free(path->step_flags);

    if (path->header) {
        free(path->header);
    }

    free(path);

}

void path_increase_id(Path * path)
{
    path->id++;
}

//We can only set the limit to new paths!
void path_set_limit(int limit, Path * p){
    assert(limit <= p->max_length);
    assert(p->length == 0);
    p->max_virtual_length = limit;
}

int path_get_limit(Path * p){
    return p->max_virtual_length;
}

void path_reset(Path * path)
{
    if (path == NULL) {
        fprintf(stderr, "[path_reset] The path is null");
        exit(-1);
    }
    while (path->length > 0) {
        path_remove_last(path);
    }

    path->seq[0] = '\0';
    path->labels[0] = Undefined;
    path->in_nodes_count = 0;
    path->out_nodes_count = 0;
    path->length = 0;
    path->new_nodes = 0;
    //  path->depth = -1;
    path->max_virtual_length = path->max_length;
    flags_action_clear_flags(&(path->flags));
    flags_action_set_flag(PRINT_FIRST, &(path->flags));
    flags_action_set_flag(NEW_PATH, &(path->flags));
    path_clean_stop_reason(path);

    int i;

    for(i = 0; i < path->in_nodes_capacity; i++){
        path->in_nodes[i] = 0;
    }

}

/**
 * To get a natural order of two steps.
 * At the moment we are using just the address to compare, is
 * cheaper than comparing the kmers. If the nodes are the same,
 * we use the orientations to distinguish them, and if Use with care.
 */
int path_step_compare(pathStep * a, pathStep * b){
    if (a->node < b->node) {
        return -1;
    }else if (a->node > b->node) {
        return 1;
    }else if (a->orientation < b->orientation) {
        return -1;
    }else if (a->orientation > b->orientation) {
        return 1;
    }else if (a->label < b->label) {
        return -1;
    }else if (a->label > b->label) {
        return 1;
    }else {
        return 0;
    }

}

boolean path_step_has_unvisited_edge_all_colours(pathStep * step){
    Edges e = db_node_get_edges_for_orientation_all_colours(step->node, step->orientation);
    Edges visited = step->flags & PATH_STEP_ALL_VISITED; //We relay on being using the 4 least significan bit.
    return e != visited;
}



Nucleotide path_step_get_unvisited_edge_all_colours(pathStep * step){
    Edges e = db_node_get_edges_for_orientation_all_colours(step->node, step->orientation);
    Edges visited = step->flags & PATH_STEP_ALL_VISITED; //We relay on being using the 4 least significan bit.
    Nucleotide n, nucleotide = Undefined;
    Edges edges = ~visited & e;
    for (n = 0; n < 4 && nucleotide == Undefined; n++) {
        if ((edges & 1) == 1) {
            nucleotide = n;
        }
        edges >>= 1;
    }

    return nucleotide;
}

boolean path_step_equals(pathStep * step, pathStep * other)
{
    return step->node == other->node
    && step->orientation == other->orientation
    && step->label == other->label;
}

// Compare steps, without comparing labels
boolean path_step_equals_without_label(pathStep * step, pathStep * other)
{
    return step->node == other->node
    && step->orientation == other->orientation;
}

void path_modfy_last_label(Nucleotide n, Path * p){
    int last_index = path_get_length(p)-1;
    char last_c = binary_nucleotide_to_char(n);
    p->labels[last_index] = n;
    p->seq[last_index] = last_c;
    p->step_flags[last_index] |= binary_nucleotide_to_edge(n);

}

boolean unlabelled_path_step_equals(pathStep * step, pathStep * other)
{
    return step->node == other->node
    && step->orientation == other->orientation;
}


boolean path_contains_step(pathStep * step, Path * path)
{
    pathStep tmp;
    path_step_initialise(&tmp);
    int i;
    int len = path_get_length(path);
    for (i = 0; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals(step, &tmp)) {
            return true;
        }
    }
    return false;
}

boolean path_contains_step_without_label(pathStep * step, Path * path)
{
    pathStep tmp;
    path_step_initialise(&tmp);
    int i;
    int len = path_get_length(path);
    for (i = 0; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals_without_label(step, &tmp)) {
            return true;
        }
    }
    return false;
}


boolean path_is_first_equals_to_last(Path * path)
{
    if(path->length < 2)
        return false;

    pathStep first, last;
    path_step_initialise(&first);
    path_step_initialise(&last);
    path_get_last_step(&last, path);
    path_get_step_at_index(0, &first, path);
    return path_step_equals_without_label(&first, &last);
}

boolean path_contains_step_from_index(int index, pathStep * step, Path * path)
{
    pathStep tmp;
    path_step_initialise(&tmp);
    int i;
    int len = path_get_length(path);
    for (i = index; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals(step, &tmp)) {
            return true;
        }
    }
    return false;
}

boolean path_contains_step_from_index_without_label(int index, pathStep * step, Path * path)
{
    pathStep tmp;
    path_step_initialise(&tmp);
    int i;
    int len = path_get_length(path);
    for (i = index; i < len; i++) {
        path_get_step_at_index(i, &tmp, path);
        if (path_step_equals_without_label(step, &tmp)) {
            return true;
        }
    }
    return false;
}

boolean path_has_space(Path * path)
{
    if (path->length >= path->max_virtual_length) {
        return false;
    }
    return true;

}

void path_add_in_step(pathStep * step, Path * path)
{
    assert(path != NULL);
    assert(step != NULL);
    assert(step->node != NULL);
    assert(path->in_nodes != NULL);
    assert(path->in_nodes_count < path->in_nodes_capacity);

    path->in_nodes[path->in_nodes_count++] = path->length;

    if (path->in_nodes_count == path->in_nodes_capacity) {
        int new_capacity = path->in_nodes_capacity + PATH_IN_NODES_CAPACITY_INCREASE;
        int * in_tmp_ptr = realloc(path->in_nodes, new_capacity * sizeof(int));
        int i;
        if (in_tmp_ptr == NULL){
            fprintf(stderr, "out of memory while increasing the number  to %d of in paths  in path: \n",
                    new_capacity);
            path_to_fasta(path, stderr);
            exit(1);
        }

        path->in_nodes = in_tmp_ptr;

        for (i=path->in_nodes_count; i<new_capacity; i++) {
            path->in_nodes[i] = 0;
        }

        path->in_nodes_capacity = new_capacity;
    }

    assert(path->in_nodes_count <= path->in_nodes_capacity);
}

boolean path_has_in_step(pathStep * step, Path * path){
    pathStep tmp_step;
    path_step_initialise(&tmp_step);

    boolean found = false;
    int i;
    for(i = 0; i < path->in_nodes_count && !found; i++ ){
        path_get_step_at_index(path->in_nodes[i], &tmp_step, path);
        if(path_step_equals_without_label(&tmp_step, step)){
            found = true;
        }
   }
   return found;
}

// I don't know what the function path_has_in_step is doing, but I think it should do this.
boolean path_has_in_step2(pathStep * step, Path * path)
{
    pathStep tmp_step;
    path_step_initialise(&tmp_step);
    for(int i = 0; i < path->length; ++i)
    {
       path_get_step_at_index(i, &tmp_step, path);
       if(path_step_equals_without_label(&tmp_step, step))
       {
            return true;
       }
    }
    return false;
}

int path_index_of_last_in_node(Path * path){
    int i = path->in_nodes_count;
    if( i !=  0){
        return path->in_nodes[i-1];
    }
    return -1;

}

int path_get_first_in_node_after(int pos, Path * path){
    assert(pos >= 0);
    int i, arr_pos = -1;
    for (i = path->in_nodes_count-1; i > 0; i--) {
        if(path->in_nodes[i] > pos){
            arr_pos = i;
        }
    }
    return arr_pos != -1? path->in_nodes[arr_pos]: -1;
}


boolean path_add_node(pathStep * step, Path * path)
{
    assert(path!=NULL);
    assert(step != NULL);
    assert(step->node != NULL);

    if (step->node == NULL) {
        fprintf(stderr, "[path_add_node] The node is null\n");
        //exit(-1);
        return false;	//If you try to add a null node, you dont add it. Currently never happens because the assert...
    }
    if (path->length == path->max_virtual_length) {
        //log_and_screen_printf("[path_add_node] The node_%qd has reached the limit[%d]\n",
        //	     path->id, path->max_virtual_length);
        //path_to_fasta(path, <#FILE *fout#>)
        return false;
    }

    dBNode *node = step->node;
    if (DEBUG) {
        char tmp_seq[path->kmer_size + 1];
        tmp_seq[path->kmer_size] = '\0';
        printf
        ("[path_add_node] Node %i in path(%lli): %s add %c (%s).\n",
         path->length, path->id,
         binary_kmer_to_seq(element_get_kmer(node), path->kmer_size,
                            tmp_seq),
         binary_nucleotide_to_char(step->label),
         step->orientation == reverse ? "reverse" : "forward");
    }
    Orientation orientation = step->orientation;
    Nucleotide nucleotide = step->label;
    pathStep first;
    path_step_initialise(&first);
    if(path->length > 0){
        path_get_step_at_index(0, &first, path);
        if(path_step_equals_without_label(&first, step)){ //TODO this is not enough for a cycle as the cycle doesn't need to close on first node (it is also disregarding the orientation)
            flags_action_set_flag(IS_CYCLE, &(path->flags));
        }
        if(path_has_in_step2(step, path)){
            flags_action_set_flag(IS_CYCLE, &(path->flags));
        }


    }
    int edges_in = db_node_edges_count_all_colours(step->node,  opposite_orientation(step->orientation));
    if(edges_in > 1){
        path_add_in_step(step, path);
    }
    if(!db_node_check_flag_visited(step->node)){
        path->new_nodes++;
    }

    int edges_out = db_node_edges_count_all_colours(step->node, step->orientation);
    if(edges_out > 1){
        path->out_nodes_count++;
    }

    /*flags_action_set_flag(FIND_AGAIN, &(path->flags)); */
    flags_action_unset_flag(NEW_PATH | EMPTY_PATH, &(path->flags));

    path->nodes[path->length] = node;
    path->orientations[path->length] = orientation;
    path->labels[path->length] = nucleotide;
    if(nucleotide != Undefined){
        path->step_flags[path->length] = step->flags | binary_nucleotide_to_edge(nucleotide);//We add the current nucleotide to the flags. The step may come with already marked nucleotides.
    }
    path->seq[path->length] = nucleotide == Undefined ? '\0'
    : binary_nucleotide_to_char(nucleotide);

    path->seq[path->length + 1] = '\0';

    path->length++;

    return true;

}

void path_step_print(pathStep * step, int kmer_size, FILE * f)
{
    char tmp_seq[kmer_size + 1];
    tmp_seq[kmer_size] = '\0';
    fprintf(f, "pathStep: %s %c (%s).",
            binary_kmer_to_seq(element_get_kmer(step->node), kmer_size,
                               tmp_seq),
            binary_nucleotide_to_char(step->label),
            step->orientation == reverse ? "reverse" : "forward");
}

Nucleotide path_last_nucleotide(Path * path)
{
    return (path->length <= 0) ? Undefined : path->labels[path->length - 1];
}

Flags path_last_flags(Path * path)
{
    return (path->length <= 0) ? 0 : path->step_flags[path->length - 1];
}

static void compute_label(dBNode * node, Orientation o, char *label)
{
    int i = 0;
    Nucleotide n;
    for (n = Adenine; n < Undefined; n++) {
        if (db_node_edge_exist(node, n, o)) {
            label[i] = binary_nucleotide_to_char(n);
            i++;
        }
    }

    label[i] = '\0';
}

static void check_print_first(Path * path)
{
    boolean ignoreFirst = false;
#ifndef SHORT_FLAGS
    Flags fs = db_node_get_flags(path->nodes[0], PRINT_FORWARD | PRINT_REVERSE);
    // Include the first node if:
    // - first node orientation is forward and PRINT_FORWARD is set
    // - first node orientation is reverse and PRINT_REVERSE is set
    ignoreFirst = (path->orientations[0] == forward) ? ((fs == PRINT_FORWARD) ? false : true) : ((fs == PRINT_REVERSE) ? false : true);
#endif
    if (!ignoreFirst) {
        flags_action_set_flag(PRINT_FIRST, &(path->flags));
    } else {
        flags_action_unset_flag(PRINT_FIRST, &(path->flags));
    }
}


/*
 * Returns a path array with all the perfect paths from the path. If the path is a perfect path, it returns an array containing a copy of the path.
 * WARNING! dont forget to use the free from buffer functions, and destroy the path array
 */
PathArray * path_split_in_perfect_paths(Path * p){
    assert(p != NULL);
    assert(path_get_length(p) > 0);
    PathArray * pa =    path_array_new(1);
    int path_length = path_get_length(p);
    int i, count_fwd, count_rev;
    Path * current_path = path_get_buffer_path();
    path_array_add_path(current_path, pa);
    pathStep tmp;
    path_step_initialise(&tmp);
    i = 0;
    path_get_step_at_index(i, &tmp, p);
    count_rev = db_node_edges_count_all_colours(tmp.node, reverse);
    path_add_node(&tmp, current_path);
    if(count_rev > 1){
        current_path = path_get_buffer_path();
        path_add_node(&tmp, current_path);
        path_array_add_path(current_path, pa);
    }
    for(i = 1; i < path_length; i++){
        path_get_step_at_index(i, &tmp, p);
        count_fwd =db_node_edges_count_all_colours(tmp.node, forward);
        count_rev =db_node_edges_count_all_colours(tmp.node, reverse);
        path_add_node(&tmp, current_path);
        if(count_fwd > 1 || count_rev > 1 ){
            current_path = path_get_buffer_path();
            path_add_node(&tmp, current_path);
            path_array_add_path(current_path, pa);
        }
    }
    return pa;
}
//Merges the first array into the second one, and frees the original array, but it doesnt destroy the paths inside. Careful with the leaks
void path_array_merge(PathArray ** from, PathArray * to){
    PathArray * from_p = *from;
    int paths = path_array_get_number_of_paths(from_p);
    int i;
    for(i = 0; i < paths; i++){
        path_array_add_path(path_array_get(i, from_p), to);
    }
    path_array_destroy_struct(from);


}

// Create a new path consisting of all paths in the array merged into one (if possible).
// Caller is responsible for freeing memory. Does not destroy original path array.
Path* path_array_merge_to_path(PathArray* pa, boolean reverse_array_order, dBGraph* db_graph)
{
    //log_printf("[path_array_merge_to_path] Merging %i paths.\n", pa->number_of_paths);
    int path_length = 0;
    for(int i = 0; i < pa->number_of_paths; i++)
    {
        path_length += pa->paths[i]->length;
    }
    Path* merged_path = path_new(path_length, pa->kmer_size);
    
    for(int i = 0; i < pa->number_of_paths; i++)
    {
        int index = reverse_array_order ? pa->number_of_paths - 1 -i : i;
        Path* current_path = pa->paths[index];
                
        for(int j = 0; j < current_path->length; j++)
        {
            pathStep path_step;
            path_step.node = current_path->nodes[j];
            path_step.flags = current_path->step_flags[j];
            path_step.label = current_path->labels[j];
            path_step.orientation = current_path->orientations[j];
            
            if(j == current_path->length - 1 && i < pa->number_of_paths - 1)
            {
                assert(path_step.label == Undefined);             
                int next_index = reverse_array_order ? index - 1 : index + 1;
                dBNode* next_node = pa->paths[next_index]->nodes[0];
                Orientation next_orientation = pa->paths[next_index]->orientations[0];
/*
                char end_kmer[db_graph->kmer_size + 1];
                binary_kmer_to_seq(element_get_kmer(path_step.node), db_graph->kmer_size, end_kmer);
                end_kmer[db_graph->kmer_size] = '\0';
                        
                char start_kmer[db_graph->kmer_size + 1];
                binary_kmer_to_seq(element_get_kmer(next_node), db_graph->kmer_size, start_kmer);
                start_kmer[db_graph->kmer_size] = '\0';
                log_printf("[path_array_merge_to_path] End node %s.\n", end_kmer);
                log_printf("[path_array_merge_to_path] Start node %s.\n", start_kmer);
*/
                if(next_node == path_step.node && next_orientation == path_step.orientation)
                {
                    // continue without adding path
                    continue;
                }
                                     
                Nucleotide rev_edge;
                Orientation orientation;
                for(int nucleotide = 0; nucleotide < 4; nucleotide++)
                {
                    if(db_graph_get_next_node(  path_step.node, path_step.orientation, 
                                                &orientation, nucleotide, &rev_edge, db_graph) == next_node)
                    {
                        path_step.label = nucleotide;
                        break;
                    }
                }
                if(path_step.label == Undefined)
                {
                    log_printf("[path_array_merge_to_path] Warning: Could not merge paths, are these really consecutive? Ignoring...\n");
                    path_destroy(merged_path);
                    return NULL;
                }
            }
            
            if(merged_path->length >0)
            {
                if(path_step.node == merged_path->nodes[merged_path->length - 1] && 
                   path_step.orientation == merged_path->orientations[merged_path->length - 1])
                {
                    continue;
                }
            }
            path_add_node(&path_step, merged_path);
        }
    }
    
    return merged_path;
}

Path* path_array_get_last_path(PathArray* pa)
{
    return pa->paths[pa->number_of_paths - 1];
}

void path_array_remove_last_path(PathArray* pa)
{
    Path* last_path = pa->paths[pa->number_of_paths - 1];
    if(last_path != NULL)
    {
        path_destroy(last_path);
        pa->paths[pa->number_of_paths - 1] = NULL;
    }
    else
    {
        printf("[path_array_remove_last_path] Warning: Attempted to remove last path from path_array"
                " but path was NULL!\n");
    }
    pa->number_of_paths--;
}

int path_array_get_total_size(PathArray* pa)
{
    int total_size = 0;
    for(int i = 0; i < pa->number_of_paths; i++)
    {
        total_size += pa->paths[i]->length;
    }
    return total_size;
}

void path_to_fasta_debug(Path * path, FILE * fout)
{
    short kmer_size = path->kmer_size;
    //int length = path->length;

    // Sanity checking
    if (path == NULL) {
        printf("[path_to_fasta] trying to print a null Path\n");
        exit(-1);
    }

    if (fout == NULL) {
        printf("[path_to_fasta] trying to print to a null FILE\n");
        exit(-1);
    }

    if (DEBUG) {
        printf("[path_to_fasta] About to print a path\n");
    }
    // Set value of PRINT_FIRST (whether first node included)
    //check_print_first(path);

    // Get coverage statistics from path
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

    // Get orientation of first and last node
    Orientation fst_orientation;
    //if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
    fst_orientation = path->orientations[0];
    //} else {
    //	fst_orientation = path->orientations[1];
    //}
    Orientation lst_orientation = path->orientations[path->length];

    // Get the first node - this will be nodes[0] if PRINT_FIRST is
    // specified, or nodes[1] otherwise.
    dBNode *fst_node;
    if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
        if (path->length == 0) {
            printf("[path_to_fasta] Trying to print an empty path[1]!\n");
            return;
        }
        fst_node = path->nodes[0];
    } else {
        if (path->length < 2) {
            fprintf(stderr, "[path_to_fasta] Trying to print an empty path[2]!\n");
            return;
        }
        fst_node = path->nodes[1];
    }

    // Get the last node
    dBNode *lst_node = path->nodes[path->length - 1];

    // Make a set of labels for first and last nodes which list the
    // acceptable forward and reverse path labels
    char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
    compute_label(fst_node, forward, fst_f);
    compute_label(fst_node, reverse, fst_r);
    compute_label(lst_node, forward, lst_f);
    compute_label(lst_node, reverse, lst_r);

    // Places to store first and last kmer sequence
    char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
    fst_seq[kmer_size] = '\0';
    lst_seq[kmer_size] = '\0';
    BinaryKmer tmp_kmer;

    // Get the first kmer sequence (or complement)
    BinaryKmer fst_kmer;
    binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
    if (fst_orientation == reverse) {
        binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
    if (DEBUG) {
        printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }

    // Get the last kmer sequence (or complement)
    BinaryKmer lst_kmer;
    binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
    if (lst_orientation == reverse) {
        binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);

    // Output to file
    /*fprintf(fout,
     ">node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s\n",
     path->id,
     (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length +
     kmer_size : length + kmer_size - 1), avg_coverage,
     min_coverage, max_coverage,
     element_get_coverage_all_colours(fst_node), fst_seq,
     (fst_orientation == forward ? fst_r : fst_f),
     (fst_orientation == forward ? fst_f : fst_r),
     element_get_coverage_all_colours(lst_node), lst_seq,
     (lst_orientation == forward ? lst_r : lst_f),
     (lst_orientation == forward ? lst_f : lst_r));
     */
    binary_kmer_to_seq(&fst_kmer, flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? kmer_size : kmer_size - 1, fst_seq);

    int i, current;
    for(i = 0, current = 1; i < path->kmer_size; i++, current++) {
        fprintf(fout, "%c", fst_seq[i]);
        if(current % PATH_FASTA_LINE == 0) {
            fprintf(fout, "\n");
        }
    }

    size_t len = strlen(path->seq);
    for(i = 0; i < len; i++, current++){
        if (path->step_flags[i] & PRINT_LABEL_AS_N) {
            fprintf(fout, "N");
        } else {
            fprintf(fout, "%c",  path->seq[i]);
        }
        if(current % PATH_FASTA_LINE == 0){
            fprintf(fout, "\n");
        }
    }
    fprintf(fout, "\n");
    fflush(fout);
}

void output_seq_with_line_breaks(char* seq, FILE* fout, int* current)
{
    int i;

    for (i=0; i <strlen(seq); i++) {
        fprintf(fout, "%c", seq[i]);
        if (*current % PATH_FASTA_LINE == 0) {
            fprintf(fout, "\n");
        }
        *current = * current + 1;
    }
}
void output_seq_without_line_breaks(char* seq, FILE* fout)
{
    int i;

    for (i=0; i <strlen(seq); i++) {
        fprintf(fout, "%c", seq[i]);
    }
}

void * initalise_gfa_stats(gfa_stats * gfa, int max_length){
  gfa->H_count=0;
  gfa->S_count=0;
  gfa->L_count=0;
  gfa->P_count=0;
  gfa->current_S_line=0;
  gfa->overlap=0;
  gfa->max_length=max_length;
  gfa->gap_or_comma[0]=' ';
  gfa->orient[0]='3';
  gfa->P_line = calloc(gfa->max_length, sizeof(char));
  gfa->P_line_overlap = calloc(gfa->max_length, sizeof(char));

  return 0;
}

void * destroy_gfa_stats(gfa_stats * gfa){
  if (gfa){
    free(gfa->P_line);
    free(gfa->P_line_overlap);
    free(gfa);
  }
  return 0;
}

void post_polymorph_L_lines(FILE * file_gfa, gfa_stats * gfa){
  while (gfa->S_count > gfa->current_S_line ){
    output_L_line(file_gfa, gfa);
    gfa->current_S_line++;
  }
}

void check_orient_gfa(gfa_stats * gfa, Orientation orientation){
  if(orientation==forward){
    gfa->orient[0]='+';
  }
  else if(orientation==reverse){
    gfa->orient[0]='-';
  }
  else{
    gfa->orient[0]='3';
  }
}

void output_S_line(FILE * file_gfa, gfa_stats * gfa, char* seq){
  fprintf(file_gfa, "\nS %qd_%d ", gfa->H_count, gfa->S_count);
  if(seq!=NULL){
    output_seq_without_line_breaks(seq, file_gfa);
  }
}

void output_L_line(FILE * file_gfa, gfa_stats * gfa){
  fprintf(file_gfa, "\nL %qd_%d %c %qd_%d %c %dM", gfa->H_count, gfa->current_S_line, gfa->orient[0], gfa->H_count, gfa->S_count, gfa->orient[0], gfa->overlap);
}

void add_to_P_line(gfa_stats * gfa){
  //check for sizes of P_line's (not greater than max size of string)
  char text[100]; // HACK FOR NOW
  char overlap[100]; // HACK FOR NOW

  sprintf(text, "%c%qd_%d%c", gfa->gap_or_comma[0], gfa->H_count, gfa->S_count, gfa->orient[0]);
  sprintf(overlap, "%c%dM", gfa->gap_or_comma[0], gfa->overlap);
  if (strlen(gfa->P_line)>(gfa->max_length-strlen(text)) || strlen(gfa->P_line_overlap)>(gfa->max_length-strlen(overlap))){
    // output error? output P line and start a new one?
    log_and_screen_printf("ERROR: gfa P line grew too large\n");
  }
  else{
    if(strlen(gfa->P_line)>1){
      // add to P_line_overlap if this is NOT the first entry
      if(strlen(gfa->P_line_overlap)==0){ // HACK - runs through this for every P entry, but only relevant for first
        // if first P_overlap entry
        sprintf(overlap, " %dM", gfa->overlap);
      }
      strcat(gfa->P_line_overlap, overlap);
    }

    // always add to P_line
    strcat(gfa->P_line, text);
  }
}

void output_P_line(FILE * file_gfa, gfa_stats * gfa){
  if(strlen(gfa->P_line)>0){
    gfa->S_count++;
    fprintf(file_gfa, "\nP %qd_%d", gfa->H_count, gfa->S_count);
    output_seq_without_line_breaks(gfa->P_line, file_gfa);
    output_seq_without_line_breaks(gfa->P_line_overlap, file_gfa);
  }
}


//boolean output_polymorphism(Path* path, int* path_pos, dBGraph* graph, FILE* file_fastg, FILE* file_gfa, int* current, gfa_stats * gfa_count){
//    return true;
//}


boolean output_polymorphism(Path* path, int* path_pos, dBGraph* graph, FILE* file_fastg, FILE* file_gfa, int* current, gfa_stats * gfa_count)
{
    Path* paths[4];
    dBNode* current_node = path->nodes[*path_pos];
    Orientation orientation = path->orientations[*path_pos];
    int chosen_edge = path->labels[*path_pos];
    int max_path_length=(path->kmer_size)*2;
    pathStep current_step, reverse_step, next_step;
    path_step_initialise(&current_step);
    path_step_initialise(&reverse_step);
    path_step_initialise(&next_step);
    int j;
    int p = 1;
    int differ_pos = -1;
    int count = 0;
    char tempseq[max_path_length];
    boolean keep_going = true;

    // Get sub path to next branch point
    current_step.node = current_node;
    current_step.label = chosen_edge;
    current_step.orientation = orientation;
    current_step.flags = 0;
    paths[chosen_edge] = path_new(max_path_length, graph->kmer_size);
    db_graph_get_perfect_path_with_first_edge_all_colours(&current_step, &db_node_action_do_nothing, paths[chosen_edge], graph);
    log_printf("Chosen edge: %d length %d seq %s\n", chosen_edge, paths[chosen_edge]->length, paths[chosen_edge]->seq);

    // Find paths that end up at the same point as the chosen path
    void check_edge(Nucleotide nucleotide) {
      // ? path_new called shortly after this, not sure what this is for
        //if (nucleotide != chosen_edge) {
            // initialise path
        //    paths[nucleotide] = 0;
        //}

        if (nucleotide != chosen_edge) {
            // memory not allocated at this point yet?
            paths[nucleotide] = 0;
            if (db_node_edge_exist_any_colour(current_node, nucleotide, orientation)) {
                current_step.node = current_node;
                current_step.label = nucleotide;
                current_step.orientation = orientation;
                current_step.flags = 0;
                db_graph_get_next_step(&current_step, &next_step, &reverse_step, graph);

                paths[nucleotide] = path_new(max_path_length, graph->kmer_size);
                db_graph_get_perfect_path_with_first_edge_all_colours(&current_step, &db_node_action_do_nothing, paths[nucleotide], graph);
                // count++;
                log_printf("Got path %d length %d seq %s\n", nucleotide, paths[nucleotide]->length, paths[nucleotide]->seq);

                if (paths[nucleotide]->nodes[paths[nucleotide]->length-1] == paths[chosen_edge]->nodes[paths[chosen_edge]->length - 1]) {
									  char seq[1024];
									  binary_kmer_to_seq(&(paths[chosen_edge]->nodes[paths[chosen_edge]->length - 1]->kmer), graph->kmer_size, seq);
                    log_printf("Got matching path at end node %s\n", seq);
                    count++;
                } else {
                    path_destroy(paths[nucleotide]);
                    paths[nucleotide] = 0;
                    log_printf("Destroyed non-matching path\n");
                }
            }
        }
    }

    count = 0;
    nucleotide_iterator(&check_edge);

    // If no other edges in this orientation, then return and continue
    if (count == 0) {
        return false;
    }

    log_printf("Paths to compare:\n");
    for (j=0; j<4; j++) {
			  if (paths[j] != 0) {
					    if (j == chosen_edge) {
								  log_printf("%dE\t", j);
							} else {
								log_printf("%d\t", j);
							}
							log_printf("%s\n",paths[j]->seq);
				}
		}

    // Now to work out the difference between the paths.
    // Start at end and find point paths differ
    p=1;
    while (keep_going && (differ_pos == -1)) {
        for (j=0; j<4; j++) {
            if (paths[j] != 0) {
                if (p > paths[j]->length) {
                    keep_going = false;
                    differ_pos = p;
                    break;
                } else {
                    if (j != chosen_edge) {
                        if (paths[j]->seq[paths[j]->length - p] != paths[chosen_edge]->seq[paths[chosen_edge]->length - p]) {
                            differ_pos = p;
                        }
                    }
                }
            }
        }
        p++;
    }
    p = differ_pos;

    // If we haven't found the point they differ, then something has gone badly wrong.
    if (differ_pos == -1) {
        printf("Error: Something went wrong calculating difference!");
        exit(1);
    }


    // check through paths - is this a snp, an indel?
    count = 0;
    int insert_size=0;
    int best_path_cov=0;
    int best_path_length=0;
    int max_coverage_nucleotide=0;
    int indel_flag = 0;

    for (j=0; j<4; j++) {
        if (paths[j] != 0) {
            log_printf("Paths:%i\n", j);
            strncpy(tempseq, paths[j]->seq, paths[j]->length - differ_pos + 1);
            tempseq[paths[j]->length - differ_pos + 1] = 0;
            if (strlen(tempseq) > 0) {
                count++;

                // is this the highest coverage path?
                double avg_coverage=0;
                uint32_t min_coverage=0;
                uint32_t max_coverage=0;
                path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, paths[j]);
                if(avg_coverage>=best_path_cov){
                  if(strlen(tempseq)>best_path_length){
                    max_coverage_nucleotide=j;
                    best_path_length=strlen(tempseq);
                  }
                }

                // find the longest length path in the bubble
                if (strlen(tempseq) > insert_size) {
                  insert_size=strlen(tempseq);
                }
            }
            else{
              // indel - length is zero, but still meets same point. Can only be one of these, impossible for more.
              indel_flag = 1;
            }
        }
    }

    if (indel_flag){
      // for an indel, add an L (link) directly between segments before and after bubble.
      int S_count=gfa_count->S_count;
      gfa_count->S_count=S_count+count+1;
      if(file_gfa!=NULL){
        output_L_line(file_gfa, gfa_count);
      }
      gfa_count->S_count=S_count;
    }

    // Output the highest coverage alternative state BEFORE the brackets
    strncpy(tempseq, paths[max_coverage_nucleotide]->seq, paths[max_coverage_nucleotide]->length - differ_pos + 1);
    if (strlen(tempseq) < 1){
      printf("Error: Something went wrong looking at bubble path lengths!");
      exit(1);
    }
    else{
      output_seq_with_line_breaks(tempseq, file_fastg, current);
      if(file_gfa!=NULL){
        //if a polymorphism has previously been reached, and the perfect path
        //  after it output to gfa, then the L lines must be completed as well.

        if((file_gfa!=NULL)&&(gfa_count->S_count>gfa_count->current_S_line)){
          post_polymorph_L_lines(file_gfa, gfa_count);
        }

        // 'S' line at this branching point will overlap with immediately subsequent 'S' lines
        gfa_count->current_S_line=gfa_count->S_count;
        gfa_count->S_count++;
        output_S_line(file_gfa, gfa_count, tempseq);
        output_L_line(file_gfa, gfa_count);
        check_orient_gfa(gfa_count, paths[max_coverage_nucleotide]->orientations[0]);
        add_to_P_line(gfa_count);
      }
    }

    // Now output difference in square brackets
    output_seq_with_line_breaks("[", file_fastg, current);

    if (insert_size==1){
      output_seq_with_line_breaks("1:alt:allele|", file_fastg, current);
    }
    else if (insert_size>1){
      output_seq_with_line_breaks("1:alt|", file_fastg, current);
    }
    else{
      // shouldn't happen - insert_size<1?
    }



    for (j=0; j<4; j++) {
        if (paths[j] != 0) {
            strncpy(tempseq, paths[j]->seq, paths[j]->length - differ_pos + 1);
            tempseq[paths[j]->length - differ_pos + 1] = 0;

            if (strlen(tempseq) > 0) {
                output_seq_with_line_breaks(",", file_fastg, current);
                output_seq_with_line_breaks(tempseq, file_fastg, current);

                if((file_gfa!=NULL)&&(j!=max_coverage_nucleotide)){
                  // want to include whole of kmer here instead, and overlap?
                  gfa_count->S_count++;
                  output_S_line(file_gfa, gfa_count, tempseq);
                  output_L_line(file_gfa, gfa_count);
                }
                count++;
            }
        }
    }

    output_seq_with_line_breaks("]", file_fastg, current);

    // outout remaining sequence after the branching section
    strcpy(tempseq, paths[chosen_edge]->seq + (paths[chosen_edge]->length - differ_pos + 1));

    output_seq_with_line_breaks(tempseq, file_fastg, current);

    if(file_gfa!=NULL){
      gfa_count->S_count++;
      gfa_count->current_S_line++;

      //fprintf(file_gfa, "S %d ",gfa_count->S_count);
      output_S_line(file_gfa, gfa_count, NULL);
      check_orient_gfa(gfa_count, paths[chosen_edge]->orientations[0]);
      add_to_P_line(gfa_count);
      output_seq_without_line_breaks(tempseq, file_gfa);
      //fprintf(file_gfa, "\n");
    }

    // Update path pos
    *path_pos = *path_pos + paths[chosen_edge]->length - 1;

    // Destroy paths
    for (j=0; j<4; j++) {
        if (paths[j] != 0) path_destroy(paths[j]);
    }

    // Go through path, making Edge array, one entry for each polymorphism
    //

    // typedef struct {
    //     Nucleotide labels[MAX_LENGTH] ; String of labels
    //     boolean is_highest_coverage;
    // } PolyQueueItem;
    // for (i=0; i<path->length; i++) {
    //     if (node has POLYMORPHISM set) {
    //          Make a list of valid bubble edges in "edges"
    //          if (queue is empty) {
    //              create new items in queue and mark which is highest coverage
    //          } else {
    //              for each item in queue {
    //                  remove item from queue;
    //                  for each label that is part of edges {
    //                      add label to item and add back to queue
    //                      if (item->is_highest_coverage) {
    //                          if (label is not highest coverage) {
    //                              item->is_highest_coverage = false;
    //                          }
    //                      }
    //                  }
    //              }
    //          }
    //     }
    // }

    // Then to output
    // Loop through each item in queue
    // Output path, using the bases in the item

    return true;
}

void path_to_fastg_gfa(Path * path, FILE * file_fastg, FILE * file_gfa, HashTable* graph)
{
    short kmer_size = path->kmer_size;
    int length = path->length;
    int max_length = path->max_length;
    int path_pos = 0;
    gfa_stats * gfa_count = malloc(sizeof(gfa_stats));

    initalise_gfa_stats(gfa_count, max_length);
    gfa_count->H_count=path->id;



    // Sanity checking
    if (path == NULL) {
        fprintf(stderr,	"[path_to_fasta] trying to print a null Path\n");
        exit(-1);
    }

    if (file_fastg == NULL) {
        fprintf(stderr,	"[path_to_fasta] trying to print to a null FILE\n");
        exit(-1);
    }

    // file_gfa may be NULL, no sense in checking

    if (DEBUG) {
        printf("[path_to_fasta] About to print a path\n");
    }

/*
    if (length == max_length) {
        log_and_screen_printf("contig length equals max length [%i] for node_%i\n", max_length, path->id);
    }
*/

    // Set value of PRINT_FIRST (whether first node included)
    //check_print_first(path);

    // Get coverage statistics from path
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);


    // Get orientation of first and last node
    Orientation fst_orientation;
    fst_orientation = path->orientations[0];
    Orientation lst_orientation = path->orientations[path->length];

    // Get the first node - this will be nodes[0] if PRINT_FIRST is
    // specified, or nodes[1] otherwise.
    dBNode *fst_node;
    if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
        if (path->length == 0) {
            fprintf(stderr, "[path_to_fasta] Trying to print an empty path[1]!\n");
            return;
        }
        fst_node = path->nodes[0];
    } else {
        if (path->length < 2) {
            fprintf(stderr,	"[path_to_fasta] Trying to print an empty path[2]!\n");
            return;
        }
        fst_node = path->nodes[1];
    }

    // Get the last node
    dBNode *lst_node = path->nodes[path->length - 1];

    // Make a set of labels for first and last nodes which list the
    // acceptable forward and reverse path labels
    char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
    compute_label(fst_node, forward, fst_f);
    compute_label(fst_node, reverse, fst_r);
    compute_label(lst_node, forward, lst_f);
    compute_label(lst_node, reverse, lst_r);

    // Places to store first and last kmer sequence
    char fst_seq[kmer_size + 1];
    //char lst_seq[kmer_size + 1];
    fst_seq[kmer_size] = '\0';
    //lst_seq[kmer_size] = '\0';
    BinaryKmer tmp_kmer;

    // Get the first kmer sequence (or complement)
    BinaryKmer fst_kmer;
    binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
    if (fst_orientation == reverse) {
        binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
    if (DEBUG) {
        printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }

    // Output to file
    fprintf(file_fastg,
            ">node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_r:%s fst_f:%s lst_coverage:%i lst_r:%s lst_f:%s\n",
            path->id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length + kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));


    // generate the initial kmer in the path
    binary_kmer_to_seq(&fst_kmer, flags_check_for_flag(PRINT_FIRST,	&(path->flags)) ? kmer_size : kmer_size - 1, fst_seq);

    int i, current= 1;

    for(i = 0, current = 1 ; i < path->kmer_size; i++, current++) {
        fprintf(file_fastg, "%c", fst_seq[i]);
        if (current % PATH_FASTA_LINE == 0) {
            fprintf(file_fastg, "\n");
        }
    }

    // Output to gfa file, start of first segment
    if(file_gfa!=NULL){
      gfa_count->S_count++;
      gfa_count->current_S_line=gfa_count->S_count;
      output_S_line(file_gfa, gfa_count, fst_seq);
      check_orient_gfa(gfa_count, path->orientations[0]);
      add_to_P_line(gfa_count);
      gfa_count->gap_or_comma[0]=',';
    }


    /*if(file_fastg!=NULL){
      for(i = 0, current = 1 ; i < path->kmer_size; i++, current++) {
          fprintf(file_fastg, "%c", fst_seq[i]);
          // never have a line break
      }
    }*/

    path_pos = 0;
    //gfa_count->new_gfa_S = false;
    while (path_pos < strlen(path->seq)) {
        boolean skip_this = false;

        if (path->nodes[path_pos]->flags & POLYMORPHISM) {
            char seq[1024];
						binary_kmer_to_seq(&(path->nodes[path_pos]->kmer), graph->kmer_size, seq);
            log_printf("Found polymorphism to output at node %s\n", seq);

            /*if(file_fastg!=NULL){
              for(i = 0, current = 1 ; i < path->kmer_size; i++, current++) {
                  fprintf(file_fastg, "%c", fst_seq[i]);
                  // never have a line break
              }
            }*/

            skip_this = output_polymorphism(path, &path_pos, graph, file_fastg, file_gfa, &current, gfa_count);
        }

        if (skip_this == false) {
            fprintf(file_fastg, "%c",  path->seq[path_pos]);

            if(file_gfa!=NULL){
              fprintf(file_gfa, "%c",  path->seq[path_pos]);
            }

            if (current % PATH_FASTA_LINE == 0) {
                fprintf(file_fastg, "\n");
            }

            path_pos++;
            current++;
        }

    }

    // catch any remaining L lines
    if(file_gfa!=NULL){
      post_polymorph_L_lines(file_gfa, gfa_count);
      output_P_line(file_gfa, gfa_count);
    }

    fprintf(file_fastg, "\n");
    fflush(file_fastg);

    if(file_gfa!=NULL){
      fprintf(file_gfa, "\n");
      fflush(file_gfa);
    }

    destroy_gfa_stats(gfa_count);
}

void path_to_fasta_with_statistics(Path * path, FILE * fout, double avg_coverage, uint32_t min_coverage, uint32_t max_coverage)
{
    short kmer_size = path->kmer_size;
    int length = path->length;
    //int max_length = path->max_length;

    // Sanity checking
    if (path == NULL) {
        fprintf(stderr,	"[path_to_fasta] trying to print a null Path\n");
        exit(-1);
    }

    if (fout == NULL) {
        fprintf(stderr,	"[path_to_fasta] trying to print to a null FILE\n");
        exit(-1);
    }

    if (DEBUG) {
        printf("[path_to_fasta] About to print a path\n");
    }

/*
    if (length == max_length) {
        log_and_screen_printf("contig length equals max length [%i] for node_%i\n", max_length, path->id);
    }
*/

    // Get orientation of first and last node
    Orientation fst_orientation;

    fst_orientation = path->orientations[0];

    Orientation lst_orientation = path->orientations[path->length - 1];

    // Get the first node - this will be nodes[0] if PRINT_FIRST is
    // specified, or nodes[1] otherwise.
    dBNode *fst_node;
    if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
        if (path->length == 0) {
            fprintf(stderr, "[path_to_fasta] Trying to print an empty path[1]!\n");
            return;
        }
        fst_node = path->nodes[0];
    } else {
        if (path->length < 2) {
            fprintf(stderr,	"[path_to_fasta] Trying to print an empty path[2]!\n");
            return;
        }
        fst_node = path->nodes[1];
    }

    // Get the last node
    dBNode *lst_node = path->nodes[path->length - 1];

    // Make a set of labels for first and last nodes which list the
    // acceptable forward and reverse path labels
    char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
    compute_label(fst_node, forward, fst_f);
    compute_label(fst_node, reverse, fst_r);
    compute_label(lst_node, forward, lst_f);
    compute_label(lst_node, reverse, lst_r);

    // Places to store first and last kmer sequence
    char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
    fst_seq[kmer_size] = '\0';
    lst_seq[kmer_size] = '\0';
    BinaryKmer tmp_kmer;

    // Get the first kmer sequence (or complement)
    BinaryKmer fst_kmer;
    binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
    if (fst_orientation == reverse) {
        binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);
    if (DEBUG) {
        printf("[path_to_fasta] First kmer: %s\n", fst_seq);
    }

    // Get the last kmer sequence (or complement)
    BinaryKmer lst_kmer;
    binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
    if (lst_orientation == reverse) {
        binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);
    
    char* path_id;
    if(path->subpath_id == 0)
    {
        asprintf(&path_id, "node_%qd", path->id);
    }
    else
    {
        asprintf(&path_id, "node_%qd.%hi", path->id, path->subpath_id);
    }

    // Output to file
    fprintf(fout,
            ">%s length:%i average_coverage:%.2f min_coverage:%u max_coverage:%u fst_coverage:%u fst_r:%s fst_f:%s lst_coverage:%u lst_r:%s lst_f:%s\n",
            path_id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length + kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));

    binary_kmer_to_seq(&fst_kmer, flags_check_for_flag(PRINT_FIRST,	&(path->flags)) ? kmer_size : kmer_size - 1, fst_seq);
    
    free(path_id);
    
    int i, current= 1;

    for(i = 0, current = 1 ; i < path->kmer_size; i++, current++) {
        fprintf(fout, "%c", fst_seq[i]);
        if(current % PATH_FASTA_LINE == 0) {
            fprintf(fout, "\n");
        }
    }

    size_t len = strlen(path->seq);
    for(i = 0; i < len; i++, current++) {
        if (path->step_flags[i] & PRINT_LABEL_AS_N) {
            fprintf(fout, "*");
        } else {
            fprintf(fout, "%c",  path->seq[i]);
        }

        if(current % PATH_FASTA_LINE == 0) {
            fprintf(fout, "\n");
        }
    }
    fprintf(fout, "\n");

    fflush(fout);    
}

void path_to_fasta(Path * path, FILE * fout)
{
    // Get coverage statistics from path
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);
    path_to_fasta_with_statistics(path, fout, avg_coverage, min_coverage, max_coverage);
}

// Cloned from path_to_fasta with changes for file colour coding
// NOTE: Has not been updated for PRINT_LABEL_AS_N flag
void path_to_fasta_colour(Path * path, FILE * fout, char *id)
{
    short kmer_size = path->kmer_size;

    // Sanity checking
    if (path == NULL) {
        fprintf(stderr,	"[path_to_fasta_colour] trying to print a null Path\n");
        exit(-1);
    }

    if (fout == NULL) {
        fprintf(stderr,	"[path_to_fasta_colour] trying to print to a null FILE\n");
        exit(-1);
    }

    if (DEBUG) {
        printf("[path_to_fasta_colour] About to print a path\n");
    }
    // Set value of PRINT_FIRST (whether first node included)
    //check_print_first(path);

    // Get coverage statistics from path
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

    // Get orientation of first and last node
    Orientation fst_orientation;
    if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
        fst_orientation = path->orientations[0];
    } else {
        fst_orientation = path->orientations[1];
    }
    Orientation lst_orientation = path->orientations[path->length];

    // Get the first node - this will be nodes[0] if PRINT_FIRST is
    // specified, or nodes[1] otherwise.
    dBNode *fst_node;
    if (flags_check_for_flag(PRINT_FIRST, &(path->flags))) {
        if (path->length == 0) {
            fprintf(stderr,	"[path_to_fasta_colour] Trying to print an empty path[1]!\n");
            return;
        }
        fst_node = path->nodes[0];
    } else {
        if (path->length < 2) {
            fprintf(stderr,	"[path_to_fasta_colour] Trying to print an empty path[2]!\n");
            return;
        }
        fst_node = path->nodes[1];
    }

    // Get the last node
    dBNode *lst_node = path->nodes[path->length - 1];

    // Make a set of labels for first and last nodes which list the
    // acceptable forward and reverse path labels
    char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
    compute_label(fst_node, forward, fst_f);
    compute_label(fst_node, reverse, fst_r);
    compute_label(lst_node, forward, lst_f);
    compute_label(lst_node, reverse, lst_r);

    // Places to store first and last kmer sequence
    char fst_seq[kmer_size + 1], lst_seq[kmer_size + 1];
    fst_seq[kmer_size] = '\0';
    lst_seq[kmer_size] = '\0';
    BinaryKmer tmp_kmer;

    // Get the first kmer sequence (or complement)
    BinaryKmer fst_kmer;
    binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
    if (fst_orientation == reverse) {
        binary_kmer_reverse_complement(&fst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);

    if (DEBUG) {
        printf("[path_to_fasta_colour] First kmer: %s\n", fst_seq);
    }

    // Get the last kmer sequence (or complement)
    BinaryKmer lst_kmer;
    binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));
    if (lst_orientation == reverse) {
        binary_kmer_reverse_complement(&lst_kmer, kmer_size, &tmp_kmer);
        binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
    }
    binary_kmer_to_seq(&lst_kmer, kmer_size, lst_seq);

    // Output to file
    fprintf(fout,
            ">%s length:%i %s average_coverage:%.2f min_coverage:%u max_coverage:%u fst_coverage:%u fst_r:%s fst_f:%s lst_coverage:%u lst_r:%s lst_f:%s\n",
            id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? ((int)strlen(path->seq) + kmer_size) : ((int)strlen(path->seq) + kmer_size - 1)),
            path->header,
            avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));

    fprintf(fout, "%s", binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq));
    fprintf(fout, "%s\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? path->seq : path->seq + 1);

    fflush(fout);

    log_and_screen_printf("path_to_fasta: PRINT_FIRST %i path->seq length %i path->length %i\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 1:0, strlen(path->seq), path->length);
}

void path_to_coverage(Path * path, FILE * fout)
{
    if (path == NULL) {
        fprintf(stderr,	"[path_to_coverage] trying to print a null Path");
        exit(-1);
    }

    if (fout == NULL) {
        fprintf(stderr,	"[path_to_coverage] trying to print to a null FILE");
        exit(-1);
    }

    check_print_first(path);

    int i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;

    fprintf(fout, ">node_%qd \n", path->id);

    for (; i < path->kmer_size - 1; i++) {
        fprintf(fout, "%u ", element_get_coverage_all_colours(path->nodes[flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1]));
    }

    for (i = 0; i < path->length; i++) {
        fprintf(fout, "%u ", element_get_coverage_all_colours(path->nodes[i]));
    }

    fprintf(fout, "\n");
}

// Clone of path_to_coverage for handling colour coded files
void path_to_coverage_colour(Path * path, FILE * fout, char *id, short colour)
{
    short kmer_size = path->kmer_size;

    // Sanity checking
    if (path == NULL) {
        fprintf(stderr,	"[path_to_coverage_colour] trying to print a null Path");
        exit(-1);
    }

    if (fout == NULL) {
        fprintf(stderr,	"[path_to_coverage_colour] trying to print to a null FILE");
        exit(-1);
    }
    // Set value of PRINT_FIRST (whether first node included)
    //check_print_first(path);
    int i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;

    // If colour 0, write the label
    if (colour == 0) {
        fprintf(fout, ">%s\n", id);
    }
    // Write first kmer coverage
    for (; i < kmer_size - 1; i++) {
        fprintf(fout, "%u ", element_get_coverage_by_colour(path->nodes[flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 :	1], colour));
    }

    // Write path coverage
    i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;
    for (; i < path->length; i++) {
        fprintf(fout, "%u ", element_get_coverage_by_colour(path->nodes[i], colour));
    }

    fprintf(fout, "\n");

    log_and_screen_printf("path_to_coverage: PRINT_FIRST %i path->length %i\n", flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 1:0, path->length);
}

void path_iterator_from_index(int index, void (*step_action) (pathStep * step), Path * path)
{
    int i = 0;
    pathStep ps;
    path_step_initialise(&ps);

    if (index >= path->length) {
        printf("[path_iterator_from_index] Error: index (%d) greater than path length (%d)\n", index, path->length);
        exit(1);
    }

    for (i = index; i < path->length; i++) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        step_action(&ps);
    }
}

void path_inner_iterator(void (*step_action) (pathStep * step), Path * path)
{

    int i;
    pathStep ps;
    path_step_initialise(&ps);

    for (i = 1; i < path->length - 1; i++) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        step_action(&ps);
    }
}

void path_iterator(void (*step_action) (pathStep * step), Path * path)
{
    int i = 0;
    pathStep ps;
    path_step_initialise(&ps);
    for (i = 0; i < path->length; i++) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        step_action(&ps);
    }
}

void path_iterator_with_args(void (*step_action) (pathStep * , void *),void * args,  Path * path)
{
    int i = 0;
    pathStep ps;
    path_step_initialise(&ps);
    for (i = 0; i < path->length; i++) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        step_action(&ps, args);
    }
}

boolean path_is_singleton(int length, Path * path){
    boolean sing = false;
    double avg_cov;
    uint32_t min_cov, max_cov;
    pathStep first;
    pathStep last;
    path_step_initialise(&first);
    path_step_initialise(&last);
    path_get_step_at_index(0, &first, path);
    path_get_last_step(&last, path);

    if(path_get_length(path) < length){
        if(path_is_blunt(forward, path ) && path_is_blunt(reverse, path)){
            sing = true;
        }else if(path_get_length(path) == 1 && (path_is_blunt(forward, path ) || path_is_blunt(reverse, path))){
            sing = true;
        }else if(first.node == last.node){
            sing = true;
        }else{
            path_get_statistics(&avg_cov, &min_cov, &max_cov, path);
            if (avg_cov  < 2) {
                sing = true;
            }
        }

    }

    return sing;
}


boolean path_is_repetitive(double graph_cov, Path * p)
{
    double avg_cov = 0;
    uint32_t min_cov = 0;
    uint32_t max_cov = 0;
    boolean rep= false;
    if(path_has_stop_reason(FIRST, PATH_FLAG_DIVERGING_PATHS, p) &&path_has_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, p) ){
        if(path_get_length(p)  < p->kmer_size * 2){
            rep = true;
        }else{
            path_get_statistics(&avg_cov, &min_cov, &max_cov, p);
            if(avg_cov > (graph_cov * 3) ){
                rep = true;
            }
        }

    }

    if(path_get_length(p)  < p->kmer_size * 2){
        if(avg_cov > (graph_cov * 3) ){
            rep = true;
        }
    }

    return rep;

}

void path_iterator_reverse(void (*step_action) (pathStep * step), Path * path)
{
    int i = 0;
    pathStep ps;
    path_step_initialise(&ps);
    for (i = path->length - 1; i >= 0; i--) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        ps.flags = 0;
        step_action(&ps);
    }
}

void path_iterator_with_index(void (*step_action) (int index, pathStep * step), Path * path)
{
    int i = 0;
    pathStep ps;
    path_step_initialise(&ps);
    for (i = 0; i < path->length; i++) {
        ps.node = path->nodes[i];
        ps.label = path->labels[i];
        ps.orientation = path->orientations[i];
        ps.flags = 0;
        step_action(i,&ps);
    }
}



/**
 *
 *Gets the index of the step in the path, or -1 if it
 * is not found.
 */
int path_index_of_step(pathStep * step, Path * path)
{

    //TODO fill it.
    return -1;
}

pathStep *path_get_last_step_reverse(pathStep * step, Path * path)
{
    return path_get_step_reverse(step, path, path->length - 1);
}

// Reverse a step at a given index
pathStep *path_get_step_reverse(pathStep * step, Path * path, int index)
{
    assert (path != NULL) ;
    assert(step!=NULL);
    //assert(path->length != 1);
    assert(path->length >=1);
    assert(index < path->length);

    step->node = path->nodes[index];
    step->orientation = opposite_orientation(path->orientations[index]);
    step->flags = path->step_flags[index];

    if (path->length == 1) {
        step->label = Undefined;
    } else {
        if (index == 0) {
            step->label = Undefined;
        } else {
            BinaryKmer second_to_last_kmer;
            binary_kmer_assignment_operator(second_to_last_kmer, path->nodes[index - 1]->kmer);

            if (path->orientations[index - 1] == reverse) {
                step->label = binary_kmer_get_last_nucleotide(&second_to_last_kmer);
            } else {
                BinaryKmer tmp_kmer;
                BinaryKmer *rev_kmer = binary_kmer_reverse_complement(&second_to_last_kmer, path->kmer_size, &tmp_kmer);
                step->label = binary_kmer_get_last_nucleotide(rev_kmer);
            }
        }
    }
    step->path = path;
    return step;
}

pathStep *path_get_step_at_index(int index, pathStep * step, Path * path)
{
    assert(path != NULL);
    assert(step != NULL);
    assert(index < path->length);
    if (path == NULL) {
        fprintf(stderr,
                "[path_get_step_at_index] passing a  null Path");
        exit(-1);
    }
    if (step == NULL) {
        fprintf(stderr,
                "[path_get_step_at_index] passing a  pathStep ");
        exit(-1);
    }
    if (index >= path->length) {
        fprintf(stderr,
                "[path_get_step_at_index] The queried index (%d) is greater than the path length (%d)",
                index, path->length);
        exit(-1);
    }
    step->node = path->nodes[index];
    step->orientation = path->orientations[index];
    step->label = path->labels[index];
    step->flags = path->step_flags[index];
    step->path = path;
    return step;
}

char *path_get_seq(char *tmp, Path * path)
{
    strcpy(tmp, path->seq);
    return tmp;
}

void path_graphviz_open_header(FILE * f){
    fprintf(f, "graph{\n"
            "graph[page=\"8.5,11\",size=\"7.5,7\",ratio=fill,center=1];\n"
            "fontsize=18;\n"
            "node [shape = circle];\n"
            "edge [dir=none];\n");


}

void path_graphviz_close_header(FILE * f){
    fprintf(f, "}\n");

}

void path_graphviz_line(FILE * f, Path * p){

    if(p == NULL || p->length == 0)
        return;

    pathStep  first ;
    path_step_initialise(&first);
    pathStep last;
    path_step_initialise(&last);

    path_get_step_at_index(0, &first, p);
    path_get_last_step(&last, p);

    short kmer_size = p->kmer_size;
    char  seq_f[kmer_size + 1];
    char  seq_l[kmer_size + 1];
    BinaryKmer bk;
    Key k = &bk;
    BinaryKmer * tmp = element_get_kmer(first.node) ;
    binary_kmer_to_seq(element_get_key(tmp ,kmer_size, k), kmer_size, seq_f);

    tmp = element_get_kmer(last.node) ;
    binary_kmer_to_seq(element_get_key(tmp ,kmer_size, k), kmer_size, seq_l);
    double size = ((double)p->length/(100000));
    size =  size;
    //double weight  = fabs(log(1/size));
    fprintf(f, "%lld [ shape=point, color=red];\n", p->id);
    fprintf(f, "%s [ shape=point, color=blue];\n",seq_f);
    fprintf(f, "%s [ shape=point, color=blue];\n",seq_l);
    fprintf(f, "%s -- %lld [len=\"%f\"];\n", seq_f, p->id, size);
    fprintf(f, "%lld -- %s [len=\"%f\"];\n", p->id, seq_l, size);
}

void path_print_contig_with_details(FILE * fout, Path * p){

    pathStep  ps;
    path_step_initialise(&ps);
    int i;
    fprintf(fout, ">node_%qd \n", p->id);
    for(i=0; i<path_get_length(p); i++){
        path_get_step_at_index(i, &ps, p);
        //   fprintf(fout, );
    }

}

void path_get_statistics_between_points(double *avg_coverage, uint32_t *min_coverage, uint32_t *max_coverage, Path * path, int start, int end)
{
    /**
    * TODO: validate if the path is empty... think about singletons....
    */
    assert(end <= path->length);
    assert(start >= 0);
    assert(start < end);
    int i = start;
    if(i == 0)
    {
        i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;
    }
    *max_coverage = 0;
    *min_coverage = UINT_MAX;
    uint32_t sum_coverage = 0;

    for (; i < end; i++) {	//Calculate the return values for the current path.

      uint32_t coverage = element_get_coverage_all_colours(path->nodes[i]);
      sum_coverage += coverage;
      *max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
      *min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;

    }
    int length = end-start;
    *avg_coverage = (double)sum_coverage / length;

    if (*min_coverage == UINT_MAX) {
      *min_coverage = 0;
    }  
}
void path_get_statistics(double *avg_coverage, uint32_t *min_coverage, uint32_t *max_coverage, Path * path)
{
    path_get_statistics_between_points(avg_coverage, min_coverage, max_coverage, path, 0, path->length);
}

void path_get_coverage_standard_deviation(double* standard_deviation, double avg_coverage, Path* path)
{
    if(path->length == 1)
    {
        *standard_deviation = 0;
        return;
    }
    int i = 0;
    i = flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? 0 : 1;
    
    int sum_deviations = 0;

    for (; i < path->length; i++) 
    {
      uint32_t coverage = element_get_coverage_all_colours(path->nodes[i]);
      double diff = (avg_coverage - coverage) * (avg_coverage - coverage);
      sum_deviations += diff;
    }

    *standard_deviation = sqrt(sum_deviations/(path->length - 1));
}

dBNode *path_last_node(Path * path)
{
    return path->length > 0 ? path->nodes[path->length - 1] : NULL;
}

int path_get_length(Path * path){
    assert(path!=NULL);
    return path->length;
}

pathStep *path_get_last_step(pathStep * ps, Path * path)
{
    assert(path->length != 0);
    ps->node = path_last_node(path);
    ps->orientation = path_last_orientation(path);
    ps->label = path_last_nucleotide(path);
    ps->flags = path_last_flags(path);
    return ps;
}

boolean path_is_empty(Path * path)
{
    //return flags_check_for_flag(EMPTY_PATH, &(path->flags));
    assert(path != NULL);
    return path->length == 0;
}

/**
 * Returns if a path is blunt.
 * If orientation is forward, it checks if the last node is blunt
 * If the orientation is reverse, it checks if the first node is blunt.
 */
boolean path_is_blunt(Orientation o, Path * p){
    if (path_get_length(p) > 0) {
        pathStep ps;
        path_step_initialise(&ps);

        if(o == forward){
            path_get_last_step(&ps, p);
        }else if(o == reverse){
            path_get_step_reverse(&ps, p, 0);
        }else{
            fprintf(stderr, "Invalid orientation at path_is_blunt\n");
            exit(-1);
        }

        return db_node_is_blunt_end(ps.node, ps.orientation);
    }else{
        return  false;
    }

}

void path_add_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
    if(o == FIRST){
        flags_action_set_flag(f, &(path->stop_reasons_first));
    }else if(o == LAST){
        flags_action_set_flag(f, &(path->stop_reasons_last));
    }
}

boolean path_has_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
    return o == FIRST?flags_check_for_flag(f, &(path->stop_reasons_first)): flags_check_for_flag(f, &(path->stop_reasons_last));
}

boolean path_has_any_stop_reason(PathEnd o, Flags f, Path * path){
    assert(o == FIRST || o == LAST);
    return o == FIRST?flags_check_for_any_flag(f, &(path->stop_reasons_first)): flags_check_for_flag(f, &(path->stop_reasons_last));
}

void path_clean_stop_reason(Path * path){

    path->stop_reasons_first = 0;
    path->stop_reasons_last = 0;
    //flags_action_clear_flags(&(path->stop_reasons_first));
    //flags_action_clear_flags(&(path->stop_reasons_last));
}



boolean path_is_cycle(Path * path)
{

    /*if(path->len < 2){
     return false;
     }
     if(!flags_check_for_flag(IS_CYCLE, &(path->flags))){
     int len = path->length;
     int i,j;
     pathStep psi, psj;

     for(i = 0; i < len; i++){
     path_get_step_at_index(i, &psi, path);

     for(j = i+i; j < len; j++){
     path_get_step_at_index(j, &psj, path);

     if(path_step_equals_without_label(&psi, &psj)){
					i = len;
					j = len;
					flags_action_set_flag(IS_CYCLE, &(path->flags));
					//if (DEBUG) {
     fprintf(stdout,
     "[path_add_node] The node_%qd is cycle!\n",
     path->id);
					//}
     }
     }
     }

     }*/
    return flags_check_for_flag(IS_CYCLE, &(path->flags));
}

Orientation path_last_orientation(Path * path)
{
    if(path->length == 0)
        return undefined;
    return path->orientations[path->length - 1];
}

int path_free_spaces(Path * path)
{
    return path->max_length - path->length;
}

void path_remove_last(Path * path)
{
    if (path->length == 0) {
        printf
        ("[path_remove_last] Trying to remove node from the node_%qd which is empty\n",
         path->id);
        flags_action_set_flag(ERROR_PATH, &(path->flags));
        flags_action_unset_flag(FIND_AGAIN, &(path->flags));
        assert(path->length > 0);
        return;		//TODO: use the flags to tell the caller the invalid state
    }

    if (DEBUG) {
        char tmp_seq[path->kmer_size + 1];
        tmp_seq[path->kmer_size] = 0;	//'\0';
        printf
        ("[path_remove_last] Removing node %i in path(%lld): %s %c\n",
         path->length - 1, path->id,
         binary_kmer_to_seq(element_get_kmer(path_last_node(path)),
                            path->kmer_size, tmp_seq),
         path->seq[path->length - 1] ==
         0 ? 'N' : path->seq[path->length - 1]);
        //printNode(current_node, db_graph->kmer_size);
    }

    if(db_node_edges_count_all_colours(path_last_node(path), path_last_orientation(path)) > 1){
        path->out_nodes_count--;
        if (path->out_nodes_count < 0) {
            path->out_nodes_count = 0;
            printf("Warning: path->out_nodes_count < 0\n");
        }
        //assert(path->out_nodes_count >=0);
    }

    path->length--;

    path->seq[path->length] = '\0';
    path->nodes[path->length] = NULL;
    path->orientations[path->length] = 0;
    path->labels[path->length] = 0;
    path->step_flags[path->length] = 0;
    if ((path->in_nodes > 0) && (path->in_nodes_count > 0)) {
        if (path->in_nodes[path->in_nodes_count-1] == path->length){
            path->in_nodes[path->in_nodes_count-1] = 0;
            path->in_nodes_count--;
        }
    }
    flags_action_set_flag(FIND_AGAIN, &(path->flags));
    if (path->length == 0) {
        //flags_action_clear_flags(&(path->flags));
        //flags_action_set_flag(PRINT_FIRST, &(path->flags));
        //flags_action_set_flag(NEW_PATH, &(path->flags));
        flags_action_set_flag(EMPTY_PATH, &(path->flags));
        flags_action_unset_flag(FIND_AGAIN, &(path->flags));
    }
}

void path_do_nothing(Path * p)
{
}

boolean path_to_retry(Path * path)
{

    if (DEBUG) {
        printf("[path_to_retry]Flags: %x\n", path->flags);
    }
    return flags_check_for_any_flag(FIND_AGAIN | NEW_PATH, &(path->flags))
    && !flags_check_for_any_flag(ERROR_PATH, &(path->flags));
}

void path_step_assign(pathStep * to, pathStep * from)
{
    to->label = from->label;
    to->node = from->node;
    to->orientation = from->orientation;
    to->path = from->path;
    to->flags = from->flags;
}

void path_step_initialise(pathStep * step){
    step->node = NULL;
    step->orientation = 0;  // NOTE check this. Valid initialisation?
    step->label = 0;
    flags_action_clear_flags(&step->flags);
    step->path = NULL;
}

PathArray *path_array_new(short number_of_paths)
{
    PathArray *pa = calloc(1, sizeof(PathArray));
#ifdef THREADS
    pthread_mutex_init(&pa->mutex, NULL);
#endif
    if (pa == NULL) {
        fprintf(stderr,
                "[path_array_new] Not enough memory to allocate the PathArray");
        exit(-1);
    }
    pa->number_of_paths = 0;
    pa->capacity = number_of_paths;
    pa->paths = calloc(number_of_paths, sizeof(Path *));
    if (pa->paths == NULL) {
        fprintf(stderr,
                "[path_array_new] Not enough memory to allocate the Paths for the PathArray ");
    }
    return pa;
}

int path_array_get_number_of_paths(PathArray * pa){
    return pa->number_of_paths;
}

void path_array_destroy(PathArray * pa)
{

    //TODO: design a cleaver lock.
    while (pa->number_of_paths > 0) {
        pa->number_of_paths--;
        path_destroy(pa->paths[pa->number_of_paths]);

    }
    free(pa->paths);
    // pa->number_of_paths=0; // not sure if this needs to be freed?
    free(pa);
}

void path_array_destroy_struct (PathArray ** pa){
    free((*pa)->paths);
    free((*pa));
    (*pa)=NULL;
}

//WARNING: This is not thread safe, the method calling it must be.
static void  path_array_double_capacity(PathArray * pa){

    int new_capacity = pa->capacity * 2;
    Path ** new_array = realloc(pa->paths, new_capacity * sizeof(Path * ));
    if(new_array == NULL){
        fprintf(stderr, "[path_array_double_capacity] Unable to double the size of the PathArray (size %d)", new_capacity);
        assert(new_array != NULL);
        exit(-1);
    }
    pa->capacity = new_capacity;
    pa->paths = new_array;



}

boolean path_array_add_path(Path * p, PathArray * pa)
{
    if (pa == NULL) {
        fprintf(stderr, "[path_array_add_path] PathArray is null");
        exit(-1);
    }

    if (p == NULL) {
        fprintf(stderr, "[path_array_add_path] Path is null");
        exit(-1);
    }

    if (pa->number_of_paths >= pa->capacity) {
        path_array_double_capacity(pa);
    }
    pa->paths[pa->number_of_paths] = p;
    pa->number_of_paths++;
    return true;
}

void path_action_clear_flags(Path * node)
{
    flags_action_clear_flags(&(node->flags));
}

void path_action_set_flag(Path * node, Flags f)
{
    flags_action_set_flag(f, &(node->flags));

}

void path_action_unset_flag(Path * node, Flags f)
{
    flags_action_unset_flag(f, &(node->flags));
}

Flags path_get_flags(Path * node, Flags f)
{
    return node->flags & f;
}

boolean path_check_for_flag(Path * node, Flags flag)
{
    return flags_check_for_flag(flag, &(node->flags));

}

boolean path_check_for_any_flag(Path * path, Flags flag)
{
    return flags_check_for_any_flag(flag, &(path->flags));
}

int path_percentage_new_nodes(Path * path){
    return (100 * path->new_nodes)/path->length;

}

void path_reverse(Path * source, Path * destination)
{
    pathStep new_step;
    path_step_initialise(&new_step);
    int i;

    for (i = source->length - 1; i >= 0; i--) {
        path_get_step_reverse(&new_step, source, i);
        new_step.flags = new_step.flags & PATH_STEP_MASK_VISITED; //Because we are reversing, we cant keep track of all the paths which were already visited in this walk
        path_add_node(&new_step, destination);
        // Update step flags - could be more elegant
        //destination->step_flags[destination->length-1] = source->step_flags[i];
    }

}

int path_get_index_of_last_in_node(Path * p){
    int count = p->in_nodes_count;
    int index = 0;
    if (count) {
        index = p->in_nodes[count-1];
    }
    return index;

}

boolean path_append(Path * destination, Path * source){
    int i;
    pathStep new_step;
    pathStep first_step;
    pathStep last_step;
    path_step_initialise(&new_step);
    path_step_initialise(&first_step);
    path_step_initialise(&last_step);
    boolean success = true;

    if (source->length == 0) {
        fprintf(stderr, "[path_append] The source path is empty!\n");
        //exit(-1);
        assert(source->length != 0);
        return false;
    }

    if ((destination->length > 0)
        &&
        (unlabelled_path_step_equals
         (path_get_step_at_index(0, &first_step, source),
          path_get_last_step(&last_step, destination)))) {
             path_remove_last(destination);
             if (DEBUG) {
                 printf("[path_append] Removing last step.\n");
             }
         } else {
             if (DEBUG) {
                 printf("[path_append] No need to remove last step.\n");
             }
         }

    for (i = 0; i < source->length; i++) {
        new_step.node = source->nodes[i];
        new_step.orientation = source->orientations[i];
        new_step.label = source->labels[i];
        new_step.flags = source->step_flags[i];
        if (!path_add_node(&new_step, destination)) {
            log_printf("[path_append] Could not append paths successfully.\n");
            success = false;
            break;
        } else {
            // Update step flags - could be more elegant
            //destination->step_flags[destination->length - 1] = source->step_flags[i];
            //destination->in_nodes[destination->length - 1] = source->in_nodes[i];
        }
    }

    return success;
}

int path_get_nodes_count(Path * path){
    return path->length;
}

int path_get_edges_count(Path * path){
    return path->length - 1;
}

boolean paths_equal(Path * path_a, Path * path_b){
    int i;
    boolean paths_equal = true;

    if ((!path_a) || (!path_b))
        return false;

    if (path_a->length != path_b->length)
        return false;

    for (i = 0; i < path_a->length; i++) {
        if ((path_a->nodes[i] != path_b->nodes[i]) ||
            (path_a->labels[i] != path_b->labels[i]) ||
            (path_a->orientations[i] != path_b->orientations[i])) {
            paths_equal = false;
            break;
        }
    }

    return paths_equal;
}

static PathArray *path_buffers = NULL;
//static PathArray *short_path_buffers = NULL;

/**
 * Makes a copy of the "from" path into the "to" path.
 * WARNING! The to path is cleared before the copy.
 */
void path_copy(Path * to, Path * from)
{
    path_reset(to);
    to->max_virtual_length = from->max_virtual_length;
    to->flags = from->flags;
    to->stop_reasons_first = from->stop_reasons_first;
    to->stop_reasons_last = from->stop_reasons_last;
    //to->depth = from->depth;
    if (from->length > 0) {
        path_append(to, from);
    }

}

#ifdef THREADS
static pthread_once_t path_array_init_once = PTHREAD_ONCE_INIT;
#endif
static void path_buffers_init(){
    path_buffers = path_array_new(MAX_PATH_BUFFERS);
}

void path_array_initialise_buffers(short kmer_size)
{
#ifdef THREADS
    pthread_once(&path_array_init_once, path_buffers_init);
#else
    if (path_buffers == NULL) {
        path_buffers_init();
    }
#endif
    if (path_buffers->kmer_size == 0) {
        path_buffers->kmer_size = kmer_size;

        //		path_buffers_short = path_array_new(MAX_PATH_BUFFERS);
        //		for (i = 0; i < MAX_PATH_BUFFERS; i++) {
        //			tmp = path_new(MAX_PATH_LENGTH, kmer_size);
        //			tmp->id = i;
        //			path_array_add_path(tmp, path_buffers);

        //		}
    }
}

void path_array_destroy_buffers(){
    Path *tmp;
    int i;
    if (path_buffers != NULL){
        for (i = 0; i < path_buffers->number_of_paths; i++) {
            tmp = path_array_get(i, path_buffers);
            path_destroy(tmp);
        }
        free(path_buffers);
        path_buffers = NULL;
    }
}

PathArray *path_array_get_from_buffer_with_size(short size)
{
    PathArray *pa = path_array_new(size);
    int i;
    Path * tmp = NULL;
    for (i = 0; i < size; i++) {
        tmp = path_get_buffer_path();
        path_array_add_path(tmp, pa);
        tmp->id = i;
    }
    if (tmp != NULL) {
        pa->kmer_size= tmp->kmer_size;

    }
    return pa;
}

Path * path_array_get(int path, PathArray *pa){
    assert(path < pa->number_of_paths);
    return pa->paths[path];
}

void path_array_free_from_buffer(PathArray * pa)
{
    while (pa->number_of_paths > 0) {
        pa->number_of_paths--;
        path_free_buffer_path(pa->paths[pa->number_of_paths]);

    }
    free(pa->paths);
    free(pa);

}

void path_array_to_fasta(FILE * f, PathArray * pa){
    int i;
    for (i = 0; i<pa->number_of_paths; i++) {
        if (path_get_length(path_array_get(i,pa)) > 0) {
            path_to_fasta(path_array_get(i, pa), f);
        }

    }
}

Path *path_get_buffer_path()
{
    //TODO make this thread safe
    assert(path_buffers !=NULL);

#ifdef THREADS
    pthread_mutex_lock(&path_buffers->mutex);
#endif

    Path *tmp = NULL, *found = NULL;
    int i;
    for (i = 0; i < path_buffers->number_of_paths && found == NULL; i++) {
        tmp = path_buffers->paths[i];
        if (!tmp->used) {
            found = tmp;

        }
    }
    tmp = NULL;
    assert(i < MAX_PATH_BUFFERS); //TODO: make this a growing array.
    if(found == NULL){

        tmp = path_new(MAX_PATH_LENGTH, path_buffers->kmer_size);
        tmp->id = i;
        path_array_add_path(tmp, path_buffers);
        found = path_buffers->paths[i];

        //This is the old logic, that was in the initializer.
        //		path_buffers_short = path_array_new(MAX_PATH_BUFFERS);
        //		for (i = 0; i < MAX_PATH_BUFFERS; i++) {
        //			tmp = path_new(MAX_PATH_LENGTH, kmer_size);
        //			tmp->id = i;
        //			path_array_add_path(tmp, path_buffers);

        //		}
    }

    assert(found != NULL);
    found->used = true;
#ifdef THREADS
    pthread_mutex_unlock(&path_buffers->mutex);
#endif
    return found;
}

void path_free_buffer_path(Path * path)
{
#ifdef THREADS
    pthread_mutex_lock(&path_buffers->mutex);
#endif
    path_reset(path);
    path->used = false;

#ifdef THREADS
    pthread_mutex_unlock(&path_buffers->mutex);
#endif
}

void path_step_mark_as_uncertain(int i, Path * path, boolean as_n) {
    if (as_n) {
        path->step_flags[i] |= PRINT_LABEL_AS_N;
    } else {
        path->step_flags[i] |= PRINT_LABEL_LOWERCASE;
    }
}

boolean is_step_marked_as_uncertain(int i, Path * path) {
    boolean r = (path->step_flags[i] & (PRINT_LABEL_AS_N | PRINT_LABEL_LOWERCASE)) > 0;

    return r;
}

static void step_mark_visited(pathStep * ps) {
    db_node_action_set_flag_visited(ps->node);
}

void path_mark_as_visited(Path* path) {
    path_iterator(&step_mark_visited, path);
}

// TO DO: This needs to be done properley!
void path_pairs_to_fasta(PathArray* pa, int distances[], FILE* fout) {
    int i, j;
    int total_length = 0;
    int kmer_size;

    // Sanity checking
    if (pa == NULL) {
        fprintf(stderr,
                "[path_pairs_to_fasta] trying to print a null Path\n");
        exit(-1);
    }

    if (pa->number_of_paths < 2) {
        fprintf(stderr,
                "[path_pairs_to_fasta] trying to print less than one path\n");
        exit(-1);
    }

    if (fout == NULL) {
        fprintf(stderr,
                "[path_pairs_to_fasta] trying to print to a null FILE\n");
        exit(-1);
    }

    kmer_size = pa->paths[0]->kmer_size;

    for (i=0; i<pa->number_of_paths; i++) {
        total_length += strlen(pa->paths[i]->seq) + kmer_size;
        if (i < (pa->number_of_paths-1)) {
            total_length += distances[i];
        }
    }

    int current = 1;
    fprintf(fout, ">rpnode_%qd length:%i\n", pa->paths[0]->id, total_length);
    for (i=0; i<pa->number_of_paths; i++) {
        Path *path = pa->paths[i];
        dBNode* fst_node = path->nodes[0];
        BinaryKmer fst_kmer;
        char fst_seq[kmer_size+1];
        binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));
        binary_kmer_to_seq(&fst_kmer, kmer_size, fst_seq);

        // Print first kmer
        for (j=0; j<kmer_size; j++, current++) {
            fprintf(fout, "%c", fst_seq[j]);
            if(current % PATH_FASTA_LINE == 0){
                fprintf(fout, "\n");
            }
        }

        // Print rest
        for (j=0; j<strlen(path->seq); j++, current++) {
            if (path->step_flags[j] & PRINT_LABEL_AS_N) {
                fprintf(fout, "N");
            } else {
                fprintf(fout, "%c",  path->seq[j]);
            }
            if(current % PATH_FASTA_LINE == 0){
                fprintf(fout, "\n");
            }
        }
        fprintf(fout, "\n");

        // Print Ns
        if (i < (pa->number_of_paths-1)) {
            for (j=0; j<distances[i]; j++, current++) {
                if(current % PATH_FASTA_LINE == 0){
                    fprintf(fout, "\n");
                }
                fprintf(fout, "N");
            }
        }
    }
}

void path_counts_reset(PathCounts * pc)
{
    pc->blunt_ends = 0;
    pc->converging_paths = 0;
    pc->diverging_paths = 0;
    pc->is_double_y = 0;
    pc->is_cycle = 0;
    pc->longer_than_buffer = 0;

    pc-> minimum_double_y = 0;
    pc-> total_double_y_lenght = 0;
}

void path_counts_print_and_log(PathCounts * pc)
{
    log_and_screen_printf("Blunt Ends\t%'lld\n", pc->blunt_ends);
    log_and_screen_printf("Converging paths\t%'lld\n", pc->converging_paths);
    log_and_screen_printf("Diverging paths\t%'lld\n", pc->diverging_paths);
    log_and_screen_printf("Is cycle\t%'lld\n", pc->is_cycle);
    log_and_screen_printf("Is double y\t%'lld\n", pc->is_double_y);
    log_and_screen_printf("Longer tha buffer\t%'lld\n", pc->longer_than_buffer);
}

void path_counts_add(Path * p, PathCounts * pc)
{
    if(path_has_stop_reason(LAST, PATH_FLAG_STOP_BLUNT_END,  p)){
        pc->blunt_ends++;
    }
    if(path_has_stop_reason(LAST, PATH_FLAG_CONVERGING_PATHS, p)){
        pc->converging_paths++;
    }
    if(path_has_stop_reason(LAST, PATH_FLAG_DIVERGING_PATHS, p)){
        pc->diverging_paths++;
    }

    if(path_has_stop_reason(LAST, PATH_FLAG_IS_CYCLE, p)){
        pc->is_cycle++;
    }

    if(path_has_stop_reason(LAST, PATH_FLAG_IS_DOUBLE_Y, p)){
        pc->is_double_y++;
    }

    if(path_has_stop_reason(FIRST, PATH_FLAG_STOP_BLUNT_END,  p)){
        pc->blunt_ends++;
    }
    if(path_has_stop_reason(FIRST, PATH_FLAG_CONVERGING_PATHS, p)){
        pc->converging_paths++;
    }
    if(path_has_stop_reason(FIRST, PATH_FLAG_DIVERGING_PATHS, p)){
        pc->diverging_paths++;
    }

    if(path_has_stop_reason(FIRST, PATH_FLAG_IS_CYCLE, p)){
        pc->is_cycle++;
    }

    if(path_has_stop_reason(FIRST, PATH_FLAG_IS_DOUBLE_Y, p)){
        pc->is_double_y++;
    }
}

boolean path_copy_subpath(Path* dest_path, const Path* source_path, int start, int end)
{
    assert(source_path != NULL);
    assert(end - start > 0);
    assert(start >= 0);
    assert(end <= source_path->length);
    if(dest_path == NULL)
    {
        // some error message
        printf("[path_copy_subpath] Error: Dest path NULL. ");
        return false;
    }
    int length = end - start;
    if(length > dest_path->max_length)
    {
        printf("[path_copy_subpath] Error: Dest path too small to copy subpath into. ");
        return false;
    }
    
    path_reset(dest_path);  
    for(int i = start; i < end; i++)
    {
        pathStep next_step;
        next_step.node = source_path->nodes[i];
        next_step.orientation = source_path->orientations[i];
        next_step.label = source_path->labels[i];
        next_step.flags = source_path->step_flags[i];
        if(!path_add_node(&next_step, dest_path))
        {
            char kmer_string[source_path->kmer_size + 1];
            kmer_string[source_path->kmer_size] = '\0';
            BinaryKmer kmer; 
            binary_kmer_assignment_operator(kmer, *element_get_kmer(next_step.node));
            binary_kmer_to_seq(&kmer, source_path->kmer_size, kmer_string);
            printf("[path_copy_subpath] Error: Could not copy subpath. Attempting to add node with kmer %s failed.", kmer_string);
            return false;          
        }
    }
    
    for(int j = 0; j < dest_path->length-1; j++)
    {
        if(dest_path->seq[j] == '\0')
        {
            printf("[path_copy_subpath] Warning: path with irregular sequence. Sequence: %s, length: %i.\n", dest_path->seq, dest_path->length);
        }
    }
    dest_path->flags = source_path->flags;
    dest_path->id = source_path->id;
    return true;
}

void write_fastg_alt(const char* sequence1, const char* sequence2, FILE* file_fastg)
{
    int length = strlen(sequence1);
    fprintf(file_fastg, "[%i:alt|%s,%s]", length, sequence1, sequence2);
}

void path_mark_path_with_flag(Path* path, Flags f)
{
    for(int i = 0; i < path->length; ++i)
    {
        flags_action_set_flag(f, &path->nodes[i]->flags);
    }   
}

void path_unmark_path_with_flag(Path* path, Flags f)
{
    for(int i = 0; i < path->length; ++i)
    {
        flags_action_unset_flag(f, &path->nodes[i]->flags);
    }   
}

void path_to_gfa2_and_fastg(Path* path, dBGraph* graph, FILE* file_gfa, FILE* file_fastg)
{
    // Sanity checking
    if (path == NULL) {
        fprintf(stderr,	"[path_to_gfa_and_fastg] trying to print a null Path\n");
        exit(-1);
    }   

    assert(file_gfa != NULL);
    assert(file_fastg != NULL);
    
    // Clear graph of CURRENT_PATH_FORWARD/REVERSE flags
    hash_table_traverse_no_progress_bar(&db_node_action_unset_flag_current_path, graph);
  
    // write fastg header  
    short kmer_size = path->kmer_size;
    int length = path->length;
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);

    // Get orientation of first and last node
    Orientation fst_orientation;
    fst_orientation = path->orientations[0];
    Orientation lst_orientation = path->orientations[path->length];

    // Get the first node - this will be nodes[0] if PRINT_FIRST is
    // specified, or nodes[1] otherwise.
    dBNode *fst_node;
    int start_pos = 0;
    if(flags_check_for_flag(PRINT_FIRST, &(path->flags)))
    {
        if(path->length == 0) 
        {
            fprintf(stderr, "[path_to_fasta] Trying to print an empty path[1]!\n");
            return;
        }
        fst_node = path->nodes[0];
    } 
    else
    {
        if(path->length < 2)
        {
            fprintf(stderr,	"[path_to_fasta] Trying to print an empty path[2]!\n");
            return;
        }
        fst_node = path->nodes[1];
        start_pos = 1;
    }

    // Get the last node
    int end_pos = path->length;
    dBNode *lst_node = path->nodes[path->length - 1];

    // Make a set of labels for first and last nodes which list the
    // acceptable forward and reverse path labels
    char fst_f[5], fst_r[5], lst_f[5], lst_r[5];
    compute_label(fst_node, forward, fst_f);
    compute_label(fst_node, reverse, fst_r);
    compute_label(lst_node, forward, lst_f);
    compute_label(lst_node, reverse, lst_r);
    
    // Output to file
    fprintf(file_fastg,
            "\n>node_%qd length:%i average_coverage:%.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_r:%s fst_f:%s lst_coverage:%i lst_r:%s lst_f:%s\n",
            path->id,
            (flags_check_for_flag(PRINT_FIRST, &(path->flags)) ? length + kmer_size : length + kmer_size - 1), avg_coverage,
            min_coverage,
            max_coverage,
            element_get_coverage_all_colours(fst_node),
            (fst_orientation == forward ? fst_r : fst_f),
            (fst_orientation == forward ? fst_f : fst_r),
            element_get_coverage_all_colours(lst_node),
            (lst_orientation == forward ? lst_r : lst_f),
            (lst_orientation == forward ? lst_f : lst_r));
    
    gfa_file_wrapper file_wrapper;
    file_wrapper.m_file = file_gfa;
    file_wrapper.m_segment_count = 0;
    file_wrapper.m_path_id = path->id;
    
    //start the recursion!

    fastg_recursion_level = 0;
    gfa_segment_array* in_segments = NULL;
    boolean skip_first = false;
    do
    {
        out_struct out = write_paths_between_nodes(path, start_pos, end_pos, graph, in_segments, skip_first, &file_wrapper, file_fastg);

        if(in_segments)
        {
            gfa_segment_array_destroy(in_segments);
            in_segments = NULL;
        }
        
        start_pos = out.m_new_start_pos;
        skip_first = out.m_skip_first;
        in_segments = out.m_segments;
        
    }
    while(skip_first || start_pos < end_pos - 1);
    
    // destroy in segments  
    if(in_segments)
    {
        gfa_segment_array_destroy(in_segments);
    }    
 
}

/*------------------------------------------------------------------------------------------------------------------------------------------*
 * Function: write_paths_between_nodes                                                                                           
 * Purpose: Recursive function to write the gfa and fastg entries for section of path between to points, including polymorphisms.
 * Given a path in the graph (the "main path"), this function iterates through each node in the path from start_pos to end_pos.
 * If there is a branch at a node, then all paths from this node are explored.
 * If any of these paths overlap with the main path, then this is written as an alternative sequence in the fastg and gfa files,
 * and the algorithm is repeated for each subpath (from the branch, to where they join). Then, the algorithm is repeated from
 * the join to the end of the main path.
 * For fastg, only one level of polymorphism is written to the file. For gfa, there is no limit. (i.e. a snp in an allele would appear in
 * the gfa file, but only the allele would appear in the fastg file.)
 * Parameters: - Path* path (const)                             - Pointer to the path to write a gfa and fastg entry for.
 *             - int start_pos                                  - The position in the path to start from.
 *             - int end_pos                                    - The position in the path to end at.
 *             - dBGraph* graph (const)                         - The graph that contains the path.
 *             - gfa_segment_array* in_segments                 - Pointer to array of path segments from previous call that join at
 *                                                                start_pos.
 *             - boolean skip_first                             - Whether or not to skip the first node in the sequence. Used for indels.
 *             - gfa_file_wrapper* file_gfa                     - Pointer to wrapper for FILE object representing the gfa output file.
 *             - FILE* file_fastg                               - Pointer to FILE object representing the fastg output file.
 * Returns: gfa_segment_array giving all paths that finish at endpos. Usually there is just one, but in the case of a polymorphism there
 *          could be many. 
 *------------------------------------------------------------------------------------------------------------------------------------------*/
out_struct write_paths_between_nodes(  Path* path, 
                                int start_pos, 
                                int end_pos, dBGraph* graph, 
                                gfa_segment_array* in_segments,
                                boolean skip_first,
                                gfa_file_wrapper* file_gfa, 
                                FILE* file_fastg)
{
    //log_printf("\n--Write paths between nodes for sequence %s, between %i and %i, skip-first %i --\n", path->seq, start_pos, end_pos, skip_first);
    fastg_recursion_level++;
    //log_printf("Recursion level: %i\n", fastg_recursion_level);
    
    int current_pos = start_pos;
    subpath subpaths[4];
    for(int i = 0; i < 4; i++)
    {
        subpaths[i].m_join_pos = -1;
        subpaths[i].m_length = -1;
        subpaths[i].m_nucleotide = i;
        subpaths[i].m_path = NULL;
    }
    boolean polymorphism = false;
    int polymorphism_end_pos = -1;
    int polymorphism_start_pos = -1;
    
    int main_subpath = -1;
    
    if(skip_first)
    {
        current_pos += 1;
    }
    
    // Go through all the nodes in the path (except the last one) to check for branches.
    while(current_pos < end_pos)
    {       
        dBNode* current_node = path->nodes[current_pos];
        Orientation current_orientation = path->orientations[current_pos];
        int edges = db_node_edges_count(current_node, current_orientation);
        if(fastg_recursion_level <= MAX_GFA_RECURSIONS && edges > 1 && current_pos != end_pos - 1)
        {
            //log_printf("Potential polymorphism at position %i\n", current_pos);
            // reset the subpaths
            for(int i = 0; i < 4; i++)
            {
                subpaths[i].m_join_pos = -1;
                subpaths[i].m_length = -1;
                subpaths[i].m_nucleotide = i;
                subpaths[i].m_path = NULL;
            }
            
            // this node is a candidate for a polymorphism
            // check every edge leaving this node
            for(int i = 0; i < 4; i++)
            {
                if(db_node_edge_exist_any_colour(current_node, i, current_orientation))
                {
                    //log_printf("Branch at edge %i\n", i);

                    pathStep current_step;
                    current_step.node = current_node;
                    current_step.label = i;
                    current_step.orientation = current_orientation;
                    
                    pathStep next_step, rev_step;
                    db_graph_get_next_step(&current_step, &next_step, &rev_step, graph);
                    if(element_get_coverage_all_colours(next_step.node) >= graph->path_coverage_minimum)
                    {
                        if(next_step.node == path->nodes[current_pos + 1])
                        {
                            //log_printf("Found main path.\n"); 
                            main_subpath = i;

                            subpaths[i].m_join_pos = 0;
                            subpaths[i].m_nucleotide = i;
                            subpaths[i].m_length = current_pos - start_pos;
                        }                
                        else
                        {                
                            //log_printf("Creating path %i\n", i);
                            pathStep join_step = db_graph_search_for_bubble(path, &current_step, &subpaths[i].m_path, graph);
                            if(join_step.node != NULL)
                            {
                                assert(subpaths[i].m_path != NULL);
                                assert(subpaths[i].m_path->length > 0);
                                for(int j = current_pos; j < path->length; ++j)
                                {
                                    if(path->nodes[j] == join_step.node && path->orientations[j] == join_step.orientation)
                                    {
                                        polymorphism = true;                             
                                        subpaths[i].m_nucleotide = i;
                                        subpaths[i].m_join_pos = j;
                                        subpaths[i].m_length = subpaths[i].m_path->length;
                                        //log_printf("Path %i has length %i\n", i, subpaths[i].m_length );      
                                        break;
                                    }                   
                                }                                        
                            }
                            if(!polymorphism)
                            {
                                log_printf("[write_paths_between_nodes] Warning: Found alternate path but could not find join node. Ignoring...\n");
                                if(subpaths[i].m_path)
                                {
                                    path_destroy(subpaths[i].m_path);
                                    subpaths[i].m_path = NULL;
                                }
                                subpaths[i].m_join_pos = -1;
                                subpaths[i].m_length = -1;
                            }
                        }
                    }
                }
            }
                        
           assert(main_subpath != -1);
                                         
            if(polymorphism)
            {
                break;
            }
        }
        current_pos++;      
    }
   
    char first_kmer_string[path->kmer_size + 1];
    if(in_segments == NULL)
    {
        // write the first k nucleotides
        first_kmer_string[path->kmer_size] = '\0';
        BinaryKmer kmer; 
        binary_kmer_assignment_operator(kmer, *element_get_kmer(path->nodes[0]));
        if(path->orientations[0] == reverse)
        {
            BinaryKmer reverse_kmer;
            binary_kmer_reverse_complement(&kmer, path->kmer_size, &reverse_kmer);
            binary_kmer_to_seq(&reverse_kmer, path->kmer_size, first_kmer_string);
        }
        else
        {
            binary_kmer_to_seq(&kmer, path->kmer_size, first_kmer_string);
        }
    }
    
    // construct the segment from start_pos up to current_pos
    int sequence_length = current_pos - start_pos;
    
    uint32_t path_coverage_max = path->nodes[start_pos]->coverage[0];
    uint32_t path_coverage_min = path->nodes[start_pos]->coverage[0];
    double path_coverage_avg = (double)path->nodes[start_pos]->coverage[0];
    
    if(fastg_recursion_level == 1)
    {
        if(start_pos < path->length - 1)
        {
            int start_pos_for_coverage = 0;
            if(in_segments != NULL)
            {
                start_pos_for_coverage = start_pos;
                if(skip_first)
                {
                    start_pos_for_coverage +=1;
                }
            }
            int end_pos_for_coverage = start_pos_for_coverage == current_pos ? current_pos + 1 : current_pos;
            assert(start_pos_for_coverage < end_pos_for_coverage);
            path_get_statistics_between_points(&path_coverage_avg, &path_coverage_min, &path_coverage_max, path, start_pos_for_coverage, end_pos_for_coverage);
        }
    }
    else
    {
        if(start_pos == 0)
        {
            // this is an "alt" route, so we want the coverage over the whole path
            if(path->length > 2)
            {
                int start_pos_for_coverage = skip_first ? 1 : 0;
                assert(start_pos_for_coverage < path->length - 1);
                path_get_statistics_between_points(&path_coverage_avg, &path_coverage_min, &path_coverage_max, path, start_pos_for_coverage, path->length - 1);
            }
        }
        else
        {
            // this is a section of the main path.
            if(start_pos < path->length - 1)
            {
                int start_pos_for_coverage = start_pos;
                if(skip_first && start_pos < end_pos)
                {
                    start_pos_for_coverage += 1;
                }
                int end_pos_for_coverage = current_pos == end_pos ? end_pos + 1 : current_pos;
                assert(start_pos_for_coverage < end_pos_for_coverage);
                path_get_statistics_between_points(&path_coverage_avg, &path_coverage_min, &path_coverage_max, path, start_pos_for_coverage, end_pos_for_coverage);
            }
        }
    }

    
    gfa_segment_array* current_segment_array = gfa_segment_array_new(1);
    assert(sequence_length >= 0);
    if(sequence_length == 0)
    {
        if(in_segments == NULL)
        {
            gfa_segment_array_append(current_segment_array, (file_gfa->m_segment_count)++, first_kmer_string, forward, path_coverage_avg);
            write_gfa_segment_array(current_segment_array, file_gfa);
        }
        else
        {
            gfa_segment_array_merge(current_segment_array, in_segments);
        }
    }
    else
    {

        char sequence[sequence_length + 1];
        strncpy(sequence, path->seq + start_pos, sequence_length);
        sequence[sequence_length] = '\0';

        if(in_segments == NULL)
        {
            char new_sequence[graph->kmer_size + sequence_length + 1];
            strcpy(new_sequence, first_kmer_string);
            strcat(new_sequence, sequence);
            gfa_segment_array_append(current_segment_array, (file_gfa->m_segment_count)++, new_sequence, forward, path_coverage_avg);          
        }
        else
        {
            gfa_segment_array_append(current_segment_array, (file_gfa->m_segment_count)++, sequence, forward, path_coverage_avg);
        }
        write_gfa_segment_array(current_segment_array, file_gfa);

        if(in_segments)
        {
            gfa_segment* current_segment = gfa_segment_array_get(current_segment_array, 0);
            for(int i = 0; i < in_segments->m_length; i++)
            {
                write_gfa_edge(gfa_segment_array_get(in_segments, i), current_segment, file_gfa);
            }
        }

        // write to fastg
        if(fastg_recursion_level == 1)
        {
            if(in_segments == NULL)
            {
                fprintf(file_fastg, "%s", first_kmer_string);
            }
            fprintf(file_fastg, "%s",  sequence);
        }
    }
    
    if(!polymorphism)
    {
        //log_printf("No polymorphism\n");
        fastg_recursion_level--;
        for(int i = 0; i < 4; i++)
        {
            assert(subpaths[i].m_path == NULL);
        }
        out_struct out;
        out.m_new_start_pos = current_pos;
        out.m_segments = current_segment_array;
        out.m_skip_first = false;
        return out;
    }
    else
    {
        // for each subpath, backtrack until we find the first nucleotide where they differ.
/*
        for(int i = 0; i < 4; ++i)
        {
            if(subpaths[i].m_path)
            {
                int end_pos = subpaths[i].m_join_pos;
                int main_length = end_pos - current_pos;
                int overlap_between_paths = 0;
                int alt_length = subpaths[i].m_path->length;
                int min_path_length = alt_length < main_length ? alt_length : main_length;
                // back track along the paths to get to the first nucleotide where they differ
                while(overlap_between_paths < min_path_length && 
                      (path->seq[end_pos - 1 - overlap_between_paths] == 
                        subpaths[i].m_path->seq[alt_length - 1 - overlap_between_paths]))
                {
                    overlap_between_paths++;
                }
                subpaths[i].m_join_pos -= overlap_between_paths;
                subpaths[i].m_length = subpaths[i].m_path->length - overlap_between_paths;
                assert(subpaths[i].m_length >= 0);
            }
        } 
*/

        // sort the subpaths by length (ascending)
        int compare(const void* a, const void* b)
        {
            int x = ((subpath *)a)->m_join_pos; 
            int y = ((subpath *)b)->m_join_pos;
            if(x == y)
            {
                x = ((subpath *)a)->m_length;
                y = ((subpath *)b)->m_length; 
            }
            return x - y;
        }
        qsort((void*)subpaths, 4, sizeof(subpath), compare);
                
/*
        for(int i = 0; i < 4; i++)
        {
            log_printf("Path %i of join pos %i, length %i\n", i, subpaths[i].m_join_pos, subpaths[i].m_length);
        }
*/

        polymorphism_end_pos = subpaths[3].m_join_pos;                
        polymorphism_start_pos = current_pos;

        //log_printf("Found joining position %i, main path %i, current pos %i\n", polymorphism_end_pos, main_subpath, current_pos);
        assert(polymorphism_end_pos >= polymorphism_start_pos);
       
        if(fastg_recursion_level == 1)
        {
            int main_seq_length = polymorphism_end_pos - polymorphism_start_pos;
            assert(main_seq_length >= 0);
            
            char main_seq[main_seq_length + 1];
            strncpy(main_seq, path->seq + polymorphism_start_pos, main_seq_length);
            main_seq[main_seq_length] = '\0';
            
            int alt_seq_length = subpaths[3].m_length;
            assert(alt_seq_length >= 0);
            assert(subpaths[3].m_path);
            //assert(strlen(subpaths[3].m_path->seq) <= alt_seq_length);
            char alt_seq[alt_seq_length + 1];
            strncpy(alt_seq, subpaths[3].m_path->seq, alt_seq_length);
            alt_seq[alt_seq_length] = '\0';    

            fprintf(file_fastg, "%s", main_seq);
            write_fastg_alt(main_seq, alt_seq, file_fastg);
        }
             
        // deal with the polymorphisms     
        gfa_segment_array* out_segments = gfa_segment_array_new(1);
        gfa_segment_array_merge(out_segments, current_segment_array);
        
        assert(polymorphism_start_pos >= 0);
        int main_start_pos = polymorphism_start_pos;
        boolean first_path = true;
        int prev_end_pos = 0;
        for(int i = 0; i < 4; i++)
        {
            if(subpaths[i].m_path && subpaths[i].m_join_pos >= 0)
            {
                //log_printf("------------------------------\n");
                //log_printf("Starting recursion for path %i\n", subpaths[i].m_nucleotide);
                assert(subpaths[i].m_join_pos >= 0);
                
                //should these segments be the start of the next iteration?
                int main_end_pos = subpaths[i].m_join_pos;
                boolean new_end_pos = main_end_pos > prev_end_pos;
                
                gfa_segment_array* main_segments = NULL;
                gfa_segment_array* alt_segments = NULL;
                            
                //log_printf("End pos: %i\n", main_end_pos);
                
                // do the subpath            
                if(main_end_pos >= 0)
                {
                    //log_printf("Alt Path... length %i\n", subpaths[i].m_path->length);
                    out_struct alt_out = write_paths_between_nodes(subpaths[i].m_path, 0, subpaths[i].m_length, graph, current_segment_array, true, file_gfa, file_fastg);
                    alt_segments = alt_out.m_segments;
                }

                // Do the main path:
                if(new_end_pos)
                {
                    int length = main_end_pos - main_start_pos;
                    if(length > 0)
                    {
                        //log_printf("Main Path...\n");
                        out_struct main_out = write_paths_between_nodes(path, main_start_pos, main_end_pos, graph, out_segments, true, file_gfa, file_fastg);
                        main_segments = main_out.m_segments;
                    }
                    else
                    {
                        //log_printf("Indel: no segment\n");
                        if(first_path)
                        {
                            main_segments = gfa_segment_array_new(1);
                            gfa_segment_array_merge(main_segments, current_segment_array);
                        }
                    }
                }
                
                //get ready for the next iteration
                if(new_end_pos)
                {
                    gfa_segment_array_destroy(out_segments);
                    out_segments = gfa_segment_array_new(2);
                    gfa_segment_array_merge(out_segments, main_segments);
                    gfa_segment_array_merge(out_segments, alt_segments);
                }
                else
                {
                    gfa_segment_array_merge(out_segments, alt_segments);
                }
                gfa_segment_array_destroy(main_segments);
                gfa_segment_array_destroy(alt_segments);               
                
                main_start_pos = main_end_pos;
                first_path = false;
                prev_end_pos = main_end_pos;
            }
        }
        
        // All return segments finish at the same place, end_node, so we continue from there.      
        //int new_start_pos = polymorphism_end_pos - overlap_between_paths;
        int new_start_pos = polymorphism_end_pos;
        boolean include_last_step = false;
        if(new_start_pos == current_pos)
        {
            // Backtracked to current pos, so this polymorphism was an indel.
            // Start from the next node and add the previous step to the first sequence.
            new_start_pos += 1;
            include_last_step = true;
        }
        
        if(current_segment_array)
        {
            gfa_segment_array_destroy(current_segment_array);
        }
        
        out_struct out;
        out.m_new_start_pos = new_start_pos;
        out.m_segments = out_segments;
        out.m_skip_first = include_last_step;
        
        fastg_recursion_level--;
        
        return out;
    }
}

PathArray* path_split_at_min_coverages(Path* path, int min_coverage)
{
    assert(path != NULL);
    assert(path_get_length(path) > 0);
    PathArray* pa = path_array_new(2);
    
    int start = -1;
    short subpath_id = 1;
    for(int i = 0; i < path->length; ++i)
    {
        dBNode* current_node = path->nodes[i];
        uint32_t current_coverage = element_get_coverage_all_colours(current_node);
        
        if(start < 0 && current_coverage >= min_coverage)
        {
            start = i;
        }
        else if(current_coverage < min_coverage)
        {
            if(start >= 0)
            {
                int length = i - start;
                if(length > 1)
                {
                    Path* new_subpath = path_new(length, path->kmer_size);
                    new_subpath->id = path->id;
                    new_subpath->subpath_id = subpath_id++;
                    path_copy_subpath(new_subpath, path, start, i);
                    path_array_add_path(new_subpath, pa);
                }
            }
            start = -1;                  
        }
    }
    
    // add the last path
    if(start > 0 && start < path->length)
    {
        int length = path->length - start;
        Path* new_subpath = path_new(length, path->kmer_size);
        new_subpath->id = path->id;
        new_subpath->subpath_id = subpath_id;
        path_copy_subpath(new_subpath, path, start, path->length);
        path_array_add_path(new_subpath, pa);
    }
    
    return pa;
}

void path_to_GFA_sequence(Path* path, gfa_file_wrapper* file_gfa, boolean include_first_kmer)
{
    assert(file_gfa != NULL);
    long int string_length = path->length + 1;
    short kmer_size = path->kmer_size;
    int start = 0;
    if(include_first_kmer)
    {
        string_length += kmer_size;
    }
    char sequence[string_length];
    if(include_first_kmer)
    {
        BinaryKmer kmer;
        char first_kmer_string[kmer_size + 1];
        first_kmer_string[kmer_size] = '\0';
        binary_kmer_assignment_operator(kmer, *element_get_kmer(path->nodes[0]));
        if(path->orientations[0] == reverse)
        {
            BinaryKmer reverse_kmer;
            binary_kmer_reverse_complement(&kmer, kmer_size, &reverse_kmer);
            binary_kmer_to_seq(&reverse_kmer, kmer_size, first_kmer_string);
        }
        else
        {
            binary_kmer_to_seq(&kmer, kmer_size, first_kmer_string);
        }
        strncpy(sequence, first_kmer_string, kmer_size);
        start += kmer_size;
    }
    
    strcpy(sequence + start, path->seq);
    sequence[string_length - 1] = '\0';
    
    uint32_t min_coverage;
    uint32_t max_coverage;
    double avg_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, path);
    
    fprintf(file_gfa->m_file, "S\tp%llds%i\t%lu\t%s\tCV:f:%f\n", 
                    file_gfa->m_path_id,
                    file_gfa->m_segment_count,
                    string_length - 1, 
                    sequence,
                    avg_coverage);  
}


void write_path_GFA_nodes(Path* path, dBGraph* graph, gfa_file_wrapper* file_gfa)
{
    // write the path as a GFA2 sequence
    path_to_GFA_sequence(path, file_gfa, true);
    log_printf("[write_path_GFA_nodes] Writing GFA for sequence %s\n", path->seq);
    
    for(int i = 0; i < path->length; i++)
    {
        dBNode* current_node = path->nodes[i];
        Orientation current_orientation = path->orientations[i];
        int edges = db_node_edges_count(current_node, current_orientation);
        if(edges > 1)
        {
            log_printf("[write_path_GFA_nodes] Found potential polymorphism\n");
            // check every edge leaving this node
            for(int e = 0; e < 4; e++)
            {
                if(db_node_edge_exist_any_colour(current_node, e, current_orientation))
                {
                    pathStep current_step;
                    current_step.node = current_node;
                    current_step.label = e;
                    current_step.orientation = current_orientation;
                    
                    pathStep next_step, rev_step;
                    db_graph_get_next_step(&current_step, &next_step, &rev_step, graph);
                    if( next_step.node != path->nodes[i+1] && 
                        element_get_coverage_all_colours(next_step.node) >= graph->path_coverage_minimum)
                    {
                        Path* alt_path = NULL;
                        log_printf("[write_path_GFA_nodes] Looking for alternate route...\n");
                        pathStep join_step = db_graph_search_for_bubble2(path, &current_step, &alt_path, graph);
                        if(join_step.node != NULL)
                        {
                            int join_pos = i;
                            for(int j = i; j < path->length; j++)
                            {
                                if(path->nodes[j] == join_step.node && path->orientations[j] == join_step.orientation)
                                {
                                    join_pos = j;
                                    break;
                                }
                            }
                            // write out GFA sequence
                            assert(alt_path != NULL);
                            file_gfa->m_segment_count++;
                            path_to_GFA_sequence(alt_path, file_gfa, false);
                            // write edge
                            int main_start = i + graph->kmer_size;
                            int main_end = join_pos + graph->kmer_size;
                            fprintf(file_gfa->m_file, "E\tp%llds0_p%llds%i\tp%llds0+\tp%llds%i+\t%i\t%i\t0\t%i$\t*\n",
                                file_gfa->m_path_id, file_gfa->m_path_id, file_gfa->m_segment_count,
                                file_gfa->m_path_id, file_gfa->m_path_id, file_gfa->m_segment_count,
                                    main_start, main_end, alt_path->length);   
                        }
                        else
                        {
                            log_printf("[write_path_GFA_nodes] Could not find path...\n");
                        }
                        if(alt_path)
                        {
                            path_destroy(alt_path);
                        }
                    }
                }
            }
        }
    }
}