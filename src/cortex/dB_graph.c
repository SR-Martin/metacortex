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

#include <structs.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <structs.h>

#include <global.h>
#include <flags.h>
#include <binary_kmer.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <cleaning.h>
#include <file_reader.h>
#include <path.h>
#include <perfect_path.h>
#include <logger.h>
#include <node_queue.h>

/**
 * Returns the next node, and the reverse to come back.
 * This function DO NOT  tells which is the next node
 * label to traverse, since it depends on the algorithm
 * to decide which one to use. It, however, gives the
 * reverse lable to be able to traverse back.
 *
 */
pathStep *db_graph_get_next_step(pathStep * current_step, pathStep * next_step, pathStep * rev_step, dBGraph * db_graph)
{
    assert(current_step != NULL);
    assert(next_step != NULL);
    assert(rev_step != NULL);

    assert(current_step->node != NULL);
    assert(current_step->label != Undefined);

    next_step->node = NULL;
    rev_step->node = current_step->node;
    next_step->flags = 0;
    rev_step->flags = current_step->flags & PATH_STEP_MASK_VISITED;

    BinaryKmer local_copy_of_kmer;
    binary_kmer_assignment_operator(local_copy_of_kmer,	current_step->node->kmer);

    BinaryKmer tmp_kmer;
    //dBNode * next_node = NULL;

    // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
    BinaryKmer *rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);

    if (current_step->orientation == reverse) {
        rev_step->label = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
        binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
    } else {//TODO: This could be avoided by reversing just the first nucleotide, not requiring to reverse always.
        rev_step->label = binary_kmer_get_last_nucleotide(rev_kmer);
    }

    binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);

    //get node from table
    next_step->node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
    rev_step->node = next_step->node;

    if (next_step->node != NULL) {
        next_step->orientation = db_node_get_orientation(&local_copy_of_kmer, next_step->node, db_graph->kmer_size);
        rev_step->orientation = opposite_orientation(next_step->orientation);
    }
//#ifdef __DEBUG
    else {
    //		if (DEBUG) {
    //
    char tmpseq[db_graph->kmer_size + 1];
    printf("[db_graph_get_next_step] Cannot find %s so get a NULL node\n", binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpseq));
    //Commented by ricardo, to reduce the log as for the traversing
    //      algorithm relays on having this as null
    //		}
    }
//#endif

    next_step->label = Undefined;

    return next_step;
}

pathStep *db_graph_get_next_step_with_reverse(pathStep * current_step, pathStep * next_step, pathStep * rev_step, dBGraph * db_graph)
{

	assert(current_step != NULL);
	assert(next_step != NULL);
	assert(rev_step != NULL);

	assert(current_step->node != NULL);
	assert(current_step->label != Undefined);

	next_step->node = NULL;
	rev_step->node = current_step->node;
    next_step->flags = 0;
    rev_step->flags = current_step->flags & PATH_STEP_MASK_VISITED;

	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer,	current_step->node->kmer);

	BinaryKmer tmp_kmer;
	//dBNode * next_node = NULL;

	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer =
    binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);

	if (current_step->orientation == reverse) {
		rev_step->label = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
		binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
	} else {//TODO: This could be avoided by reversing just the first nucleotide, not requiring to reverse always.
		rev_step->label = binary_kmer_get_last_nucleotide(rev_kmer);
	}

	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, current_step->label, db_graph->kmer_size);

	//get node from table
	next_step->node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);
	rev_step->node = next_step->node;

	if (next_step->node != NULL) {
		next_step->orientation = db_node_get_orientation(&local_copy_of_kmer, next_step->node, db_graph->kmer_size);
		rev_step->orientation = opposite_orientation(next_step->orientation);
	}
    //#ifdef __DEBUG
	else {
        //		if (DEBUG) {
        //
        char tmpseq[db_graph->kmer_size + 1];
        printf("[db_graph_get_next_step] Cannot find %s so get a NULL node\n",
               binary_kmer_to_seq(&tmp_kmer, db_graph->kmer_size, tmpseq));
        //Commented by ricardo, to reduce the log as for the traversing
        //      algorithm relays on having this as null
        //		}
	}
    //#endif

	next_step->label = Undefined;
	return next_step;
}

/**
 * This function resets the flags of the whole graph. A future improvement
 * should be to do this in a traditional look, using a single binary operation.
 * That will reduce the number of jumps, and should be paralellizable and vectorizable.
 */
void db_graph_reset_flags(dBGraph * db_graph)
{
	if (DEBUG) {
		printf("[db_graph_identify_branches] Cleaning flags...\n");
	}
	printf("Cleaning flags...\n");
#ifdef THREADS
	hash_table_threaded_traverse(&db_node_action_clear_flags, db_graph);
#else
	hash_table_traverse(&db_node_action_clear_flags, db_graph);
#endif
	printf("\n");
	fflush(stdout);
}

int db_graph_get_perfect_path_with_first_edge(pathStep * first_step, void (*node_action) (dBNode* node), Path * path, dBGraph * db_graph)
{

	dBNode *node = first_step->node;

	//Nucleotide fst_nucleotide = first_step->label;
	//printf("First node\n");
	//printNode(first_step->node, db_graph->kmer_size);
	dBNode *current_node = node;
    //Nucleotide nucleotide;
    Nucleotide nucleotide2;
	int length = 0;
	char tmp_seq[db_graph->kmer_size + 1];
	tmp_seq[db_graph->kmer_size] = '\0';

	//sanity check
	if (node == NULL) {
		printf
        ("[db_graph_get_perfect_path_with_first_edge] db_graph_get_perfect_path: can't pass a null node\n");
		exit(1);
	}

	path_reset(path);

	if (DEBUG) {
		printf
        ("[db_graph_get_perfect_path_with_first_edge] Node %i in path: %s\n",
         path->length,
         binary_kmer_to_seq(element_get_kmer(current_node),
                            db_graph->kmer_size, tmp_seq));
	}
	//first edge defined
	//nucleotide = fst_nucleotide;
	boolean added = true;

	pathStep current_step, next_step, rev_step;
	path_step_assign(&next_step, first_step);
	do {
		path_step_assign(&current_step, &next_step);
		/*current_step.node = next_step.node;
         current_step.orientation = next_step.orientation;
         current_step.label = next_step.label; */
		//printNode(current_step.node, db_graph->kmer_size);
		added = path_add_node(&current_step, path);
		if (added) {

			node_action(current_step.node);
			//                      db_graph_get_next_node(pathStep * current_step,
			//                                      pathStep * next_step, pathStep * rev_step, dBGraph * db_graph);
			db_graph_get_next_step(&current_step, &next_step,
                                   &rev_step, db_graph);

			//sanity check
			if (next_step.node == NULL) {
				fprintf(stderr,
                        "[db_graph_get_perfect_path_with_first_edge] dB_graph: didnt find node in hash table: %s %c %s\n",
                        binary_kmer_to_seq(element_get_kmer
                                           (current_step.node),
                                           db_graph->kmer_size,
                                           tmp_seq),
                        binary_nucleotide_to_char
                        (current_step.label),
                        current_step.orientation ==
                        forward ? "forward" : "reverse");
				exit(1);
			}

		}
	}
	while (added && !path_is_cycle(path) &&	//loop
	       db_node_has_precisely_one_edge(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&	//multiple entries
	       db_node_has_precisely_one_edge(next_step.node, next_step.orientation, &next_step.label));	//has one next edge only
	if (current_step.node != NULL) {
		added = path_add_node(&next_step, path);
		if (added) {
			node_action(next_step.node);
		}
	}
	if (DEBUG) {
		if (next_step.node == NULL) {
			printf("\n[db_graph_get_perfect_path_with_first_edge] Next node is null! \n");
		} else {
			printf("\n[db_graph_get_perfect_path_with_first_edge] Last node in path: %s %i length: %i\n",
                   binary_kmer_to_seq(element_get_kmer(next_step.node), db_graph->kmer_size, tmp_seq),
                   db_node_get_edges(next_step.node),	//next_node->edges,
			       path->length);
		}
	}

	return length;

}

// An all colours version of db_graph_get_perfect_path_with_first_edge
int db_graph_get_perfect_path_with_first_edge_all_colours(pathStep * first_step, void (*node_action) (dBNode * node), Path * path, dBGraph * db_graph)
{
    dBNode *node = first_step->node;

    //Nucleotide fst_nucleotide = first_step->label;
    dBNode *current_node = node;
    //Nucleotide nucleotide;
    Nucleotide nucleotide2;
    char tmp_seq[db_graph->kmer_size + 1];
    tmp_seq[db_graph->kmer_size] = '\0';

    //sanity check
    if (node == NULL) {
        printf("[db_graph_get_perfect_path_with_first_edge_all_colours] db_graph_get_perfect_path: can't pass a null node\n");
        exit(1);
    }

    path_reset(path);

    if (DEBUG) {
        printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Node %i in path: %s\n",
        path->length, binary_kmer_to_seq(element_get_kmer(current_node), db_graph->kmer_size, tmp_seq));
    }
    //first edge defined
    //nucleotide = fst_nucleotide;
    boolean added = true;

    pathStep current_step, next_step, rev_step;
    path_step_assign(&next_step, first_step);
    do {
        path_step_assign(&current_step, &next_step);
        added = path_add_node(&current_step, path);
        if (added) {
            node_action(current_step.node);
            //db_graph_get_next_node(pathStep * current_step,
            //pathStep * next_step, pathStep * rev_step, dBGraph * db_graph);
            db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);

            //sanity check
            if (next_step.node == NULL) {
                fprintf(stderr, "[db_graph_get_perfect_path_with_first_edge_all_colours] dB_graph: didnt find node in hash table: %s %c %s\n",
                binary_kmer_to_seq(element_get_kmer(current_step.node), db_graph->kmer_size, tmp_seq),
                binary_nucleotide_to_char(current_step.label),
                current_step.orientation == forward ? "forward" : "reverse");
                exit(1);
            }
        } else {
            // If not added (we ran out of space), then remove the last node and make it's label undefined...
            pathStep ps;
            path_get_last_step(&ps, path);
            ps.label = Undefined;
            path_remove_last(path);
            path_add_node(&ps, path);
            if (DEBUG) 
            {
                printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Trying to correct last label failed.\n");
            }
        }
    }
    while (added && !path_is_cycle(path) &&	//loop
           db_node_has_precisely_one_edge_all_colours(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&	//multiple entries
           db_node_has_precisely_one_edge_all_colours(next_step.node, next_step.orientation, &next_step.label));	//has one next edge only

    if (current_step.node != NULL) {
        if (  !(db_node_has_precisely_one_edge_all_colours(next_step.node, opposite_orientation(next_step.orientation), &nucleotide2) &&
                db_node_has_precisely_one_edge_all_colours(next_step.node, next_step.orientation, &next_step.label)))
        {
            next_step.label = Undefined;
            if (DEBUG) 
            {
                printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Changing last label to 'N'.\n");
            }
        }
        added = path_add_node(&next_step, path);
        if (added) {
            node_action(next_step.node);
        }
    }

    if (DEBUG) {
        if (next_step.node == NULL) {
            printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Next node is null! \n");
        } else {
            printf("[db_graph_get_perfect_path_with_first_edge_all_colours] Last node in path: %s %i length: %i\n",
            binary_kmer_to_seq(element_get_kmer(next_step.node), db_graph->kmer_size, tmp_seq), db_node_get_edges_all_colours(next_step.node), path->length);
        }
    }

    return path_get_edges_count(path);

}


void db_graph_print_status(dBGraph * db_graph)
{
	log_and_screen_printf("dBGraph:\n unique kmers: %'lld\n", db_graph->unique_kmers);
	log_and_screen_printf(" Capacity: %'lld \n", (db_graph->bucket_size * db_graph->number_buckets));
	float cap = (float)db_graph->unique_kmers /(float) (db_graph->bucket_size * db_graph->number_buckets) * 100;
	log_and_screen_printf(" Occupied: %f%%\n", cap);
}


int db_graph_get_perfect_path(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path)
{

	//sanity check
	if (node == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null node\n");
		exit(1);
	}
	if (path == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null path\n");
		exit(1);
	}
	if (db_graph == NULL) {
		printf("[db_graph_get_perfect_path] can't pass a null db_graph\n");
		exit(1);
	}
	perfect_path_get_path(node, orientation, node_action, db_graph, path);

	return path_get_edges_count(path);
}



// limit is the max length
// min_coverage, max_coverage and avg_coveragte refer to the internal nodes
int db_graph_supernode(dBNode * node, void (*node_action) (dBNode * node), Path * path, dBGraph * db_graph)
{
	return db_graph_get_perfect_path(node, undefined, node_action, db_graph, path);

}

void printNode(dBNode * dbn, short int kmerSize)
{
	if (dbn != NULL) {
		char seq1[kmerSize + 1];
		printf("MemAdd:\%p\n", dbn);
		printf("Node:\t%s \n",
		       binary_kmer_to_seq(&(dbn->kmer), kmerSize, seq1));
		printf("Flags:\t%x \n", dbn->flags);
		printf("Edges:\t%x \n", db_node_get_edges(dbn));
        int i;
        for (i = 0; i < NUMBER_OF_COLOURS; i++) {
            printf("Col %d:\t%d \n", i, dbn->coverage[i]);
        }
	} else {
		printf("NULL node");
	}
	return;
}


void db_graph_add_node_action(WalkingFunctions * wf, void (*node_action)(dBNode * node)){

    assert(wf);

    if (node_action == NULL) {
        return;
    }

    if(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->node_callbacks.callback[wf->node_callbacks.used] = NULL;

    wf->node_callbacks.callback[wf->node_callbacks.used++] = (void (*) ) node_action; //Make sure the action is not in the array yet.


}

void db_graph_add_path_callback(WalkingFunctions * wf, void (*path_callback)(Path * path)){

    assert(wf != NULL);
    if (path_callback == NULL) {
        return;
    }

    if(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->path_callbacks.args[wf->path_callbacks.used] = NULL;
    wf->path_callbacks.callback[wf->path_callbacks.used++] = (void (*) ) path_callback;

}

boolean db_graph_remove_path_callback(WalkingFunctions * wf, void * funct){
    int used_orig = wf->path_callbacks.used;
    int removed = 0;
    int i;

    for (i = 0; i < used_orig; i++) {
        if(wf->path_callbacks.callback[i] == funct) {
            removed++;
        }
        if (i + removed < MAX_STACKED_FUNCTIONS) {
            wf->path_callbacks.args[i] = wf->path_callbacks.args[i + removed];
            wf->path_callbacks.callback[i] = wf->path_callbacks.callback[i + removed];

        }else{
            wf->path_callbacks.args[i] = NULL;
            wf->path_callbacks.callback[i] = NULL;
        }
    }
    wf->path_callbacks.used -= removed;
    assert( wf->path_callbacks.used >= 0);
    return  removed > 0;
}

void db_graph_add_step_action(WalkingFunctions * wf, void (*step_action)(pathStep * ps)){

    assert(wf != NULL);
    if (step_action == NULL) {
        return;
    }

    if(wf->step_actions.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->step_actions.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->step_actions.callback[wf->step_actions.used] = NULL;
    wf->step_actions.callback[wf->step_actions.used++] = (void (*) ) step_action; //Make sure the action is not in the array yet.

}

void db_graph_add_node_action_with_args(WalkingFunctions * wf, void (*node_action)(dBNode * node, void * arg), void * args){

    assert(wf);

    if (node_action == NULL) {
        return;
    }

    if(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->node_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->node_callbacks.callback[wf->node_callbacks.used] = args;
    wf->node_callbacks.callback[wf->node_callbacks.used++] = node_action; //Make sure the action is not in the array yet.


}

void db_graph_add_path_callback_with_args(WalkingFunctions * wf, void (*path_callback)(Path * path, void * arg), void * args){

    assert(wf != NULL);
    if (path_callback == NULL) {
        return;
    }

    if(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->path_callbacks.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->path_callbacks.args[wf->path_callbacks.used] = args;
    wf->path_callbacks.callback[wf->path_callbacks.used++] = path_callback;

}

void db_graph_add_step_action_with_args(WalkingFunctions * wf, void (*step_action)(pathStep * ps, void * arg), void * args){

    assert(wf != NULL);
    if (step_action == NULL) {
        return;
    }

    if(wf->step_actions.used >= MAX_STACKED_FUNCTIONS){
        fprintf(stderr, "Trying to overload more functions than we have capacity.");
        assert(wf->step_actions.used >= MAX_STACKED_FUNCTIONS);
        exit(-1);
    }

    wf->step_actions.callback[wf->step_actions.used] = args;
    wf->step_actions.callback[wf->step_actions.used++] = step_action; //Make sure the action is not in the array yet.

}


static void execute_path_callbacks(Path * p, PathCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }
}


static void execute_path_step_callbacks(pathStep * p, PathStepActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(p);
        }else{
            f(p, callbacks->args[i]);
        }
    }
}

static void execute_node_callbacks(dBNode * n, NodeActionCallbackArray * callbacks){
    int i;
    for (i = 0; i < callbacks->used; i++) {
        void  (* f)() = callbacks->callback[i];
        if (callbacks->args[i] == NULL) {
            f(n);
        }else{
            f(n, callbacks->args[i]);
        }
    }
}

// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip(dBNode * node, int limit, void (*node_action) (dBNode * node), dBGraph * db_graph)
{
	int length_tip = 0;

	length_tip = db_graph_db_node_clip_tip_with_orientation(node, forward, limit, node_action, db_graph);

	if (length_tip == 0) {
        //printf("PLEASE NOTE: I do get here!\n");
		length_tip = db_graph_db_node_clip_tip_with_orientation(node, reverse, limit, node_action, db_graph);
	}

	return length_tip;
}

void db_graph_print_supernodes(char *filename, int max_length, boolean with_coverages, dBGraph * db_graph)
{
	FILE *fout;
	FILE *fout_cov;
	fout = fopen(filename, "w");

	if (with_coverages) {
		char filename_cov[strlen(filename) + 10];
		sprintf(filename_cov, "%s_cov", filename);
		fout_cov = fopen(filename_cov, "w");
	}

	int count_nodes = 0;

	Path *path = path_new(max_length, db_graph->kmer_size);
	if (!path) {
		fprintf(stderr, "\n[db_graph_print_supernodes] Can't get memory for new path.\n\n");
		exit(-1);
	}

	long long count_kmers = 0;
	long long count_sing = 0;

	void print_supernode(dBNode * node) {

		count_kmers++;
		if (db_node_check_flag_visited(node) == false) {
			int length = db_graph_supernode(node,
                                            &db_node_action_set_flag_visited,
                                            path, db_graph);

			if (length > 0) {
				if (with_coverages) {
					path_to_coverage(path, fout_cov);
				}

				path_to_fasta(path, fout);
				path_increase_id(path);
				if (length == max_length) {
					printf("contig length equals max length [%i] for node_%i\n", max_length, count_nodes);
				}
				count_nodes++;
				path_reset(path);

			}
		}
	}

	hash_table_traverse(&print_supernode, db_graph);
	log_and_screen_printf("%'qd nodes visited [%'qd singletons]\n", count_kmers, count_sing);

	path_destroy(path);
	fclose(fout);
	if (with_coverages) {
		fclose(fout_cov);
	}
}

typedef struct {
    long long tmp_cov[NUMBER_OF_COLOURS];
    dBGraph * db_graph;
} db_graph_stats_args;

static void get_cov(dBNode * node, void * dbsa){
    db_graph_stats_args * args = (db_graph_stats_args * ) dbsa;
    int i;
    boolean all_common;
    long long curr_cov;

    if(db_node_check_flag_not_pruned(node)){
        all_common = true;
        for(i = 0; i < NUMBER_OF_COLOURS; i++){
            curr_cov = element_get_coverage_by_colour(node, i);
            if(curr_cov){
                args->tmp_cov[i] += curr_cov;
                args->db_graph->colour_kmers[i]++;
            }else{
                all_common = false;
            }
        }
        if(all_common){
            args->db_graph->common_kmers_in_all_colours ++;
        }
    }
}


void db_graph_calculate_stats(dBGraph * db_graph){
    //TODO: This should be easy to multithread...

    if(db_graph->calculated == true){
        return;
    }
    double avg_cov = 0;
    //	long long tmp_cov[NUMBER_OF_COLOURS];
    short i;
    //   long long curr_cov;
    db_graph_stats_args ** tmp_args = calloc(1, sizeof(db_graph_stats_args *));
    tmp_args[0] = calloc(1, sizeof(db_graph_stats_args));

    for(i = 0; i < NUMBER_OF_COLOURS; i++){

        tmp_args[0]->tmp_cov[i] = 0;
        db_graph->colour_kmers[i] = 0;
    }
    tmp_args[0]->db_graph = db_graph;

		log_and_screen_printf("Calculating graph stats...\n");

    db_graph->common_kmers_in_all_colours = 0;

    hash_table_traverse_with_args(&get_cov,(void **) tmp_args, db_graph);
    for(i = 0; i < NUMBER_OF_COLOURS; i++){
	    avg_cov = (double)tmp_args[0]->tmp_cov[i] / (double)db_graph->colour_kmers[i];
	    db_graph->average_coverage[i] = avg_cov;

    }
    db_graph->calculated = true;
    free(tmp_args[0]);
    free(tmp_args);

}

//This just gets the average colour of the first colour
double db_graph_get_average_coverage(dBGraph * db_graph){
	if (db_graph->calculated == false) {//We already have calculated the stats, there is no need
        db_graph_calculate_stats(db_graph);
    }

	return db_graph->average_coverage[0];
}

double db_graph_get_average_coverage_by_colour(short colour, dBGraph * db_graph){
	if (db_graph->calculated == false) {//We already have calculated the stats, there is no need
        db_graph_calculate_stats(db_graph);
    }

	return db_graph->average_coverage[colour];
}


void db_graph_print_coverage(FILE * out, dBGraph * db_graph)
{
	long long max_cov = 0;
	long long *count;
	long long tmp_cov;
	fprintf(stdout, "Getting all the coverages...\n");
	void f_max_cov(dBNode * node) {
		tmp_cov = element_get_coverage_all_colours(node);
		if (max_cov < tmp_cov) {
			max_cov = tmp_cov;
		}
	}

	hash_table_traverse(&f_max_cov, db_graph);
	count = calloc(max_cov, sizeof(long long));

	if(count == NULL){
		fprintf(stderr, "[db_graph_print_coverage] WARNING! Unable to allocate memory to count the coverages (%lli)\n", tmp_cov);
		return;
	}

	void sum_cov(dBNode * node) {
		tmp_cov = element_get_coverage_all_colours(node);
		count[tmp_cov - 1]++;
	}
	hash_table_traverse(&sum_cov, db_graph);
	fprintf(out, "Kmer\tCoverage\n");

	for (tmp_cov = 0; tmp_cov < max_cov; tmp_cov++) {
		fprintf(out, "%lld\t%lld\n", tmp_cov + 1, count[tmp_cov]);
	}

	fprintf(stdout, "done \n");
	free(count);
}

void db_graph_print_kmer_coverage(FILE * out, dBGraph * db_graph)
{

	char *tmp_kmer = calloc(db_graph->kmer_size + 1, sizeof(char));

	fprintf(stdout, "Printing kmer and coverage..");
	fprintf(out, "Kmer\tCoverage\n");

	void print_cov(dBNode * node) {
		//BinaryKmer * bk = (BinaryKmer*)node->kmer;
		binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_kmer);
		fprintf(out, "%s\t%d\n", tmp_kmer, element_get_coverage_all_colours(node));
	}

	hash_table_traverse(&print_cov, db_graph);
	fprintf(stdout, "done \n");

	free(tmp_kmer);
}

int db_graph_clip_tips(int threshold, dBGraph * db_graph)
{
	db_graph_reset_flags(db_graph);
	long long tips_clipped = 0;
	long long nodes_clipped = 0;

	void local_node_prune(dBNode * no) {
		//db_node_action_set_flag_pruned(no);
        cleaning_prune_db_node(no, db_graph);
		nodes_clipped++;
	}

	void clip_tips(dBNode * node) {
		if (!db_node_check_for_any_flag(node, PRUNED | VISITED)) {
            if (db_graph_db_node_clip_tip(node, threshold, &local_node_prune, db_graph) > 0) {
                tips_clipped ++;
            }
		}
	}

#ifdef THREADS
	hash_table_threaded_traverse(&clip_tips, db_graph);
#else
	hash_table_traverse(&clip_tips, db_graph);
#endif
	log_and_screen_printf("clip_tip removed nodes %'lld\n", nodes_clipped);

	return tips_clipped;
}


// Should combine with below function
void db_graph_dump_binary(char *filename, boolean(*condition) (dBNode * node), dBGraph * db_graph)
{
	FILE *fout;		//binary output
	fout = fopen(filename, "w");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}

	int mean_read_len=0;//Mario - you can plumb this in if you want. See cortex_var/core/graph_info.c, and how the GraphInfo object is used in my file_reader * cprtex_var.c,
	// for how I did it.
	long long total_seq=0;
	print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);


	long long count = 0;
	//routine to dump graph
	void print_node_binary(dBNode * node) {
		if (condition(node)) {
			count++;
			db_node_print_binary(fout, node, db_graph->kmer_size);
		}
	}

	hash_table_traverse(&print_node_binary, db_graph);
	fclose(fout);

	log_and_screen_printf("%'lld kmers dumped\n", count);
}

// Should combine with above function
void db_graph_dump_binary_by_colour(char *filename, boolean(*condition) (dBNode * node), short colour, dBGraph * db_graph)
{
	FILE *fout;		//binary output
	fout = fopen(filename, "w");
	if (fout == NULL) {
		fprintf(stderr, "cannot open %s", filename);
		exit(1);
	}

	int mean_read_len=0;//Mario - you can plumb this in if you want. See cortex_var/core/graph_info.c, and how the GraphInfo object is used in my file_reader * cprtex_var.c,
	// for how I did it.
	long long total_seq=0;
	print_binary_signature(fout,db_graph->kmer_size,1, mean_read_len, total_seq);


	long long count = 0;
	//routine to dump graph
	void print_node_binary(dBNode * node) {
		if ((condition(node)) && (node->coverage[colour] > 0)) {
			count++;
			db_node_print_binary_by_colour(fout, node, colour, db_graph->kmer_size);
		}
	}

	hash_table_traverse(&print_node_binary, db_graph);
	fclose(fout);

	log_and_screen_printf("%qd kmers dumped\n", count);
}





int db_graph_db_node_clip_tip_with_orientation(dBNode * node, Orientation orientation, int limit, void (*node_action) (dBNode * node), dBGraph * db_graph)
{
	Nucleotide nucleotide, reverse_nucleotide;
	int length = 0;
	int i;
	dBNode *nodes[limit];
	Orientation next_orientation;
	dBNode *next_node;
	char seq[db_graph->kmer_size + 1];

	//starting in a blunt end also prevents full loops
	if (db_node_is_blunt_end_all_colours(node, opposite_orientation(orientation))) {
		boolean join_main_trunk = false;

		while (db_node_has_precisely_one_edge_all_colours(node, orientation, &nucleotide)) {
			nodes[length] = node;

			next_node = db_graph_get_next_node(node, orientation, &next_orientation, nucleotide, &reverse_nucleotide, db_graph);

			if (next_node == NULL) {
				printf("dB_graph_db_node_clip_tip_with_orientation: didnt find node in hash table: %s\n",
                       binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, seq));
				exit(1);
			}

			length++;

			if (length > limit) {
				break;
			}
			//want to stop when we join another trunk
			if (!db_node_has_precisely_one_edge_all_colours(next_node, opposite_orientation(next_orientation), &nucleotide) || !db_node_has_precisely_one_edge_all_colours(next_node, next_orientation, &nucleotide)) {
				join_main_trunk = true;
				break;
			}
			//keep track of node

			node = next_node;
			orientation = next_orientation;
		}

		if (!join_main_trunk) {
			length = 0;
		} else {	//clear edges and mark nodes as pruned
			for (i = 0; i < length; i++) {

				if (DEBUG) {
					printf("CLIPPING node: %s\n", binary_kmer_to_seq(element_get_kmer(nodes[i]), db_graph->kmer_size, seq));
				}

				node_action(nodes[i]);
                //perhaps we want to move this inside the node action?
				//(we already did). db_node_reset_edges_all_colours(nodes[i]);
			}

			if (DEBUG) {
				printf("RESET %c BACK\n", binary_nucleotide_to_char(reverse_nucleotide));
			}
			db_node_reset_edge_all_colours(next_node, opposite_orientation(next_orientation), reverse_nucleotide);//Why is this necessary?
		}

	}

	return length;
}


//it doesn't check that it is a valid arrow -- it just assumes the arrow is fine
dBNode *db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, Orientation * next_orientation,
                               Nucleotide edge, Nucleotide * reverse_edge, dBGraph * db_graph)
{
	if (edge == Undefined) {
		return NULL;	//If it is undefined, the next node is null...
	}

	BinaryKmer local_copy_of_kmer;
	binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);

	BinaryKmer tmp_kmer;
	dBNode *next_node = NULL;

	// after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
	BinaryKmer *rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer);

	if (current_orientation == reverse) {
		*reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
		binary_kmer_assignment_operator(local_copy_of_kmer, *rev_kmer);
	} else {
		*reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
	}

	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);

	//get node from table
	next_node = hash_table_find(element_get_key(&local_copy_of_kmer, db_graph->kmer_size, &tmp_kmer), db_graph);

	if (next_node != NULL) {
		*next_orientation = db_node_get_orientation(&local_copy_of_kmer, next_node, db_graph->kmer_size);
	}

	return next_node;
}

int db_graph_generic_walk(pathStep * first_step, Path * path, WalkingFunctions * functions, dBGraph * db_graph)
{
    //dBNode * node = first_step->node;
    //sanity check
    if (first_step == NULL) {
        printf("[db_graph_generic_walk] can't pass a null node\n");
        assert(0);
        exit(1);
    }

    if (path == NULL) {
        printf("[db_graph_generic_walk] can't pass a null path\n");
        assert(0);
        exit(1);
    }

    if (functions == NULL) {
        printf("[db_graph_generic_walk] can't pass  null functions\n");
        assert(0);
        exit(1);
    }

    if (db_graph == NULL) {
        printf("[db_graph_generic_walk] can't pass  null db_graph\n");
        assert(0);
        exit(1);
    }

    pathStep current_step, next_step, rev_step;
    path_step_assign(&next_step, first_step);
    next_step.path = path; //This way, on all the assignments we keep the pointer to the path that was sent to the function originally.

    functions->get_starting_step(&next_step, db_graph);
    //boolean  try = true;
    int count = 0;
    boolean walked;
    boolean added;

    do {
        walked = false;
        do {
            added = path_add_node(&next_step, path);
            path_step_assign(&current_step, &next_step);
            //try = false;
            functions->pre_step_action(&current_step);

            added = false; // redundant?
            if (current_step.label != Undefined) {
                functions->get_next_step(&current_step, &next_step, &rev_step, db_graph);
            }

            if (added == false) {
                // Do something?
            }

            functions->step_action(&current_step);
            execute_path_step_callbacks(&current_step,&functions->step_actions);
            execute_node_callbacks(current_step.node, &functions->node_callbacks);
        }
        while (functions->continue_traversing(&current_step, &next_step, &rev_step, path,db_graph));

        if (path->length > 0) {
                walked = true;
        }

        if (walked) {
            count++;
            functions->output_callback(path);
            execute_path_callbacks(path, &functions->path_callbacks);
        }
        do {
            pathStep ps;

            path_get_last_step(&ps, path);
            if(ps.node != NULL){
                functions->post_step_action(&ps);
                path_remove_last(path);
            }
        }
        while (functions->continue_backwards(path, db_graph));

        if(path_get_length(path) > 0){//This is to enable the backtracking.
            path_get_last_step(&current_step, path);//The continue_backwards should fix the last step
            functions->get_next_step(&current_step,&next_step, &rev_step,db_graph);
        }
    }
    while (path_get_length(path) > 0);
    return count;
}

void db_graph_write_graphviz_file(char *filename, dBGraph * db_graph)
{
	FILE *fp;

	void print_graphviz_colours(dBNode * node) {
            if (node != NULL && element_get_coverage_all_colours(node) > 0) 
            {
                BinaryKmer tmp, co;
                short kmer_size = db_graph->kmer_size;
                binary_kmer_assignment_operator(tmp, node->kmer);
                binary_kmer_reverse_complement(&tmp, kmer_size, &co);
                char seq[kmer_size + 1], seqNext[kmer_size + 1], seq1[kmer_size + 1];
                binary_kmer_to_seq(&tmp, kmer_size, seq1);
                char *print = db_node_check_for_any_flag(node, STARTING_FORWARD  | BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | END_NODE_FORWARD | END_NODE_REVERSE | X_NODE) ? "ellipse" : "ellipse";
                char *node_colour = "black";

                if (db_node_check_for_any_flag(node, BRANCH_NODE_FORWARD | BRANCH_NODE_REVERSE | X_NODE)) {
                         print = "circle";
                }

                char reverse_seq1[kmer_size + 1];
                char label_string[512];

                seq_reverse_complement(seq1, kmer_size, reverse_seq1);
                //sprintf(label_string, "\"%s\\n%s\\nC0=%d C1=%d\"", seq1, reverse_seq1, node->coverage[0], node->coverage[1]);

                sprintf(label_string,
                        "<<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\">"
                        "<tr><td><font color=\"blue\">%s</font></td></tr>"
                        "<tr><td><font color=\"red\">%s</font></td></tr>"
                        "<tr><td><font color=\"black\">%d</font> </td></tr>"
                        "</table>>", seq1, reverse_seq1,
                        node->coverage[0]
                        );

                fprintf(fp, "%s [label=%s, shape=%s, color=%s]\n", seq1, label_string, print, node_colour);

                binary_kmer_to_seq(&tmp, kmer_size, seq);
                binary_kmer_left_shift(&tmp, 2, kmer_size);
                Orientation o = forward;
                Edges all_edges = db_node_get_edges_for_orientation_all_colours(node, o);
                Edges ea[NUMBER_OF_COLOURS];
                int i;

                if (all_edges < 0) {
                    // Do something???
                }

                for (i = 0; i < NUMBER_OF_COLOURS; i++) {
                        ea[i] = db_node_get_edges_for_orientation_by_colour(node, o, i);
                }

                void hasEdge(Nucleotide n) {
                    BinaryKmer bk;
                    Key k = &bk;
                    boolean colour0edge = (((ea[0] >> n) & 1) == 1);
                    boolean colour1edge = false;

                    if (colour0edge || colour1edge) {
                            char *labelcolour = (colour0edge && colour1edge) ? "black":(colour0edge ? "green":"orange");
                            binary_kmer_modify_base(&tmp, n, kmer_size, 0);
                            binary_kmer_to_seq(element_get_key(&tmp, kmer_size, k), kmer_size, seqNext);
                            fprintf(fp, "%s -> %s [ label=<<font color=\"%s\">%c</font>>, color=\"%s\"];\n", seq, seqNext, labelcolour, binary_nucleotide_to_char(n), o == forward ? "blue" : "red");
                    }
                }

                nucleotide_iterator(&hasEdge);
                binary_kmer_assignment_operator(tmp, co);
                binary_kmer_left_shift(&tmp, 2, kmer_size);

                o = reverse;
                all_edges = db_node_get_edges_for_orientation_all_colours(node, o);
                for (i = 0; i < NUMBER_OF_COLOURS; i++) {
                    ea[i] = db_node_get_edges_for_orientation_by_colour(node, o, i);
                }
                nucleotide_iterator(&hasEdge);
            }
	}

	//if (db_graph->unique_kmers < 450) {
		if (filename) 
                {
                    fp = fopen(filename, "w");
		} 
                else 
                {
                    fp = stderr;
		}
		if (fp) 
                {
			fprintf(fp, "digraph finite_state_machine {\n");
			fprintf(fp, "\trankdir=LR;\n");
			fprintf(fp, "\tsize=\"12,8\";\n");
			fprintf(fp, "\tfontsize=18;\n");
			fprintf(fp, "\tnode [shape = circle];\n");
			hash_table_traverse(&print_graphviz_colours, db_graph);
			fprintf(fp, "}\n");
			fclose(fp);
                        printf("\nWritten graphviz file %s.\n\n", filename);
		}
	//} else {
	//	printf("\nNote: Too many kmers to output graphviz file.\n\n");
	//}
}

// Build array of neighbouring nodes
int db_graph_get_neighbouring_nodes_all_colours(dBNode* start, Orientation orientation, dBNode* neighbours[], Nucleotide labels[], dBGraph * db_graph) {
    int num_neighbours = 0;
    Nucleotide nucleotide;
    for (nucleotide = Adenine; nucleotide < Undefined;
         nucleotide++) {
        if (db_node_edge_exist_any_colour(start, nucleotide, orientation)) {
            pathStep current_step, next_step, rev_step;
            current_step.node = start;
            current_step.label = nucleotide;
            current_step.orientation = orientation;
            current_step.flags = 0;
            db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);

            if (next_step.node != NULL) {
                labels[num_neighbours] = nucleotide;
                neighbours[num_neighbours] = next_step.node;
                num_neighbours++;
            }
        }
    }


    return num_neighbours;
}

Nucleotide db_graph_get_best_next_step_nucleotide(dBNode * from, dBNode * previous,  Orientation orientation, dBGraph * db_graph ){
    assert(from != NULL);
    assert(previous != NULL);
    Nucleotide best = Undefined;

    //TODO: Make this with the coverage.

    return best;
}

pathStep get_path_to_junction(pathStep* first_step, Path* new_path, dBGraph* db_graph, Queue* node_queue)
{ 
    // don't add first step
    pathStep return_step;
    return_step.node = NULL;
    return_step.label = Undefined;
    
    Nucleotide current_reverse;
    Orientation current_orientation;
    dBNode* current_node = db_graph_get_next_node(first_step->node, first_step->orientation, &current_orientation, first_step->label, 
                                                &current_reverse, db_graph);
    
    int num_edges = db_node_edges_count_all_colours(current_node, current_orientation);
    while(num_edges < 2)
    {
        if(current_orientation == forward)
        {
            if(flags_check_for_flag(VISITED_FORWARD, &(current_node->flags)))
            {
                //log_printf("[get_path_to_junction] Already visited this node. Returning NULL.\n");
                return return_step;
            }
            if(flags_check_for_flag(CURRENT_PATH_FORWARD, &(current_node->flags)))
            {
                break;
            }
            flags_action_set_flag(VISITED_FORWARD, &current_node->flags);
        }
        else
        {
            if(flags_check_for_flag(VISITED_REVERSE, &(current_node->flags)))
            {
                //log_printf("[get_path_to_junction] Already visited this node. Returning NULL.\n");
                return return_step;
            }
            if(flags_check_for_flag(CURRENT_PATH_REVERSE, &(current_node->flags)))
            {
                break;
            }
            flags_action_set_flag(VISITED_REVERSE, &current_node->flags);
        }
        queue_push(node_queue, current_node);
                
        if(num_edges == 0)
        {
/*
            char kmer_string[db_graph->kmer_size + 1];
            BinaryKmer kmer; 
            binary_kmer_assignment_operator(kmer, *element_get_kmer(current_node));
            binary_kmer_to_seq(&kmer, db_graph->kmer_size, kmer_string); 
            kmer_string[db_graph->kmer_size] = '\0';
            log_printf("[get_path_to_junction] Found dead end at node %s. Returning NULL.\n", kmer_string);
*/
            return return_step;
        }
        
        pathStep next_step;
        next_step.node = current_node;
        next_step.orientation = current_orientation;
        next_step.label = Undefined;
        Edges edges = db_node_get_edges_for_orientation_all_colours(current_node, current_orientation);
        for(int n = 0; n < 4; n++) 
        {
            if ((edges & 1) == 1) 
            {
                next_step.label = n;
                break;
            }
            edges >>= 1;
        }
        
        assert(next_step.label != Undefined);
  
        boolean added = path_add_node(&next_step, new_path);    
        if(!added)
        {
            log_printf("[get_path_to_junction] Max path length reached, could not add node.\n");
            return return_step;
        }
                  
        current_node = db_graph_get_next_node(next_step.node, next_step.orientation, &current_orientation, next_step.label, 
                                                &current_reverse, db_graph);  
        num_edges = db_node_edges_count_all_colours(current_node, current_orientation);
    }
    
    //log_printf("[get_path_to_junction] Returning junction node.\n");
    return_step.node = current_node;
    return_step.orientation = current_orientation;
    path_add_node(&return_step, new_path);
    
    // mark as visited
    db_node_action_set_flag_visited_with_orientation(current_node, current_orientation);
    queue_push(node_queue, current_node);
    
    return return_step;
}

#define BUBBLE_QUEUE_SIZE 200000 // 10000000
#define MAX_EXPLORE_BUBBLE_LENGTH 2000000
#define MAX_EXPLORE_BUBBLE_JUNCTIONS 2000
// breadth first search
pathStep db_graph_search_for_bubble2(Path* main_path, pathStep* first_step, Path** new_path_ptr, dBGraph* db_graph)
{
    assert(new_path_ptr != NULL && *new_path_ptr == NULL);
    
    dBNode* node = first_step->node;
    Orientation orientation = first_step->orientation;
    pathStep join_step;
    join_step.node = NULL;
    join_step.orientation = undefined;
    join_step.label = undefined;
    
    boolean joined = false;
     
    Queue* step_queue = queue_new(MAX_EXPLORE_BUBBLE_JUNCTIONS, sizeof(pathStep*));
    Queue* node_queue_for_flags = node_queue_new(BUBBLE_QUEUE_SIZE);
    pathStep* new_step = malloc(sizeof(pathStep));
    new_step->node = first_step->node;
    new_step->orientation = first_step->orientation;
    queue_push_step(step_queue, new_step);
   
    typedef struct
    {
        dBNode* child_node;
        dBNode* parent_node;
        Path* path;
    } parent_child;
    Queue* parent_child_list = queue_new(MAX_EXPLORE_BUBBLE_JUNCTIONS, sizeof(parent_child*));
    
    if(!step_queue || !node_queue_for_flags || !parent_child_list)
    {
        log_and_screen_printf("[db_graph_search_for_bubble2] Couldn't get memory for queue. Exiting...\n");
        exit(-1);
    }
   
    // mark the path after starting node
    boolean mark = false;
    for(int i = 0; i < main_path->length; i++)
    {
        if(mark)
        {
            if(main_path->orientations[i] == forward)
            {
                flags_action_set_flag(CURRENT_PATH_FORWARD, &main_path->nodes[i]->flags);
            }
            else
            {
                flags_action_set_flag(CURRENT_PATH_REVERSE, &main_path->nodes[i]->flags);
            }
        }
        
        if(main_path->nodes[i] == node && main_path->orientations[i] == orientation)
        {
            mark = true;
        }
    }
    
    // function to walk perfect path from node/orientation/n
    // adds node at end of path to queue if unvisited
    void walk_if_exists(Nucleotide n) 
    {
        if(node == first_step->node && n != first_step->label)
        {
            return;
        }
     
        if (db_node_edge_exist_any_colour(node, n, orientation)) 
        {
            // Get first node along this edge and check we've not already visited it...
            Orientation next_orientation;
            Nucleotide reverse_nucleotide;
            dBNode* next_node;
            next_node = db_graph_get_next_node(node, orientation, &next_orientation, n, &reverse_nucleotide, db_graph);

            // If not already visited the first node, walk it...
            if (!db_node_check_flag_visited_with_orientation(next_node, next_orientation)) 
            {
                pathStep new_first_step;

                // Get path
                new_first_step.node = node;
                new_first_step.orientation = orientation;
                new_first_step.label = n;
                new_first_step.flags = 0;
                Path* new_path = path_new(MAX_EXPLORE_BUBBLE_LENGTH, db_graph->kmer_size);
                if (!new_path) 
                {
                    log_and_screen_printf("[db_graph_search_for_bubble2] ERROR: Not enough memory to allocate new path.\n");
                    exit(-1);
                }

                db_graph_get_perfect_path_with_first_edge_all_colours(&new_first_step, &db_node_action_do_nothing, new_path, db_graph);
                // Now go through all nodes, look for one that joins the main path
                for (int i=1; i<new_path->length; i++) 
                {
                    dBNode* current_node = new_path->nodes[i];
                    Orientation current_orientation = new_path->orientations[i];
                    db_node_action_set_flag_visited_with_orientation(current_node, current_orientation);
                    queue_push(node_queue_for_flags, current_node);
                    if( (current_orientation == forward && flags_check_for_flag(CURRENT_PATH_FORWARD, &(current_node->flags))) ||
                        (current_orientation == reverse && flags_check_for_flag(CURRENT_PATH_REVERSE, &(current_node->flags))) )
                    {
                        join_step.node = current_node;
                        join_step.orientation = current_orientation;
                        joined = true;
                        
                        parent_child* pc = malloc(sizeof(parent_child));
                        pc->child_node = current_node;
                        pc->parent_node = node;
                        pc->path = path_new(i+1, db_graph->kmer_size);
                        path_copy_subpath(pc->path, new_path, 0, i+1);
                        path_destroy(new_path);
                        queue_push(parent_child_list, pc);
                        
                        return; 
                    }                  
                }
                
                // Add end node to list of nodes to visit
                assert(!joined);             
                dBNode* end_node = new_path->nodes[new_path->length-1];
                Orientation end_orientation = new_path->orientations[new_path->length-1];
                if (!db_node_check_flag_visited_with_orientation(end_node, end_orientation)) 
                {
                    if (!db_node_is_blunt_end_all_colours(end_node, end_orientation)) 
                    {
                        if(element_get_coverage_all_colours(end_node) > db_graph->path_coverage_minimum)
                        {
                            pathStep* next_step = malloc(sizeof(pathStep));
                            next_step->node = end_node;
                            next_step->orientation = end_orientation;
                            if (queue_push_step(step_queue, next_step) == NULL) 
                            {
                                log_and_screen_printf("Queue too large. Ending.\n");
                                exit(1);
                            }

                            parent_child* pc = malloc(sizeof(parent_child));
                            pc->child_node = end_node;
                            pc->parent_node = node;
                            pc->path = new_path;
                            queue_push(parent_child_list, pc);
                            return;
                        }
                    }
                }                
                path_destroy(new_path);              
            }
        }
    }
    
    while (step_queue->number_of_items > 0 && !joined) 
    {
        // Take top node from list
        pathStep* next_step = queue_pop(step_queue);
        node = next_step->node;
        orientation = next_step->orientation;
        // Look at all paths out from here
        nucleotide_iterator(&walk_if_exists);
        
        free(next_step);
    }

    // free remaining pathSteps from step_queue
    while (step_queue->number_of_items > 0)
    {
        pathStep* next_step = queue_pop(step_queue);
        free(next_step);
    }
    queue_free(step_queue);
    
    PathArray* pa = path_array_new(10);
    if(joined)
    {
        // Now backtrace to create bubble path
        dBNode* current_child = join_step.node;
        assert(parent_child_list->number_of_items > 0);
        do
        {
            for(int i = 0; i < parent_child_list->number_of_items; i++)
            {
                parent_child* pc = parent_child_list->items[i];
                if(pc->child_node == current_child)
                {
                    path_array_add_path(pc->path, pa);
                    current_child = pc->parent_node;
                    pc->path = NULL;
                    break;
                }
            }
        }
        while(current_child != first_step->node);

        Path* new_path = path_array_merge_to_path(pa, true, db_graph);
        *new_path_ptr = new_path;
        
/*
        assert(join_step.node != NULL);
        assert(new_path->nodes[new_path->length - 1] == join_step.node);
        assert(new_path->orientations[new_path->length - 1] == join_step.orientation);
        log_printf("[db_graph_search_for_bubble2] Found bubble with sequence %s\n", new_path->seq);
        log_printf("[db_graph_search_for_bubble2] and length %i\n", new_path->length);
*/
    }
    
    // unmark the path
    for(int i = 0; i < main_path->length; i++)
    {
        flags_action_unset_flag(CURRENT_PATH_FORWARD, &main_path->nodes[i]->flags);
        flags_action_unset_flag(CURRENT_PATH_REVERSE, &main_path->nodes[i]->flags);
    }
    
    // unmark the VISITED_FORWARD/REVERSE
    while(node_queue_for_flags->number_of_items > 0) 
    {
        dBNode* queue_node = (dBNode*)queue_pop(node_queue_for_flags);
        db_node_action_unset_flag_visited_forward_reverse(queue_node);
    }
    queue_free(node_queue_for_flags);
    
    // free all the paths
    for(int i = 0; i < parent_child_list->number_of_items; i++)
    {
        parent_child* pc = parent_child_list->items[i];
        if(pc->path)
        {
            path_destroy(pc->path);
        }
    }
    queue_free(parent_child_list);
    path_array_destroy(pa);
    return join_step;
}

// iterative depth first search choosing path with max coverage at each step.
pathStep db_graph_search_for_bubble(Path* main_path, pathStep* first_step, Path** new_path_ptr, dBGraph* db_graph)
{
    assert(new_path_ptr != NULL && *new_path_ptr == NULL);
    
    typedef struct
    {
        Nucleotide base;
        uint32_t coverage;
    } base_coverage_pair;

    int compare(const void* a, const void* b)
    {
        int x = ((base_coverage_pair *)a)->coverage; 
        int y = ((base_coverage_pair *)b)->coverage;
        return y - x;
    }
    
    double avg_coverage;
    uint32_t min_coverage;
    uint32_t max_coverage;
    path_get_statistics(&avg_coverage, &min_coverage, &max_coverage, main_path);
    
    int max_path_size = main_path->length;
    
    // mark the main path from the first step
    boolean mark = false;
    for(int i = 0; i < main_path->length; i++)
    {
        if(mark)
        {
            if(main_path->orientations[i] == forward)
            {
                flags_action_set_flag(CURRENT_PATH_FORWARD, &main_path->nodes[i]->flags);
            }
            else
            {
                flags_action_set_flag(CURRENT_PATH_REVERSE, &main_path->nodes[i]->flags);
            }
        }
        if(main_path->nodes[i] == first_step->node)
        {
            mark = true;
        }
    }
    assert(mark);
    // setup array etc.
    PathArray* path_array = path_array_new(2);
    Queue* step_queue = queue_new(200000, sizeof(pathStep*));
    Queue* node_queue_for_flags = node_queue_new(BUBBLE_QUEUE_SIZE);
    
    // This step is freed once it has been popped
    pathStep* new_step = malloc(sizeof(pathStep));
    //TODO: consider writing a copy_step function
    new_step->flags = first_step->flags;
    new_step->label = first_step->label;
    new_step->node = first_step->node;
    new_step->orientation = first_step->orientation;
    new_step->path = first_step->path;
    
    queue_push_step(step_queue, new_step);
    
    // the step to return
    pathStep join_step;
    join_step.node = NULL;
    join_step.orientation = undefined;
    join_step.label = undefined;
       
    boolean joined_path = false;
    boolean add_first_step = true;
    
    while(step_queue->number_of_items > 0)
    {
        pathStep* next_step = queue_pop_step(step_queue);

        // Remove paths from array that don't join this one. This is a depth first search
        // so we will get back to the path that ended in the node in next_step eventually.
        boolean removed = false;
        while(path_array->number_of_paths > 0)
        {
            Path* last_path = path_array_get_last_path(path_array);
            assert(last_path->length > 0);
            if(last_path->nodes[last_path->length - 1] == next_step->node)
            {
                break;
            }
            else
            {
                path_array_remove_last_path(path_array);
                removed = true;
            }
        }
        
        if(path_array->number_of_paths == 0 && removed)
        {
            // we removed all the paths from the array!
            log_printf("[db_graph_search_for_bubble] Path array emptied. Step stack: %i.\n", step_queue->number_of_items);
            free(next_step);
            break;
        }
       
        Path* new_path = path_new(max_path_size, db_graph->kmer_size);
        if(add_first_step)
        {
            path_add_node(first_step, new_path);
            add_first_step = false;
        }
        pathStep junction_step = get_path_to_junction(next_step, new_path, db_graph, node_queue_for_flags);
        if(junction_step.node == NULL)
        {
            // found dead end or revisited node. Go back to top of loop.
            path_destroy(new_path);
            free(next_step); 
            continue;
        }
        if(path_array_get_total_size(path_array) + new_path->length > max_path_size)
        {
            log_printf("[db_graph_search_for_bubble] Warning: Maximum path array total length reached. Ending search for alternative route.\n");
            path_destroy(new_path);
            free(next_step);
            break;
        }
        
        dBNode* junction_node = junction_step.node;
        Orientation junction_orientation = junction_step.orientation;
                
        // Check to see if junction_node is back on the main path
        // if so, break from this loop
        if( (junction_orientation == forward && flags_check_for_flag(CURRENT_PATH_FORWARD, &(junction_node->flags))) ||
            (junction_orientation == reverse && flags_check_for_flag(CURRENT_PATH_REVERSE, &(junction_node->flags))) )
        {

            //add path to path array
            path_array_add_path(new_path, path_array);
            joined_path = true;
            free(next_step);
            break;
        }
       
        // Otherwise order the edges by coverage and add them to the queue (actually a stack) to be explored
        int num_edges = db_node_edges_count_all_colours(junction_node, junction_orientation);
        assert(num_edges > 1);

        // get the coverage for each edge
        base_coverage_pair bcp_array[4];
        for(Nucleotide base = 0; base < 4; base++) 
        {
            bcp_array[base].base = base;
            bcp_array[base].coverage = 0;
            if(db_node_edge_exist_any_colour(junction_node, base, junction_orientation))
            {
                Orientation next_orientation;
                Nucleotide reverse_edge;
                dBNode* next_node = db_graph_get_next_node(junction_node, junction_orientation, &next_orientation, base, &reverse_edge, db_graph);
                if(db_node_check_flag_visited_with_orientation(next_node, next_orientation))
                {
                    continue;
                }
                if(next_node == junction_node && next_orientation == junction_orientation)
                {
                    continue;
                }
                bcp_array[base].coverage = element_get_coverage_all_colours(next_node);                   
            }      
        }
        
        // sort via coverage (highest first)
        qsort((void*)bcp_array, 4, sizeof(base_coverage_pair), compare);

        //TODO: Check that the total length of the path array isn't too big!
        boolean added = false;
        // reverse iterate so that the smallest is added to the stack first
        for(int i = 3; i >= 0; i--) 
        {
            if(bcp_array[i].coverage >= db_graph->path_coverage_minimum)
            {
                Nucleotide base = bcp_array[i].base;

                pathStep* step_for_label = malloc(sizeof(pathStep));
                step_for_label->node = junction_node;
                step_for_label->orientation = junction_orientation;
                step_for_label->label = base;
                if(queue_push_step(step_queue, step_for_label))
                {
                    added = true;
                }
                else
                {
                    log_printf("[db_graph_search_for_bubble] Warning: Could not add step to queue.\n");
                    free(step_for_label);
                    break;
                }
            }
        }

        if(added)
        {
            assert(new_path->length > 0);
            path_array_add_path(new_path, path_array);
        }
        else
        {
            path_destroy(new_path);
        }
        
        free(next_step);
    }
        
    
    if(joined_path)
    {
        Path* alt_path = path_array_merge_to_path(path_array, false, db_graph);
        if(alt_path)
        {
            int end = -1;
            for(int i = alt_path->length - 1; i >= 0; i--)
            {
                if((alt_path->orientations[i] == forward && flags_check_for_flag(CURRENT_PATH_FORWARD, &(alt_path->nodes[i]->flags))) ||
                   (alt_path->orientations[i] == reverse && flags_check_for_flag(CURRENT_PATH_REVERSE, &(alt_path->nodes[i]->flags))) )
                {
                    end = i;
                    break;
                }
            }

            if(end > 0)
            {
                Path* new_path = path_new(end, alt_path->kmer_size);

                path_copy_subpath(new_path, alt_path, 0, end);

                join_step.node = alt_path->nodes[end];
                join_step.orientation = alt_path->orientations[end];            

                 *new_path_ptr = new_path;

                 assert(join_step.node != NULL);
            }
            else
            {
                log_printf("[db_graph_search_for_bubble] Warning: Found alternative path but could not find join. Ignoring...\n");
            }
            path_destroy(alt_path);
        }
    }
    
    // free used memory
    path_array_destroy(path_array);
    queue_free(step_queue);
    
    // unmark the VISITED_FORWARD/REVERSE
    while(node_queue_for_flags->number_of_items > 0) 
    {
        dBNode* queue_node = (dBNode*)queue_pop(node_queue_for_flags);
        db_node_action_unset_flag_visited_forward_reverse(queue_node);
    }
    queue_free(node_queue_for_flags);
    
    // unmark the path
    for(int i = 0; i < main_path->length; i++)
    {
        flags_action_unset_flag(CURRENT_PATH_FORWARD, &main_path->nodes[i]->flags);
        flags_action_unset_flag(CURRENT_PATH_REVERSE, &main_path->nodes[i]->flags);
    }
    
    return join_step;
}

#ifdef DEBUG_CLEANUP
void db_graph_cleanup_graph(dBGraph * db_graph)
{
	Orientation orientation;

	void clean_node(dBNode * node) {
		void check_edge(Nucleotide nucleotide) {
			// Check if edge exists
			if (db_node_edge_exist_any_colour(node, nucleotide, orientation)) {
				// Now see if kmer it leads to exists
				pathStep current_step, next_step, rev_step;
				current_step.node = node;
				current_step.label = nucleotide;
				current_step.orientation = orientation;
				db_graph_get_next_step(&current_step, &next_step, &rev_step, db_graph);
				if (next_step.node == NULL) {
					printf("Missing node\n");
					db_node_reset_edge_all_colours(node, orientation, nucleotide);
				}
			}
		}

		orientation = forward;
		nucleotide_iterator(&check_edge);

		orientation = reverse;
		nucleotide_iterator(&check_edge);
	}

	printf("\nTidying up graph...\n");
	hash_table_traverse(&clean_node, db_graph);
}
#endif
