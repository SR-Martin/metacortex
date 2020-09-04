/*----------------------------------------------------------------------*
 * File:    nodequeue.h                                                 *
 * Purpose: A queue for nodes!                                          *
 * Author:  Richard Leggett                                             *
 * History: 12-Aug-10: RML: Created                                     *
 *----------------------------------------------------------------------*/

#ifndef NODE_QUEUE_H
#define NODE_QUEUE_H

/*----------------------------------------------------------------------*
 * Structures                                                           *
 *----------------------------------------------------------------------*/
typedef struct {
	dBNode* node;
	int depth;
} QueueItem;

typedef struct {
	int max_size;
	int number_of_items;
	int head;
	int tail;
	void** items;
        boolean warn_queue_full;
} Queue;

/*----------------------------------------------------------------------*
 * General purpose queue functions                                      *
 *----------------------------------------------------------------------*/
Queue* queue_new(int n, size_t item_size);
void* queue_push(Queue* q, void* item);
void* queue_pop(Queue* q);
void queue_free(Queue *q);

/*----------------------------------------------------------------------*
 * Queue wrappers for nodes                                             *
 *----------------------------------------------------------------------*/
Queue* node_queue_new(int n);
QueueItem* queue_push_node(Queue* q, dBNode* n, int d);
dBNode* queue_pop_node(Queue* q, int* d);

/*----------------------------------------------------------------------*
 * Queue wrappers for pathSteps                                         *
 *----------------------------------------------------------------------*/
pathStep* queue_push_step(Queue* q, pathStep* step);
pathStep* queue_pop_step(Queue* q);

#endif // NODE_QUEUE_H