/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  This code is based on R. Sedgewick "Algorithms in C    **/
/**  Parts 1-4", 3rd Edition, Addison-Wesley (1998).        **/
/**  *****************************************************  **/
/**  Abstract data type implementation for skip list        **/
/**  symbol-table whose nodes have a counter that is        **/
/**  incremented each time an item of an existing node is   **/
/**  searched.                                              **/
/**  The skip list data structure was developed by Pugh in  **/
/**  1990. It is an ordered linked list where each node     **/
/**  contains a variable number of links (set randomly).    **/
/**  There by, during a search it skips through large       **/
/**  portions of the list at a time due to the extra links  **/
/**  in the nodes. This characteristics provides a search,  **/
/**  insertion and removal of O(log N).                     **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Item.h"
#include "Utils.h"
#include "ST.h"

#define lgNmax 30 /* maximum number of levels */

/* Each node has an item, an array of links with length 'sz' */
/* and a counter 'cont' that is incremented each time the    */
/* item is searched.                                         */
typedef struct STnode *link;
struct STnode{ Item item; link *next; int sz; unsigned long cont; };

static link head;
static int N; /* number of items in the list */
static int lgN; /* actual number of levels */
static link maxItem; /* element with higher frequency */


/* ********************************************************* */
/* Creates a new 'STnode' with 'item' content and 'k' links, */
/* and returns its pointer.                                  */
static link NEW (Item item, int k)
{
   int i;
   link x;

   x = UTILmalloc (sizeof *x); /* allocates the 'STnode' */
   x->next = UTILmalloc (k * sizeof (link)); /* array of links */

   x->item = item; /* content */
   x->sz = k; /* number of links */
   for (i = 0; i < k; i++)
      x->next[i] = NULL; /* initializes the links with NULL pointer */
   x->cont = 1; /* initializes the counter */

   /* If 'maxItem' has not yet been initialized... */
   if (maxItem == NULLitem)
      maxItem = x; /* initializes the maxItem */

   return x;

} /* NEW */


/* ********************************************************* */
/* Initializes the skip list. Initializes the item's counter */
/* 'N' and the actual number of levels 'lgN' with 0, creates */
/* the head node with 'NULLitem' content and with 'lgNmax'   */
/* links and initializes the pointer to the element with     */
/* higher frequency with 'NULLitem'.                         */
void STinit ()
{
   N = 0; /* item's counter */
   lgN = 0; /* actual number of levels */
   head = NEW (NULLitem, lgNmax); /* creates the head node */
   maxItem = NULLitem; /* initializes the maxItem pointer */

} /* STinit */


/* ********************************************************* */
/* Generates a pseudo-random number 'i' (number of links of  */
/* a new node which will be inserted) with probability       */
/* '1/j^i' (where here 'j' is set equal 2).                  */
/* Obs.: for large skip lists, the memory space occupied by  */
/* the extra pointers becomes non-trivial and, in this case, */
/* the user can change the algorithm below by replacing      */
/* 'j = 2' by 'j = 3' and 'j = j * 2' by 'j = j * 3'.        */
static int randX ()
{
   int i, j, t;

   /* Generates a pseudo-random integer in [0,RAND_MAX]. */
   t = rand ();

   /* Probability: 1/2, 1/2^2, 1/2^3, ... */
   for (i = 1, j = 2; i < lgNmax; i++, j = j * 2)
      if (t > RAND_MAX / j)
	 break;

   if (i > lgN)
      lgN = i; /* updates the actual number of levels */

   return i;

} /* randX */


/* ********************************************************* */
/* Recursive function that inserts 'x' after 't' at level    */
/* 'k', preserving the lexicographic order.                  */
static void insertR (link t, link x, int k)
{
   Key v;

   /* Gets 'x' key. */
   v = key(x->item);
   if ((t->next[k] == NULL) || less (v, key (t->next[k]->item))) {
      if (k < x->sz) { /* dancing links... */
	 x->next[k] = t->next[k];
	 t->next[k] = x;
      }
      if (k == 0) /* the insertion ended */
	 return;
      insertR (t, x, k-1); /* tail recursion */
      return;
   }

   /* 'x' is after 't->next[k]' */
   insertR (t->next[k], x, k);

} /* insertR */


/* ********************************************************* */
/* Adds a new item (envelope function to be exported).       */
void STinsert (Item item)
{
   /* The insertion starts from the 'head'   */
   /* and at the actual highest level 'lgN'. */
   insertR (head, NEW (item, randX()), lgN);
   N++; /* one more item in the list */

} /* STinsert */


/* ********************************************************* */
/* Recursive function that searches an item with key 'v'     */
/* after 't' (taking into account that the list is in        */
/* lexicographic order) at level 'k'. The recursion proceeds */
/* as follows: it moves to the next node in the list on      */
/* level 'k' if its key is smaller than the search key 'v'   */
/* or down to level 'k-1' if its key is not smaller.         */
static Item searchR (link t, Key v, int k)
{
   if (t->next[k] == NULL) { /* end of the list of level 'k' */
      if (k == 0) /* not found */
	 return NULLitem;
      return searchR (t, v, k-1); /* down to level 'k-1' */
   }
   if (eq (v, key (t->next[k]->item))) { /* it was found */
      t->next[k]->cont++; /* increments the item counter */
      /* Checks if the item counter exceeded the actual maxItem. */
      if (maxItem->cont < t->next[k]->cont)
      	 maxItem = t->next[k];
      return t->next[k]->item;
   }
   if (less (v, key (t->next[k]->item))) { /* 'v' is smaller */
      if (k == 0) /* not found */
	 return NULLitem;
      return searchR (t, v, k-1); /* down to level 'k-1' */
   }
   return searchR (t->next[k], v, k); /* searches from the next item */

} /* searchR */


/* ********************************************************* */
/* Searches an item with a given key 'v' (envelope function  */
/* to be exported).                                          */
Item STsearch (Key v)
{
   /* The search begins from the 'head' and */
   /* at the actual highest level 'lgN'.    */
   return searchR (head, v, lgN);

} /* STsearch */


/* ********************************************************* */
/* Returns a pointer to the highest score element of the     */
/* list.                                                     */
/* Obs.: it is only to be used to correct the 'maxItem'      */
/* pointer (which may occur, for example, at 'deleteR'       */
/* function).                                                */
static link searchMaxItem ()
{
   link newMax, t;

   t = head;
   if (t->next[0] == NULL) /* empty list */
      return NULLitem;

   newMax = t->next[0];
   t = t->next[0];
   while (t->next[0] != NULL) { /* scans the list */
      t = t->next[0];
      if (t->cont > newMax->cont)
	 newMax = t;
   }

   return newMax;

} /* searchMaxItem */

/* ********************************************************* */
/* Recursive function that removes the item with key 'v'. It */
/* proceeds quite similar to the 'searchR' function. If the  */
/* item found has a counter greater than 1 this function     */
/* only decrements its counter, otherwise it is necessary to */
/* unlink it at each level that we find a link to it and     */
/* free the entire node when it reaches the bottom level. If */
/* the item is the highest score item, the 'maxItem' pointer */
/* have to be fixed by calling 'searchMaxItem' function.     */
static void deleteR (link t, Key v, int k)
{
   link x = t->next[k];
   if (t->next[k] == NULL) { /* end of the list of level 'k' */
      if (k > 0)
	 deleteR(t, v, k-1); /* down to level 'k-1' */
   }
   else if (eq (v, key (t->next[k]->item))) {
      if (t->next[k]->cont > 1) { /* only decrements the item's counter */
	 t->next[k]->cont--;
	 if (maxItem == t->next[k])
	    maxItem = searchMaxItem (); /* fixes maxItem pointer */
	 return;
      }
      t->next[k] = x->next[k]; /* unlink at level k */
      if (k == 0) { /* reached the bottom level */
	 if (maxItem == x) {
	    free (x->item); /* frees the node's item */
	    free (x->next); /* frees the node's vector of links */
	    free (x); /* frees the node */
	    maxItem = searchMaxItem (); /* fixes maxItem pointer */
	 }
	 else {
	    free (x->item); /* frees the node's item */
	    free (x->next); /* frees the node's vector of links */
	    free (x); /* frees the node */
	 }
	 return;
      }
      deleteR (t, v, k-1); /* down to level 'k-1' */
   }
   else if (less (v, key (t->next[k]->item))) { /* 'v' is smaller */
      if (k > 0)
	 deleteR (t, v, k-1); /* down to level 'k-1' */
   }
   else deleteR (t->next[k], v, k); /* 'v' is greater */

} /* deleteR */


/* ********************************************************* */
/* Removes an item with key 'v' (envelope function to be     */
/* exported).                                                */
void STdelete (Key v)
{
   /* It starts looking for the item to be deleted from */
   /* the 'head' and at the actual highest level 'lgN'. */
   deleteR (head, v, lgN);
   N--; /* one less item in the list */

} /* STdelete */


/* ********************************************************* */
/* Returns the 'k'-th smallest item or returns 'NULLitem' if */
/* the list is has less then 'k' items.                      */
Item STselect (int k)
{
   int i;
   link t = head;

   for (i = 0; i < k && t->next[0] != NULL; i++)
      t = t->next[0];

   /* If there is less than 'k' items in the list. */
   if (i != k)
      return NULLitem;

   return t->item;

} /* STselect */


/* ********************************************************* */
/* Visit the items in the order of their keys (calling a     */
/* procedure passed as an argument for each item).           */
void STsort (FILE *std, void (*visit)(FILE *std, Item))
{
   link t = head;

   while (t->next[0] != NULL) {
      t = t->next[0];
      visit(std, t->item); /* call the procedure */
   }

} /* STsort */


/* ********************************************************* */
/* Return the quantity of different items.                   */
int STcount ()
{
   return N;

} /* STcount */


/* ********************************************************* */
/* Returns the key of the highest score item.                */
Key STmaxItem ()
{
   return (key(maxItem->item));

} /* STmaxItem */


/* ********************************************************* */
/* Prints at 'std' the key of the highest score item.        */
void STshowMaxItem (FILE *std)
{
   ITEMshow (std, key(maxItem->item));

} /* STshowMaxItem */


/* ********************************************************* */
/* Returns the score of the highest score item.              */
unsigned long STmaxCont ()
{
   return maxItem->cont;

} /* STmaxCont */


/* ********************************************************* */
/* Returns the sum of the scores of all items.               */
unsigned long STtotalCount ()
{
   unsigned long count = 0;
   link t = head;

   while (t->next[0] != NULL) {
      t = t->next[0];
      count += t->cont;
   }

   return count;

} /* STtotalCount */


/* ********************************************************* */
/* Frees memory of all nodes (destroys the symbol-table).    */
void STfree ()
{
   link t;

   while (head->next[0] != NULL) {
      t = head->next[0];
      head->next[0] = t->next[0];
      free (t->item); /* frees the node's item */
      free (t->next); /* frees the node's vector of links */
      free (t); /* frees the node */
   }

   free (head); /* finally, frees the head */

} /* STfree */

