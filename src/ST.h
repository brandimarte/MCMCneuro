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
/**  Abstract data type interface for symbol-table whose    **/
/**  items have a counter that is incremented each time     **/
/**  the item is searched.                                  **/
/**  *****************************************************  **/

/* Initializes. */
void STinit ();

/* Adds a new item. */
void STinsert (Item);

/* Searches an item with a given key. */
Item STsearch (Key);

/* Removes an item. */
void STdelete (Item);

/* Returns the "int"-th smallest item. */
Item STselect (int);

/* Visit the items in the order of their keys (calling */
/* a procedure passed as an argument for each item).   */
void STsort (FILE *std, void (*visit)(FILE *std, Item));

/* Return the quantity of different items. */
int STcount ();

/* Returns the key of the highest score item. */
Key STmaxItem ();

/* Prints at 'std' the key of the highest score item. */
void STshowMaxItem (FILE *std);

/* Returns the score of the highest score item. */
unsigned long STmaxCont ();

/* Returns the sum of the scores of all items. */
unsigned long STtotalCount ();

/* Frees memory (destroys the symbol-table). */
void STfree ();
