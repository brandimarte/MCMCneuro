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
/**  Abstract data type interface for graph items, which    **/
/**  are represented by strings containing '0's and '1's.   **/
/**  Therefore, the associated operations for an item are   **/
/**  those for strings and a few others.                    **/
/**  One of the most common representations for a graphs    **/
/**  is the adjacent matrix, i.e., a square matrix where    **/
/**  an element 'a_ij = 1' indicates that there is an edge  **/
/**  from vertex 'i' to vertex 'j' and 'a_ij = 0'           **/
/**  indicates that these vertices are not connected.       **/
/**  However, considering that the items are undirected     **/
/**  graphs, the adjacent matrix is symmetric and it is     **/
/**  necessary only to know the elements above (or below)   **/
/**  the main diagonal. Thus, these items (undirected       **/
/**  graphs) are strings that represents the elements       **/
/**  above the main diagonal. The way of indexing a string  **/
/**  V to match the adjacent matrix AM is given by:         **/
/**     for a_ij in AM:  if i < j,  V[i+1+(j+1)*(j-2)/2]    **/
/**                      otherwise, V[j+1+(i+1)*(i-2)/2]    **/
/**  *****************************************************  **/

typedef char *Item;
typedef char *Key;

#define key(A) (A)
#define eq(A, B) (strcmp(A,B)==0)
#define less(A, B) (strcmp(A,B)<0)
#define size(A) (strlen(A))
#define copy(A, B) (strcpy(A,B))
#define comp(A, B) (strcmp(A,B))
#define NULLitem (NULL)
#define idxAdj(i, j) (((i)<(j)) ? (i+1+(j+1)*(j-2)/2) : (j+1+(i+1)*(i-2)/2))

/* Reads a key from standard input. */
int ITEMscan (Item *);

/* Prints at 'std' the key of an item. */
void ITEMshow (FILE *std, Item);

/* Changes the 'idx' element from an item. */
void ITEMgenerator (char *new, int idx);

/* Picks a random element from an item and returns its index */
/* if the element is '0' or the negative value of the index. */
int ITEMrandIdx (char *item);
