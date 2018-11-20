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
/**  Abstract data type implementation for graph items,     **/
/**  which are represented by strings containing '0's and   **/
/**  '1's.                                                  **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Item.h"
#include "Utils.h"


/* ********************************************************* */
/* Reads a key from standard input and copies it into 'x'.   */
int ITEMscan (char **x)
{
   int t;
   char buf[1000];

   t = scanf ("%s", buf);
   *x = UTILmalloc (size (buf) * sizeof (char));
   (*x)[0] = '\0';
   copy (*x, buf);

   return t;

} /* ITEMscan */


/* ********************************************************* */
/* Prints at 'std' the key of an item 'x'.                   */
void ITEMshow (FILE *std, char *x)
{
   fprintf (std, "%s\n", x);

} /* ITEMshow */


/* ********************************************************* */
/* Performs the unary operator "not" (logical negation) on   */
/* an element index '|idx|-1' from an item 'new'.            */
/* If 'idx > 0' indicates that the element 'new[idx-1]' is   */
/* '0' and then changes it to '1'. If 'idx < 0' indicates    */
/* that the element 'new[-idx-1]' is '1' and then changes it */
/* to '0'.                                                   */
void ITEMgenerator (char *new, int idx)
{
   if (idx > 0)
      new[idx-1] = '1';
   else
      new[-idx-1] = '0';

} /* ITEMgenerator */


/* ********************************************************* */
/* Generates a pseudo-random index. Receives the length      */
/* 'delta' of a partition, generates a pseudo-random number  */
/* 'u' uniformly distributed in [0,1] and checks to which    */
/* partition 'part' this 'u' belongs (that's the index).     */
static int geraIdx (double delta)
{
   int part;
   double u;
   
   /* Generates u ~ Unif[0,1]. */
   u = 1.0 * rand() / RAND_MAX;

   /* Checks where is 'u'. */
   for (part = 1; delta * part < 1.0; part++)
      if (u < delta * part)
	 return (part - 1);
   return (part - 1);

} /* geraIdx */


/* ********************************************************* */
/* Picks a random element from an item and returns its index */
/* plus 1 ('idx + 1') if the element is '0' or the negative  */
/* value 'idx - 1'. The addition (or subtraction) of '1' is  */
/* necessary to avoid mistake when the index is zero.        */
int ITEMrandIdx (char *item)
{
   int idx;
   double delta;

   /* Computes the length of a partition of [0,1]. */
   delta = 1.0 / size(item);

   /* Generates a pseudo-random index. */
   idx = geraIdx (delta);

   if (item[idx] == '0')
      return (idx + 1);

   return (-idx - 1);

} /* ITEMrandIdx */

