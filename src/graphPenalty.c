/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  This is a (client) program that estimates for each     **/
/**  mouse the graph that best represents the observed      **/
/**  data in the first and third parts of the experiment    **/
/**  (before and after the mouse got in touch with          **/
/**  geometric objects) considering different values for    **/
/**  the penalty constant and considering 3 different       **/
/**  methods for computing the posterior probability:       **/
/**    1. <Xi|Xj>=1 if Xi=1 and Xj=1                        **/
/**    2. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj)       **/
/**    3. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=-1 (if Xi!=Xj)      **/
/**  This is done via  Markov Chain Monte Carlo method.     **/
/**  The output files are sent to the chosen output folder  **/
/**  (default is 'out/outPenalty'). For example:            **/
/**    - the file 'outputM6HPp1Met1.dat' contains general   **/
/**    information about the run for mouse 6 (M6) in the    **/
/**    hippocampus region (HP) in the first part of the     **/
/**    experiment (p1) and where the posterior probability  **/
/**    were computed with method 1 (Met1).                  **/
/**    - the file 'penal1M6HPp1Met1.dat' contains the       **/
/**    penalty constant versus the logarithm of the         **/
/**    non-normalized posterior probability (with penalty   **/
/**    included).                                           **/
/**    - the file 'penal2M6HPp1Met1.dat' contains the       **/
/**    penalty constant versus the logarithm of the         **/
/**    non-normalized posterior probability (without        **/
/**    penalty).                                            **/
/**    - the file 'penal3M6HPp1Met1.dat' contains the       **/
/**    penalty constant versus the empirical probability    **/
/**    (obtained from the Monte Carlo).                     **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include "Neuro.h"

int main (int nargs, char *arg[])
{

   /* Checks if the input were typed correctly. */
   if (nargs != 7) {
      fprintf (stderr, "\n\n Wrong number of arguments!\n");
      fprintf (stderr, "\n Use: ./graphPenalty"); /* arg[0] */
      fprintf (stderr, " [directory with data paths]"); /* arg[1] */
      fprintf (stderr, " [directory for output files]"); /* arg[2] */
      fprintf (stderr, " [available memory]"); /* arg[3] */
      fprintf (stderr, " [fixed # of MC steps option: 0 or 1]"); /* arg[4] */
      fprintf (stderr, " [brain region option]"); /* arg[5] */
      fprintf (stderr, " [chosen mouse]\n\n"); /* arg[6] */
      exit (EXIT_FAILURE);
   }

   /* Sets available memory for graphs storage. */
   NEUROsetMem (arg[3]);

   /* Sets the option for a fixed number of MC steps. */
   NEUROsetMCsteps (arg[4]);

   /* Sets the brain region chosen for study. */
   NEUROsetRegion (arg[5]);

   /* Sets the mouse chosen for study. */
   NEUROsetMouse (arg[6]);

   /* Brain region specified by the user. */
   NEUROpenalAnalysis (arg[1], arg[2]);

   return 0;
   
} /* main */

