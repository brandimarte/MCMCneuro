/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  This is a (client) program that computes a Markov      **/
/**  Chain Monte Carlo on graphs representing neuronal      **/
/**  connectivity to estimate the most representative       **/
/**  graph for each set of observed neurons. Given a        **/
/**  penalty constant (chosen form a previous analysis      **/
/**  with 'graphPenalty' program) and a specific method     **/
/**  for computing the posterior probability                **/
/**    1. <Xi|Xj>=1 if Xi=1 and Xj=1                        **/
/**    2. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj)       **/
/**    3. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=-1 (if Xi!=Xj)      **/
/**  it estimates the graph that best represents each       **/
/**  observed brain region of each observed mouse.          **/
/**  The output files are sent to the chosen output folder  **/
/**  (default is 'out/outBestGraph'). For example:          **/
/**    - the file 'outputHPMet1.dat' contains general       **/
/**    information about the run for all mouses where the   **/
/**    hippocampus region (HP) where observed and where     **/
/**    the posterior probability were computed with method  **/
/**    1 (Met1).                                            **/
/**    - the file 'adjM6HPp1Met1.dat' is the adjacent       **/
/**    matrix of the most probable graph of mouse 6 (M6),   **/
/**    in the hippocampus region (HP), in the first part    **/
/**    of the experiment (p1) and where the posterior       **/
/**    probability was computed with method 1 (Met1).       **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include "Neuro.h"

int main (int nargs, char *arg[])
{

   /* Checks if the input were typed correctly. */
   if (nargs != 9) {
      fprintf (stderr, "\n\n Wrong number of arguments!\n");
      fprintf (stderr, "\n Use: ./bestGraph"); /* arg[0] */
      fprintf (stderr, " [directory with data paths]"); /* arg[1] */
      fprintf (stderr, " [directory for output files]"); /* arg[2] */
      fprintf (stderr, " [available memory]"); /* arg[3] */
      fprintf (stderr, " [fixed # of MC steps option: 0 or 1]"); /* arg[4] */
      fprintf (stderr, " [brain region option]"); /* arg[5] */
      fprintf (stderr, " [chosen mouse]"); /* arg[6] */
      fprintf (stderr,
	       " [method for probability computation (0,1,2)]"); /* arg[7] */
      fprintf (stderr, " [penalty value]\n\n"); /* arg[8] */
      exit (EXIT_FAILURE);
   }

   /* Sets the available memory for graphs storage. */
   NEUROsetMem (arg[3]);

   /* Sets the option for a fixed number of MC steps. */
   NEUROsetMCsteps (arg[4]);

   /* Sets the brain region chosen for study. */
   NEUROsetRegion (arg[5]);

   /* Sets the mouse chosen for study. */
   NEUROsetMouse (arg[6]);

   /* Sets the method for posterior probability computation. */
   NEUROsetMethod (arg[7]);

   /* Sets the penalty constant. */
   NEUROsetPenal (arg[8]);

   /* Brain region specified by the user. */
   NEURObestGraph (arg[1], arg[2]);

   return 0;
   
} /* main */

