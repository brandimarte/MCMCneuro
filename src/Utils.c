/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  Implementation of alternative versions of well known   **/
/**  functions from 'stdio.h' and 'stdlib.h' to avoid code  **/
/**  repetition.                                            **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include "Utils.h"


/* ********************************************************* */
/* Allocates a block of bytes if there are enough memory.    */
/* Otherwise returns an error message and exits the program. */
void *UTILmalloc (unsigned int nbytes)
{
   void *ptr;
   ptr = malloc (nbytes);

   if (ptr == NULL) {    
      fprintf (stderr, "\n Insufficient memory.\n");
      exit(EXIT_FAILURE);
   }

   return ptr;

} /* UTILmalloc */


/* ********************************************************* */
/* Reallocates a block of bytes pointed by 'ptr' and change  */
/* its size to 'nbytes' if there are enough memory.          */
/* Otherwise returns an error message and exits the program. */
void *UTILrealloc (void *ptr1, unsigned int nbytes)
{
   void *ptr2;
   ptr2 = realloc (ptr1, nbytes);

   if (ptr2 == NULL) {    
      fprintf (stderr, "\n Insufficient memory.\n");
      exit (EXIT_FAILURE);
   }

   return ptr2;

} /* UTILrealloc */


/* ********************************************************* */
/* Opens the file named 'filename' in order to execute an    */
/* operation specified by 'mode' (operations of the 'fopen'  */
/* function from 'stdio.h') and verifies error.              */
FILE *UTILfopen (const char *filename, const char *mode)
{
   FILE *pfile;

   pfile = fopen (filename, mode);
   if (pfile == NULL) {
      fprintf (stderr, "\n Error: Unable to open the file '%s'!\n\n", filename);
      exit (EXIT_FAILURE);
   }

   return pfile;

} /* UTILfopen */


/* ********************************************************* */
/* Receives the integer 'info' returned from a call to the   */
/* 'fscanf' function to read data from the file named        */
/* 'filename' and checks if 'info = EOF', which indicates    */
/* that an input failure happened before any data could be   */
/* successfully read.                                        */
void UTILcheckFscan (int info, const char *filename)
{
   if (info == EOF) {
      fprintf (stderr, "\n Error: Unable to read the file '%s'!\n\n", filename);
      exit (EXIT_FAILURE);
   }
   
} /* UTILcheckFscan */


