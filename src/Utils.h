/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface for alternative versions of well known       **/
/**  functions from 'stdio.h' and 'stdlib.h' to avoid code  **/
/**  repetition.                                            **/
/**  *****************************************************  **/

/* Allocates a block of bytes if there are     */
/* enough memory, otherwise exits the program. */
void *UTILmalloc (unsigned int nbytes);

/* Change the size of a block of bytes if there */
/* are enough memory or exits the program.      */
void *UTILrealloc (void *ptr1, unsigned int nbytes);

/* Opens the file named 'filename' in order to    */
/* execute a 'mode' operation and verifies error. */
FILE *UTILfopen (const char *filename, const char *mode);

/* Verifies the returned value of a call to the 'fscanf' */
/* function to read data from the file named 'filename'. */
void UTILcheckFscan (int info, const char *filename);

