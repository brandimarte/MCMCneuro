/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface of Markov Chain Monte Carlo on graphs        **/
/**  representing neuronal connectivity.                    **/
/**  *****************************************************  **/

/* Matrix index in row-major order. */
#define idx2d(i, j, n) ((i) * (n) + j)

/* Sets the available memory for graphs storage. */
void NEUROsetMem (char *freeMem);

/* Sets the option for a fixed number of Monte Carlo steps. */
void NEUROsetMCsteps (char *steps);

/* Sets the brain region chosen for study. */
void NEUROsetRegion (char *rg);

/* Sets the mouse chosen for study. */
void NEUROsetMouse (char *mouse);

/* Sets the method for posterior probability computation. */
void NEUROsetMethod (char *method);

/* Sets the penalty constant chosen by the user. */
void NEUROsetPenal (char *penalty);

/* Estimates for each mouse the graph that best represents the  */
/* observed data in the first and third parts of the experiment */
/* for a fixed penalty value and method (1, 2 and 3) of         */
/* computing the posterior probability.                         */
void NEURObestGraph (char *dataPath, char *outPath);

/* Estimates for each mouse the graph that best represents the  */
/* observed data in the first and third parts of the experiment */
/* considering different penalty values and the different       */
/* methods (1, 2 and 3) of computing the posterior probability. */
void NEUROpenalAnalysis (char *dataPath, char *outPath);
