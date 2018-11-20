/**  *****************************************************  **/
/**       ** Finding the Most Representative Graph **       **/
/**       **    Model for Neuronal Interactions    **       **/
/**       **     via Markov Chain Monte Carlo      **       **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**  *****************************************************  **/
/**  Implementation of Markov Chain Monte Carlo on graphs   **/
/**  representing neuronal connectivity. Given the          **/
/**  experimental data, it generates a Markov Chain, on an  **/
/**  undirected graph space, whose limit distribution is    **/
/**  given by the posterior probability 'P(g|X)'. This is   **/
/**  done via Monte Carlo method with Metropolis            **/
/**  algorithm.                                             **/
/**  The generated graphs are stored in a skip list         **/
/**  symbol-table whose nodes have a counter that is        **/
/**  incremented each time an item of an existing node is   **/
/**  searched (to avoid the need to store each generated    **/
/**  graph, i.e., only distinct graphs are stored). The     **/
/**  use of skip list structure has also the advantage      **/
/**  that search, insertion and removal are O(log N).       **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Utils.h"
#include "Item.h"
#include "ST.h"
#include "Neuro.h"

#define Trange 0.01 /* time range = 0.01 s = 10 ms */
#define Tstep 100 /* consider 1 Trange window after Tstep=100*Trange */
#define Jij 1.0 /* interaction "energy" */
#define Nsteps 100000 /* # Monte Carlo steps */

/* Structure to store the considered spikes read from the    */
/* data file, where 'label' is the neuron label and 'spikes' */
/* is a vector of spikes in the time interval considered.    */
typedef struct NEUROspk spkInfo;
struct  NEUROspk { char label[5]; int *spikes; };

static long MEM; /* available memory */
static int fixSteps; /* fixed number of MC steps option (0 or 1) */
static double penal; /* penalty constant */
static int met; /* label of the method for probability computation */
static int spkRange; /* time range of considered spikes */
static int rat; /* mouse id */
static int Nneuron; /* number of neurons */
static int Nedges; /* maximum number of edges */
static int part; /* part of the experiment (1 or 3) */
static char region[6]; /* brain region */


/* ********************************************************* */
/* Sets the available memory for graphs storage.             */
void NEUROsetMem (char *freeMem)
{
   MEM = atol (freeMem);

} /* NEUROsetMem */


/* ********************************************************* */
/* Sets the option for a fixed number of Monte Carlo steps   */
/* ('1' means that a fixed number was chosen by the user and */
/* '0' represents the opposite).                             */
void NEUROsetMCsteps (char *steps)
{
   fixSteps = atoi (steps);

} /* NEUROsetMCsteps */


/* ********************************************************* */
/* Sets the brain region chosen for study.                   */
void NEUROsetRegion (char *rg)
{
   region[0] = '\0';
   copy (region, rg);

} /* NEUROsetRegion */


/* ********************************************************* */
/* Sets the mouse chosen for study.                          */
void NEUROsetMouse (char *mouse)
{
   rat = atoi (mouse);

} /* NEUROsetMouse */


/* ********************************************************* */
/* Sets the method for posterior probability computation.    */
/* If first method was chosen (method = '1'), than it is     */
/* considered: <Xi|Xj>=1 if Xi=1 and Xj=1.                   */
/* If second method was chosen (method = '2'), than it is    */
/* considered: <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj).  */
/* If third method was chosen (method = '3'), than it is     */
/* considered: <Xi|Xj>=1 if Xi=Xj=1, <Xi|Xj>=-1 if Xi!=Xj    */
/* or <Xi|Xj>=0 if Xi=Xj=0.                                  */
void NEUROsetMethod (char *method)
{
   met = atoi (method);

} /* NEUROsetMethod */


/* ********************************************************* */
/* Sets the penalty constant chosen by the user.             */
void NEUROsetPenal (char *penalty)
{
   penal = atof (penalty);

} /* NEUROsetPenal */


/* ********************************************************* */
/* Prints the adjacency matrix representation of the key     */
/* 'max' at 'std' file. The key is a graph in vector         */
/* representation and, in order to maintain the coherence    */
/* with other parts of the program, the rule is fill the     */
/* elements below the main diagonal in row-major order.      */
static void printAdjMatrix (FILE *std, Key max, spkInfo *tkt)
{
   int i, j;
   Key adjMatrix;

   /* Creates the matrix. */
   adjMatrix = UTILmalloc (Nneuron * Nneuron * sizeof (adjMatrix));

   /* Initializes the matrix with "0s". */
   for (i = 0; i < Nneuron; i++)
      for (j = 0; j < Nneuron; j++)
	 adjMatrix[idx2d(i,j,Nneuron)] = '0';

   /*  *** This is an important "trick" of the algorithm!! ***  */
   /* Copies the values of vector 'max' to the symmetric matrix */
   /* following the rule of filling the elements below the      */
   /* main diagonal in row-major order. Note that, 'i > j' and  */
   /* 'idxAdj(i,j)' will give the correct index (see Item.h).   */
   for (i = 1; i < Nneuron; i++)
      for (j = 0; j < i; j++) {
	 adjMatrix[idx2d(i,j,Nneuron)] = max[idxAdj(i,j)];
	 adjMatrix[idx2d(j,i,Nneuron)] = max[idxAdj(i,j)];
      }

   /* Prints the matrix in 'std' file with the   */
   /* electrodes labels at each column and line. */
   fprintf (std, "    ");
   for (i = 0; i < Nneuron; i++)
      fprintf (std, " %s ", tkt[i].label);
   fprintf (std, "\n");
   for (i = 0; i < Nneuron; i++) {
      fprintf (std, "%s ", tkt[i].label);
      for (j = 0; j < Nneuron; j++)
	 fprintf (std, "  %c  ", adjMatrix[idx2d(i,j,Nneuron)]);
      fprintf (std, "\n");
   }

} /* printAdjMatrix */


/* ********************************************************* */
/* Receives an output file name 'outName' and prints a lot   */
/* of relevant information about the run (the code is quite  */
/* self explanatory).                                        */
static void output (char *outName, spkInfo *tkt, unsigned long accept,
		    unsigned long steps, unsigned long maxMCsteps)
{
   int i;
   FILE *out; /* file for general output */

   out = UTILfopen (outName, "a"); /* opens the general output file */

   fprintf (out, "** Mouse %d - %s - part %d **\n", rat, region, part);
   fprintf (out, "\nNumber of neurons: %d", Nneuron);
   fprintf (out, "\nNeurons labels: ");
   for (i = 0; i < Nneuron; i++)
      fprintf (out, " %s ", tkt[i].label);
   fprintf (out, "\nMC steps: %lu", steps);
   fprintf (out, "\nMaximum allowed MC steps: %lu", maxMCsteps);
   fprintf (out, "\nTotal graphs counted: %lu", STtotalCount());
   fprintf (out, "\nPenalty constant: %.5f", penal);
   fprintf (out, "\nDistinct graphs: %d", STcount());

   /* *** This might interest you!!! *** */
   /* For those who do not believe that this program really */
   /* works please uncomment the following line code. It    */
   /* will print all distinct graphs at the output file,    */
   /* but be aware that it can be a lot of graphs.          */
   /* STsort (stdout, ITEMshow); */

   fprintf (out, "\nAccepted graphs: %lu", accept);
   fprintf (out, "\nMost representative graph counter = %lu", STmaxCont());
   fprintf (out, "\nMost representative graph probability = %.5f",
	    1.0*STmaxCont()/(steps + 1));
   fprintf (out, "\nMost representative graph (vectorial form):\n");
   STshowMaxItem (out); /* show in vectorial form */
   fprintf (out, "\nMost representative graph (adjacency matrix):\n");
   printAdjMatrix (out, STmaxItem(), tkt); /* show adjacency matrix */
   fprintf (out, "\n\n");

   fclose (out); /* closes the general output file */

} /* output */


/* ********************************************************* */
/* Receives an output file name 'outName' and prints the     */
/* adjacency matrix of the highest score graph in this file. */
static void outputAdjM (char *outName, spkInfo *tkt)
{

   FILE *out; /* file for adjacency matrix output */

   out = UTILfopen (outName, "w"); /* opens the file */
   printAdjMatrix (out, STmaxItem(), tkt); /* prints the adjacency matrix */
   fclose (out); /* closes the file */

} /* outputAdj */


/* ********************************************************* */
/* Creates the Monte Carlo starting state randomly (the      */
/* probability of having each edge is 0.5).                  */
static Key MCinit ()
{
   int i;
   double u;
   Key ini;

   /* Allocates the graph. */
   ini = UTILmalloc ((Nedges + 1) * sizeof (ini));

   /* Decides (randomly) if each possible edge exists. */
   for (i = 0; i < Nedges; i++) {
      u =  1.0 * rand() / RAND_MAX; /* pseudo-random u ~ Unif[0,1] */
      if (u < 0.5) /* probability = 0.5 */
	 ini[i] = '0'; /* no edge */
      else
	 ini[i] = '1'; /* with edge */
   }
   ini[i] = '\0'; /* terminating null character */

   return ini;

} /* MCinit */


/* ********************************************************* */
/* Given two observed data 'obs1' and 'obs2', computes their */
/* their "interaction" according to the method chosen:       */
/*    1:  <Xi|Xj>=1 if Xi=1 and Xj=1.                        */
/*    2:  <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj).       */
/*    3:  <Xi|Xj>=1 if Xi=Xj=1, <Xi|Xj>=-1 if Xi!=Xj         */
/*        or <Xi|Xj>=0 if Xi=Xj=0.                           */
static int gibbsEnMethod (int obs1, int obs2)
{
   if (met == 1) /* <Xi|Xj>=1 if Xi=1 and Xj=1. */
      return (obs1 * obs2);
   if (met == 2) { /* <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj). */
      if (obs1 == obs2)
	 return 1;
      return 0;
   }
   /* <Xi|Xj>=1 if Xi=Xj=1, <Xi|Xj>=-1 if Xi!=Xj */
   /* or <Xi|Xj>=0 if Xi=Xj=0.                   */
   if (obs1 == 1 && obs2 == 1)
      return 1;
   else if (obs1 != obs2)
      return -1;
   return 0;

} /* gibbsEnMethod */


/* ********************************************************* */
/* Computes the "interaction energies" between all possible  */
/* neighbors.                                                */
static double *gibbsEn (spkInfo *tkt)
{
   int i, j, k, g;
   double *Vij;

   /* Allocates the "interaction energy" vector. */
   Vij = UTILmalloc (Nedges * sizeof (double));

   /*  *** This is an important "trick" of the algorithm!! ***  */
   /* The index of the "interaction energy" vector must have    */
   /* the same rule of the graph index when writing the         */
   /* adjacency matrix. Here this rule is fill the elements     */
   /* below the main diagonal in row-major order.               */
   for (g = 0, i = 1; i < Nneuron; i++)
      for (j = 0; j < i; j++) {
	  /* Scan the spikes. */
	 for (Vij[g] = 0.0, k = 0; k < spkRange; k++)
	    Vij[g] += gibbsEnMethod (tkt[i].spikes[k], tkt[j].spikes[k]);
	 Vij[++g] *= Jij; /* "interaction energy" between neighbors */
      }

   return Vij;

} /* *gibbsEn */


/* ********************************************************* */
/* Receives a candidate for a new state (the changed 'edge') */
/* and returns '1' if it is accepted or '0' if not. Uses the */
/* acceptance probability proposed originally by Metropolis  */
/* et al (1954) with a symmetric sampling kernel ('Qij=Qji') */
/*           'a(i,j) = min {1, P(g_j|X)/P(g_i|X)}'           */
/* where 'g_i' is the current state, 'g_j' is the candidate  */
/* and 'P(g_i|X)' is the posterior probability. It is easy   */
/* to see that the ratio 'P(g_j|X)/P(g_i|X)}' is given by    */
/* one of the following options:                             */
/*   1. 'exp(penal*spkRange-Jij*<Xi|Xj>)' if 'g_j' has one   */
/*      less edge than 'g_i' (an edge has been removed).     */
/*   2. 'exp(Jij*<Xi|Xj>-penal*spkRange)' if 'g_j' has one   */
/*      more edge than 'g_i' (an edge has been inserted).    */
/* The information about which edge has been changed comes   */
/* from the integer 'edge': if 'edge<0' it means that an     */
/* edge has been removed, otherwise an edge has been         */
/* inserted.                                                 */
/* Obs.: the addition of '1' in the 'edge' value was         */
/* necessary to avoid mistake when the index is zero (see    */
/* 'ITEMrandIdx' function at Item.c).                        */
static int metropolis (double *gibbsVij, int edge)
{
   double u;

   /* Generates u ~ Unif[0,1]. */
   u = 1.0 * rand() / RAND_MAX;

   /* An edge has been removed. */
   if (edge < 0) {
      if (exp (penal * spkRange - gibbsVij[-edge-1]) > u)
	 return 1; /* candidate accepted */
   }
   else { /* An edge has been inserted. */
      if (exp (gibbsVij[edge-1] - penal * spkRange) > u)
	 return 1; /* candidate accepted */
   }
   return 0; /* candidate rejected */

} /* metropolis */


/* ********************************************************* */
/* Given an initial state 'gr', it computes 'Nsteps' Monte   */
/* Carlo steps without including the generated graphs into   */
/* the accepted graphs list. This "thermalization steps" are */
/* used reduce the importance of choosing the initial state. */
/* After all the steps, the resulting graph 'gr' is taken as */
/* the initial state of the Monte Carlo.                     */
static void mcThermSteps (double *gibbsVij, Key gr)
{
   int i;
   int edge;

   /* 'Nsteps' Monte Carlo steps. */
   for (i = 0; i < Nsteps; i++) {

      /* Chooses a random edge to change. */
      edge = ITEMrandIdx (gr);

      /* metropolis = 1 if the candidate is accepted. */
      if (metropolis (gibbsVij, edge) == 1) {
	 
	 /* Changes the edge. */
	 if (edge >= 0)
	    gr[edge-1] = '1'; /* inserts the edge */
	 else
	    gr[-edge-1] = '0'; /* removes the edge */

      }

   } /* for (i = 0; i ... */

} /* mcThermSteps */


/* ********************************************************* */
/* Given an initial state 'gr', it computes 'Nsteps' Monte   */
/* Carlo steps. Each accepted state is included into the     */
/* accepted graphs list (skip list symbol-table - see ST.c). */
/* Returns the number of accepted graphs.                    */
static int mcSteps (double *gibbsVij, Key *gr)
{
   int i;
   int accept; /* # of accepted graphs */
   int edge; /* index of the changed edge */
   Item item; /* skip list object */
   Key grProx; /* candidate graph */

   /* 'Nsteps' Monte Carlo steps. */
   for (accept = 0, i = 0; i < Nsteps; i++) {

      /* Chooses a random edge to change. */
      edge = ITEMrandIdx (*gr);

      /* metropolis = 1 if the candidate is accepted. */
      if (metropolis (gibbsVij, edge) == 1) {

	 accept++; /* one more graph */

	 /* Allocates memory for the new graph. */
	 grProx = UTILmalloc ((Nedges + 1 ) * sizeof (grProx));
	 copy (grProx, *gr); /* grProx = *gr*/
	 ITEMgenerator (grProx, edge); /* changes the edge */

	 /* Searches for 'grProx' at graphs list. Returns its   */
	 /* pointer and increments its counter if it was found, */
	 /* otherwise returns 'NULLitem' (see ST.c).            */
	 *gr = STsearch(grProx);

	 if (*gr == NULLitem) { /* there is no grProx at graphs list */
	    key(item) = grProx;
	    STinsert (item); /* adds to the list */
	    *gr = grProx; /* grProx is the current state */
	 }
	 else /* grProx exists already at graphs list */
	    free (grProx);

      }
      else /* candidate rejected */
	 *gr = STsearch(*gr); /* increments 'gr' counter */

   } /* for (i = 0; i ... */

   return accept;

} /* mcSteps */


/* ********************************************************* */
/* Generates a Markov Chain, on an undirected graph space,   */
/* whose limit distribution is given by the posterior        */
/* probability 'P(g|X)'. This is done via Monte Carlo method */
/* with Metropolis algorithm.                                */
static void mcmc (int type, char *outPath, spkInfo *tkt, double *logPP)
{
   int i, length;
   unsigned long steps, accept; /* MC steps counter, # accepted graphs */
   Item item; /* skip list object */
   Key gr; /* current graph */
   double *Vij; /* "interaction energy": Vij = Jij * <Xi|Xj> */
   unsigned long maxMCsteps; /* maximum MC steps */
   char *outName; /* file name for general output */

   /* Initializes variables. */
   Nedges = Nneuron * (Nneuron - 1) / 2; /* # of edges */
   gr = MCinit (); /* first graph generated randomly */

   /* File name for general output. */
   length = size (outPath);
   outName = UTILmalloc ((length + 35) * sizeof (char));
   outName[0] = '\0';

   /* Maximum Monte Carlo steps. */
   if (fixSteps) /* chosen by the user */
      maxMCsteps = MEM;
   else /* depends on the memory available */
      maxMCsteps = 3 * (MEM / (2 * (Nedges + 35)));

   /* Computes all possible "interaction energies". */
   Vij = gibbsEn (tkt);

   /* "Thermalization" steps. */
   mcThermSteps (Vij, gr);

   /* Creates and initializes the skip list. */
   key(item) = gr;
   STinit (); /* creates the skip list */
   STinsert (item); /* insert 'gr' in the skip list */
   accept = 1; /* # of accepted graphs */

   /* Monte Carlo steps. */
   for (steps = 0; steps < maxMCsteps; steps += Nsteps)
      accept += mcSteps (Vij, &gr);

   /* Computes some results. */
   if (type == 0) { /* 'penalty analysis' run */
      logPP[0] = logPP[1] = logPP[2] = 0.0;
      gr = STmaxItem();
      for (i = 0; i < Nedges; i++) {
	 /* Non-normalized log-posterior probability (with penalty). */
	 logPP[0] += (gr[i] - '0') * (Vij[i] - penal * spkRange);
	 /* Non-normalized log-posterior probability (without penalty). */
	 logPP[1] += (gr[i] - '0') * Vij[i];
      }
      /* Empirical probability (obtained from the Monte Carlo). */
      logPP[2] = 1.0*STmaxCont()/(steps + 1);

      /* Writes output data. */
      sprintf (outName, "%soutputM%d%sp%dMet%d.dat",
	       outPath, rat, region, part, met);
      output (outName, tkt, accept, steps, maxMCsteps);

      /* /\* Output of the adjacency matrix. *\/ */
      /* outName[0] = '\0'; */
      /* sprintf (outName, "%sadjM%d%sp%dMet%dPen%.5f.dat", */
      /* 	       outPath, rat, region, part, met, penal); */
      /* outputAdjM (outName, tkt); */

      /* Frees memory. */
      free (Vij);
      STfree ();
      free (outName);

   }
   else { /* 'best graph' run */
      /* Writes output data. */
      sprintf (outName, "%soutputM%d%sp%dMet%d.dat",
	       outPath, rat, region, part, met);
      output (outName, tkt, accept, steps, maxMCsteps);

      /* Output of the adjacency matrix. */
      outName[0] = '\0';
      sprintf (outName, "%sadjM%d%sp%dMet%d.dat",
	       outPath, rat, region, part, met);
      outputAdjM (outName, tkt);

      /* Frees memory. */
      free (Vij);
      STfree ();
      free (outName);
   }

} /* mcmc */


/* ********************************************************* */
/* Receives the name of a file containing the observed       */
/* spikes of a given neuron and the minimum time to be       */
/* considered. Reads the file and generates a vector of      */
/* spikes where an element is equal to '1' if there was a    */
/* spike at a 'Trange' time window after each 'Tstep' time   */
/* interval or '0' if not. This vector is then returned.     */
static int *spkRead (char *spkPath, double min)
{
   int j;
   int *tktSPK; /* vector of spikes */
   double time, aux;
   FILE *spk; /* file with spikes */

   /* Opens the file with spikes. */
   spk = UTILfopen (spkPath, "r");

   /* Allocates the vector of spikes. */
   tktSPK = UTILmalloc (spkRange * sizeof (int));

   /* Scans the file considering 1 'Trange' window after each 'Tstep'. */
   aux = 0.0;
   time = min;
   for (j = 0; j < spkRange; j++) {

      /* Reads the file until it reaches the next starting point. */
      while (aux < time)
	 UTILcheckFscan (fscanf (spk, "%lf", &aux), spkPath);

      tktSPK[j] = 0;
      while (aux < time + Trange) { /* checks if there is a spike */
	 UTILcheckFscan (fscanf (spk, "%lf", &aux), spkPath);
	 tktSPK[j] = 1;
      }

      /* Jumps the 'Tstep' time interval. */
      time += Tstep * Trange;

   } /* for (j = 0; j < ... */
   
   fclose (spk); /* closes the file */

   return tktSPK;
   
} /* *spkRead */


/* ********************************************************* */
/* Frees memory of the stored spikes data (i.e. neuron label */
/* and its spikes in the time interval considered).          */
static void spkFree (spkInfo *tkt)
{
   int i;

   for (i = 0; i < Nneuron; i++)
      free (tkt[i].spikes);
   free (tkt);
   
} /* spkFree */


/* ********************************************************* */
/* Function for penalty analysis. It calls the 'mcmc'        */
/* function (Markov Chain Monte Carlo) with different        */
/* penalty values and with the different methods (1, 2 and   */
/* 3) for computing the posterior probability. It receives a */
/* path 'outPath' for writing the output files, the penalty  */
/* interval ('ini' to 'end') to be considered and the rate   */
/* of change of the penalty. Creates 3 files containing the  */
/* penalty value versus:                                     */
/*    penal1: Non-normalized log-posterior probability       */
/*            with penalty.                                  */
/*    penal2: Non-normalized log-posterior probability       */
/*            without penalty.                               */
/*    penal3: Empirical probability obtained from the        */
/*            Monte Carlo.                                   */
static void penalMetMCMC (char *outPath, spkInfo *tkt,
			  double ini, double end, double delta)
{
   int length, cont;
   double *logPP; /* log-posterior probability */
   char *outName1, *outName2, *outName3; /* names of the output files */
   FILE *out1, *out2, *out3; /* files for output */

   /* Initializes variables. */
   logPP = UTILmalloc (3 * sizeof (double)); /* log-posterior probability */
   length = size (outPath);
   outName1 = UTILmalloc ((length + 26) * sizeof (char));
   outName2 = UTILmalloc ((length + 26) * sizeof (char));
   outName3 = UTILmalloc ((length + 26) * sizeof (char));
   outName1[0] = '\0';
   outName2[0] = '\0';
   outName3[0] = '\0';
   sprintf (outName1, "%spenal1M%d%sp%dMet%d.dat",
	    outPath, rat, region, part, met);
   sprintf (outName2, "%spenal2M%d%sp%dMet%d.dat",
	    outPath, rat, region, part, met);
   sprintf (outName3, "%spenal3M%d%sp%dMet%d.dat",
	    outPath, rat, region, part, met);

   /* Some user interaction. */
   printf ("\n  method %d", met);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Opens the output files. */
   out1 = UTILfopen (outName1, "w");
   out2 = UTILfopen (outName2, "w");
   out3 = UTILfopen (outName3, "w");

   /* Looping over penalty values. */
   for (cont = 0, penal = ini; penal < end; penal += delta) {
      /* Markov Chain Monte Carlo. */
      mcmc (0, outPath, tkt, logPP);

      /* Writes results. */
      fprintf (out1, "%.7f  %.10f\n", penal, logPP[0]);
      fflush (out1); /* print now! */
      fprintf (out2, "%.7f  %.10f\n", penal, logPP[1]);
      fflush (out2); /* print now! */
      fprintf (out3, "%.7f  %.10f\n", penal, logPP[2]);
      fflush (out3); /* print now! */

      /* Just to avoid unnecessary computation. */
      if (logPP[0] < 0.0000001) {
	 cont++;
	 if (cont > 10)
	    break;
      }

   } /* for (cont = 0, penal = ... */

   /* Closes the output files. */
   fclose (out1);
   fclose (out2);
   fclose (out3);

   /* Frees memory. */
   free (outName1);
   free (outName2);
   free (outName3);
   free (logPP);

} /* penalMetMCMC */


/* ********************************************************* */
/* Function for penalty analysis. It receives the name of a  */
/* file 'dataFile' containing summarized information about   */
/* the experimental data from the chosen brain region (like  */
/* (the mouse ID, number of neurons, start and end times of  */
/* observation, label and path to the observed spikes of     */
/* each neuron. Receives also a path 'outPath' for writing   */
/* the output files. Then for each mouse in 'dataFile' it    */
/* gets its spikes information and computes the Markov Chain */
/* Monte Carlo with different methods for computing the      */
/* posterior probability and with different penalty values   */
/* by calling the 'penalMetMCMC' function.                   */
static void neuroPenal (char *dataFile, char *outPath)
{
   int i, mouse;
   double min, max; /* 'min' and 'max' spike times in a set */
   spkInfo *tkt; /* vector with electrode label and spikes */
   char spkPath[150]; /* path to the spikes data */
   char aux1[5], aux2[150];
   FILE *summary; /* file containing a summary of data */

   /* Opens the input file with data paths. */
   summary = UTILfopen (dataFile, "r");

   /* Looping over the mouses. */
   while (fscanf (summary, "%d%d", &mouse, &Nneuron) == 2) {

      /* Reads the 'min' and 'max' values of spike time in the set. */
      UTILcheckFscan (fscanf (summary, "%lf%lf", &min, &max), dataFile);

      if (mouse == rat) { /* mouse chosen for study */

	 /* Time range considered (discounting the first and the last ~5min). */
	 min += 300.0;
	 max -= 300.0;
	 spkRange = (int) (max - min);

	 /* Computes the Monte Carlo only if time is greater than ~30 min. */
	 if (spkRange > 1000) {

	    /* Allocates a vector for neuron's label and spikes. */
	    tkt = UTILmalloc (Nneuron * sizeof (spkInfo));

	    /* Looping over the neurons. */
	    for (i = 0; i < Nneuron; i++) {
	       
	       /* Reads the neuron's label and the path to the spikes data. */
	       UTILcheckFscan (fscanf (summary, "%s%s",
				       tkt[i].label, spkPath), dataFile);
	 
	       /* Gets the spikes. */
	       tkt[i].spikes = spkRead (spkPath, min);
	    }

	    /* Computes the MCMC with different methods for computing the */
	    /* posterior probability and with different penalty values.   */

	    /* First method: <Xi|Xj>=1 if Xi=1 and Xj=1. */
	    met = 1;
	    penalMetMCMC (outPath, tkt, 0.00001, 0.1, 0.0001);

	    /* Second method: <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj). */
	    met = 2;
	    penalMetMCMC (outPath, tkt, 0.6, 1.0, 0.001);

	    /* Third method: <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=-1 (if Xi!=Xj). */
	    met = 3;
	    penalMetMCMC (outPath, tkt, 0.2, 1.0, 0.001);	 

	    /* Frees memory. */
	    spkFree (tkt);
	 }
	 else {
	    /* Looping over the neurons. */
	    for (i = 0; i < Nneuron; i++)
	       /* Reads the neuron's label and the path to the spikes data. */
	       UTILcheckFscan (fscanf (summary, "%s%s", aux1, aux2), dataFile);

	    /* Some user interaction. */
	    printf (" (insufficient data)");
	    setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
	 }
      }
      else { /* this is not the mouse chosen for study */
	 /* Looping over the neurons. */
	 for (i = 0; i < Nneuron; i++)
	    /* Reads the neuron's label and the path to the spikes data. */
	    UTILcheckFscan (fscanf (summary, "%s%s", aux1, aux2), dataFile);
      }

   } /* while */

   /* Closes the input file. */
   fclose (summary);

} /* neuroPenal */


/* ********************************************************* */
/* Function for best graph computation. It receives the name */
/* of a file 'dataFile' containing summarized information    */
/* about the experimental data from the chosen brain region  */
/* (like the mouse ID, number of neurons, start and end      */
/* times of observation, label and path to the observed      */
/* spikes of each neuron. Receives also a path 'outPath' for */
/* writing the output files. Then for each mouse in          */
/* 'dataFile' it gets its spikes information and computes    */
/* the Markov Chain Monte Carlo.                             */
static void neuro (char *dataFile, char *outPath)
{
   int i, mouse;
   double min, max; /* 'min' and 'max' spike times in a set */
   double *aux = &min;
   spkInfo *tkt; /* vector with electrode label and spikes */
   char spkPath[150]; /* path to the spikes data */
   char aux1[5], aux2[150];
   FILE *summary; /* file containing a summary of data */

   /* Opens the input file with data paths. */
   summary = UTILfopen (dataFile, "r");

   /* Looping over the mouses. */
   while (fscanf (summary, "%d%d", &mouse, &Nneuron) == 2) {

      /* Reads the 'min' and 'max' values of spike time in the set. */
      UTILcheckFscan (fscanf (summary, "%lf%lf", &min, &max), dataFile);

      if (mouse == rat) { /* mouse chosen for study */

	 /* Time range considered (discounting the first and the last ~5min). */
	 min += 300.0;
	 max -= 300.0;
	 spkRange = (int) (max - min); /* time range of considered spikes */

	 /* Computes the Monte Carlo only if time is greater than ~30 min. */
	 if (spkRange > 1000) {

	    /* Allocates a vector for neuron's label and spikes. */
	    tkt = UTILmalloc (Nneuron * sizeof (spkInfo));

	    /* Looping over the neurons. */
	    for (i = 0; i < Nneuron; i++) {

	       /* Reads the neuron's label and the path to the spikes data. */
	       UTILcheckFscan (fscanf (summary, "%s%s",
				       tkt[i].label, spkPath), dataFile);
	 
	       /* Gets the spikes. */
	       tkt[i].spikes = spkRead (spkPath, min);
	    }

	    /* Markov Chain Monte Carlo. */
	    mcmc (1, outPath, tkt, aux);

	    /* Frees memory. */
	    spkFree (tkt);
	 }
	 else {
	    /* Loop over the neurons. */
	    for (i = 0; i < Nneuron; i++)
	       /* Reads the neuron's label and the path to the spikes data. */
	       UTILcheckFscan (fscanf (summary, "%s%s", aux1, aux2), dataFile);

	    /* Some user interaction. */
	    printf (" (insufficient data)");
	    setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
	 }
      }
      else { /* this is not the mouse chosen for study */
	 /* Looping over the neurons. */
	 for (i = 0; i < Nneuron; i++)
	    /* Reads the neuron's label and the path to the spikes data. */
	    UTILcheckFscan (fscanf (summary, "%s%s", aux1, aux2), dataFile);
      }

   } /* while */

   /* Closes the input file. */
   fclose (summary);

} /* neuro */


/* ********************************************************* */
/* Given the path to the summarized data 'dataPath' and the  */
/* path 'outPath' for writing the output files, it           */
/* estimates, for the chosen mouse at the chosen region, the */
/* graph that best represents the observed data in the first */
/* and third parts of the experiment (before and after the   */
/* mouse got in touch with geometric objects) for a fixed    */
/* penalty value and method (1, 2 and 3) of computing the    */
/* posterior probability.                                    */
void NEURObestGraph (char *dataPath, char *outPath)
{
   int length;
   char *file;

   length = size (dataPath);
   file = UTILmalloc ((length + 20) * sizeof (char));

   /* First part of the experiment (before contacts). */
   printf ("\n Mouse %d - %s region - part 1", rat, region);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   file[0] = '\0';
   sprintf (file, "%s%smousesP1.dat", dataPath, region);
   part = 1; /* first part of the experiment */
   neuro (file, outPath);

   /* Third part of the experiment (after contacts). */
   printf ("\n\n Mouse %d - %s region - part 3", rat, region);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   file[0] = '\0';
   sprintf (file, "%s%smousesP3.dat", dataPath, region);
   part = 3; /* third part of the experiment */
   neuro (file, outPath);

   /* Frees memory. */
   free (file);

} /* NEURObestGraph */


/* ********************************************************* */
/* Given the path to the summarized data 'dataPath' and the  */
/* path 'outPath' for writing the output files, it           */
/* estimates, for the chosen mouse at the chosen region, the */
/* graph that best represents the observed data in the first */
/* and third parts of the experiment (before and after the   */
/* mouse got in touch with geometric objects) considering    */
/* different penalty values and the different methods (1, 2  */
/* and 3) of computing the posterior probability.            */
void NEUROpenalAnalysis (char *dataPath, char *outPath)
{
   int length;
   char *file;

   length = size (dataPath);
   file = UTILmalloc ((length + 20) * sizeof (char));

   /* First part of the experiment (before contacts). */
   printf ("\n Mouse %d - %s region - part 1", rat, region);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   file[0] = '\0';
   sprintf (file, "%s%smousesP1.dat", dataPath, region);
   part = 1; /* first part of the experiment */
   neuroPenal (file, outPath);

   /* Third part of the experiment (after contacts). */
   printf ("\n\n Mouse %d - %s region - part 3", rat, region);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   file[0] = '\0';
   sprintf (file, "%s%smousesP3.dat", dataPath, region);
   part = 3; /* third part of the experiment */
   neuroPenal (file, outPath);

   /* Frees memory. */
   free (file);

} /* NEUROpenalAnalysis */

