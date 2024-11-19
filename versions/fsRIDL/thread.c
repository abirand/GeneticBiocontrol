#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sched.h>
#include <sys/ipc.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include "parameters.h"
#include "random.h"
#include "niche.h"

#if _HABITAT
extern u64b     Habitat_logs;
extern FILE    *Habitatf;
#endif

#if _DISTANCE
extern u64b     Distance_logs;
extern FILE    *Distancef;
int             distancedata[12][_DRANGE+1];
#endif

#if _DISTANCE
extern u64b     Scentscape_logs;
extern FILE    *Scentscapef;
double          Scentdata[2][_WIDTH*_HEIGHT];
#endif

int             InocGen[]                       = _InocGen;
//int             InocCoor[_InocPatchNo][2]       = _InocCoor;  // converts the data to coordinate pairs
//int             InocIndividuals[]               = _InocIndividuals;

//#if _MONITOR
// monitoring locations post culling
//int             MonitorGen[]                    = _MonitorGen;
//int             MonitorCoor[_MonitorPatchNo][2] = _MonitorCoor;  // converts the data to coordinate pairs
//#endif

int             CullGen[]                       = _CullGen;
RndState        R;

// Mutex and condition variables for thread syncronization
extern volatile pthread_cond_t  Master,Slaves;
extern volatile pthread_mutex_t Mutex;
extern volatile int             Go[_CPUS],Done[_CPUS];
extern volatile pThreadData     Threads[_CPUS];

// Random Seed, and starting space state.  Copied by the threads
extern          u64b            Seed;
extern          Space           Start;
extern          u64b            UniqueCounter;

// Array based parameters (to be coppied to threads) 
extern          double          Xi[];

// The current and total generation numbers
extern          u64b            Gen;
extern          u64b            NGen;
extern          int             OffspringPopSize;
extern          int             Mortality;
extern          int             MonogamousMating;
extern          int             SinglePaternity;
extern          int             MultiplePaternity;
extern          int             Surviving;
extern          int             Inoculants;
extern          int             Total_Females;
extern          int             Total_Population;


// Code shared with niche.c
#include "shared.h"


/***************************************************************************************
 * Functions dealing with the creation and syncronization of threads
 ***************************************************************************************/

// Signals that the thread td has completed and then waits for a go (or exit) signal
static int SignalDone(pThreadData td) 
{
  // Signal done and wait for Go (or exit) singal
  pthread_mutex_lock((pthread_mutex_t*)&Mutex);
  Done[td->id] = 1;
  pthread_cond_signal((pthread_cond_t*)&Master);
  while( !Go[td->id] )
    pthread_cond_wait((pthread_cond_t*)&Slaves, (pthread_mutex_t*)&Mutex);
  Go[td->id] = 0;
  pthread_mutex_unlock((pthread_mutex_t*)&Mutex);

  // For now there is no exit condition
  return 1;
}


/***************************************************************************************
 * Functions dealing directly with maintenance of Individuals
 ***************************************************************************************/

// I'm commenting out these functions as they are not used:
// this will keep the compiler warnings to a minimun.
#if 0

// Returns the number of individuals used in the "next generation"
// Note: This count _will_ include wanderers if there are any.
static u64b nNextIndividuals()
{
  u64b i,m;

  for(m=i=0; i<_CPUS; i++)
    m += Threads[i]->nINext;
  return m;
}


// Returns the number of individuals used in the entire simulation.
// Note: This count _will_ include wanderers if there are any.
static u64b nIndividuals()
{
  u64b i,m;

  for(m=i=0; i<_CPUS; i++)
    m += Threads[i]->nICurrent + Threads[i]->nINext;
  return m;
}
#endif


// Fills in traits from characteristic array
static void SyncIndividual(pThreadData td, pIndividual i)
{

  u64b j,k;

  // Sync patch related traits from characteristic array
  //for(j=0;j<_k;j++) {
  //  for(i->x[j]=0.,k=_Le*j; k<_Le*(j+1); k++)
  //    i->x[j] += ((i->x0[k/LPW_E]>>(_Ae*(k%LPW_E)))&MASK_E) + ((i->x1[k/LPW_E]>>(_Ae*(k%LPW_E)))&MASK_E);
    //for(i->y[j]=0.,k=_Lp*j; k<_Lp*(j+1); k++)
      //i->y[j] += ((i->y0[k/LPW_P]>>(_Ap*(k%LPW_P)))&MASK_P) + ((i->y1[k/LPW_P]>>(_Ap*(k%LPW_P)))&MASK_P);
    // Scale traits
  //  i->x[j] /= MASK_E*(_Le<<1);
    //i->y[j] /= MASK_P*(_Lp<<1);
  //}

  // Sync and scale maker
  //for(i->k=0.,k=0; k<_Lk; k++)
  //  i->k += ((i->k0[k/LPW_K]>>(_Ak*(k%LPW_K)))&MASK_K) + ((i->k1[k/LPW_K]>>(_Ak*(k%LPW_K)))&MASK_K);
  //i->k /= MASK_K*(_Lk<<1);
  
  // Sync and scale mating like in others [0,1]
  for(i->m=0.,k=0; k<_Lm; k++)
    i->m += ((i->m0[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M) + ((i->m1[k/LPW_M]>>(_Am*(k%LPW_M)))&MASK_M);
  i->m /= MASK_M*(_Lm<<1);

  // Sync and scale mating preference
  //for(i->f=0.,k=0; k<_Lf; k++)
  //  i->f += ((i->f0[k/LPW_F]>>(_Af*(k%LPW_F)))&MASK_F) + ((i->f1[k/LPW_F]>>(_Af*(k%LPW_F)))&MASK_F);
  //i->f /= MASK_F*(_Lf<<1);

  // Sync and scale sex chromosomes (in interval [0,1], and ONLY 0 and 1)
  //for(i->s=0,k=0; k<_Ln; k++)
  //  i->s += ((i->z0[k/LPW_N]>>(_An*(k%LPW_N)))&MASK_N) + ((i->z1[k/LPW_N]>>(_An*(k%LPW_N)))&MASK_N);
  //i->s = (i->s/(MASK_N*_Ln)) - 0.5;
    
}


// Installs individual i into patch p in a safe manner
static int InstallIndividual(pPatch p, pIndividual i)
{
  if(p->ni == _MAX_INDIVIDUALS) return 0;
  p->i[p->ni++] = i;            return 1;
}


// Returns a pointer to an available individual structure.  NULL for none.
static pIndividual GetNextIndividual(pThreadData td)
{
  // Record where this individual was allocated from
  td->INext[td->nINext]->id = td->nINext;
  return td->INext[td->nINext++];
}


// Releases an Individual back to the its individual array
static void ReleaseNextIndividual(pThreadData td, pIndividual i)
{
  // Swap last in use with i and decriment
  td->INext[i->id]      = td->INext[--td->nINext];
  td->INext[td->nINext] = i;

  // Update allocated from (id) field
  td->INext[i->id]->id   = i->id;
}


// Clears all individuals by restoring the initial sate of individual allocation
static void ReleaseAllCurrentIndividuals(pThreadData td)
{
  td->nICurrent = 0;
  memcpy(td->ICurrent, td->ICurrentCleared, MAX_THREAD_POP*sizeof(pIndividual));

  // this is to shuffle the ind in the memory (to prevent potential first-come-first-serve bias)
  {
    int         i;
    pIndividual ind;

    for(i=0; i<MAX_THREAD_POP; i++) {
      ind = td->ICurrent[(MAX_THREAD_POP-1)-i];
      td->ICurrent[(MAX_THREAD_POP-1)-i] = td->ICurrent[i];
      td->ICurrent[i] = ind;
    }
  }
}

// I'm commenting out these functions as they are not used:
// this will keep the compiler warnings to a minimun.
#if 0
// Clears all individuals by restoring the initial sate of individual allocation
static void ReleaseAllNextIndividuals(pThreadData td)
{
  td->nINext = 0;
  memcpy(td->INext, td->INextCleared, MAX_THREAD_POP*sizeof(pIndividual));
}
#endif


// Copies all data in a safe manner from Individual s to d
static void CopyIndividual(pIndividual d, pIndividual s)
{
  u64b did;

  // Save id field
  did  = d->id;

  // Copy data
  memcpy(d,s,sizeof(Individual));

  // Re-link and restore id field
  LinkIndividual(d);
  d->id  = did;
}

/***************************************************************************************
 * Functions which compute probabilities, properties, or preferences of individuals
 ***************************************************************************************/

// effective population size, the number of individuals present in dispersal range (for dispersal)
static double NeDispersal(pThreadData td, pPatch p)
{
  int     j;
  double  c = 0.;

  // Transform p (which is a pointer into next) into "current" version of p
  p = &(td->Current[p->lat][p->lon]);
    
  // Go through the range of the individual
  for(j=0; j<(p->ddne); j++) {
    // loop through and add all the number of individuals in the roaming range
    c += p->dne[j]->ni;
  }
  return c;
}


// The number of offspring to be produced, this is the seed for Poisson distribution
static double w(pThreadData td, pIndividual i, pPatch p)
{
  double r, Kz;
  
  // Transform p (which is a pointer into next) into "current" version of p
  p = &(td->Current[p->lat][p->lon]);
  
  Kz = _K0;
  
  r = ( ((double)(p->ni)) / ((double)(Kz)) );
  
  return (1.0/(1.0+(((double)(_b))/2.0-1.0)*r));
}


static double CrossFlips(byte x) 
{
  /* initial 0 <--> previous cross point */
  int a = 1-(x&1), i = 7; 

  for ( ; i--; x >>= 1) a += (x&1)^((x>>1)&1);
  return (double) a;
}


//Returns 1 if crossover mask x is valid (does not split a multi-bit loci), 0 on invalid
static int CrossValid(byte x, int abits)
{
  int a,i;

  /* Move through all bits */
  for (i=a=0; x; x>>=1){
    /* Sum the number of ones */
    a += x&1;
    /* Check to see if we are done with a locus here */
    if (!((++i)%abits)){
      /* If this locus is split, return 0 */
      if (a%abits) 
	return 0;
      /* Reset a and keep searching */
      a = 0;
    }
  }

  /* Test last locus and return valid/not valid */
  return ((a%abits)? 0: 1);
}

// Fills in crosover distribution
static void CrossProb(pThreadData td, int w, int fields, int abits)
{
  int i;  double s;

  /* Allocate space for all possible crossover masks */
  td->CrossDist[w] = allocdist(256);
  // note: The C library function void *memset(void *str, int c, size_t n)
  // copies the character c (an unsigned char) to the first n characters of the string pointed to,
  // by the argument str.
  // sizeof(double) = 8 bytes
  // so this is making everything zero?
  memset(td->CrossDist[w]->p, 0, 256*sizeof(double));

  /* Consider all crossover masks */
  for (s=0.,i=0; i<256; i++){
    /* Make sure mask is valid and fill in it's prob */
    // when is it not valid???
    if (CrossValid(i,abits))
      s += td->CrossDist[w]->p[i] = pow(Xi[w],CrossFlips(i)) * pow(1.-Xi[w],fields-CrossFlips(i));
    //printf("w, Xi[w], i, s, td->CrossDist[w]->p[i]: %i, %lf, %i, %0.18llf, %0.18llf\n",w, Xi[w], i,s, td->CrossDist[w]->p[i]);
    }
    /* Call initdist to initialize this distribution */
    /*
     in random.c it says: "Note: d->p must have d->n elements which sum to s on entry to initdist.
    The elements of d->p and d->a are overwritten by the initialization process."
    */
  initdist(td->CrossDist[w], s);
}


// Initialize crossover
static void CrossInit(pThreadData td)
{
  // Initialize distributions ->e and ->p
  
  CrossProb(td, 0, FIELDS_E, _Ae);
  // Initialize distributions ->m and ->k
  CrossProb(td, 1, FIELDS_M, _Am);
  CrossProb(td, 2, FIELDS_K, _Ak);
  // Initialize distributions ->f
  CrossProb(td, 3, FIELDS_F, _Af);
  // Initialize distributions ->z
  CrossProb(td, 4, FIELDS_N, _An);
}


// Returns a pointer to the patch ox,oy spaces from x,y, or NULL if out of bounds
// This is not a generic function. It is mostly specific to NeighboringDispersePatches()
static pPatch ValidNeighborPatch(pThreadData td, u64b x, u64b y, int ox, int oy)
{

  // Width bounds (Only bound first left and last right)
  if( (x == 0)        && (ox < 0) )     return NULL;
  if( (x == WIDTH-1)  && (ox > 0) )     return NULL;
  if( ((x + ox) < 0) )                  return NULL;
  if( ((x + ox) > (WIDTH-1)) )          return NULL;

  // Height bounds
  if( (y == 0)        && (oy < 0) )     return NULL;
  if( (y == HEIGHT-1) && (oy > 0) )     return NULL;
  if( ((y + oy) < 0) )                  return NULL;
  if( ((y + oy) > (HEIGHT - 1)) )       return NULL;
  
  return &td->Next[x+ox][y+oy];
}


// Returns a pointer to the patch ox,oy spaces from x,y, or NULL if out of bounds
// This is not a generic function. It is mostly specific to NeighboringDispersePatches()
static int ValidMatePatch(pThreadData td, u64b x, u64b y, int ox, int oy)
{

  // Width bounds (Only bound first left and last right)
  if( (x == 0)        && (ox < 0) )     return 0;
  if( (x == WIDTH-1)  && (ox > 0) )     return 0;
  if( ((x + ox) < 0) )                  return 0;
  if( ((x + ox) > (WIDTH-1)) )          return 0;

  // Height bounds
  if( (y == 0)        && (oy < 0) )     return 0;
  if( (y == HEIGHT-1) && (oy > 0) )     return 0;
  if( ((y + oy) < 0) )                  return 0;
  if( ((y + oy) > (HEIGHT - 1)) )       return 0;

  return 1;
}


// Fills in the neighboring dispersal patches field of td->Next[x][y].
static void NeighboringDispersePatches(pThreadData td, u64b x, u64b y)
{

  int i,j;
  int dist=0;
  int count=0;
  
  for(td->Next[x][y].ddne=0,i=-((_DRANGE - 1.) / 2.); i<((_DRANGE - 1.) / 2.)+1.; i++){
    for(j=-((_DRANGE - 1.) / 2.); j<((_DRANGE - 1.) / 2.)+1.; j++) {
      if( (td->Next[x][y].dne[td->Next[x][y].ddne] = ValidNeighborPatch(td,x,y,i,j)) ){
        // density and distance dependent based dispersal
#if _DISTANCEDENSITYDISPERSAL
        // find distance
        if( abs(i) > abs(j) ){
          dist = abs(i);
        } else {
          dist = abs(j);
        }
      
        //double denom = 0.;
        double range = 0.;
        range = (((double)(_DRANGE)) - 1.)/2. + 1;
        
        td->Next[x][y].disdist[count] = dist; // distance data
        count++;
#endif
        td->Next[x][y].ddne++;
      }
    }
  }
}


// Swaps Current and Next generations and deps
static void SwapGenerations(pThreadData td)
{

  u64b t;  pIndividual *tp;  Space ts;  

  // Swap spaces
  ts = td->Current;
  td->Current = td->Next;
  td->Next = ts;

  // Swap Individual resources
  t = td->nICurrent;
  td->nICurrent = td->nINext;
  td->nINext = t;

  tp = td->ICurrentCleared;
  td->ICurrentCleared = td->INextCleared;
  td->INextCleared = tp;

  tp = td->ICurrent;
  td->ICurrent = td->INext;
  td->INext = tp;
}


// Initializes any data structures needed by a thread (mostly [if not all] local)
static void InitThread(pThreadData td)
//void InitThread(pThreadData td)
{
  u64b i,j,k,l;
  
  // Initialize Individuals and set up the "cleared" states
  ZALLOC(td->ICurrentData,   MAX_THREAD_POP*sizeof(Individual),  "Could not allocate IndividualData! (Cur)\n");
  ALLOC(td->ICurrent,        MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Individuals! (Cur)\n");
  ZALLOC(td->INextData,      MAX_THREAD_POP*sizeof(Individual),  "Could not allocate IndividualData! (Next)\n");
  ALLOC(td->INext,           MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Individuals! (Next)\n");
  ALLOC(td->ICurrentCleared, MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Cleared State! (Cur)\n");
  ALLOC(td->INextCleared,    MAX_THREAD_POP*sizeof(pIndividual), "Could not allocate Cleared State! (Next)\n");
  for(td->nICurrent=td->nINext=i=0; i<MAX_THREAD_POP; i++){
    td->ICurrent[i]      = td->ICurrentData+i;
    td->INext[i]         = td->INextData+i;
    LinkIndividual(td->ICurrent[i]);
    LinkIndividual(td->INext[i]);
  }
  memcpy(td->ICurrentCleared, td->ICurrent, MAX_THREAD_POP*sizeof(pIndividual));
  memcpy(td->INextCleared,    td->INext,    MAX_THREAD_POP*sizeof(pIndividual));

  // Malloc for our patch spaces
  ALLOC(td->Current,(WIDTH+2)*sizeof(Patch*),   "Could not allocate Current patch space!\n");
  for(i=0; i<(WIDTH+2); i++)
    ZALLOC(td->Current[i],HEIGHT*sizeof(Patch), "Could not allocate Current patch space!\n");
  ALLOC(td->Next,(WIDTH+2)*sizeof(Patch*),      "Could not allocate Next patch space!\n");
  for(i=0; i<(WIDTH+2); i++)
    ZALLOC(td->Next[i],HEIGHT*sizeof(Patch),    "Could not allocate Next patch space!\n");

  // Initialize the random number generator
  initrand(&R, Seed);

  
  // Figure out what patches to copy from start
#if _CPUS == 1
  // Only one thread:    Do not copy any wandering zones
  k = 1, l = WIDTH+1;
#else
  switch(td->id){
  case 0:         // Left-most thread:   Only copy right wandering zone
    k = 1, l = WIDTH+2;
    break;
  case (_CPUS-1): // Right-most thread:  Only  copy left wandering zone
    k = 0, l = WIDTH+1;
    break;
  default:        // Middle threads:     Copy both wandering zones
    k = 0, l = WIDTH+2;
  }
#endif

  // Actually copy them
  for(i=0; i<WIDTH; i++){
    // Copy into current
    memcpy(td->Current[i],    Start[i],          HEIGHT*sizeof(Patch));
    // Next matches current (no individuals)
    memcpy(td->Next[i],       td->Current[i],    HEIGHT*sizeof(Patch));
    for(j=0; j<HEIGHT; j++) td->Next[i][j].ni = 0;
  }

  // Fill in patch individuals
  for(i=0; i<WIDTH; i++){
    // *for(i=1; i<(WIDTH+1); i++) {
    // Copy individuals from start
    for(j=0; j<HEIGHT; j++){
      for(k=0; k<td->Current[i][j].ni; k++){
        //printf("\nenters the loop. start: %d, current: %d", Start[i][j].ni, td->Current[i][j].ni);
        CopyIndividual(td->Current[i][j].i[k]=GetCurrentIndividual(td),Start[i][j].i[k]);
        // Start[][]'s individuals are unsynced
        SyncIndividual(td,td->Current[i][j].i[k]);
      }
    }
  }

  // Now that patches are copied, pre-compute neighbors, but Since NeighboringPatches()
  // only fills in the Next generation, swap and process twice.
  SwapGenerations(td);
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++){
      NeighboringDispersePatches(td,i,j);
    }
  }

  SwapGenerations(td);
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++){
      NeighboringDispersePatches(td,i,j);
    }
  }

  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++) {
      td->Current[i][j].nothers=0;
    }
  }

  // calculates who is using which patch for how long
  for(i=0; i<WIDTH; i++){
    for(j=0; j<HEIGHT; j++) {
      for(k=0,l=td->Current[i][j].ni; k < l; k++) {
        TimeSpentRoam(td->Current[i][j].i[k], &td->Current[i][j]);
      }
    }
  }

  // Setup recombination  distributions
  CrossInit(td);
}

#if ENDIAN == NICHE_BIG_ENDIAN
// This macro will do an index translation that should be 
// functionally equvalent to byte-swapping.  
static int BSWAP[32] = {
 7,  6,  5,  4,  3,  2,  1,  0,
 15, 14, 13, 12, 11, 10, 9,  8,
 23, 22, 21, 20, 19, 18, 17, 16,
 31, 30, 29, 28, 27, 26, 25, 24
};
#define Endian(x) (BSWAP[x])
#else
#if ENDIAN == NICHE_LITTLE_ENDIAN
#define Endian(x) (x)
#else
#error ENDIAN _must_ be set to either NICHE_BIG_ENDIAN or NICHE_LITTLE_ENDIAN!
#endif
#endif


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineEco(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)){ t=x; x=y; y=t;}

#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_E; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[0]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_E; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[0]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineNeutral(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
  
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_N; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[4]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_N; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[4]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineMating(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
    
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_M; i++){
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[1]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_M; i++){
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[1]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#endif
}


// Handles recombination and mutation for a single destinaation chromosome
static void RecombineMarker(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }
    
#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_K; i++) {
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[2]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_K; i++) {
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[2]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#endif
}


// Handles recombination and mutation for a single destination chromosome
static void RecombineFPref(pThreadData td, byte *z, byte *x, byte *y)
{
  u64b i,r;  byte *t;
  
  // Swap sources
  if(rnd(&R, 2)) { t=x; x=y; y=t; }

#if SKIP_MUTATION
  // Move through each byte
  for(i=0,r=128; i<PLENGTH_F; i++) {
    // Obtain a crossover mask for this byte and recombine
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[3]);
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#else
  // Normal mutation
  for(i=r=0; i<PLENGTH_F; i++) {
    // Obtain a crossover mask, recombine and mutate
    r = (((r>>7)&1)-1) ^ drand(&R, td->CrossDist[3]);
    // No mutation
    z[Endian(i)] = (x[Endian(i)] & r) | (y[Endian(i)] & ~r);
  }
#endif
}


// Combines parents m and f and fills in child c with correct recombination and mutation
static void FastRecombine(pThreadData td, pIndividual c, pIndividual m, pIndividual f)
{
  
  if (m->s != 0) {
    Error("!!there is an error in FastRecombine(), m should be a male: m.sex: %d, f.sex: %d\n", m->s, f->s);
  }
  
  // Even though these are called Recombine, it only deals with
  // chromosomal recombination shuffling loci on the same 'chromosome' i.e. trait,
  // independent traits that are on different chromosomes, i.e. are independent anyway
  // there is no cross-over between loci (there is one loci per trait so not relevant)
  // but randomly allocates alleles to offspring and skews allocation in GeneDrive
  RecombineEco      (td, (byte*)(c->x0), (byte*)(m->x0), (byte*)(m->x1));
  RecombineEco      (td, (byte*)(c->x1), (byte*)(f->x0), (byte*)(f->x1));
  RecombineMarker   (td, (byte*)(c->k0), (byte*)(m->k0), (byte*)(m->k1));
  RecombineMarker   (td, (byte*)(c->k1), (byte*)(f->k0), (byte*)(f->k1));
  RecombineNeutral  (td, (byte*)(c->z0), (byte*)(f->z0), (byte*)(f->z1));
  RecombineNeutral  (td, (byte*)(c->z1), (byte*)(m->z0), (byte*)(m->z1));
  RecombineFPref    (td, (byte*)(c->f0), (byte*)(m->f0), (byte*)(m->f1));
  RecombineFPref    (td, (byte*)(c->f1), (byte*)(f->f0), (byte*)(f->f1));
  RecombineMating   (td, (byte*)(c->m0), (byte*)(m->m0), (byte*)(m->m1));
  RecombineMating   (td, (byte*)(c->m1), (byte*)(f->m0), (byte*)(f->m1));
    
}


#if _DISTANCE
void RecordDistance(double density, int distance)
{
  if( (density >= 0.0) && (density < 0.1) ) { distancedata[0][distance]++; }
  if( (density >= 0.1) && (density < 0.2) ) { distancedata[1][distance]++; }
  if( (density >= 0.2) && (density < 0.3) ) { distancedata[2][distance]++; }
  if( (density >= 0.3) && (density < 0.4) ) { distancedata[3][distance]++; }
  if( (density >= 0.4) && (density < 0.5) ) { distancedata[4][distance]++; }
  if( (density >= 0.5) && (density < 0.6) ) { distancedata[5][distance]++; }
  if( (density >= 0.6) && (density < 0.7) ) { distancedata[6][distance]++; }
  if( (density >= 0.7) && (density < 0.8) ) { distancedata[7][distance]++; }
  if( (density >= 0.8) && (density < 0.9) ) { distancedata[8][distance]++; }
  if( (density >= 0.9) && (density < 1.0) ) { distancedata[9][distance]++; }
  if( (density >= 1.0) && (density < 1.1) ) { distancedata[10][distance]++; }
  if( (density >= 1.1) )                    { distancedata[11][distance]++; }
}


void LogDistance()
{

  int     x,y;
  double  maxdis = 0.;
  
  maxdis = ( ((double)(_DRANGE)) - 1. )/2. + 1.;
  
  // Gen+1 is for viewing
  // otherwise since dispersal occurs at gen 0, it is viewed at gen 0
  // we want to see it 'after' it happened at gen 1
  
  if( !(Gen%_DISTANCE) )
    fprintf(Distancef, "gen: %llu ", Gen+1);
  for( x=0; x<12; x++){
    //i-- > 0 ;
    //for( y=0; y<((int)(maxdis)); y++){
    for( y= ((int)(maxdis)); y-- >0;){
      fprintf(Distancef, "%d ", distancedata[x][y]);
      //printf("%d ", distancedata[x][y]);
    }
  }
  fprintf(Distancef, "\n");
  fflush(Distancef);
 // Increment counter
 Distance_logs++;
}
#endif

#if _SCENTSCAPE
void LogScentscape()
{
    int     x,i,j;
    
    // Gen+1 is for viewing
    // otherwise since dispersal occurs at gen 0, it is viewed at gen 0
    // we want to see it 'after' it happened at gen 1
    
    //if( !(Gen%_SCENTSCAPE) )
    //    fprintf(Scentscapef, "gen: %llu \n", Gen);
    for( x = 0; x<2; x++){
        //two sexes
        //printf("*%d\n", x);
        for( i = 0; i<_WIDTH; i++){
            for( j = 0; j<_HEIGHT; j++){
                //printf("(%d, %d) %d - %lf\n", i, j, i*_HEIGHT+j, Scentdata[x][i*_HEIGHT+j]);
                fprintf(Scentscapef, "%llu, %d, %d, %d, %d, %lf\n", Gen, x, i, j, i*_HEIGHT+j, Scentdata[x][i*_HEIGHT+j]);
            }
        }
        //fprintf(Scentscapef, "\n");
    }
    fflush(Scentscapef);
    // Increment counter
    Scentscape_logs++;

}
#endif


// dispersal is done w.r.t. one of the neighboring patches of the mother's
// this could be modified to make the dispersal from the patch where mating took place
static void Disperse(pThreadData td, pIndividual i, pPatch p)
{
  int     distance, k, scent;
  static  distribution *d  = NULL;
  static  distribution *cd = NULL;
  double  s, sd, maxdis, x, normalizeddis;
  
  x             = 0.;
  normalizeddis = 0.;
  maxdis        = ((double)((_DRANGE - 1.) / 2.));
  
#if _DISTANCEDENSITYDISPERSAL
  
  d = allocdist(((int)(maxdis))+1);
  d->n = ((int)(maxdis))+1;
 
  cd = allocdist(p->ddne);
  cd->n = p->ddne;
  
  // Transform p (which is a pointer into next) into "current" version of p to get density
  p = &(td->Current[p->lat][p->lon]);
  
  
  // density and distance dependent based dispersal
  for( s=0., k=0; k<((int)(maxdis))+1; k++ ){
    d->p[k] = 0.0;
    normalizeddis = ((double)(k))/maxdis;
    x = (_DISPERSALCOEF*normalizeddis - _DENSITYCOEF*p->density); //NDDD
    //x = (_DISPERSALCOEF*normalizeddis - (1-_DENSITYCOEF*p->density)); // PDDD
    s += ( d->p[k] = Exp(x*x) );
  }
/*
  //***
    //if occupied, skip; otherwise assign probability
    // works but does not prevent if an ind already picked it in the same day,
    // since they are written in 'Next', and biologically hard to know as well, right?
    if( p->dne[k]->ni > 0){
      sd += ( cd->p[k] = 0. );
        //if( p->lat == 208 && p->lon == 96 ){
        //if( Gen == 588 ){
          //  printf("in ni>0, sum = %lf (%d - %lf), from (%llu, %llu) to (%llu, %llu) N: %llu\n", sd, k, cd->p[k], p->lat, p->lon, p->dne[k]->lat, p->dne[k]->lon, p->dne[k]->ni);
        //}
        
    } else {
      // 1. calculate distance
      distance = p->distance[k];
      
      // 2. check scent
      if( i->s == 0 ){
        if( p->dne[k]->malscentstrength == 0 ){
          //sd += ( cd->p[k] = 1. );
          scent = 0;
        } else {
          //sd += ( cd->p[k] = (_MAXSCENT-p->dne[k]->malscentstrength));
          scent = _MAXSCENT - (p->dne[k]->malscentstrength);
        }
      } else {
        if( p->dne[k]->femscentstrength == 0 ){
          //sd += ( cd->p[k] = 1. );
          scent = 0;
        } else {
          //sd += ( cd->p[k] = (_MAXSCENT-p->dne[k]->femscentstrength) );
          scent = _MAXSCENT - (p->dne[k]->malscentstrength);
        }
      }
      
    */
    
  //***
  if( s != 0. ){
    initdist(d,s);
    distance = drand(&R, d);
  }
  
  for( sd=0.,k=0; k<p->ddne; k++ ){
    if( distance == p->disdist[k] ){
      sd += ( cd->p[k] = 1. );
    } else {
      sd += ( cd->p[k] = 0.);
    }
  }
  
  k = 0;
  if( sd != 0. ){
    initdist(cd,sd);
    k = drand(&R, cd);
  }
  
  //if( Gen == 112 ){ printf("(%d, %d) - (%d,%d)\n", p->lat, p->lon, p->dne[k]->lat, p->dne[k]->lon); }
#if _DISTANCE
  RecordDistance(p->density, distance);
#endif
#endif

#if _RAND_DISPERSAL
  // dispersal is random to one
  // Child will randomly to one of the patches within the range of the parent (female)
  d->n = p->ddne;
  k = rnd(&R, p->ddne);

  if( abs(p->dne[k]->lat - p->lat) > abs(p->dne[k]->lon - p->lon)){
    distance = abs(p->dne[k]->lat - p->lat);
  } else {
    distance = abs(p->dne[k]->lon - p->lon);
  }

#if _DISTANCE
  RecordDistance(p->density, distance);
#endif
  
#endif
  
  freedist(d);
  freedist(cd);
  // Transform p back to "next"
  p = &(td->Next[p->lat][p->lon]);
  
  if(!InstallIndividual(p->dne[k],i)) {
    ReleaseNextIndividual(td,i);
    printf("release next\n");
    Error("!! HITTING UPPER BOUND ON MAX_INDS\n");
  }
}


// dispersal is done w.r.t. one of the neighboring patches of the mother's
// this could be modified to make the dispersal from the patch where mating took place
static void MateSearchMovement(pThreadData td, pIndividual i, pPatch p)
{
  int     distance, k;
  static  distribution *d  = NULL;
  static  distribution *cd = NULL;
  double  s, sd, maxdis, x, normalizeddis;
  
  x             = 0.;
  normalizeddis = 0.;
  maxdis  = ((double)((_DRANGE - 1.) / 2.));

  
#if _DISTANCEDENSITYDISPERSAL
  
  d = allocdist(((int)(maxdis))+1);
  d->n = ((int)(maxdis))+1;
 
  cd = allocdist(p->ddne);
  cd->n = p->ddne;
  
  // Transform p (which is a pointer into next) into "current" version of p to get density
  p = &(td->Current[p->lat][p->lon]);
  
  // density and distance dependent based dispersal
  for( s=0., k=0; k<((int)(maxdis))+1; k++ ){
    d->p[k] = 0.0;
    normalizeddis = ((double)(k))/maxdis;
    x = (_DISPERSALCOEF*normalizeddis - _DENSITYCOEF*p->density); // NDDD
    //x = (_DISPERSALCOEF*normalizeddis - (1-_DENSITYCOEF*p->density)); //PDDD
    s += ( d->p[k] = Exp(x*x) );
  }

  if( s != 0. ){
    initdist(d,s);
    distance = drand(&R, d);
  }

  for( sd=0.,k=0; k<p->ddne; k++ ){
    if( distance == p->disdist[k] ){
      sd += ( cd->p[k] = 1. );
    } else {
      sd += ( cd->p[k] = 0.);
    }
  }
  
  k = 0;
  if( sd != 0. ){
    initdist(cd,sd);
    k = drand(&R, cd);
  }
  
  //if (Gen == 52) printf("from (%d, %d) to (%d, %d)\n", p->lat, p->lon, p->dne[k]->lat, p->dne[k]->lon);
    
#if _DISTANCE
  RecordDistance(p->density, distance);
#endif
#endif

#if _RAND_DISPERSAL
  // dispersal is random to one
  // Child will randomly to one of the patches within the range of the parent (female)
  d->n = p->ddne;
  k = rnd(&R, p->ddne);
#endif

  freedist(d);
  freedist(cd);
  // Transform p back to "next"
  p = &(td->Next[p->lat][p->lon]);
  
  if(!InstallIndividual(p->dne[k],i)) {
    ReleaseNextIndividual(td,i);
    printf("release next\n");
    Error("!! HITTING UPPER BOUND ON MAX_INDS\n");
  }
}


#if _SCENTSCAPE
// dispersal is done w.r.t. one of the neighboring patches of the mother's
// this could be modified to make the dispersal from the patch where mating took place
static void ExtensiveMateSearchMovement(pThreadData td, pIndividual i, pPatch p)
{
  int     maxdistance, distance, k, x, targetlat, targetlon;
  static  distribution *cd = NULL;
  double  sd, scent, adjustedscent;
  
  //printf("in ExtensiveMateSearch()\n");
  
  maxdistance = ((double)((_DRANGE - 1.) / 2.));
  
  cd = allocdist(p->ddne);
  cd->n = p->ddne;
  
  // Transform p (which is a pointer into next) into "current" version of p to get density
  p = &(td->Current[p->lat][p->lon]);
  
  for( sd=0.,k=0; k<p->ddne; k++ ){
      // find the patches with maximum distance on the edge
      if( maxdistance == p->disdist[k] ){
          // find the relevant locations on Scentdata
          targetlat = p->dne[k]->lat;
          targetlon = p->dne[k]->lon;
          // based on opposite of sex i->s
          if( i->s == 0 ){ scent = Scentdata[1][targetlat * _HEIGHT + targetlon]; }
          if( i->s == 1 ){ scent = Scentdata[0][targetlat * _HEIGHT + targetlon]; }
          sd += ( cd->p[k] = (scent*scent*scent*scent) );
          //printf("from (%d, %d) to (%d, %d), s: %d, %d/%d, scent: %lf, sd: %lf\n", p->lat, p->lon, targetlat, targetlon, i->s, k+1, p->ddne, scent, sd);
      } else {
        sd += ( cd->p[k] = 0.);
      }
  }
  
  x = 0;
  if( sd != 0. ){
    initdist(cd,sd);
    x = drand(&R, cd);
  } else {
    x = rnd(&R, k);
  }
  
  if( abs(p->lat - p->dne[x]->lat) > abs(p->lon - p->dne[x]->lon) ){
      distance = abs(p->lat - p->dne[x]->lat);
  } else {
      distance = abs(p->lon - p->dne[x]->lon);
  }
  
  //printf("picked: %d (%d, %d) - distance: %d\n", x, p->dne[x]->lat, p->dne[x]->lon, distance);

  // HERE calculate distance to record
  RecordDistance(p->density, distance);

  freedist(cd);
  // Transform p back to "next"
  p = &(td->Next[p->lat][p->lon]);
  
  if(!InstallIndividual(p->dne[x],i)) {
    ReleaseNextIndividual(td,i);
    printf("release next\n");
    Error("!! HITTING UPPER BOUND ON MAX_INDS\n");
  }
}
#endif


// mates m and f, and there is fertility selection
static void Mate(pThreadData td, pIndividual f, pIndividual m, pPatch p)
{

  pIndividual i;
  int         j, rn;
  double      nc = 0.;
  u64b        l;
  double randomno = 0.;
  int g, cullflag = 0;
  
  //printf("in Mate()\n");
  //if ( *(m->m0) == 1 || *(m->m1) == 1 ) Error("!!*This ERROR shouldn't print, check Mate(), toxicmale should get into this function\n");
  //if ( f->toxicpregnancy == 1 )         Error("!!*This ERROR shouldn't print, check Mate(), toxicpregrantfemale should get into this function\n");
  if ( m->s == 1 )                      Error("!!*This ERROR shouldn't print, check Mate() m is NOT a male\n");
  if ( f->s == 0 )                      Error("!!*This ERROR shouldn't print, check Mate() f is NOT a female\n");
  if ( *(f->m0) == 1 || *(f->m1) == 1 ) Error("!!*This ERROR shouldn't print, check Mate() f cannot be a carrier\n");
  
    nc = _b*w(td,f,p)*(1/_PROBMATING);
  rn = Poisson(&R, nc);
  //printf("off. no: %d\n", rn);
    
  for( g = 0; g < _NoCull; g++ ){
    if( Gen == CullGen[g] ) {
      cullflag = 1;
    }
  }
  //if( Gen == 112){ printf("off. nu.:%d\n", rn); }
  for(j=0; j<((int)(rn)); j++) {
      randomno = U01(&R);
      if( cullflag == 1 && randomno > _CullOffspringSurvival ){
      //printf("offspring is not produced/culling\n");
      continue;
    } else {
      //if( Gen >= 54){ printf("off. nu.: %d of %d\n", j+1, rn); }
      // Get an individual to use as a child
      i = GetNextIndividual(td);
        
      // Set the child's sex (not used)
      i->s = rnd(&R, 2);
      
      if( (*(m->m0) == 1 || *(m->m1) == 1 )){
        //printf("male: (%d, %d), offspring's s: %d\n", *(m->m0), *(m->m1), i->s);
        if(i->s == 1){
          //printf("lethal female offspring\n");
          continue;
        }
      }
    
      i->age = 0;
      i->uniqueid = UniqueCounter;
      UniqueCounter++;
            
      // Recombine and mutate
      FastRecombine(td, i, m, f);

      //if( Gen >= 54 ) printf("f: (%d/%d), m: (%d, %d), o[%d]: (%d/%d)\n", *(m->m0), *(m->m1), *(f->m0), *(f->m1), i->s, *(i->m0), *(i->m1));
      //if( *(f->m0) == 1 || *(f->m1) == 1 ){
      //  printf("f: (%d/%d), m: (%d, %d), o[%d]: (%d/%d)\n", *(m->m0), *(m->m1), *(f->m0), *(f->m1), i->s, *(i->m0), *(i->m1));
      //}
      // Sync, Handle dispersal and install child
      SyncIndividual(td,i);
      Disperse(td, i, p);
      OffspringPopSize++;
    }
  }
  MonogamousMating++;
}


static int SpermCompIndex(pIndividual m_one, pIndividual m_two)
{
  // this is to assign probabilities of siring, based on specific genotypes of m_one and m_two
  // option 1: wildtype, wildtype or drive, drive - default first male adv
  // option 2: wildtype, drive - increased first male adv, or increased drive adv if order is not imp
  // option 3: drive, wildtype - decreased first male adv, or increased drive adv if order is not imp
  // initialize the spermcomp index
  int spermindex = 0;
  
  // find genotypes of m_one and m_two based on strategies
  // if wildtype first then t-allele (option 2 - increased first male adv)
  if( (*(m_one->m0) == 0 && *(m_one->m1) == 0) && (*(m_two->m0) != 0 || *(m_two->m1) != 0) ){
    spermindex = 1;
  } else if( (*(m_one->m0) != 0 || *(m_one->m1) != 0) && (*(m_two->m0) == 0 && *(m_two->m1) == 0) ){
    // if t-allele first then wildtype (option 3 - decreased first male adv)
    spermindex = -1;
  } else {
    // both males wildtype; 2: both males are t-carrying (at least one nonzero allele);
    spermindex = 0;
  }
  return spermindex;
}


// mates m and f, and there is fertility selection
static void PolyMate(pThreadData td, pIndividual f, pIndividual m_one, pIndividual m_two, pPatch p)
{
  
  pIndividual i;
  pIndividual m;  // this is the siring male
  int         j, rn, paternityindex;
  double      nc = 0.;
  u64b        l;
  double      randomno = 0.;
  int         g, cullflag = 0;
   
  //printf("*in PolyMate()\n");
  if ( m_one->s == 1 )                  Error("!!*This ERROR shouldn't print, check PolyMate() m_one is NOT a male\n");
  if ( m_two->s == 1 )                  Error("!!*This ERROR shouldn't print, check PolyMate() m_two is NOT a male\n");
  if ( f->s == 0 )                      Error("!!*This ERROR shouldn't print, check PolyMate() f is NOT a female\n");
  if ( *(f->m0) == 1 || *(f->m1) == 1 ) Error("!!*This ERROR shouldn't print, check PolyMate() f cannot be a carrier\n");
    //if ( *(m_one->m0) == 1 || *(m_one->m1) == 1 ) Error("!!*This ERROR shouldn't print, check PolyMate(), toxic male should get into this function\n");
  //if ( *(m_two->m0) == 1 || *(m_two->m1) == 1 ) Error("!!*This ERROR shouldn't print, check PolyMate(), toxic male should get into this function\n");
  //if ( f->toxicpregnancy == 1 )         Error("!!*This ERROR shouldn't print, check PolyMate(), toxicpregnant should get into this function\n");
    
  nc = _b*w(td,f,p)*(1/_PROBMATING);
  rn = Poisson(&R, nc);
  //printf("off. no (poly): %d\n", rn);
    
  for( g = 0; g < _NoCull; g++ ){
    if( Gen == CullGen[g] ) {
      cullflag = 1;
    }
  }
    
  paternityindex = 0;
  //if( Gen >= 54){ printf("off. nu.:%d\n", rn); }
  for(j=0; j<((int)(rn)); j++){
      randomno = U01(&R);
      if( cullflag == 1 && randomno > _CullOffspringSurvival ){
          //printf("offspring is not produced/culling\n");
          continue;
      } else {
          //if( Gen >= 54){ printf("off. nu.: %d of %d\n", j+1, rn); }
          // Get an individual to use as a child
          i = GetNextIndividual(td);
        
          // Set the child's sex
          i->s = rnd(&R, 2);
          i->age = 0;
          i->uniqueid = UniqueCounter;
          UniqueCounter++;
            
          // pick the father here amon n=2 males and initialize to m and continue as before
          int spermindex = 0;
          int n = 2;
          int siringmaleinteger = 3;
          double siringprobsum = 0.0;
            
          // initialize the siring male probability distribution
          distribution *mdist;
          mdist = allocdist(n);
          mdist->n = 2;
            
          spermindex = SpermCompIndex(m_one,m_two);
            
          siringprobsum += (mdist->p[0] = _FIRSTSPERMADVANTAGE + (spermindex * _SPERMCOMPETITIONCOEF));
          siringprobsum += (mdist->p[1] = (1 - _FIRSTSPERMADVANTAGE) - (spermindex * _SPERMCOMPETITIONCOEF));
            
          initdist(mdist,siringprobsum);
          siringmaleinteger = rnd(&R, 2);
            
          if( siringmaleinteger == 3 ){
              Error("!!something's wrong in PolyMate(), not picking a male, why?\n");
          }
          if( siringmaleinteger == 0 ){
              m = m_one;
          }
          if( siringmaleinteger == 1 ){
              m = m_two;
          }
           
          if( (*(m->m0) == 1 || *(m->m1) == 1 )){
            //printf("male: (%d, %d), offspring's s: %d\n", *(m->m0), *(m->m1), i->s);
            if(i->s == 1){
              //printf("lethal female offspring\n");
              continue;
            }
          }
          
          // use this index to calculate multiple paternity if all the father are same,
          // paternity index==0 if if sired by first male; or paternity index==numberofoffpsring if sired by the latter.
          // if not equal, then multiple paternity
          paternityindex += siringmaleinteger;
            
          freedist(mdist);
          FastRecombine(td, i, m, f);
            
          //if( Gen >= 54 ) printf("f: (%d/%d), m: (%d, %d), o[%d]: (%d/%d)\n", *(m->m0), *(m->m1), *(f->m0), *(f->m1), i->s, *(i->m0), *(i->m1));
          //if( *(f->m0) == 1 || *(f->m1) == 1 ){
          //  printf("f: (%d/%d), m: (%d, %d), o[%d]: (%d/%d)\n", *(m->m0), *(m->m1), *(f->m0), *(f->m1), i->s, *(i->m0), *(i->m1));
          //}
          // Sync, Handle dispersal and install child
          SyncIndividual(td,i);
          Disperse(td, i, p);
          OffspringPopSize++;
      }
  }
  
  if( paternityindex != 0 && paternityindex != rn ){
      MultiplePaternity++;
  } else {
      SinglePaternity++;
  }
}


// survival of an individual
static void Survive(pThreadData td, pPatch p, pIndividual m)
{
  pIndividual i;
  u64b        l;
  
  // copy all details
  i = GetNextIndividual(td);
  for(l=0; l<T_LENGTH; l++){
    i->d[l] = m->d[l];
  }
  i->s = m->s;
  
  i->uniqueid = m->uniqueid;
  i->age = (m->age+1); //increment age, here since copying to Next
  //if (i->s == 1){
  //  i->toxicpregnancy = m->toxicpregnancy;
  //} else {
  //  i->toxicpregnancy = 0;
  //}
  SyncIndividual(td,i);
  MateSearchMovement(td, i, p);
}

#if _SCENTSCAPE
// survival of an individual
static void ExtensiveSearchSurvive(pThreadData td, pPatch p, pIndividual m)
{
  pIndividual i;
  u64b        l;
  
  // copy all details
  i = GetNextIndividual(td);
  for(l=0; l<T_LENGTH; l++){
    i->d[l] = m->d[l];
  }
  i->s = m->s;
  
  i->uniqueid = m->uniqueid;
  i->age = (m->age+1); //increment age, here since copying to Next
  //if (i->s == 1){
  //  i->toxicpregnancy = m->toxicpregnancy;
  //} else {
  //  i->toxicpregnancy = 0;
  //}
  SyncIndividual(td,i);
  ExtensiveMateSearchMovement(td, i, p);
}
#endif


#if _INOC
// mates m and f, and there is fertility selection

static void Inoculate(pThreadData td, pPatch p, int inocindividuals)
{
  
  pIndividual i;
  int         k,j;
  u64b        l,c;
  
  for( j=0; j<inocindividuals; j++ ) {
      i = GetNextIndividual(td);
      i->s = 0;
      for(l=0; l<T_LENGTH; l++) i->d[l] = 0;

      // this is trangenic genotype m1m1 (homozygous transgenic)
      for(c=0; c<_Lm; c++)    SetLocus(i->m0, c, 1, 1);
      for(c=0; c<_Lm; c++)    SetLocus(i->m1, c, 1, 1);
            
#if _MALEDRIVE
      // here _Ln is used as a sex chromosome; heteroz (0,1) are males; homoz. (1,1) are females.
      // gene-drive could be male only, or two sex
      // male only drive
      //i->s = 0;
      //********
      i->s = 0;
#else
      i->s = rnd(&R, 2);
#endif
      i->age = 0;
                        
      SyncIndividual(td,i);
      MateSearchMovement(td, i, p);
      
      //printf("inoc i[%d]: (%d/%d)\n", i->s, *(i->m0), *(i->m1));
  }
}
#endif


// Does all the real work involved in the simulation.
void* THREAD(void *arg)
{

  ThreadData    td;
  u64b          h,i,j,k,l;
  int           g,b,matingflagingeneration;
  distribution *fdist;
  distribution *pdist;
  pIndividual   secondmatingmale;
  pIndividual   matingfemale;
  pIndividual   matingmale;
  pIndividual   female;
  pPatch        matingpatch;
  double        matingrandomno;
    
  // Initialize this thread's portion of the simulation
  td.id = *((int*)arg);
  Threads[td.id] = &td;
  InitThread(&td);
  
  // two identical distributions, in fdist, after the first male is chosen in polyandry
  // its prob. is reduced to 0 so that it is not picked again
  fdist = allocdist(_MAX_INDIVIDUALS*_DRANGE*_DRANGE);
  pdist = allocdist(_MAX_INDIVIDUALS*_DRANGE*_DRANGE);

  while(SignalDone(&td)){

#if _DISTANCE
    for( i=0; i<12; i++){
      for( j=0; j<((int)(_DRANGE+1)); j++){
        distancedata[i][j] = 0;
      }
    }
#endif

#if _SCENTSCAPE
    //printf("%d\n", Gen);
    for( i=0; i<WIDTH; i++){
        for( j=0; j<HEIGHT; j++){
            Scentdata[0][i*HEIGHT + j] = 0.;
            //printf("%lf", Scentdata[0][i*HEIGHT + j]);
            Scentdata[1][i*HEIGHT + j] = 0.;
            //printf("%lf", Scentdata[1][i*HEIGHT + j]);
        }
    }
#endif
  
      
#if _INOC
    // this is for monitoring!
    // writing the monitoring info to Next
    for(i=0; i<WIDTH; i++){
        for(j=0; j<HEIGHT; j++){
            for( g = 0; g < _NoInoc; g++ ) {
                if( Gen == InocGen[g] ) {
                    //printf("inoc gen: %d\n", Gen);
                    // monitoring and inoculation patches are the SAME ("logistically easy")
                    if( td.Current[i][j].inocpatch == 1 ) {
                        //N's should be calculated in 'next' since this is after all mating survival is done and written to next
                        //printf("%d- (%d %d), N_next: %d\n",b,i,j,td.Next[i][j].ni);
#if _MONITOR
                        if( !td.Current[i][j].ni ) continue;
#if _TOXICMALEDETECTION
                        int toxicmaleno = 0;
                        for(h=0, l=td.Current[i][j].ni; h<l; h++){
                            if( *(td.Current[i][j].i[h]->m0) == 1 || *(td.Current[i][j].i[h]->m1) == 1 ) toxicmaleno++; // count toxic males
                        }
                        //printf("tox.m no: %d - N: %d (diff: %d)\n", toxicmaleno, td.Next[i][j].ni, td.Next[i][j].ni - toxicmaleno);
                        if( td.Current[i][j].ni - toxicmaleno >= 1 ){
                            td.Next[i][j].monitorgen    = 1; // in the next because it will be swapped
                            td.Current[i][j].monitorgen = 0; // let's clean for swap
                            //printf("*(tox. male detection) gen: %d, at (%d, %d)\n", Gen, i,j);
                        }
                      
#else
                        // should not be here uncless there are individuals
                        td.Next[i][j].monitorgen    = 1; // in the next because it will be swapped
                        td.Current[i][j].monitorgen = 0; // let's clean for swap
                        //printf("*(no tox. m. detection) gen: %d, at (%d, %d)\n", Gen, i,j);
#endif

#else
                        // monitor gens and inco gens are same hence it is turned on below even if there is no monitoring
                        td.Next[i][j].monitorgen    = 1; // in the next because it will be swapped
                        td.Current[i][j].monitorgen = 0; // let's clean for swap
                        //printf("*(no monitoring) gen: %d, at (%d, %d)\n", Gen, i,j);
#endif
                      }
                  }
              }
          }
      }
#endif
      

    OffspringPopSize        = 0;
    Mortality               = 0;
    MonogamousMating        = 0;
    SinglePaternity         = 0;
    MultiplePaternity       = 0;
    Surviving               = 0;
    
    matingflagingeneration  = 0;
    
    int numberoffemales     = 0;
    Total_Females    = 0;
    Total_Population = 0;
      
      for(i=0; i<WIDTH; i++){
          for(j=0; j<HEIGHT; j++){
              
              int numberofmales       = 0;
              int mating              = 0;
              int success             = 0;
              int maleinteger         = 0;
              int nothers             = 0;
              double maleprobsum      = 0.;
              matingrandomno          = 0.;
              
#if _DISTANCEDENSITYDISPERSAL
              //td.Current[i][j].density = ( ((double)(NeDispersal(&td,&td.Current[i][j]))) / ((double)(td.Current[i][j].ddne * _K0)) );
              td.Current[i][j].density = ((double)(td.Current[i][j].ni)) / ((double)(_INIT_INDIVIDUALS));
#endif
               
              td.Current[i][j].nf = 0;
              td.Current[i][j].nm = 0;
              
              l=td.Current[i][j].ni;
              Total_Population += l;
              
              for(h=0 ; h<l; h++){          // loops through individuals in a patch
                  
                  //if( *(td.Current[i][j].i[h]->m0) == 1 || *(td.Current[i][j].i[h]->m1) == 1 )
                  //  printf("%llu %llu\n", *(td.Current[i][j].i[h]->m0), *(td.Current[i][j].i[h]->m1));
                  
                  //printf("\n\n[%d, %d] %d/%d, pointer: %llu, s:%d, tp: [%d]\n", i, j, h+1, l, td.Current[i][j].i[h], td.Current[i][j].i[h]->s, td.Current[i][j].i[h]->toxicpregnancy);
                  
                  if( td.Current[i][j].i[h]->s == 0 ){
                      td.Current[i][j].nm++; //increment the number of males
                      continue; // Skip any males we encounter
                  }
                  //printf("here\n");
                  td.Current[i][j].nf++; //increment the number of females
                  Total_Females++;
                  
                  matingrandomno = U01(&R);
                  //printf("%llf/%llf\n", matingrandomno, _PROBMATING);
                  if( matingrandomno > _PROBMATING )                continue; // Skip if random no is greater than prob
                  //if( td.Current[i][j].i[h]->toxicpregnancy == 1 )  continue; // Skip any females with toxic pregnancy (shouldn't be here!)
                  
                  //printf("here\n");
                  matingfemale = td.Current[i][j].i[h];             // the female is fertile, find a mate:
                  success = 0;
                  numberoffemales++;
                  
                  //printf("**found mating female [%d] with m:(%d, %d), (at %d, %d), tp: %d \n", td.Current[i][j].i[h]->s, *(td.Current[i][j].i[h]->m0), *(td.Current[i][j].i[h]->m1), i, j, td.Current[i][j].i[h]->toxicpregnancy);
                  
                  
                  //--------------------------------------------------------------------------------------
                  // SEARCH 1
                  // searching at the 'home' patch
                  int n;
                  numberofmales   = 0;
                  maleprobsum     = 0.0;
                  
                  for(n=0; n<td.Current[i][j].ni; n++ ){
                      //printf("here\n");
                      if( td.Current[i][j].i[n]->s == 0 ){            // if individual is a "male"
                          numberofmales++;
                          //printf("%d\n", td.Current[i][j].i[n]->s);
                          double matechoosingprob = 0.;
                          matechoosingprob = 1.;
                          maleprobsum += ( fdist->p[n] = matechoosingprob );
                          fdist->p[n] = matechoosingprob;   //two identical distributions pdist and fdist
                          pdist->p[n] = matechoosingprob;   //two identical distributions pdist and fdist
                      } else {
                          // if not male
                          fdist->p[n] = 0.0;  //two identical distributions pdist and fdist
                          pdist->p[n] = 0.0;  //two identical distributions pdist and fdist
                      }
                  }
                  
                  if( maleprobsum != 0. ){
                      maleinteger = 0;
                      fdist->n = td.Current[i][j].ni;
                      initdist(fdist,maleprobsum);
                      maleinteger = drand(&R, fdist);
                      //printf("male integer: %d\n", maleinteger);
                      matingmale = td.Current[i][j].i[maleinteger];
                      success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                  }
                  
                  
                  
                  //printf("here\n");
                  //if( success == 1 && ((*(matingmale->m0) == 1 && *(matingmale->m1) == 0) ||
                  //   (*(matingmale->m0) == 0 && *(matingmale->m1) == 1)) ){
                      //matingfemale->toxicpregnancy = 1;
                      //printf("female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                  //    continue;
                  //}
                  
                  
                  double polymatingprob     = 0.;
                  int    secondmaleinteger  = 0;
                  int    polymating         = 0;
                  int    polysuccess        = 0;
                  
                  polymatingprob = U01(&R);
                  
                  // if more than 1 male, the first male is chosen with success and prob is satisfied, then:
                  if( numberofmales > 1 && success == 1 && polymatingprob <= _POLYANDRYPROB) {
                      pdist->p[maleinteger] = 0.0; //initialize the first mates prob to zero so not chosen again
                      pdist->n = td.Current[i][j].ni;
                      
                      initdist(pdist, maleprobsum-1);
                      
                      secondmaleinteger = drand(&R, pdist);
                      if( (secondmaleinteger - maleinteger) == 0){
                          Error("!!warning! picking the same male twice at (%d, %d): %d - %d\n", i, j, maleinteger, secondmaleinteger);
                      }
                      secondmatingmale = td.Current[i][j].i[secondmaleinteger];
                      //printf("2nd mating male (%llu): (%d/%d) \n", secondmatingmale, *(secondmatingmale->m0), *(secondmatingmale->m1));
                      
                      //if( (*(secondmatingmale->m0) == 1 && *(secondmatingmale->m1) == 0) ||
                      //   (*(secondmatingmale->m0) == 0 && *(secondmatingmale->m1) == 1) ){
                          //matingfemale->toxicpregnancy = 1;
                          //printf("female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, secondmatingmale, *(secondmatingmale->m0),*(secondmatingmale->m1));
                      //    continue;
                      //}
                      polysuccess = 1; // if 2nd male is not sterile then second mate choice is success
                  }
                  
                  if( polysuccess == 1 ) {
                      PolyMate(&td,matingfemale,matingmale,secondmatingmale,&td.Next[i][j]);
                      if( matingflagingeneration == 0 ) matingflagingeneration = 1;
                      //printf("S1-PM\n");
                      //if( matingfemale->toxicpregnancy == 1 )
                      //    printf("PMate-female (%llu) toxic pregnancy [%d]; 1st male: (%llu) %d/%d; 2nd male: (%llu) (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1), secondmatingmale, *(secondmatingmale->m0),*(secondmatingmale->m1));
                      polymating++;
                  }
                   
                  //printf("here\n");
                  if( polysuccess == 0 ){
                      if( success != 0 ){
                          //printf("1 mating female (%d, %d) ind's sex: %d, tp: %d \n", i, j, td.Current[i][j].i[h]->s, td.Current[i][j].i[h]->toxicpregnancy);
                          Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                          if( matingflagingeneration == 0 ) matingflagingeneration = 1;
                          //printf("S1-M\n");
                          mating++;
                      }
                  }
                  
                  //--------------------------------------------------------------------------------------
                  // SEARCH 2
                  // if no males in 'home' patch, search among the roamers at the 'home' patch
                  int     numberofmalesroam   = 0;
                  double  maleprobsumroam     = 0.;
                  //printf("S2\n");
                  //printf("here2\n");
                  //printf("maleprobsum: %d, numberofmales: %d\n", maleprobsum, numberofmales);
                  success = 0;
                  //printf("here3\n");
                  if( !maleprobsum || !numberofmales ){
                      int n;
                      nothers = td.Current[i][j].nothers;
                      //printf("(%d, %d) NOTHERS: %d\n", i,j, nothers);
                      //printf("here4\n");
                      if( !nothers) continue;
                      //printf("here5\n");
                      for(n=0; n<nothers; n++ ){
                          //printf("s: %d (%d/%d)\n",td.Current[i][j].Other[n].ind->s, n+1,nothers);
                          if( td.Current[i][j].Other[n].ind->s == 0 ){
                              numberofmalesroam++;
                              double matechoosingprob = 0.;
                              matechoosingprob = 1.;
                              maleprobsumroam += ( fdist->p[n] = matechoosingprob );
                          } else {
                              // if not male
                              fdist->p[n] = 0.0;
                          }
                      }
                      //printf("here6\n");
                      if( maleprobsumroam != 0. && numberofmalesroam !=0 ){
                          maleinteger = 0;
                          fdist->n = nothers;
                          initdist(fdist,maleprobsumroam);
                          maleinteger = drand(&R, fdist);
                          matingmale = td.Current[i][j].Other[maleinteger].ind;
                          success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                      }
                      //printf("here7\n");
                      //printf("S2 mating male (%llu): %d/%d \n", matingmale, *(matingmale->m0), *(matingmale->m1));
                      
                      //if( success == 1 && ((*(matingmale->m0) == 1 && *(matingmale->m1) == 0) ||
                      //   (*(matingmale->m0) == 0 && *(matingmale->m1) == 1)) ){
                          //matingfemale->toxicpregnancy = 1;
                          //printf("female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                      //    continue;
                      //}
                      
                      if( success != 0 ){
                          //printf("2 mating female (%d, %d) ind's sex: %d, tp: %d \n", i, j, td.Current[i][j].i[h]->s, td.Current[i][j].i[h]->toxicpregnancy);
                          Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                          if( matingflagingeneration == 0 ) matingflagingeneration = 1;
                          //printf("S2-M\n");
                          //if( matingfemale->toxicpregnancy == 1 )
                              //printf("S2-female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                          mating++;
                      }
                  }
                  
                  
                  //--------------------------------------------------------------------------------------
                  // SEARCH 3
                  // if no males in 'home' patch (+ no roamers)
                  // but this search 3 is futile, if there were roamers in home patch, she would
                  // have intersected with their central patch while she was roaming
                  int distmatingsuccess;
                  //printf("S3\n");
                  if( !numberofmales && !success ){
                      
                      //int mr; u64b w; double m;
                      int matedist,nPatches,totalPatches,matei,matej,v;
                      distmatingsuccess = 0;
                      // skips matedist=0, should have already found it if it was there.
                      // if no males in 'home' patch, search among the roamers at the 'home' patch
                      
                      //for( matedist=1; matedist<((_RRANGE-1)/2+1); matedist++ ){
                      for( matedist=1; matedist<((_DRANGE-1)/2+1); matedist++ ){
                          
                          if( distmatingsuccess == 1 ) continue;
                          
                          // for a given distance initiate an array for coordinates to randomly pick from
                          totalPatches = (matedist*8);
                          int coordArray[_DRANGE*_DRANGE][2];
                          int tempArray[1][2];
                          int randomPIndex;
                          nPatches = 0;
                          
                          for( matei=-matedist; matei<(matedist+1); matei++ ) {
                              for( matej=-matedist; matej<(matedist+1); matej++ ) {
                                  if( (!matei) && (!matej) )                continue;
                                  if( !ValidMatePatch(&td,i,j,matei,matej)) continue;
                                  coordArray[nPatches][0] = i+matei;
                                  coordArray[nPatches][1] = j+matej;
                                  nPatches++;
                              }
                          }
                          
                          // or randomly rearrange the array and go through systematically
                          for( v = 0; v < (nPatches); v++ ) {
                              tempArray[0][0] = coordArray[v][0];
                              tempArray[0][1] = coordArray[v][1];
                              randomPIndex = rnd(&R, nPatches);
                              coordArray[v][0] = coordArray[randomPIndex][0];
                              coordArray[v][1] = coordArray[randomPIndex][1];
                              coordArray[randomPIndex][0] = tempArray[0][0];
                              coordArray[randomPIndex][1] = tempArray[0][1];
                          }
                          
                          int mg;
                          int rlat = 0, rlon = 0, success = 0, maleintr = 0;
                          // loop through as many times as there are patches until a male is found
                          
                          for(mg=0; mg < (nPatches); mg++ ){
                              // random (lat, long)
                              //printf("searching for N in S3: %d/%d \n ", mg+1,nPatches);
                              rlat = coordArray[mg][0];
                              rlon = coordArray[mg][1];
                              //printf("(%d, %d) N:%d \n", rlat,rlon, td.Current[rlat][rlon].ni);
                              
                              if( !td.Current[rlat][rlon].ni) continue;         // if no ind. here in the random patch skip it
                              matingpatch = &td.Current[rlat][rlon];
                              
                              // there is someone in the patch, loop through individuals
                              int gn;
                              numberofmales = 0;
                              maleprobsum   = 0.0;
                              for(gn=0; gn<matingpatch->ni; gn++ ){
                                  if( matingpatch->i[gn]->s == MALE ){
                                      numberofmales++;
                                      
                                      double matechoosingprob = 0.;
                                      matechoosingprob = 1.;
                                      maleprobsumroam += ( fdist->p[gn] = matechoosingprob );
                                  } else {
                                      // if not male
                                      fdist->p[gn] = 0.0;
                                  }
                              } // closes the gn loop
                              
                              
                              if( maleprobsum != 0. && distmatingsuccess != 1 ){
                                  maleinteger = 0;
                                  fdist->n = matingpatch->ni;
                                  initdist(fdist,maleprobsum);
                                  maleinteger = drand(&R, fdist);
                                  matingmale = matingpatch->i[maleinteger];
                                  success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                              }
                              
                              //printf("S3 mating male (%llu): %d/%d \n", matingmale, *(matingmale->m0), *(matingmale->m1));
                              
                              //if( success == 1 && ((*(matingmale->m0) == 1 && *(matingmale->m1) == 0) ||
                              //   (*(matingmale->m0) == 0 && *(matingmale->m1) == 1)) ){
                                  //matingfemale->toxicpregnancy = 1;
                                  //printf("female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                              //    continue;
                              //}
                              
                              if( success != 0 && distmatingsuccess != 1){
                                  //printf("3 mating female (%d, %d) ind's sex: %d, tp: %d \n", i, j, td.Current[i][j].i[h]->s, td.Current[i][j].i[h]->toxicpregnancy);
                                  Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                                  if( matingflagingeneration == 0 ) matingflagingeneration = 1;
                                  //printf("S3-M\n");
                                  //if( matingfemale->toxicpregnancy == 1 )
                                  //    printf("S3-female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                                  distmatingsuccess = 1;
                                  mating++;
                                  break;
                              }
                          }
                           
                          //printf("finished search 3\n");
                          
                          //--------------------------------------------------------------------------------------
                          // SEARCH 4
                          // check the roaming males as well here
                          int roammate = 0;
                          //printf("S4 starts\n");
                          if( !distmatingsuccess && roammate == 0){
                              for(mg=0; mg < (nPatches); mg++ ){
                                  
                                  //printf("patches: %d/%d, roammate: %d\n", mg+1, nPatches, roammate);
                                  if( roammate == 1 ) break;
                                  rlat = coordArray[mg][0];
                                  rlon = coordArray[mg][1];
                                  //printf("(%d, %d) N:%d Nothers: %d\n", rlat,rlon, td.Current[rlat][rlon].ni, td.Current[rlat][rlon].nothers);
                                  
                                  //if( !td.Current[rlat][rlon].ni) continue;         // if no ind. here in the random patch skip it
                                  matingpatch = &td.Current[rlat][rlon];
                              
                                  int     n, nothers;
                                  int     numberofmalesroam = 0;
                                  double  maleprobsumroam   = 0.;
                                  
                                  nothers = matingpatch->nothers;
                                  if( !nothers ) continue;
                                  for(n=0; n<nothers; n++ ){
                                      if( matingpatch->Other[n].ind->s == MALE ){
                                          numberofmalesroam++;
                                          
                                          double matechoosingprob = 0.;
                                          matechoosingprob = 1.;
                                          maleprobsumroam += ( fdist->p[n] = matechoosingprob );
                                      } else {
                                          // if not male
                                          fdist->p[n] = 0.0;
                                      }
                                  } // closes the n loop
                                  
                                  if( !numberofmalesroam || !maleprobsumroam ) continue; // if no males, go to the next patch
                                  
                                  if( maleprobsumroam != 0. ){ // if there are males then mate...
                                      maleinteger = 0;
                                      fdist->n = nothers;
                                      initdist(fdist,maleprobsumroam);
                                      maleinteger = drand(&R, fdist);
                                      matingmale = matingpatch->Other[maleinteger].ind;
                                      success = 1;  // initialize the success to 1, so that we don't force mating with maleintr = 0;
                                      roammate = 1; // even though not mated, flag before break
                                  }
                                  
                                  // skip t-allele homozygous male if picked, the female loses a mating opportunity
                                  //if( (*(matingmale->m0) == 1 && *(matingmale->m1) == 0) ||
                                  //   (*(matingmale->m0) == 0 && *(matingmale->m1) == 1) ){
                                      //matingfemale->toxicpregnancy = 1;
                                      //printf("female (%llu) toxic pregnancy [%d], male (%llu) is (%d/%d)\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                                  //    break;
                                  //}
                                  
                                  //if( matingfemale->toxicpregnancy == 1 ) break;
                                  //printf("S4-mating male (%llu): (%d/%d) \n", matingmale, *(matingmale->m0), *(matingmale->m1));
                                  
                                  
                                  if( success != 0 ){
                                      //printf("4 mating female (%d, %d) ind's sex: %d, tp: %d \n", i, j, matingfemale->s, matingfemale->toxicpregnancy);
                                      Mate(&td,matingfemale,matingmale,&td.Next[i][j]);
                                      roammate = 1;
                                      if( matingflagingeneration == 0 ) matingflagingeneration = 1;
                                      //printf("S4-M\n");
                                      //if( matingfemale->toxicpregnancy == 1 )
                                          //printf("S4-female (%llu) toxic pregnancy [%d], male (%llu) is %d/%d\n", matingfemale, matingfemale->toxicpregnancy, matingmale, *(matingmale->m0),*(matingmale->m1));
                                      distmatingsuccess = 1;
                                      mating++;
                                      break;
                                  }
                                  
                                  if( roammate == 1 ) break;
                              } // closes search 4
                          } // closes mg loop
                      } // closes mate dist loop
                      
                  } // closes if no males in 'home' patch (+roamers)
              } // closes the female search loop (going through all N in patch)
              if( td.Current[i][j].nf + td.Current[i][j].nm != td.Current[i][j].ni)
                  Error("! %d + %d != %d [%d, %d]\n", td.Current[i][j].nf, td.Current[i][j].nm, td.Current[i][j].ni, i, j);
          } // closes j loop
      } // closes j loop
       
      
      //**********************
#if _SCENTSCAPE
      
      // write number of females (.nf) and males (.nm) in each patch (add in niche.h and write while going though females above?)
      // go through landscape
      int ii, jj, distance;
      
      for(i=0; i<WIDTH; i++){
          for(j=0; j<HEIGHT; j++){
              //printf("(%d, %d) %d\n", i, j, i*HEIGHT+j);
              //if( td.Current[i][j].ni > 0 ){
              //    Scentdata[1][i*HEIGHT+j] = ((double)(td.Current[i][j].nf)); // add female data
              //    Scentdata[0][i*HEIGHT+j] = ((double)(td.Current[i][j].nm)); // add male data
                  //printf("here\n");
              //} else {
                  double sumf = 0.;
                  double summ = 0.;
                  int count = 0;
                  for(ii=0; ii<WIDTH; ii++){
                      for(jj=0; jj<HEIGHT; jj++){
                          if( i == ii && j == jj ) continue;
                          if( abs(i-ii) >= abs(j-jj) ){
                              distance = abs(i-ii);
                          } else {
                              distance = abs(j-jj);
                          }
                          //printf("(%d, %d) -> (%d, %d), distance: %d, F: %d, M: %d\n", i,j, ii, jj, distance, td.Current[i][j].nf, td.Current[i][j].nm );
                          sumf += ((double)(td.Current[ii][jj].nf)) / ((double)(distance));
                          summ += ((double)(td.Current[ii][jj].nm)) / ((double)(distance));
                          //count++;
                          //printf("sf: %lf, sm: %lf, distance: %d\n", sumf, summ, distance);
                      }
                  }
                  Scentdata[1][i*HEIGHT+j] = sumf;
                  Scentdata[0][i*HEIGHT+j] = summ;
                  //printf("here\n");
              //}
          }
      }
      
      LogScentscape();
#endif
    //**********************
                         
      for(i=0; i<WIDTH; i++){
          for(j=0; j<HEIGHT; j++){
              for(h=0, l=td.Current[i][j].ni; h<l; h++){
                  //printf("checking survival (%d, %d) %d of %d age: %d\n", i,j,h+1,l, td.Current[i][j].i[h]->age);
                  if( !l ) continue;                                          // skip empty patches
                  if( td.Current[i][j].i[h]->age > _MAXAGE ) continue;        // Skip any inds older than maximum age
                  //if( td.Current[i][j].i[h]->s == 1 && td.Current[i][j].i[h]->toxicpregnancy == 1 ) continue; // Skip any toxic females
          
                  double randomno = 0.;
                  double survivalprobability = 0.;
                  randomno = U01(&R);
                  int g, cullflag = 0;

                  for( g = 0; g < _NoCull; g++ ){
                      if( Gen == CullGen[g] ) {
                          cullflag = 1;
                      }
                  }

                  if( cullflag == 1 ){
                      survivalprobability = _SURVIVALPROBABILITY * _CullSurvival;
                      //printf("sur: %lf\n", survivalprobability);
                  } else {
                      survivalprobability = _SURVIVALPROBABILITY;
                  }
                  //printf("sur: %lf\n", survivalprobability);
                  if( randomno < survivalprobability ) {
#if _SCENTSCAPE
                      if( matingflagingeneration == 1 ){
#endif
                          Survive(&td,&td.Next[i][j],td.Current[i][j].i[h]);
#if _SCENTSCAPE
                      }
#endif
#if _SCENTSCAPE
                      if( matingflagingeneration == 0 ){
                          ExtensiveSearchSurvive(&td,&td.Next[i][j],td.Current[i][j].i[h]);
                      }
#endif
                      Surviving++;
                  } else {
                      Mortality++;
                  }
                  
              }
          } // closes j loop
      } //closes i loop

#if _DISTANCEDENSITYDISPERSAL
     // after swapping, next is current, so I can check for survival here
     for(i=0; i<WIDTH; i++){
         for(j=0; j<HEIGHT; j++){
             // critical for survival calculations!
             //td.Current[i][j].density = ( ((double)(NeDispersal(&td,&td.Current[i][j]))) / ((double)(td.Current[i][j].ddne * _K0)) );
             //td.Next[i][j].density = ( ((double)(NeDispersal(&td,&td.Next[i][j]))) / ((double)(td.Next[i][j].ddne * _K0)) );
             td.Current[i][j].density = ((double)(td.Current[i][j].ni)) / ((double)(_INIT_INDIVIDUALS));
         }
     }
#endif

      
#if _INOC
      int inoculateindividuals = 0;
      for(i=0; i<WIDTH; i++){
          for(j=0; j<HEIGHT; j++){
              for( g = 0; g < _NoInoc; g++ ) {
                  if( Gen == InocGen[g] ) {
                      if( td.Next[i][j].monitorgen == 1 ){
                          
                          //inoculateindividuals  = InocIndividuals[g];
                          //Inoculants += inoculateindividuals;
                          //Inoculate(&td,&td.Next[i][j], inoculateindividuals);
                          Inoculants += _InocIndividuals;
                          Inoculate(&td,&td.Next[i][j], _InocIndividuals);
                          //printf("inoc for next gen: %d, at (%d, %d), Ni: %d\n", Gen, i,j, _InocIndividuals);
                      }
                  }
              }
          }
      }
#endif
      
/*
#if _INOC
      for(i=0; i<WIDTH; i++){
          for(j=0; j<HEIGHT; j++){

          int inoculateflag        = 0;
          int inoculateindividuals = 0;
          
          for( g = 0; g < _NoInoc; g++ ) {
              if( Gen == InocGen[g] ) {
                  //printf("inoc gen: %d\n", Gen);
                  if( td.Current[i][j].inocpatch == 1 ) {
                      //N's should be calculated in 'next' since this is after all mating survival is done and written to next
                      //printf("%d- (%d %d), N_next: %d\n",b,i,j,td.Next[i][j].ni);
#if _MONITOR
                      if( !td.Next[i][j].ni ) continue;
#if _TOXICMALEDETECTION
                      int toxicmaleno = 0;
                      for(h=0, l=td.Next[i][j].ni; h<l; h++){
                          if( *(td.Next[i][j].i[h]->m0) == 1 || *(td.Next[i][j].i[h]->m1) == 1 ) toxicmaleno++; // count toxic males
                      }
                      //printf("tox.m no: %d - N: %d (diff: %d)\n", toxicmaleno, td.Next[i][j].ni, td.Next[i][j].ni - toxicmaleno);
                      if( td.Next[i][j].ni - toxicmaleno >= 1 ){
                          inoculateflag         = 1;
                          inoculateindividuals  = InocIndividuals[g];
                          //printf("(tox. male detection) gen: %d, at (%d, %d), Ni: %d\n", Gen, i,j, inoculateindividuals);
                          Inoculants += inoculateindividuals;
                      }
                      
#else
                      // should not be here uncless there are individuals
                      inoculateflag         = 1;
                      inoculateindividuals  = InocIndividuals[g];
                      //printf("(no tox. m. detection) gen: %d, at (%d, %d), Ni: %d\n", Gen, i,j, inoculateindividuals);
                      Inoculants += inoculateindividuals;
#endif

#else
                      inoculateflag         = 1;
                      inoculateindividuals  = InocIndividuals[g];
                      Inoculants += inoculateindividuals;
                      //printf("(no monitoring) gen: %d, at (%d, %d), Ni: %d\n", Gen, i,j, inoculateindividuals);
#endif
                  }
              }
          }
          if(inoculateflag == 1){
              Inoculate(&td,&td.Next[i][j], inoculateindividuals);
          }
          inoculateflag = 0;
      }
  }
#endif
 */
      
  for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
          td.Current[i][j].ni = 0;
      }
  }

  // Free and swap to generations to make all next gen current.
  ReleaseAllCurrentIndividuals(&td);
  SwapGenerations(&td);
#if _DISTANCE
  LogDistance();
#endif

    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        td.Current[i][j].nothers=0;
        td.Current[i][j].monitorgen = 0;
      }
    }
    
    // calculates who is using which patch for how long
    for(i=0; i<WIDTH; i++){
      for(j=0; j<HEIGHT; j++){
        for(k=0,l=td.Current[i][j].ni; k < l; k++) {
            TimeSpentRoam(td.Current[i][j].i[k], &td.Current[i][j]);
        }
      }
    }

  } // while signal done
  return NULL;
} // void main
