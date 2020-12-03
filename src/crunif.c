
#include <R.h>
#include <Rmath.h>

/* call before generating any random variates */
void F77_SUB(rndstart)(void) { GetRNGstate(); }

/* call after done generating all random variates */
void F77_SUB(rndend)(void) { PutRNGstate(); }

/* call to generate one uniform random variate */
double F77_SUB(frunif)(void)
{
  return runif(0.0, 1.0);
}
