#ifndef OWT_H
#define OWT_H

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <math.h>

typedef struct Criteria {
  char* title;
  int num_pitches;
  double* ideal_intervals;
  double* interval_weights; 
  double* key_weights; 
  double repeat_factor;
} Criteria;

typedef struct OWTResults {
  double chisq;
  gsl_vector* optimal_tuning;
} OWTResults;

OWTResults owt_optimize_temperament(Criteria* criteria); 

gsl_vector* owt_interval_error(int num_pitches, int interval_stride, double interval_ideal,
    double octave, gsl_vector* tuning);

#endif
