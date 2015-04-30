#ifndef OWT_H
#define OWT_H

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <math.h>

typedef struct OWTCriteria {
  char* title;
  int num_pitches;
  double* ideal_intervals;
  double* interval_weights; 
  double* key_weights; 
  double repeat_factor;
} OWTCriteria;

typedef struct OWTResults {
  double chisq;
  gsl_vector* optimal_tuning;
} OWTResults;

typedef struct OWTSearchParams {
  OWTCriteria* criteria;
  gsl_vector* tuning;
} OWTSearchParams;

OWTResults owt_optimize_temperament(OWTCriteria* criteria); 
void owt_populate_source_matrix(gsl_matrix* X, OWTCriteria* criteria);
void owt_populate_ideal_interval_vector(gsl_vector* y, OWTCriteria* criteria);
void owt_populate_weights_vector(gsl_vector* w, OWTCriteria* criteria);
double owt_deviation_from_tuning(const gsl_vector *x, void *params);
void owt_criteria_memcpy(OWTCriteria* dest, OWTCriteria* src);

gsl_vector* owt_interval_error(int num_pitches, int interval_stride, double interval_ideal,
    double octave, gsl_vector* tuning);

#endif
