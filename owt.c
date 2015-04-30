#include <gsl/gsl_blas.h>
#include "owt.h"

OWTResults owt_optimize_temperament(OWTCriteria* criteria) {
  OWTResults results;

  // matrix information
  gsl_matrix *X, *cov;
  gsl_vector *w, *y, *c;

  // Data to be minimized (A, in Python)
  X = gsl_matrix_alloc(criteria->num_pitches * (criteria->num_pitches-1), criteria->num_pitches-1);
  // observations (B, in Python)
  y = gsl_vector_alloc(criteria->num_pitches * (criteria->num_pitches-1));
  // weights (W) - a vector of each interval weight multiplied by the key weight
  w = gsl_vector_alloc(criteria->num_pitches * (criteria->num_pitches-1));
  // best-fit results
  results.optimal_tuning = gsl_vector_alloc(criteria->num_pitches-1);
  // covariance results, or (X^T * W * X)^{-1}
  cov = gsl_matrix_alloc(criteria->num_pitches-1, criteria->num_pitches-1);

  owt_populate_source_matrix(X, criteria);
  owt_populate_ideal_interval_vector(y, criteria);
  owt_populate_weights_vector(w, criteria);

  gsl_multifit_linear_workspace *work = 
    gsl_multifit_linear_alloc(criteria->num_pitches * (criteria->num_pitches-1), 
        criteria->num_pitches-1);
  gsl_multifit_wlinear(X, w, y, results.optimal_tuning, cov, &results.chisq, work);
  gsl_multifit_linear_free(work);

  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(w);
  gsl_vector_free(y);

  return results;
}

void owt_populate_source_matrix(gsl_matrix* X, OWTCriteria* criteria) {
  int row = 0;
  for (int i = 0; i < (criteria->num_pitches - 1); i++) {
    for (int j = 0; j < criteria->num_pitches; j++) {
      int ed = i + j;
      if (j > 0) {
        gsl_matrix_set(X, row, j - 1, -1.0);
      }
      if (ed < (criteria->num_pitches - 1)) {
        gsl_matrix_set(X, row, ed, 1.0);
      } else if (ed > (criteria->num_pitches - 1)) {
        ed = ed % criteria->num_pitches;
        gsl_matrix_set(X, row, ed, 1.0);
      }
      row++;
    }
  }
}

void owt_populate_ideal_interval_vector(gsl_vector* y, OWTCriteria* criteria) {
  int row = 0;
  for (int i = 0; i < (criteria->num_pitches - 1); i++) {
    for (int j = 0; j < criteria->num_pitches; j++) {
      int ed = i + j;
      gsl_vector_set(y, row, criteria->ideal_intervals[i]);
      if (ed >= (criteria->num_pitches - 1)) {
        gsl_vector_set(y, row, gsl_vector_get(y, row) - criteria->repeat_factor);
      }
      row++;
    }
  }
}

double owt_deviation_from_tuning(const gsl_vector *x, void *params) {
  OWTSearchParams* search_params = (OWTSearchParams *)params;
  // copy from the vector x into the key_weights
  for (int i = 0; i < (search_params->criteria->num_pitches); i++) {
    search_params->criteria->key_weights[i] = gsl_vector_get(x, i);
  }
  // get the optimal tuning for those criteria
  OWTResults results = owt_optimize_temperament(search_params->criteria);
  gsl_vector* error = gsl_vector_alloc(search_params->criteria->num_pitches - 1);
  // calculate signed error
  gsl_vector_memcpy(error, results.optimal_tuning);
  gsl_vector_sub(error, search_params->tuning);
  // take absolute value and return it
  double accum = 0.0;
  for (int i = 0; i < search_params->criteria->num_pitches - 1; i++) {
    accum += fabs(gsl_vector_get(error, i));
  }
  gsl_vector_free(error);
  gsl_vector_free(results.optimal_tuning);
  return accum;
}

void owt_criteria_memcpy(OWTCriteria* dest, OWTCriteria* src) {
  strcpy(dest->title, src->title);
  dest->num_pitches = src->num_pitches;

  dest->ideal_intervals = malloc(sizeof(double) * (src->num_pitches - 1));
  memcpy(dest->ideal_intervals, src->ideal_intervals, sizeof(double) * (src->num_pitches - 1));

  dest->interval_weights = malloc(sizeof(double) * (src->num_pitches - 1));
  memcpy(dest->interval_weights, src->interval_weights, sizeof(double) * (src->num_pitches - 1));

  dest->key_weights = malloc(sizeof(double) * src->num_pitches);
  memcpy(dest->key_weights, src->interval_weights, sizeof(double) * src->num_pitches);

  dest->repeat_factor = src->repeat_factor;
}

/* void owt_populate_ideal_interval_combination_vectors(gsl_vector* y, OWTCriteria* criteria_one, OWTCriteria* criteria_two) { */
/*   gsl_vector* ideal_intervals = gsl_vector_calloc(criteria_one->num_pitches * (criteria_one->num_pitches - 1)); */
/*   gsl_vector_memcpy(ideal_intervals, criteria_one->ideal_intervals); */
/*  */
/*   for (unsigned index = 0; index < (criteria_one->num_pitches - 1); index++) { */
/*     // If the two options are identical, don't do anything */
/*     if (criteria_one->ideal_intervals[index] == criteria_two->ideal_intervals[index]) { */
/*       continue; */
/*     } */
/*     // Otherwise, try every combination of the two */
/*     for (unsigned combo = 0; combo < pow(2, criteria_one->num_pitches), combo++) { */
/*     } */
/*   } */
/*  */
/*   int row = 0; */
/*   for (int i = 0; i < criteria->num_pitches - 1; i++) { */
/*     for (int j = 0; j < criteria->num_pitches; j++) { */
/*       int ed = i + j; */
/*       // Set the ideal interval of the row */
/*       if (ed >= (criteria->num_pitches - 1)) { */
/*         gsl_vector_set(y, row, gsl_vector_get(y, row) - criteria->repeat_factor); */
/*       } */
/*       row++; */
/*     } */
/*   } */
/* } */

void normalize_vector(gsl_vector* vec, double min, double max) {
  double old_min = gsl_vector_min(vec);
  double old_max = gsl_vector_max(vec);
  double range = old_max - old_min;
  gsl_vector_add_constant(vec, min - old_min);
  gsl_vector_scale(vec, 1.0 / range);
}

void owt_populate_weights_vector(gsl_vector* w, OWTCriteria* criteria) {
  for (int i = 0; i < criteria->num_pitches-1; i++) {
    for (int j = 0; j < criteria->num_pitches; j++) {
      gsl_vector_set(w, i * criteria->num_pitches + j,
          criteria->key_weights[j] * criteria->interval_weights[i]);
    }
  }
  double scale = 1.0 / gsl_blas_dasum(w);
  gsl_vector_scale(w, scale);
}

gsl_vector* owt_interval_error(int num_pitches, int interval_stride, 
    double interval_ideal, double octave, gsl_vector* tuning) 
{
  gsl_vector *errors, *intervals, *old_tuning;
  errors = gsl_vector_alloc(num_pitches);
  intervals = gsl_vector_alloc(num_pitches);
  // The octave below
  int inversion = interval_stride - num_pitches;

  old_tuning = tuning;
  tuning = gsl_vector_calloc(num_pitches + 1);
  for (int i = 0; i < num_pitches-1; i++) 
    gsl_vector_set(tuning, i+1, gsl_vector_get(old_tuning, i));

  for (int i = 0; i < num_pitches; i++) {
    double interval, tonic;
    tonic = gsl_vector_get(tuning, i);
    if (i + interval_stride >= num_pitches) {
      interval = gsl_vector_get(tuning, i + inversion) + octave;
    } else {
      interval = gsl_vector_get(tuning, i + interval_stride);
    }
    // Set the intervals vector to the various thirds
    gsl_vector_set(intervals, i, interval - tonic);
  }
  // Set the errors vector to the ideal interval
  gsl_vector_set_all(errors, interval_ideal);

  // Subtract the actual intervals from the expected intervals
  gsl_vector_sub(errors, intervals);

  for (int i = 0; i < num_pitches; i++) {
    double error = gsl_vector_get(errors, i);
    gsl_vector_set(errors, i, fabs(error));
  }

  gsl_vector_free(intervals);
  return errors;
}

// next: use multidimensional minimization to minimize the vector of weights
// from a given tuning system.

