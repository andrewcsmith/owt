#include "owt.h"

OWTResults owt_optimize_temperament(Criteria* criteria) {
  OWTResults results;

  // matrix information
  gsl_matrix *X, *cov;
  gsl_vector *w, *y, *c;

  // Data to be minimized (A, in Python)
  X = gsl_matrix_alloc(criteria->num_pitches * (criteria->num_pitches-1), criteria->num_pitches-1);
  // observations (B, in Python)
  y = gsl_vector_alloc(criteria->num_pitches * (criteria->num_pitches-1));

  // TODO: use the gsl vector/matrix operations for all of this
  // Assign the portions of the X matrix and y vector
  int row = 0;
  for (int i = 1; i < criteria->num_pitches; i++) {
    for (int j = 0; j < criteria->num_pitches; j++) {
      // I'm not sure what this means
      int ed = i + j;

      // The y vector should begin as contiguous blocks of each ideal interval
      gsl_vector_set(y, row, criteria->ideal_intervals[i-1]);

      // Set all but the first column to -1.0 initially?
      if (j > 0) {
        gsl_matrix_set(X, row, j-1, -1.0);
      }

      if (ed < criteria->num_pitches) {
        // set the element of the current row to 1.0
        gsl_matrix_set(X, row, ed-1, 1.0);
      } else if (ed == criteria->num_pitches) {
        // the "ideal" is the inversion. decrease it by an octave.
        gsl_vector_set(y, row, gsl_vector_get(y, row) - criteria->repeat_factor);
      } else {
        ed = ed % criteria->num_pitches;
        // next two lines are same as above (should refactor...)
        gsl_matrix_set(X, row, ed-1, 1.0);
        gsl_vector_set(y, row, gsl_vector_get(y, row) - criteria->repeat_factor);
      }
      row++;
    }
  }

  // weights (W) - a vector of each interval weight multiplied by the key weight
  w = gsl_vector_alloc(criteria->num_pitches * (criteria->num_pitches-1));
  // set the weights vector, TODO: use the gsl vector library
  for (int i = 0; i < criteria->num_pitches-1; i++) {
    for (int j = 0; j < criteria->num_pitches; j++) {
      gsl_vector_set(w, i*criteria->num_pitches + j, 
          criteria->key_weights[j] * criteria->interval_weights[i]);
    }
  }

  // best-fit results
  results.optimal_tuning = gsl_vector_alloc(criteria->num_pitches-1);
  // covariance results, or (X^T * W * X)^{-1}
  cov = gsl_matrix_alloc(criteria->num_pitches-1, criteria->num_pitches-1);

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

gsl_vector* owt_interval_error(int num_pitches, int interval_stride, double interval_ideal,
    double octave, gsl_vector* tuning) 
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

