#include "owt.h"

gsl_vector* owt_optimize_temperament(int num_pitches, double* ideal_intervals,
    double octave, double* interval_weights, double* key_weights, double* chisq) 
{
  // matrix information
  gsl_matrix *X, *cov;
  gsl_vector *w, *y, *c;

  // Data to be minimized (A, in Python)
  X = gsl_matrix_alloc(num_pitches * (num_pitches-1), num_pitches-1);
  // observations (B, in Python)
  y = gsl_vector_alloc(num_pitches * (num_pitches-1));

  // TODO: use the gsl vector/matrix operations for all of this
  // Assign the portions of the X matrix and y vector
  int row = 0;
  for (int i = 1; i < num_pitches; i++) {
    for (int j = 0; j < num_pitches; j++) {
      // I'm not sure what this means
      int ed = i + j;

      // The y vector should begin as contiguous blocks of each ideal interval
      gsl_vector_set(y, row, ideal_intervals[i-1]);

      // Set all but the first column to -1.0 initially?
      if (j > 0) {
        gsl_matrix_set(X, row, j-1, -1.0);
      }

      if (ed < num_pitches) {
        // set the element of the current row to 1.0
        gsl_matrix_set(X, row, ed-1, 1.0);
      } else if (ed == num_pitches) {
        // the "ideal" is the inversion. decrease it by an octave.
        gsl_vector_set(y, row, gsl_vector_get(y, row) - octave);
      } else {
        ed = ed % num_pitches;
        // next two lines are same as above (should refactor...)
        gsl_matrix_set(X, row, ed-1, 1.0);
        gsl_vector_set(y, row, gsl_vector_get(y, row) - octave);
      }
      row++;
    }
  }

  // weights (W) - a vector of each interval weight multiplied by the key weight
  w = gsl_vector_alloc(num_pitches * (num_pitches-1));
  // set the weights vector, TODO: use the gsl vector library
  for (int i = 0; i < num_pitches-1; i++) {
    for (int j = 0; j < num_pitches; j++) {
      gsl_vector_set(w, i*num_pitches + j, key_weights[j] * interval_weights[i]);
    }
  }

  // Weights vector is made up of contiguous blocks of each interval in each
  // key. For a repeat factor containing n notes, the first n elements are the
  // weight of interval 1 multiplied by the weight of each successive key
  // weight. This cycles through the entire matrix. 
  //
  // In the following representation, the rows are each interval, and the
  // columns are each key.
  /*
  printf("== weights vector:\n");
  for (int i = 0; i < num_pitches-1; i++) {
    for (int j = 0; j < num_pitches; j++) {
      printf("%0.3f\t", gsl_vector_get(w, i*num_pitches + j));
    }
    printf("\n");
  }
  */

  // best-fit results
  c = gsl_vector_alloc(num_pitches-1);
  // covariance results, or (X^T * W * X)^{-1}
  cov = gsl_matrix_alloc(num_pitches-1,num_pitches-1);

  gsl_multifit_linear_workspace *work = 
    gsl_multifit_linear_alloc(num_pitches * (num_pitches-1), num_pitches-1);
  gsl_multifit_wlinear(X, w, y, c, cov, chisq, work);
  gsl_multifit_linear_free(work);

  /*
  printf("\n== covariance matrix:\n");
  for (int i = 0; i < num_pitches-1; i++) {
    for (int j = 0; j < num_pitches-1; j++) {
      printf("%0.5f\t", gsl_matrix_get(cov, i, j));
    }
    printf("\n");
  }
  */

  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(w);
  gsl_vector_free(y);

  return c;
}

double owt_interval_std(int num_pitches, int interval_stride, double interval_ideal, 
    double octave, gsl_vector* tuning)
{
  gsl_vector *errors, *intervals;
  errors = gsl_vector_alloc(num_pitches);
  intervals = gsl_vector_alloc(num_pitches);
  // The octave below
  int inversion = interval_stride - num_pitches;

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
  printf("\nerrors:\n");
  for (int i = 0; i < num_pitches; i++)
    printf("%0.2f\t", gsl_vector_get(errors, i));
  return gsl_stats_sd(errors, 1, num_pitches);
}

double owt_interval_error_max(int num_pitches, int interval_stride, double interval_ideal,
		double octave, gsl_vector* tuning)
{
  gsl_vector *errors, *intervals;
  errors = gsl_vector_alloc(num_pitches);
  intervals = gsl_vector_alloc(num_pitches);
  // The octave below
  int inversion = interval_stride - num_pitches;

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
  return abs(gsl_vector_get(errors, gsl_blas_idamax(errors)));
}
// next: use multidimensional minimization to minimize the vector of weights
// from a given tuning system.

