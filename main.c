#include <stdio.h>
#include <gsl/gsl_multifit.h>

gsl_vector* optimize_temperament(int num_pitches, double* ideal_intervals, double octave,
                                 double* interval_weights, double* key_weights) 
{
  // matrix information
  gsl_matrix *X, *cov;
  gsl_vector *w, *y, *c;
  double chisq;

  // Data to be minimized (A)
  X = gsl_matrix_alloc(num_pitches * (num_pitches-1), num_pitches-1);
  // observations (B)
  y = gsl_vector_alloc(num_pitches * (num_pitches-1));

  // Assign the portions of the X matrix and y vector
  int row = 0;
  for (int i = 1; i < num_pitches; i++) {
    for (int j = 0; j < num_pitches; j++) {
      int ed = i + j;
      int st = ed - i;

      gsl_vector_set(y, row, ideal_intervals[i-1]);
      if (st > 0) {
        gsl_matrix_set(X, row, st-1, -1.0);
      }
      if (ed < num_pitches) {
        gsl_matrix_set(X, row, ed-1, 1.0);
      } else if (ed == num_pitches) {
        gsl_vector_set(y, row, gsl_vector_get(y, row) - octave);
      } else {
        ed = ed % num_pitches;
        gsl_matrix_set(X, row, ed-1, 1.0);
        gsl_vector_set(y, row, gsl_vector_get(y, row) - octave);
      }
      row++;
    }
  }

  // weights (W) - a vector of each interval weight multiplied by the key weight
  w = gsl_vector_alloc(num_pitches * (num_pitches-1));
  // set the weights vector
  for (int i = 0; i < num_pitches-1; i++) {
    for (int j = 0; j < num_pitches; j++) {
      gsl_vector_set(w, i*num_pitches + j, key_weights[j] * interval_weights[i]);
    }
  }

  // best-fit results
  c = gsl_vector_alloc(num_pitches-1);
  // covariance results, or (X^T * W * X)^{-1}
  cov = gsl_matrix_alloc(num_pitches-1,num_pitches-1);

  gsl_multifit_linear_workspace *work = 
    gsl_multifit_linear_alloc(num_pitches * (num_pitches-1), num_pitches-1);
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  gsl_vector_free(w);
  gsl_vector_free(y);

  return c;
}

int main(void) {
  int num_pitches = 12;
  double ideal_intervals[11] = {100, 200, 300, 386, 498, 600, 702, 814, 900, 1000, 1100};
  double octave = 1200.0;

  double interval_weights[11] = {1, 1, 1, 150, 1600, 1, 1600, 150, 1, 1, 1};
  double key_weights[12] = {50, 60, 7, 1, 7, 60, 50, 60, 7, 1, 7, 80};

  gsl_vector* c;

  c = optimize_temperament(num_pitches, ideal_intervals, octave, interval_weights, key_weights);
  
  printf("Best fit tuning: \n");
  for (int i = 0; i < num_pitches-1; i++) {
    printf("\t%.2f\n", gsl_vector_get(c, i));
  }

  return 0;
}

