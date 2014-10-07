#include <stdio.h>
#include <gsl/gsl_multifit.h>

int main(void) {
  // matrix information
  gsl_matrix *X, *cov;
  gsl_vector *w, *y, *c;
  double chisq;

  int num_pitches = 12;
  double idealIntervals[11] = {100, 200, 300, 386, 498, 600, 702, 814, 900, 1000, 1100};
  double octave = 1200.0;

  double intervalWeights[11] = {1, 1, 1, 150, 1600, 1, 1600, 150, 1, 1, 1};
  double keyWeights[12] = {50, 60, 7, 1, 7, 60, 50, 60, 7, 1, 7, 80};
  
  // Data to be minimized (A)
  X = gsl_matrix_alloc(num_pitches * (num_pitches-1), num_pitches-1);
  // observations (B)
  y = gsl_vector_alloc(num_pitches * (num_pitches-1));

  int row = 1;
  for (int i = 1; i < num_pitches; i++) {
    for (int j = 0; j < num_pitches; j++) {
      int ed = i + j;
      int st = ed - i;

      gsl_vector_set(y, row-1, idealIntervals[i-1]);
      if (st > 0) {
        gsl_matrix_set(X, row-1, st-1, -1.0);
      }
      if (ed < num_pitches) {
        gsl_matrix_set(X, row-1, ed-1, 1.0);
      } else if (ed == num_pitches) {
        gsl_vector_set(y, row-1, gsl_vector_get(y, row-1) - octave);
      } else {
        ed = ed % num_pitches;
        gsl_matrix_set(X, row-1, ed-1, 1.0);
        gsl_vector_set(y, row-1, gsl_vector_get(y, row-1) - octave);
      }
      row++;
    }
  }

  // weights (W) - a vector of each interval weight multiplied by the key weight
  w = gsl_vector_alloc(num_pitches * (num_pitches-1));
  // set the weights vector
  for (int i = 0; i < num_pitches-1; i++) {
    for (int j = 0; j < num_pitches; j++) {
      gsl_vector_set(w, i*num_pitches + j, keyWeights[j] * intervalWeights[i]);
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

  printf("Best fit tuning: \n");
  for (int i = 0; i < num_pitches-1; i++) {
    printf("\t%.2f\n", gsl_vector_get(c, i));
  }

  return 0;
}

