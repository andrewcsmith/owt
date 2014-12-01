#include <stdio.h>
#include "owt.h"

int main(void) {
  int num_pitches = 12;
  double octave = 1200.0;

  // Intervals specified for septimal OWT2
  double ideal_intervals[11] = 
    {100, 204, 267, 386, 498, 600, 702, 800, 900, 969, 1100}; 
  double interval_weights[11] = 
    {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 1, 0.001, 0.001, 0.001, 0.001}; 
  double key_weights[12] = 
    {1, 0.001, 0.001, 0.001, 0.001, 1, 0.001, 1, 0.001, 0.001, 0.001, 0.001};

  double chisq;

  gsl_vector* c;

  c = owt_optimize_temperament(num_pitches, ideal_intervals, octave,
      interval_weights, key_weights, &chisq);
  
  printf("\nChi-squared results: %f\n", chisq);
  printf("Best fit tuning: \n");
  for (int i = 0; i < num_pitches-1; i++) {
    printf("%.2f\t", gsl_vector_get(c, i));
  }
  printf("\nStandard deviation: %0.2f\n", owt_interval_std(num_pitches-1, 4, 386.0, octave, c));
  printf("\nWorst fifth: %0.2f\n", owt_interval_error_max(num_pitches-1, 7, 702.0, octave, c));
  printf("\n");

  gsl_vector_free(c);

  return 0;
}

