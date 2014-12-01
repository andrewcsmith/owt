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

  double chisq, third_error_sd, fifth_error_max;

  gsl_vector* c, *third_error, *fifth_error;


  c = owt_optimize_temperament(num_pitches, ideal_intervals, octave,
      interval_weights, key_weights, &chisq);

  third_error = owt_interval_error(num_pitches, 4, 386.0, octave, c);
  double third_error_data[num_pitches];

  printf("\nThirds error:\n");
  for (int i = 0; i < num_pitches; i++) {
    third_error_data[i] = gsl_vector_get(third_error, i);
    printf("%0.2f\t", third_error_data[i]);

  }

  third_error_sd = gsl_stats_sd(third_error_data, 1, num_pitches);

  fifth_error = owt_interval_error(num_pitches, 7, 702.0, octave, c);
  fifth_error_max = fabs(gsl_vector_get(fifth_error, gsl_blas_idamax(fifth_error)));
  printf("\nFifths error:\n");
  for (int i = 0; i < num_pitches; i++)
    printf("%0.2f\t", gsl_vector_get(fifth_error, i));
  
  printf("\nChi-squared results: %f\n", chisq);
  printf("Best fit tuning: \n");
  for (int i = 0; i < num_pitches-1; i++) {
    printf("%.2f\t", gsl_vector_get(c, i));
  }
  printf("\nStandard deviation: %0.2f\n", third_error_sd);
  printf("\nWorst fifth: %0.2f\n", fifth_error_max);
  printf("\n");

  gsl_vector_free(c);
  gsl_vector_free(third_error);
  gsl_vector_free(fifth_error);

  return 0;
}

