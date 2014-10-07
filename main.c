#include <stdio.h>
#include "owt.h"

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

  gsl_vector_free(c);

  return 0;
}

