#include <gsl/gsl_multifit.h>

gsl_vector* optimize_temperament(int num_pitches, double* ideal_intervals, double octave,
                                 double* interval_weights, double* key_weights); 
