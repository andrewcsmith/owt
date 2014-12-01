#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <math.h>

gsl_vector* owt_optimize_temperament(int num_pitches, double* ideal_intervals,
    double octave, double* interval_weights, double* key_weights, double*
    chisq); 
double owt_interval_std(int num_pitches, int interval_stride, double interval_ideal, 
    double octave, gsl_vector* tuning);
double owt_interval_error_max(int num_pitches, int interval_stride, double interval_ideal,
		double octave, gsl_vector* tuning);
