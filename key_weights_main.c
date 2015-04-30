#include <stdio.h>
#include <string.h>
#include <mach/mach_time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include "owt.h"
#include "owt_parser.h"

int main(int argc, char** argv) {
  char* path_to_file = malloc(1024);
  char* path_to_output_file = malloc(1024);

  /* Read in command-line arguments: 
   * location_of_file */
  for (int arg = 1; arg < argc; arg++) {
    if (!strcmp(argv[arg], "-i") && argv[arg+1]) 
      strcpy(path_to_file, argv[arg+1]);
    if (!strcmp(argv[arg], "-o") && argv[arg+1]) 
      strcpy(path_to_output_file, argv[arg+1]);
  }

  OWTCriteria criteria = owt_parse(path_to_file);

  if (criteria.num_pitches <= 0) {
    printf("Error occured with owt_parse. Whoops.\n");
    exit(-1);
  }

  printf("title: %s\n", criteria.title);
  printf("num_pitches: %i\n", criteria.num_pitches);
  printf("Ideal Intervals: ");
  for (int i = 0; i < (criteria.num_pitches - 1); i++)
    printf("%1.1f ", criteria.ideal_intervals[i]);
  printf("\ninterval_weights: ");
  for (int i = 0; i < (criteria.num_pitches - 1); i++)
    printf("%1.3f ", criteria.interval_weights[i]);
  printf("\nkey_weights: ");
  for (int i = 0; i < (criteria.num_pitches); i++)
    printf("%1.3f ", criteria.key_weights[i]);
  printf("\n");

  double werck_iii[11] = {90.226, 192.180, 294.135, 390.225, 498.045, 588.045, 696.090, 792.181, 888.270, 996.090, 1092.180};
  gsl_vector* target_tuning = gsl_vector_alloc(criteria.num_pitches - 1);
  for (int i = 0; i < 11; i++) {
    gsl_vector_set(target_tuning, i, werck_iii[i]);
  }

  gsl_vector* tuning = gsl_vector_alloc(criteria.num_pitches);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  int iter = 0;
  int status;
  double size;

  OWTSearchParams params;
  params.criteria = &criteria;
  params.tuning = target_tuning;

  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
  /* gsl_rng_set(r, 1000); */
  gsl_rng_set(r, mach_absolute_time());
  // Initializing the key weights so they sum to 1
  x = gsl_vector_alloc(criteria.num_pitches);
  for (int i = 0; i < criteria.num_pitches; i++) {
    gsl_vector_set(x, i, gsl_rng_uniform_pos(r));
  }


  gsl_vector_fprintf(stdout, x, " %1.3f ");

  // Step by 1c. each iteration
  ss = gsl_vector_alloc(criteria.num_pitches);
  gsl_vector_set_all(ss, 0.1);

  // Initialize the func
  minex_func.n = criteria.num_pitches;
  minex_func.f = &owt_deviation_from_tuning;
  minex_func.params = (void *)&params;

  s = gsl_multimin_fminimizer_alloc(T, criteria.num_pitches);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  OWTResults results;
  results = owt_optimize_temperament(params.criteria);
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    normalize_vector(s->x, 0.0, 1.0);

    /* double range = gsl_vector_max(s->x) - gsl_vector_min(s->x); */
    /* gsl_vector_add_constant(s->x, 0.0 - gsl_vector_min(s->x)); */
    /* gsl_vector_scale(s->x, 1.0 / range); */
    
    if (status) break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-2);

    if (status == GSL_SUCCESS) {
      /* printf("Success on iteration %i!\n", iter); */
      for (int i = 0; i < (params.criteria->num_pitches); i++) {
        params.criteria->key_weights[i] = gsl_vector_get(x, i);
      }
      results = owt_optimize_temperament(params.criteria);
      /* printf("Converged to minimum at\n"); */
      /* gsl_vector_fprintf(stdout, s->x, " %1.3f "); */
      /* printf("\nChi-squared is: %f\n", results.chisq); */

      // If it's really bad, start over...
      /* if (results.chisq > 50.0) { */
      /*   #<{(| printf("Starting over!\n"); |)}># */
      /*   for (int i = 0; i < criteria.num_pitches; i++) { */
      /*     gsl_vector_set(x, i, gsl_rng_uniform_pos(r)); */
      /*   } */
      /*   gsl_multimin_fminimizer_set(s, &minex_func, x, ss); */
      /*   status = GSL_CONTINUE; */
      /* } */
    }
  } while (status == GSL_CONTINUE && iter < 100000);
  printf("Converged to minimum on iteration %i\n", iter);
  gsl_vector_fprintf(stdout, s->x, " %1.3f ");
  printf("Found tuning...\n");
  gsl_vector_fprintf(stdout, results.optimal_tuning, " %1.3f ");
  printf("\nChi-squared is: %f\n", results.chisq);

  return 0;
}
