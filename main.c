#include <stdio.h>
#include <string.h>

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

  /* printf("Path to output file is %s\n", path_to_output_file); */

  /* Things we need for loading the tuning */
  char* tuning_name = NULL;
  char* title = NULL;
  double* ideal_intervals;
  double* interval_weights;
  double* key_weights;
  double repeat_factor;

  Criteria criteria = owt_parse(path_to_file);

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

  double chisq, third_error_sd, fifth_error_max;

  gsl_vector* c, *third_error, *fifth_error;

  c = owt_optimize_temperament(criteria.num_pitches, 
      criteria.ideal_intervals, criteria.repeat_factor, 
      criteria.interval_weights, criteria.key_weights, &chisq);

  third_error = owt_interval_error(criteria.num_pitches, 4, 386.0, criteria.repeat_factor, c);
  double third_error_data[criteria.num_pitches];

  printf("\nThirds error:\n");
  for (int i = 0; i < criteria.num_pitches; i++) {
    third_error_data[i] = gsl_vector_get(third_error, i);
    printf("%0.2f\t", third_error_data[i]);

  }

  third_error_sd = gsl_stats_sd(third_error_data, 1, criteria.num_pitches);

  fifth_error = owt_interval_error(criteria.num_pitches, 7, 702.0, criteria.repeat_factor, c);
  fifth_error_max = fabs(gsl_vector_get(fifth_error, gsl_blas_idamax(fifth_error)));
  printf("\nFifths error:\n");
  for (int i = 0; i < criteria.num_pitches; i++)
    printf("%0.2f\t", gsl_vector_get(fifth_error, i));
  
  printf("\nChi-squared results: %f\n", chisq);
  printf("Best fit tuning: \n");
  for (int i = 0; i < criteria.num_pitches-1; i++) {
    printf("%.2f\t", gsl_vector_get(c, i));
  }
  printf("\nStandard deviation: %0.2f\n", third_error_sd);
  printf("\nWorst fifth: %0.2f\n", fifth_error_max);
  printf("\n");

  FILE *output_file = fopen(path_to_output_file, "w");
  if ( output_file == NULL) {
    printf("Output file cannot be written to\n");
  } else {
    fprintf(output_file, "!\n%s\n %d\n!\n", criteria.title, criteria.num_pitches);
    for (int i = 0; i < criteria.num_pitches-1; i++)
      fprintf(output_file, "%0.3f\n", gsl_vector_get(c, i));
    fprintf(output_file, "%0.3f\n", criteria.repeat_factor);
    fclose(output_file);
  }

  gsl_vector_free(c);
  gsl_vector_free(third_error);
  gsl_vector_free(fifth_error);

  return 0;
}

