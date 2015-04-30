#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "owt_parser.h"

#define MAX_LINES 256
#define MAX_LINE_LENGTH 1024

/* returns the number of pitches */
OWTCriteria owt_parse(char* path_to_file) {
  /* Read in the tuning file */
  FILE *tuning_file = fopen(path_to_file, "r");
  if ( tuning_file == NULL ) {
    printf("File %s can't be read!\n", path_to_file);
  }

  char* current_line = malloc(MAX_LINE_LENGTH);
  char* token;
  char* token_source;
  char* tokenized_line[MAX_LINES];
  int num_tokens;
  char* significant_tokens = malloc(MAX_LINE_LENGTH);

  OWTCriteria criteria;

  for (int i = 0; i < MAX_LINES; i++)
    tokenized_line[i] = malloc(MAX_LINE_LENGTH);

  /* Iterate over each line */
  while ( fgets(current_line, MAX_LINE_LENGTH, tuning_file) != NULL ) {
    token_source = current_line;
    num_tokens = 0;

    /* Then break each line into tokens */
    while ( (token = strtok(token_source, " \n")) != NULL ) {
      strcpy(tokenized_line[num_tokens], token);
      token_source = NULL;
      num_tokens++;
    }

    /* Then iterate over the tokens and parse */
    for (int i = 0; i < num_tokens; i++) {
      /* ignore comments */
      if (strcmp(tokenized_line[i], "#") == 0) break;

      if (strcmp(tokenized_line[i], "title:") == 0) {
        criteria.title = malloc(256);
        while (++i < num_tokens) {
          strcat(criteria.title, tokenized_line[i]);
          strcat(criteria.title, " ");
        }
        break;
      } 
      
      if ((strcmp(tokenized_line[i], "num_pitches:") == 0) && (++i < num_tokens)) {
        /* i is now one larger */
        criteria.num_pitches = atoi(tokenized_line[i]);
        /* num_pitches means we should also load in the ideals and weights */ 
        criteria.ideal_intervals = malloc(sizeof(double) * (criteria.num_pitches-1));
        criteria.interval_weights = malloc(sizeof(double) * (criteria.num_pitches-1));
        criteria.key_weights = malloc(sizeof(double) * criteria.num_pitches);
        break;
      } 

      if (strcmp(tokenized_line[i], "repeat_factor:") == 0) {
        if (++i < num_tokens) {
          criteria.repeat_factor = atof(tokenized_line[i]);
        }
      }
      
      if (strcmp(tokenized_line[i], "ideal_intervals:") == 0) {
        int current_interval = 0;
        while (++i < num_tokens) {
          criteria.ideal_intervals[current_interval++] = atof(tokenized_line[i]);
        }
        break;
      }

      if (strcmp(tokenized_line[i], "interval_weights:") == 0) {
        int current_interval = 0;
        while (++i < num_tokens) {
          criteria.interval_weights[current_interval++] = atof(tokenized_line[i]);
        }
        break;
      }

      if (strcmp(tokenized_line[i], "key_weights:") == 0) {
        int current_interval = 0;
        while (++i < num_tokens) {
          criteria.key_weights[current_interval++] = atof(tokenized_line[i]);
        }
        break;
      }
    }
  }
  return criteria;
}
