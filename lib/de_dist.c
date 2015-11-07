#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "methods.h"

/* #define MASKED_DATABASE "assets/masked_db.fa" */
#define MASKED_DATABASE "/Users/moorer/projects/ZetaHunter/assets/masked_db.fa"

int main(int argc, char *argv[])
{
  FILE *db_f;

  int total_seqs = 0;
  int seq_len = 0;

  char *seq;
  char **seqs;
  char **bases_at_posn;

  int baseno = 0;
  int current_seq = 0;

  int c = 0;
  int i = 0;
  int j = 0;
  int wind_i = 0;

  double obs_evol_dist = 0.0;
  double database_overall_evol_dist = 0.0;
  double fitting_coef = 1.0;

  struct MismatchInfo *minfo;

  int pos = 0;
  int seq_num = 0;
  int num_a = 0;
  int num_c = 0;
  int num_t = 0;
  int num_g = 0;
  double max = 0;
  double total_counts = 0;
  double freqs[4];
  char this_base;
  double variability;
  struct Windows *dwindows;
  double total = 0.0;
  int window_size = 300;
  int win_i = 0;

  double de_total = 0.0;
  double subtraction = 0.0;
  double de = 0.0;



  /* TODO test for empyt file, test for file exists */
  db_f = fopen(MASKED_DATABASE, "r");
  assert(db_f != NULL);

  /* get num seqs and seq length */
  if (fscanf(db_f, "%d\t%d\n", &total_seqs, &seq_len) == 2) {
    seqs          = malloc(total_seqs * sizeof(char *));
    seq           = calloc((seq_len + 1), sizeof(char));

    bases_at_posn = malloc(seq_len * sizeof(char *));

    /* intialize posns, one more for '\0', just to make it work as a
       string too TODO prolly dont need this */
    for (i = 0; i < seq_len; i++) {
      bases_at_posn[i] = calloc((total_seqs + 1), sizeof(char));
    }

    /* initialize seqs */
    for (i = 0; i < total_seqs; i++) {
      seqs[i] = calloc((seq_len + 1), sizeof(char));
    }
  }

  double variabilities[seq_len];

  /* read the infile */
  while ((c = getc(db_f)) != EOF) {
    if (c == '\n') {
      seq[baseno] = '\0';

      strcpy(seqs[current_seq++], seq);

      baseno = 0;
    } else {
      seq[baseno++] = c;
    }
  }

  assert (current_seq == total_seqs);

  /* i is posn, j is seq */
  for (i = 0; i < seq_len; i++) {
    for (j = 0; j < total_seqs; j++) {
      /* TODO check this */
      bases_at_posn[i][j] = seqs[j][i];
    }
  }

  for (pos = 0; pos < seq_len; pos++) {
    num_a = 0;
    num_c = 0;
    num_t = 0;
    num_g = 0;
    for (seq_num = 0; seq_num < total_seqs; seq_num++) {
      this_base = bases_at_posn[pos][seq_num];

      /* TODO make sure seqs have only ACTG and capital  */
      switch(this_base) {
      case 'A':
        num_a += 1;
        break;
      case 'C':
        num_c += 1;
        break;
      case 'T':
        num_t += 1;
        break;
      case 'G':
        num_g += 1;
        break;
      case '-':
        ;
        break;
      default:
        fprintf(stderr, "WARNING: base %c was not counted\n", this_base);
      }
    }

    /* get max value */
    double freqs[] = { num_a / (double) total_seqs,
                       num_c / (double) total_seqs,
                       num_t / (double) total_seqs,
                       num_g / (double) total_seqs };

    max = -1;
    for (i = 0; i < 4; i++) {
      if (freqs[i] > max) {
        max = freqs[i];
      }
    }
    assert (max != -1);

    variability = (1 - ((max - 0.25) / 0.75)) * 100;
    variabilities[pos] = variability;
  }
  /* variabilities holds the positional variability */
  dwindows = double_windows(variabilities, seq_len);
  /* TODO make WINDOW_SIZE A GLOBAL */

  double windowed_avg_probs[dwindows->num_windows];
  double exp_perc_diffs[dwindows->num_windows];

  for (i = 0; i < dwindows->num_windows; i++) {
    windowed_avg_probs[i] = arr_mean(dwindows->dwindows[i], window_size);
  }

  database_overall_evol_dist = arr_mean(windowed_avg_probs,
                                        dwindows->num_windows);

  printf("database_overall_evol_dist: %.3f\n",
         database_overall_evol_dist);


  for (i = 0; i < total_seqs; i++) {
    for (j = i + 1; j < total_seqs; j++) {
      minfo = windowed_str_mismatch(seqs[i], seqs[j]);
      obs_evol_dist = arr_mean(minfo->mismatches, minfo->num_windows);
      fitting_coef = obs_evol_dist / database_overall_evol_dist;

      de_total = 0;
      for (win_i = 0; win_i < dwindows->num_windows; win_i++) {
        de_total += pow(minfo->mismatches[win_i] -
                        (fitting_coef * windowed_avg_probs[win_i]),
                        2);
      }
      de = sqrt(de_total / dwindows->num_windows);

      printf("%d %d %f %f\n",
             i,
             j,
             obs_evol_dist,
             de);

      minfo_free(minfo);
    }
  }

  /* Clean up */
  for (i = 0; i < seq_len; i++) {
    free(bases_at_posn[i]);
  }
  free(bases_at_posn);

  for (i = 0; i < total_seqs; i++) {
    free(seqs[i]);
  }
  free(seqs);

  free(seq);
  fclose(db_f);
  windows_free(dwindows);
  return 0;
}
