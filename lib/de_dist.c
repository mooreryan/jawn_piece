#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* #define MASKED_DATABASE "assets/masked_db.fa" */
#define MASKED_DATABASE "test.fa"

int main(int argc, char *argv[])
{
  FILE *db_f;

  int total_seqs = 0;
  int seq_len = 0;

  char *seq;
  char **seqs;

  int baseno = 0;
  int current_seq = 0;

  int c = 0;
  int i = 0;
  int j = 0;

  db_f = fopen(MASKED_DATABASE, "r");
  assert(db_f != NULL);

  if (fscanf(db_f, "%d\t%d\n", &total_seqs, &seq_len) == 2) {
    seqs = malloc(total_seqs * sizeof(char *));
    seq  = calloc((seq_len + 1), sizeof(char));

    for (i = 0; i < total_seqs; i++) {
      seqs[i] = calloc((seq_len + 1), sizeof(char));
    }
  }

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

  for (i = 0; i < total_seqs; i++) {
    for (j = i + 1; j < total_seqs; j++) {
      /* seqs[i], seqs[j] */
    }
  }

  /* Clean up */
  for (i = 0; i < total_seqs; i++) {
    free(seqs[i]);
  }
  free(seqs);
  free(seq);
  fclose(db_f);

  return 0;
}
