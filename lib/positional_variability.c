#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* see https://www.securecoding.cert.org/confluence/display/c/
   INT05-C.+Do+not+use+input+functions+to+convert+character+data+
   if+they+cannot+handle+all+possible+inputs */

/* TODO zintail problem, most of the posistions have . or - and not
   bases in the alignment */

int main(int argc, char *argv[])
{
  int i = 0;
  int c = 0;
  int max = -1;

  int seqno = 0;
  int posno = 0;

  int num_seqs = 0;
  int seq_len = 0;

  char *seq;
  char **seqs;
  char **bases_at_pos;

  /* TODO assumes ascii */
  int base_counts[128];
  /* zero the base counts */
  for (i = 0; i < 128; i++) {
    base_counts[i] = 0;
  }

  if (scanf("%d\t%d\n", &num_seqs, &seq_len) == 2) {
    seq = calloc((seq_len + 1), sizeof(char));

    /* initialize seqs */
    seqs = malloc(num_seqs * sizeof(char *));
    for (i = 0; i < num_seqs; i++) {
      seqs[i] = calloc((seq_len + 1), sizeof(char));
    }

    /* init bases_at_pos, used to store the transpose of seqs */
    bases_at_pos = malloc(seq_len * sizeof(char *));
    for (i = 0; i < seq_len; i++) {
      bases_at_pos[i] = calloc(num_seqs, sizeof(char));
    }
  }

  /* read rest of infile */
  while ((c = getchar()) != EOF) {
    if (c == '\n') {
      seq[posno] = '\0';

      strcpy(seqs[seqno++], seq);

      posno = 0;
    } else {
      seq[posno++] = c;
    }
  }

  /* TODO can i remove this loop? */
  /* transpose the seqs */
  for (posno = 0; posno < seq_len; posno++) {
    for (seqno = 0; seqno < num_seqs; seqno++) {
      bases_at_pos[posno][seqno] = seqs[seqno][posno];
    }
  }

  /* count the base freqs per position */
  for (posno = 0; posno < seq_len; posno++) {

    /* count this posistion */
    for (seqno = 0; seqno < num_seqs; seqno++) {
      base_counts[(int) bases_at_pos[posno][seqno]] += 1;
    }

    /* find max base */
    max = -1;
    for (i = 0; i < 128; i++) {
      if (base_counts[i] > max) {
        max = base_counts[i];
      }
    }

    assert(max != -1);

    /* print positional variability */
    printf("%.2f\n", (1 - max / (double) num_seqs) * 100);

    /* reset the base counts */
    for (i = 0; i < 128; i++) {
      base_counts[i] = 0;
    }
  }

  /* clean up */
  free(seq);

  for (posno = 0; posno < seq_len; posno++) {
    free(bases_at_pos[posno]);
  }
  free(bases_at_pos);

  for (seqno = 0; seqno < num_seqs; seqno++) {
    free(seqs[seqno]);
  }
  free(seqs);

  return 0;
}
