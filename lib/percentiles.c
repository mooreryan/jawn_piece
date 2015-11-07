#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define MAX_LEN 1000

int cmp (const void * a, const void * b)
{
  double f = *((double*) a);
  double s = *((double*) b);

  return (f > s) - (f < s);
}

int eq(double d)
{
  return (d - -1) < 0.00001;
}

int main(int argc, char *argv[])
{
  printf("pgroup\tcount\tdist\t50\t75\t90\t95\t99\n");

  char num_str[MAX_LEN];
  int total = 0;
  double dist = 0.0;
  double de = 0.0;
  int i = 0;
  int j = 0;
  int c = 0;

  double one = -1;
  double two = -1;
  double three = -1;
  double four = -1;
  double five = -1;

  int charno = 0;
  int lineno = 0;
  int group_des_i = 0;

  double *percentiles, *dists, *des, *des_percs;
  int *pgroups;

  double *group_des;

  int total_comparisons = 0;

  /* the first line is the total number of comparisons */
  while ((c = getchar()) != '\n') {
    num_str[charno++] = c;
  }
  num_str[charno] = '\0';
  total = atoi(num_str);

  percentiles = malloc(total * sizeof(double));
  assert(percentiles != NULL);

  dists = malloc(total * sizeof(double));
  assert(dists != NULL);

  des = malloc(total * sizeof(double));
  assert(des != NULL);

  /* TODO this can be smaller than total */
  group_des = malloc(total * sizeof(double));
  assert(des != NULL);

  pgroups = malloc(total * sizeof(int));
  assert(pgroups != NULL);

  /* set up the dists and des arrays, the input is already sorted */
  while (scanf("%lf\t%lf", &dist, &de) == 2) {
    percentiles[lineno] = lineno / (double) total * 100;
    pgroups[lineno] = (int) floor(percentiles[lineno]);
    dists[lineno] = dist;
    des[lineno++] = de;
  }

  for (i = 0; i < total; i++) {
    /* new perc group */
    if (i != 0 && pgroups[i-1] != pgroups[i]) {

      /* calculate the de_group percentiles */
      qsort(group_des, group_des_i, sizeof(double), cmp);

      des_percs = malloc(group_des_i * sizeof(double));
      for (j = 0; j < group_des_i; j++) {
        des_percs[j] = j / (double) group_des_i * 100;
        switch ((int) floor(des_percs[j])) {
        case 50:
          if (eq(one)) { one = group_des[j]; }
          break;
        case 75:
          if (eq(two)) { two = group_des[j]; }
          break;
        case 90:
          if (eq(three)) { three = group_des[j]; }
          break;
        case 95:
          if (eq(four)) { four = group_des[j]; }
          break;
        case 99:
          if (eq(five)) { five = group_des[j]; }
          break;
        }
      }
      free(des_percs);

      assert(one != -1);
      assert(two != -1);
      assert(three != -1);
      assert(four != -1);
      assert(five != -1);

      printf("%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
             pgroups[i-1],
             group_des_i,
             dists[i-1],
             one,
             two,
             three,
             four,
             five);
      /* end calculate the de_group percentiles */

      /* reset stuff */
      total_comparisons += group_des_i;
      group_des_i = 0;
      group_des[group_des_i++] = des[i];
      one = -1;
      two = -1;
      three = -1;
      four = -1;
      five = -1;
    } else if (i == total - 1) { /* the last one */
      /* calculate the de_group percentiles */
      qsort(group_des, group_des_i, sizeof(double), cmp);

      des_percs = malloc(group_des_i * sizeof(double));
      for (j = 0; j < group_des_i; j++) {
        des_percs[j] = j / (double) group_des_i * 100;
        switch ((int) floor(des_percs[j])) {
        case 50:
          if (one == -1) { one = group_des[j]; }
          break;
        case 75:
          if (two == -1) { two = group_des[j]; }
          break;
        case 90:
          if (three == -1) { three = group_des[j]; }
          break;
        case 95:
          if (four == -1) { four = group_des[j]; }
          break;
        case 99:
          if (five == -1) { five = group_des[j]; }
          break;
        }
      }
      free(des_percs);

      assert(one != -1);
      assert(two != -1);
      assert(three != -1);
      assert(four != -1);
      assert(five != -1);

      printf("%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
             pgroups[i-1],
             group_des_i,
             dists[i-1],
             one,
             two,
             three,
             four,
             five);
      /* end calculate the de_group percentiles */

      total_comparisons += group_des_i;
    } else {
      group_des[group_des_i++] = des[i];
    }
  }

  /* TODO im missing one comparison */
  fprintf(stderr, "total comparisons: %d\n", total_comparisons);

  free(percentiles);
  free(dists);
  free(des);
  free(group_des);
  free(pgroups);

  return 0;
}
