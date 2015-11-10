#include <stdio.h>

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "methods.h"

double get_de(double *observed, double *expected, int len)
{
  int i = 0;
  double total = 0.0;
  double subt = 0.0;
  for (i = 0; i < len; i++) {
    subt = observed[i] - expected[i];
    total += pow(subt, 2);
  }

  return sqrt(total / len);
}

int iupac_match(char c1, char c2)
{
  c1 = toupper(c1);
  c2 = toupper(c2);
  if (c1 == 'U') { c1 = 'T'; }
  if (c2 == 'U') { c2 = 'T'; }

  if (c1 == c2) {
    return 1;
  } else {
    switch (c1) {
    case 'A':
      if (c2 == 'R' ||
          c2 == 'W' ||
          c2 == 'M' ||
          c2 == 'D' ||
          c2 == 'H' ||
          c2 == 'V' ||
          c2 == 'N') {
        return 1;
      } else {
        return 0;
      }
    case 'C':
      if (c2 == 'Y' ||
          c2 == 'S' ||
          c2 == 'M' ||
          c2 == 'B' ||
          c2 == 'H' ||
          c2 == 'V' ||
          c2 == 'N') {
        return 1;
      } else {
        return 0;
      }
    case 'G':
      if (c2 == 'Y' ||
          c2 == 'S' ||
          c2 == 'M' ||
          c2 == 'B' ||
          c2 == 'H' ||
          c2 == 'V' ||
          c2 == 'N') {
        return 1;
      } else {
        return 0;
      }
    case 'T':
      if (c2 == 'Y' ||
          c2 == 'W' ||
          c2 == 'K' ||
          c2 == 'B' ||
          c2 == 'D' ||
          c2 == 'H' ||
          c2 == 'N') {
        return 1;
      } else {
        return 0;
      }
    case 'R':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'G') {
        return 1;
      } else {
        return 0;
      }
    case 'Y':
      if (c2 == 'N' ||
          c2 == 'C' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'S':
      if (c2 == 'N' ||
          c2 == 'G' ||
          c2 == 'C') {
        return 1;
      } else {
        return 0;
      }
    case 'W':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'K':
      if (c2 == 'N' ||
          c2 == 'G' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'M':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'C') {
        return 1;
      } else {
        return 0;
      }
    case 'B':
      if (c2 == 'N' ||
          c2 == 'C' ||
          c2 == 'G' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'D':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'G' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'H':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'C' ||
          c2 == 'T') {
        return 1;
      } else {
        return 0;
      }
    case 'V':
      if (c2 == 'N' ||
          c2 == 'A' ||
          c2 == 'C' ||
          c2 == 'G') {
        return 1;
      } else {
        return 0;
      }
    case 'N':
      if (c2 == 'A' ||
          c2 == 'C' ||
          c2 == 'T' ||
          c2 == 'G' ||
          c2 == 'R' ||
          c2 == 'Y' ||
          c2 == 'S' ||
          c2 == 'W' ||
          c2 == 'K' ||
          c2 == 'M' ||
          c2 == 'B' ||
          c2 == 'D' ||
          c2 == 'H' ||
          c2 == 'V') {
        return 1;
      } else {
        return 0;
      }
    default:
      return 0;
    }
  }
}

/* TODO doesn't ensure strings are the same length */
double perc_mismatch(char *str1, char *str2)
{
  int num_mismatches = 0;
  int i = 0;
  int posns_not_counted = 0;

  for (i = 0; str1[i] != '\0'; i++) {
    if (str1[i] == '.' || str2[i] == '.') {
      posns_not_counted += 1;
    } else if (iupac_match(str1[i], str2[i]) == 0) {
      num_mismatches += 1;
    }
  }

  return num_mismatches / (double) (i - posns_not_counted) * 100;
}

int num_windows(char *str)
{
  int i = 0;

  for (i = 0; str[i] != '\0'; i++) { ; }

  return (int) ceil((i - WINDOW_SIZE + 1) / (double) WINDOW_STEP);
}

int get_num_windows(int window_size, int window_step, int arr_len)
{
  return (int) ceil((arr_len - window_size + 1) / (double) window_step);
}

/* TODO this has a lot of duplication with windows function */
struct Windows *double_windows(double *arr, int len)
{
  int i = 0;
  int j = 0;
  int window_idx = 0;
  int window_num = 0;
  int max = 0;
  int start = 0;
  double posn = 0.0;

  int num_windows =
    get_num_windows(WINDOW_SIZE, WINDOW_STEP, len);

  struct Windows *windows = windows_new(num_windows);

  /* this max is inclusive */
  max = len - WINDOW_SIZE;

  for (i = 0; i <= max; i += WINDOW_STEP) { /* each window start */
    start = i;
    posn = ((start+1) + (start+1 + WINDOW_SIZE-1)) / 2.0;

    windows->xposns[window_num] = posn;

    for (j = start; j < start + WINDOW_SIZE; j++) {
      windows->dwindows[window_num][window_idx++] = arr[j];
    }
    window_num++;
    window_idx = 0;
  }

  return windows;
}

struct Windows *windows(char *str)
{
  int len = 0;
  int i = 0;
  int j = 0;
  int window_idx = 0;
  int window_num = 0;
  int max = 0;
  int start = 0;
  double posn = 0.0;

  char *window = calloc(WINDOW_SIZE+1, sizeof(char));

  /* TODO replace with strlen */
  /* get str len */
  for (i = 0; str[i] != '\0'; i++) { ; }
  len = i;

  int num_windows =
    get_num_windows(WINDOW_SIZE, WINDOW_STEP, len);

  struct Windows *windows = windows_new(num_windows);

  /* this max is inclusive */
  max = len - WINDOW_SIZE;

  for (i = 0; i <= max; i += WINDOW_STEP) { /* each window start */
    start = i;
    posn = ((start+1) + (start+1 + WINDOW_SIZE-1)) / 2.0;

    windows->xposns[window_num] = posn;

    for (j = start; j < start + WINDOW_SIZE; j++) {
      window[window_idx++] = str[j];
    }
    window[window_idx] = '\0';

    strcpy(windows->windows[window_num++], window);
    window_idx = 0;
  }

  free(window);

  return windows;
}

struct Windows *windows_new(int num_windows)
{
  int i = 0;

  struct Windows *windows = malloc(sizeof(struct Windows));
  assert(windows != NULL);

  windows->num_windows = num_windows;

  windows->xposns = malloc(num_windows * sizeof(double));
  assert(windows->xposns != NULL);

  windows->windows = malloc(num_windows * sizeof(char *));
  assert(windows->windows != NULL);

  windows->dwindows = malloc(num_windows * sizeof(double *));
  assert(windows->dwindows != NULL);

  for (i = 0; i < num_windows; i++) {
    windows->windows[i] = calloc(WINDOW_SIZE+1, sizeof(char));
    assert(windows->windows[i] != NULL);

    windows->dwindows[i] = calloc(WINDOW_SIZE+1, sizeof(double));
    assert(windows->dwindows[i] != NULL);
  }

  return windows;
}

void windows_free(struct Windows *windows)
{
  int i = 0;

  assert(windows != NULL);

  assert(windows->xposns != NULL);
  free(windows->xposns);

  assert(windows->windows != NULL);
  assert(windows->dwindows != NULL);

  for (i = 0; i < windows->num_windows; i++) {
    assert(windows->windows[i] != NULL);
    free(windows->windows[i]);

    assert(windows->dwindows[i] != NULL);
    free(windows->dwindows[i]);
  }

  free(windows->windows);
  free(windows->dwindows);

  free(windows);
}

void print_array(double *arr, int len)
{
  int i = 0;
  for (i = 0; i < len; i++ ) {
    printf("%f ", arr[i]);
  }

  printf("\n");

}


struct TwoWindows *twds_new(struct Windows *w1, struct Windows *w2)
{
  struct TwoWindows *twds = malloc(sizeof(struct TwoWindows));
  assert(twds != NULL);

  assert(w1->num_windows == w2->num_windows);
  twds->num_windows = w1->num_windows;

  twds->xposns = w1->xposns;

  twds->mismatches = malloc(twds->num_windows * sizeof(double));

  assert(w1 != NULL);
  twds->w1 = w1;

  assert(w2 != NULL);
  twds->w2 = w2;

  return twds;
}

void twds_free(struct TwoWindows *twds)
{
  assert(twds != NULL);

  assert(twds->w1 != NULL);
  windows_free(twds->w1);

  assert(twds->w2 != NULL);
  windows_free(twds->w2);

  assert(twds->mismatches != NULL);
  free(twds->mismatches);

  free(twds);
}

void twds_set_perc_mismatches(struct TwoWindows *twds)
{
  int which_window = 0;
  int num = twds->num_windows;
  double mismatch = 0.0;

  char *w1, *w2;

  for (which_window = 0; which_window < num; which_window++) {
    w1 = twds->w1->windows[which_window];
    w2 = twds->w2->windows[which_window];

    mismatch = perc_mismatch(w1, w2);
    twds->mismatches[which_window] = mismatch;
  }
}

/* frees the twds afterwards */
struct MismatchInfo *minfo_new(struct TwoWindows *twds)
{
  int i = 0;
  assert(twds != NULL);

  assert(twds->xposns != NULL);
  assert(twds->mismatches != NULL);

  struct MismatchInfo *minfo = malloc(sizeof(struct MismatchInfo));

  minfo->xposns = malloc(twds->num_windows * sizeof(double));
  assert(minfo->xposns != NULL);
  for (i = 0; i < twds->num_windows; i++) {
    minfo->xposns[i] = twds->xposns[i];
  }

  minfo->mismatches = malloc(twds->num_windows * sizeof(double));
  assert(minfo->mismatches != NULL);
  for (i = 0; i < twds->num_windows; i++) {
    minfo->mismatches[i] = twds->mismatches[i];
  }

  minfo->num_windows = twds->num_windows;

  twds_free(twds);

  return minfo;
}

void minfo_free(struct MismatchInfo *minfo)
{
  assert(minfo->xposns != NULL);
  free(minfo->xposns);

  assert(minfo->mismatches != NULL);
  free(minfo->mismatches);

  free(minfo);
}

struct MismatchInfo *windowed_str_mismatch(char *str1, char *str2)
{
  struct Windows *w1, *w2;
  struct TwoWindows *twds;
  struct MismatchInfo *minfo;

  w1 = windows(str1);
  w2 = windows(str2);

  assert(w1->num_windows == w2->num_windows);

  twds = twds_new(w1, w2);

  twds_set_perc_mismatches(twds);

  minfo = minfo_new(twds);

  return minfo;
}

double arr_mean(double *arr, int len)
{
  double total = 0.0;
  int i = 0;

  for (i = 0; i < len; i++) {
    total += arr[i];
  }

  return total / len;
}

/* int main() */

/* { */
/*   int i = 0; */

/*   struct Windows *the_windows; */
/*   int num = num_windows("1234567890"); */

/*   the_windows = windows("1234567890"); */

/*   for (i = 0; i < num; i++) { */
/*     printf("%f %s\n", the_windows->xposns[i], the_windows->windows[i]); */
/*   } */

/*   struct MismatchInfo *minfo = */
/*     windowed_str_mismatch("AAAACCCCTTTTGGGG", "AAACCCCTTTTGGGGA"); */

/*   assert(minfo != NULL); */

/*   for (i = 0; i < minfo->num_windows; i++) { */
/*     printf("%d %f %f\n", i, minfo->xposns[i], minfo->mismatches[i]); */
/*   } */

/*   /\* print_array(minfo->xposns, minfo->num_windows); *\/ */
/*   /\* print_array(minfo->mismatches, minfo->num_windows); *\/ */

/*   minfo_free(minfo); */
/*   windows_free(the_windows); */
/*   return 0; */
/* } */
