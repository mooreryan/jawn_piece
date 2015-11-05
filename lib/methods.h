int iupac_match(char c1, char c2);
double perc_mismatch(char *str1, char *str2);

struct Windows {
  int num_windows;
  double *xposns;
  char **windows;
};

struct Windows *windows_new(int num_windows);
void windows_free(struct Windows *num_windows);
struct Windows *windows(char *str);

struct TwoWindows {
  int num_windows;
  struct Windows *w1;
  struct Windows *w2;
  double *xposns;
  double *mismatches;
};

struct TwoWindows *twds_new(struct Windows *w1, struct Windows *w2);
void twds_free(struct TwoWindows *twds);
void twds_set_perc_mismatches(struct TwoWindows *twds);

struct MismatchInfo {
  int num_windows;
  double *xposns;
  double *mismatches;
};

struct MismatchInfo *windowed_str_mismatch(char *str1, char *str2);
struct MismatchInfo *minfo_new(struct TwoWindows *twds);
void minfo_free(struct MismatchInfo *minfo);
