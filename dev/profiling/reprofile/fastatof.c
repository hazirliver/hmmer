/* Prototype: exact fast path for atof() on plain decimal tokens.
 * Clinger's fast path: a token [+-]?digits[.digits] with at most 15
 * significant digits and at most 22 fraction digits has value m / 10^e
 * with m < 2^53 and 10^e exact in double, so one IEEE division gives the
 * correctly rounded result -- exactly what glibc's strtod returns.
 * Anything else goes to atof().
 *   gcc -O3 -o fastatof fastatof.c && ./fastatof Pfam-A.hmm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

static const double p10[23] = { 1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,
                                1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22 };
static long nfall = 0;

static inline double
fast_atof(const char *s)
{
  const char *p = s;
  int      neg = 0, nd = 0, ne = 0, sawdig = 0;
  uint64_t m = 0;

  if (*p == '-') { neg = 1; p++; } else if (*p == '+') p++;
  while (*p >= '0' && *p <= '9') { if (m || *p != '0') { m = m*10 + (*p - '0'); nd++; } sawdig = 1; p++; if (nd > 15) goto fallback; }
  if (*p == '.') {
    p++;
    while (*p >= '0' && *p <= '9') { if (m || *p != '0') { m = m*10 + (*p - '0'); nd++; } ne++; sawdig = 1; p++; if (nd > 15 || ne > 22) goto fallback; }
  }
  if (!sawdig || *p != '\0') goto fallback;   /* whole token must be a plain decimal */
  { double d = (double) m / p10[ne]; return neg ? -d : d; }
 fallback:
  nfall++;
  return atof(s);
}

static int check(const char *tok, long *nbad)
{
  double a = atof(tok), b = fast_atof(tok);
  if (memcmp(&a, &b, sizeof a) != 0) { if (*nbad < 10) printf("MISMATCH '%s' atof=%a fast=%a\n", tok, a, b); (*nbad)++; return 1; }
  return 0;
}

int
main(int argc, char **argv)
{
  long nbad = 0, ntok = 0;
  /* 1. every whitespace token of the given files */
  for (int f = 1; f < argc; f++) {
    FILE *fp = fopen(argv[f], "r"); char tok[256];
    if (!fp) { perror(argv[f]); return 1; }
    while (fscanf(fp, "%255s", tok) == 1) { check(tok, &nbad); ntok++; }
    fclose(fp);
  }
  printf("file tokens: %ld checked, %ld mismatches, %ld fallbacks\n", ntok, nbad, nfall);

  /* 2. random decimal strings, incl. hard-rounding cases (long digit strings) */
  srand(7); long nr = 0; nfall = 0;
  char buf[64];
  for (long i = 0; i < 50000000; i++) {
    int ni = rand() % 6, nf = rand() % 18, k = 0;
    if (rand() % 4 == 0) buf[k++] = '-';
    for (int j = 0; j < ni; j++) buf[k++] = '0' + rand() % 10;
    if (nf || ni == 0) { buf[k++] = '.'; for (int j = 0; j < nf; j++) buf[k++] = '0' + rand() % 10; }
    buf[k] = '\0';
    check(buf, &nbad); nr++;
  }
  printf("random tokens: %ld checked, total mismatches %ld, %ld fallbacks\n", nr, nbad, nfall);

  /* 3. speed on typical HMM values */
  const char *v[] = { "2.35312", "0.00000", "3.46071", "1.38629", "0.97875", "4.20134", "12.0021", "0.69315" };
  double acc = 0; clock_t t0 = clock();
  for (long i = 0; i < 20000000; i++) acc += atof(v[i & 7]);
  double ta = (double)(clock()-t0)/CLOCKS_PER_SEC; t0 = clock();
  for (long i = 0; i < 20000000; i++) acc += fast_atof(v[i & 7]);
  double tb = (double)(clock()-t0)/CLOCKS_PER_SEC;
  printf("20M conversions: atof %.2fs (%.0f ns each), fast %.2fs (%.0f ns each)  [%g]\n", ta, ta*50, tb, tb*50, acc);
  return nbad != 0;
}
