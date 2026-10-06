/* Prototype: exact replacement for (float) log((double) x), x a float.
 * Fast path: fdlibm-style double log (error < 1 ulp of double), then a
 * check that the double result is not near a float rounding midpoint;
 * otherwise (or for x <= 0, subnormal, inf, nan, x == 1) call libm.
 * Exhaustive test over all 2^32 float bit patterns.
 *   gcc -O3 -o fastlog fastlog.c -lm && ./fastlog
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

static const double
  ln2_hi = 6.93147180369123816490e-01, /* 3fe62e42 fee00000 */
  ln2_lo = 1.90821492927058770002e-10, /* 3dea39ef 35793c76 */
  Lg1 = 6.666666666666735130e-01,
  Lg2 = 3.999999999940941908e-01,
  Lg3 = 2.857142874366239149e-01,
  Lg4 = 2.222219843214978396e-01,
  Lg5 = 1.818357216161805012e-01,
  Lg6 = 1.531383769920937332e-01,
  Lg7 = 1.479819860511658591e-01;

static long nfallback = 0;

static inline float
fast_logf_exact(float x)
{
  uint32_t ix; memcpy(&ix, &x, 4);
  /* fast path only for normal, positive, finite x != 1 */
  if (ix - 0x00800000u >= 0x7f800000u - 0x00800000u || ix == 0x3f800000u) { nfallback++; return (float) log((double) x); }

  int      k  = (int) (ix >> 23) - 127;
  uint32_t mb = (ix & 0x007fffffu) | 0x3f800000u;   /* m in [1,2) */
  float    mf; memcpy(&mf, &mb, 4);
  double   m  = mf;
  if (m > 1.4142135623730951) { m *= 0.5; k += 1; } /* m in [sqrt2/2, sqrt2), exact */

  double f    = m - 1.0;                             /* exact */
  double s    = f / (2.0 + f);
  double z    = s * s;
  double w    = z * z;
  double t1   = w * (Lg2 + w * (Lg4 + w * Lg6));
  double t2   = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
  double R    = t2 + t1;
  double hfsq = 0.5 * f * f;
  double dk   = (double) k;
  double y    = dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);

  /* float rounding drops the low 29 bits of the double mantissa; it is
   * safe when those bits are far from the halfway pattern 0x10000000 */
  uint64_t uy; memcpy(&uy, &y, 8);
  uint32_t low = (uint32_t) (uy & 0x1fffffffu);
  uint32_t d   = low > 0x10000000u ? low - 0x10000000u : 0x10000000u - low;
  if (d <= 4096) { nfallback++; return (float) log((double) x); }
  return (float) y;
}

int
main(void)
{
  uint64_t n = 0, nbad = 0;
  clock_t t0 = clock();
  for (uint64_t b = 0; b <= 0xffffffffull; b++)
    {
      uint32_t ib = (uint32_t) b; float x; memcpy(&x, &ib, 4);
      float a = (float) log((double) x);
      float c = fast_logf_exact(x);
      if (memcmp(&a, &c, 4) != 0 && !(isnan(a) && isnan(c))) {
        if (nbad < 10) printf("MISMATCH x=%a ref=%a fast=%a\n", x, a, c);
        nbad++;
      }
      n++;
    }
  printf("checked %llu inputs, %llu mismatches, %ld fallbacks to libm (%.2e)  %.1fs\n",
         (unsigned long long) n, (unsigned long long) nbad, nfallback, (double) nfallback / n,
         (double) (clock() - t0) / CLOCKS_PER_SEC);

  /* speed: sum over a range of typical bias-filter row maxima */
  volatile float sink = 0;
  float acc = 0; t0 = clock();
  for (int rep = 0; rep < 20; rep++) for (uint32_t ib = 0x3c000000u; ib < 0x3f800000u; ib += 97) { float x; memcpy(&x, &ib, 4); acc += (float) log((double) x); }
  double tl = (double) (clock() - t0) / CLOCKS_PER_SEC; sink = acc; acc = 0; t0 = clock();
  for (int rep = 0; rep < 20; rep++) for (uint32_t ib = 0x3c000000u; ib < 0x3f800000u; ib += 97) { float x; memcpy(&x, &ib, 4); acc += fast_logf_exact(x); }
  double tf = (double) (clock() - t0) / CLOCKS_PER_SEC; sink += acc;
  printf("throughput loop: libm %.3fs  fast %.3fs\n", tl, tf);
  return 0;
}
