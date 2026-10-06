/* The SSV filter's DP on AVX2: 32 byte lanes instead of 16.
 *
 * This is ssvfilter.c's get_xE() with __m128i replaced by __m256i: the
 * same saturating signed-byte recurrence along each diagonal, the same
 * unsigned max, on the same scores (om->sbw: om->sbv's values striped
 * for 32 lanes, see simd.c). Read ssvfilter.c for how it works.
 *
 * Why the result equals the SSE one: the lane count only changes which
 * diagonals share a register. Every real cell (k <= M) follows the same
 * recurrence from the same begin value (-128), and padding cells (score
 * 127) saturate to -128, so until some cell overflows, both versions
 * compute the same cells and the same max. A cell can only overflow
 * from a predecessor >= 256 - bias, which both versions have already
 * counted in xE, so both report the overflow (p7_SSVFilter() then
 * doesn't use xE's value). simd_utest checks this against the SSE code.
 *
 * The one real change: a striped vector's shift by one lane has to cross
 * the two 128-bit halves (lshift_byte()).
 *
 * Compiled with -mavx2 (P7_AVX2_CFLAGS); p7_SSVFilter() calls it only
 * for profiles carrying AVX2 copies, which p7_oprofile_SetSIMD() makes
 * only on CPUs with AVX2.
 */
#include <p7_config.h>

#ifdef p7_ENABLE_AVX2

#include <immintrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"

/* 16 ymm registers, as for SSE's xmm ones (ssvfilter.c). Must stay
 * <= p7O_EXTRA_SB_AVX2 + 1. */
#define  MAX_BANDS 14

/* Shift the 32 bytes of <v> up by one lane; lane 0 gets 0. (The SSE
 * version's _mm_slli_si128(v, 1); the AVX2 byte shift works within
 * 128-bit halves, so the byte crossing them comes from a permute.) */
#define LSHIFT_BYTE(v) _mm256_alignr_epi8((v), _mm256_permute2x128_si256((v), (v), 0x08), 15)

#define STEP_SINGLE(sv)                         \
  sv   = _mm256_subs_epi8(sv, *rsc); rsc++;     \
  xEv  = _mm256_max_epu8(xEv, sv);

#define LENGTH_CHECK(label)                     \
  if (i >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1()   STEP_SINGLE(sv00)
#define STEP_BANDS_2()   STEP_BANDS_1()  STEP_SINGLE(sv01)
#define STEP_BANDS_3()   STEP_BANDS_2()  STEP_SINGLE(sv02)
#define STEP_BANDS_4()   STEP_BANDS_3()  STEP_SINGLE(sv03)
#define STEP_BANDS_5()   STEP_BANDS_4()  STEP_SINGLE(sv04)
#define STEP_BANDS_6()   STEP_BANDS_5()  STEP_SINGLE(sv05)
#define STEP_BANDS_7()   STEP_BANDS_6()  STEP_SINGLE(sv06)
#define STEP_BANDS_8()   STEP_BANDS_7()  STEP_SINGLE(sv07)
#define STEP_BANDS_9()   STEP_BANDS_8()  STEP_SINGLE(sv08)
#define STEP_BANDS_10()  STEP_BANDS_9()  STEP_SINGLE(sv09)
#define STEP_BANDS_11()  STEP_BANDS_10() STEP_SINGLE(sv10)
#define STEP_BANDS_12()  STEP_BANDS_11() STEP_SINGLE(sv11)
#define STEP_BANDS_13()  STEP_BANDS_12() STEP_SINGLE(sv12)
#define STEP_BANDS_14()  STEP_BANDS_13() STEP_SINGLE(sv13)

#define CONVERT_STEP(step, length_check, label, sv, pos)        \
  length_check(label)                                           \
  rsc = (const __m256i *) om->sbw[dsq[i]] + pos;                \
  step()                                                        \
  sv = LSHIFT_BYTE(sv);                                         \
  sv = _mm256_or_si256(sv, beginv);                             \
  i++;

#define CONVERT_1(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv00, Q - 1)
#define CONVERT_2(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv01, Q - 2)  CONVERT_1(step, LENGTH_CHECK, label)
#define CONVERT_3(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv02, Q - 3)  CONVERT_2(step, LENGTH_CHECK, label)
#define CONVERT_4(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv03, Q - 4)  CONVERT_3(step, LENGTH_CHECK, label)
#define CONVERT_5(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv04, Q - 5)  CONVERT_4(step, LENGTH_CHECK, label)
#define CONVERT_6(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv05, Q - 6)  CONVERT_5(step, LENGTH_CHECK, label)
#define CONVERT_7(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv06, Q - 7)  CONVERT_6(step, LENGTH_CHECK, label)
#define CONVERT_8(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv07, Q - 8)  CONVERT_7(step, LENGTH_CHECK, label)
#define CONVERT_9(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv08, Q - 9)  CONVERT_8(step, LENGTH_CHECK, label)
#define CONVERT_10(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv09, Q - 10) CONVERT_9(step, LENGTH_CHECK, label)
#define CONVERT_11(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv10, Q - 11) CONVERT_10(step, LENGTH_CHECK, label)
#define CONVERT_12(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv11, Q - 12) CONVERT_11(step, LENGTH_CHECK, label)
#define CONVERT_13(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv12, Q - 13) CONVERT_12(step, LENGTH_CHECK, label)
#define CONVERT_14(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv13, Q - 14) CONVERT_13(step, LENGTH_CHECK, label)

#define RESET_1()   register __m256i sv00 = beginv;
#define RESET_2()   RESET_1()  register __m256i sv01 = beginv;
#define RESET_3()   RESET_2()  register __m256i sv02 = beginv;
#define RESET_4()   RESET_3()  register __m256i sv03 = beginv;
#define RESET_5()   RESET_4()  register __m256i sv04 = beginv;
#define RESET_6()   RESET_5()  register __m256i sv05 = beginv;
#define RESET_7()   RESET_6()  register __m256i sv06 = beginv;
#define RESET_8()   RESET_7()  register __m256i sv07 = beginv;
#define RESET_9()   RESET_8()  register __m256i sv08 = beginv;
#define RESET_10()  RESET_9()  register __m256i sv09 = beginv;
#define RESET_11()  RESET_10() register __m256i sv10 = beginv;
#define RESET_12()  RESET_11() register __m256i sv11 = beginv;
#define RESET_13()  RESET_12() register __m256i sv12 = beginv;
#define RESET_14()  RESET_13() register __m256i sv13 = beginv;

/* As ssvfilter.c's CALC(), with Q = the number of 32-lane vectors */
#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQV(om->M, 32);            \
  const __m256i *rsc;                           \
                                                \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = (const __m256i *) om->sbw[dsq[i]] + i + q; \
      step()                                    \
    }                                           \
                                                \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
                                                \
 for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
   {                                            \
     for (i = 0; i < Q - w; i++)                \
       {                                        \
         rsc = (const __m256i *) om->sbw[dsq[i2 + i]] + i; \
         step()                                 \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = (const __m256i *) om->sbw[dsq[i2 + i]] + i; \
     step()                                     \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
                                                \
 return xEv;

#define BAND_FUNCTION(n)                                                                            \
static __m256i                                                                                      \
calc_band_avx2_##n(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m256i beginv, __m256i xEv) \
{                                                                                                   \
  CALC(RESET_##n, STEP_BANDS_##n, CONVERT_##n, n)                                                  \
}

BAND_FUNCTION(1)
BAND_FUNCTION(2)
BAND_FUNCTION(3)
BAND_FUNCTION(4)
BAND_FUNCTION(5)
BAND_FUNCTION(6)
BAND_FUNCTION(7)
BAND_FUNCTION(8)
BAND_FUNCTION(9)
BAND_FUNCTION(10)
BAND_FUNCTION(11)
BAND_FUNCTION(12)
BAND_FUNCTION(13)
BAND_FUNCTION(14)


/* Function:  p7_SSVFilter_xE_avx2()
 * Synopsis:  AVX2 version of the SSV filter's DP (ssvfilter.c: get_xE()).
 *
 * Purpose:   Returns the maximum diagonal score xE, in the same offset
 *            byte units as get_xE(), for sequence <dsq> of length <L>
 *            against the AVX2 copy of profile <om>'s SSV scores
 *            (<om->simd == p7_SIMD_AVX2>). p7_SSVFilter() turns it into
 *            a score.
 */
uint8_t
p7_SSVFilter_xE_avx2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m256i xEv;                       /* E state: keeps max for Mk->E as we go                     */
  __m256i beginv;                    /* begin scores                                              */
  __m128i xE16;                      /* xEv's two halves, max'ed                                  */
  int     q;                         /* counter over vectors 0..nq-1                              */
  int     Q      = p7O_NQV(om->M, 32); /* segment length: # of vectors                            */
  int     bands;                     /* the number of bands (rounds) to use                       */
  int     last_q = 0;                /* for saving the last q value to find band width            */
  int     i;                         /* counter for bands                                         */

  /* function pointers for the various number of vectors to use */
  __m256i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, __m256i, __m256i)
    = { NULL,
        calc_band_avx2_1,  calc_band_avx2_2,  calc_band_avx2_3,  calc_band_avx2_4,  calc_band_avx2_5,
        calc_band_avx2_6,  calc_band_avx2_7,  calc_band_avx2_8,  calc_band_avx2_9,  calc_band_avx2_10,
        calc_band_avx2_11, calc_band_avx2_12, calc_band_avx2_13, calc_band_avx2_14 };

  beginv = _mm256_set1_epi8(-128);
  xEv    = beginv;

  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) {
    q = (Q * (i + 1)) / bands;

    xEv = fs[q-last_q](dsq, L, om, last_q, beginv, xEv);

    last_q = q;
  }

  xE16 = _mm_max_epu8(_mm256_castsi256_si128(xEv), _mm256_extracti128_si256(xEv, 1));
  return esl_sse_hmax_epu8(xE16);
}

#endif /*p7_ENABLE_AVX2*/
