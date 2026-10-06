/* The SSV filter's DP on AVX-512BW: 64 byte lanes.
 *
 * ssvfilter_avx2.c with __m512i: the same saturating signed-byte
 * diagonal DP and unsigned max (see ssvfilter.c), on om->sbw striped
 * for 64 lanes (p7_oprofile_SetSIMD(om, p7_SIMD_AVX512)). The argument
 * for exactness is the AVX2 one (ssvfilter_avx2.c). With 32 vector
 * registers a band holds up to 28 vectors.
 *
 * Compiled with -mavx512f -mavx512bw (P7_AVX512_CFLAGS); p7_SSVFilter()
 * calls it only for profiles carrying AVX-512 copies, which
 * p7_oprofile_SetSIMD() makes only on CPUs with AVX-512BW.
 */
#include <p7_config.h>

#ifdef p7_ENABLE_AVX512

#include <immintrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"

/* Must stay <= p7O_EXTRA_SB_AVX512 + 1 */
#define  MAX_BANDS 28

/* Shift the 64 bytes of <v> up by one lane; lane 0 gets 0. The byte
 * shift works within 128-bit lanes; the bytes crossing them come from
 * <v> moved up by one 128-bit lane (valignq). */
#define LSHIFT_BYTE(v) _mm512_alignr_epi8((v), _mm512_alignr_epi64((v), _mm512_setzero_si512(), 6), 15)

#define STEP_SINGLE(sv)                         \
  sv   = _mm512_subs_epi8(sv, *rsc); rsc++;     \
  xEv  = _mm512_max_epu8(xEv, sv);

#define LENGTH_CHECK(label)                     \
  if (i >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1()   STEP_SINGLE(sv00)
#define STEP_BANDS_2()  STEP_BANDS_1() STEP_SINGLE(sv01)
#define STEP_BANDS_3()  STEP_BANDS_2() STEP_SINGLE(sv02)
#define STEP_BANDS_4()  STEP_BANDS_3() STEP_SINGLE(sv03)
#define STEP_BANDS_5()  STEP_BANDS_4() STEP_SINGLE(sv04)
#define STEP_BANDS_6()  STEP_BANDS_5() STEP_SINGLE(sv05)
#define STEP_BANDS_7()  STEP_BANDS_6() STEP_SINGLE(sv06)
#define STEP_BANDS_8()  STEP_BANDS_7() STEP_SINGLE(sv07)
#define STEP_BANDS_9()  STEP_BANDS_8() STEP_SINGLE(sv08)
#define STEP_BANDS_10()  STEP_BANDS_9() STEP_SINGLE(sv09)
#define STEP_BANDS_11()  STEP_BANDS_10() STEP_SINGLE(sv10)
#define STEP_BANDS_12()  STEP_BANDS_11() STEP_SINGLE(sv11)
#define STEP_BANDS_13()  STEP_BANDS_12() STEP_SINGLE(sv12)
#define STEP_BANDS_14()  STEP_BANDS_13() STEP_SINGLE(sv13)
#define STEP_BANDS_15()  STEP_BANDS_14() STEP_SINGLE(sv14)
#define STEP_BANDS_16()  STEP_BANDS_15() STEP_SINGLE(sv15)
#define STEP_BANDS_17()  STEP_BANDS_16() STEP_SINGLE(sv16)
#define STEP_BANDS_18()  STEP_BANDS_17() STEP_SINGLE(sv17)
#define STEP_BANDS_19()  STEP_BANDS_18() STEP_SINGLE(sv18)
#define STEP_BANDS_20()  STEP_BANDS_19() STEP_SINGLE(sv19)
#define STEP_BANDS_21()  STEP_BANDS_20() STEP_SINGLE(sv20)
#define STEP_BANDS_22()  STEP_BANDS_21() STEP_SINGLE(sv21)
#define STEP_BANDS_23()  STEP_BANDS_22() STEP_SINGLE(sv22)
#define STEP_BANDS_24()  STEP_BANDS_23() STEP_SINGLE(sv23)
#define STEP_BANDS_25()  STEP_BANDS_24() STEP_SINGLE(sv24)
#define STEP_BANDS_26()  STEP_BANDS_25() STEP_SINGLE(sv25)
#define STEP_BANDS_27()  STEP_BANDS_26() STEP_SINGLE(sv26)
#define STEP_BANDS_28()  STEP_BANDS_27() STEP_SINGLE(sv27)

#define CONVERT_STEP(step, length_check, label, sv, pos)        \
  length_check(label)                                           \
  rsc = (const __m512i *) om->sbw[dsq[i]] + pos;                \
  step()                                                        \
  sv = LSHIFT_BYTE(sv);                                         \
  sv = _mm512_or_si512(sv, beginv);                             \
  i++;

#define CONVERT_1(step, LENGTH_CHECK, label)  CONVERT_STEP(step, LENGTH_CHECK, label, sv00, Q - 1)
#define CONVERT_2(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv01, Q - 2) CONVERT_1(step, LENGTH_CHECK, label)
#define CONVERT_3(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv02, Q - 3) CONVERT_2(step, LENGTH_CHECK, label)
#define CONVERT_4(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv03, Q - 4) CONVERT_3(step, LENGTH_CHECK, label)
#define CONVERT_5(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv04, Q - 5) CONVERT_4(step, LENGTH_CHECK, label)
#define CONVERT_6(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv05, Q - 6) CONVERT_5(step, LENGTH_CHECK, label)
#define CONVERT_7(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv06, Q - 7) CONVERT_6(step, LENGTH_CHECK, label)
#define CONVERT_8(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv07, Q - 8) CONVERT_7(step, LENGTH_CHECK, label)
#define CONVERT_9(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv08, Q - 9) CONVERT_8(step, LENGTH_CHECK, label)
#define CONVERT_10(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv09, Q - 10) CONVERT_9(step, LENGTH_CHECK, label)
#define CONVERT_11(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv10, Q - 11) CONVERT_10(step, LENGTH_CHECK, label)
#define CONVERT_12(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv11, Q - 12) CONVERT_11(step, LENGTH_CHECK, label)
#define CONVERT_13(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv12, Q - 13) CONVERT_12(step, LENGTH_CHECK, label)
#define CONVERT_14(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv13, Q - 14) CONVERT_13(step, LENGTH_CHECK, label)
#define CONVERT_15(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv14, Q - 15) CONVERT_14(step, LENGTH_CHECK, label)
#define CONVERT_16(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv15, Q - 16) CONVERT_15(step, LENGTH_CHECK, label)
#define CONVERT_17(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv16, Q - 17) CONVERT_16(step, LENGTH_CHECK, label)
#define CONVERT_18(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv17, Q - 18) CONVERT_17(step, LENGTH_CHECK, label)
#define CONVERT_19(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv18, Q - 19) CONVERT_18(step, LENGTH_CHECK, label)
#define CONVERT_20(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv19, Q - 20) CONVERT_19(step, LENGTH_CHECK, label)
#define CONVERT_21(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv20, Q - 21) CONVERT_20(step, LENGTH_CHECK, label)
#define CONVERT_22(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv21, Q - 22) CONVERT_21(step, LENGTH_CHECK, label)
#define CONVERT_23(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv22, Q - 23) CONVERT_22(step, LENGTH_CHECK, label)
#define CONVERT_24(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv23, Q - 24) CONVERT_23(step, LENGTH_CHECK, label)
#define CONVERT_25(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv24, Q - 25) CONVERT_24(step, LENGTH_CHECK, label)
#define CONVERT_26(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv25, Q - 26) CONVERT_25(step, LENGTH_CHECK, label)
#define CONVERT_27(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv26, Q - 27) CONVERT_26(step, LENGTH_CHECK, label)
#define CONVERT_28(step, LENGTH_CHECK, label) CONVERT_STEP(step, LENGTH_CHECK, label, sv27, Q - 28) CONVERT_27(step, LENGTH_CHECK, label)

#define RESET_1()   register __m512i sv00 = beginv;
#define RESET_2()  RESET_1() register __m512i sv01 = beginv;
#define RESET_3()  RESET_2() register __m512i sv02 = beginv;
#define RESET_4()  RESET_3() register __m512i sv03 = beginv;
#define RESET_5()  RESET_4() register __m512i sv04 = beginv;
#define RESET_6()  RESET_5() register __m512i sv05 = beginv;
#define RESET_7()  RESET_6() register __m512i sv06 = beginv;
#define RESET_8()  RESET_7() register __m512i sv07 = beginv;
#define RESET_9()  RESET_8() register __m512i sv08 = beginv;
#define RESET_10()  RESET_9() register __m512i sv09 = beginv;
#define RESET_11()  RESET_10() register __m512i sv10 = beginv;
#define RESET_12()  RESET_11() register __m512i sv11 = beginv;
#define RESET_13()  RESET_12() register __m512i sv12 = beginv;
#define RESET_14()  RESET_13() register __m512i sv13 = beginv;
#define RESET_15()  RESET_14() register __m512i sv14 = beginv;
#define RESET_16()  RESET_15() register __m512i sv15 = beginv;
#define RESET_17()  RESET_16() register __m512i sv16 = beginv;
#define RESET_18()  RESET_17() register __m512i sv17 = beginv;
#define RESET_19()  RESET_18() register __m512i sv18 = beginv;
#define RESET_20()  RESET_19() register __m512i sv19 = beginv;
#define RESET_21()  RESET_20() register __m512i sv20 = beginv;
#define RESET_22()  RESET_21() register __m512i sv21 = beginv;
#define RESET_23()  RESET_22() register __m512i sv22 = beginv;
#define RESET_24()  RESET_23() register __m512i sv23 = beginv;
#define RESET_25()  RESET_24() register __m512i sv24 = beginv;
#define RESET_26()  RESET_25() register __m512i sv25 = beginv;
#define RESET_27()  RESET_26() register __m512i sv26 = beginv;
#define RESET_28()  RESET_27() register __m512i sv27 = beginv;

/* As ssvfilter.c's CALC(), with Q = the number of 64-lane vectors */
#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQV(om->M, 64);            \
  const __m512i *rsc;                           \
                                                \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = (const __m512i *) om->sbw[dsq[i]] + i + q; \
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
         rsc = (const __m512i *) om->sbw[dsq[i2 + i]] + i; \
         step()                                 \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = (const __m512i *) om->sbw[dsq[i2 + i]] + i; \
     step()                                     \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
                                                \
 return xEv;

#define BAND_FUNCTION(n)                                                                              \
static __m512i                                                                                        \
calc_band_avx512_##n(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i xEv) \
{                                                                                                     \
  CALC(RESET_##n, STEP_BANDS_##n, CONVERT_##n, n)                                                    \
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
BAND_FUNCTION(15)
BAND_FUNCTION(16)
BAND_FUNCTION(17)
BAND_FUNCTION(18)
BAND_FUNCTION(19)
BAND_FUNCTION(20)
BAND_FUNCTION(21)
BAND_FUNCTION(22)
BAND_FUNCTION(23)
BAND_FUNCTION(24)
BAND_FUNCTION(25)
BAND_FUNCTION(26)
BAND_FUNCTION(27)
BAND_FUNCTION(28)


/* Function:  p7_SSVFilter_xE_avx512()
 * Synopsis:  AVX-512BW version of the SSV filter's DP (ssvfilter.c: get_xE()).
 *
 * Purpose:   As p7_SSVFilter_xE_avx2(), for a profile with AVX-512
 *            copies (<om->simd == p7_SIMD_AVX512>).
 */
uint8_t
p7_SSVFilter_xE_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m512i xEv;                       /* E state: keeps max for Mk->E as we go                     */
  __m512i beginv;                    /* begin scores                                              */
  __m256i xE32;                      /* xEv's two 256-bit halves, max'ed                          */
  __m128i xE16;                      /* ... and down to 128 bits                                  */
  int     q;                         /* counter over vectors 0..nq-1                              */
  int     Q      = p7O_NQV(om->M, 64); /* segment length: # of vectors                            */
  int     bands;                     /* the number of bands (rounds) to use                       */
  int     last_q = 0;                /* for saving the last q value to find band width            */
  int     i;                         /* counter for bands                                         */

  /* function pointers for the various number of vectors to use */
  __m512i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, __m512i, __m512i)
    = { NULL,
        calc_band_avx512_1, calc_band_avx512_2, calc_band_avx512_3, calc_band_avx512_4, calc_band_avx512_5,
        calc_band_avx512_6, calc_band_avx512_7, calc_band_avx512_8, calc_band_avx512_9, calc_band_avx512_10,
        calc_band_avx512_11, calc_band_avx512_12, calc_band_avx512_13, calc_band_avx512_14, calc_band_avx512_15,
        calc_band_avx512_16, calc_band_avx512_17, calc_band_avx512_18, calc_band_avx512_19, calc_band_avx512_20,
        calc_band_avx512_21, calc_band_avx512_22, calc_band_avx512_23, calc_band_avx512_24, calc_band_avx512_25,
        calc_band_avx512_26, calc_band_avx512_27, calc_band_avx512_28 };

  beginv = _mm512_set1_epi8(-128);
  xEv    = beginv;

  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) {
    q = (Q * (i + 1)) / bands;

    xEv = fs[q-last_q](dsq, L, om, last_q, beginv, xEv);

    last_q = q;
  }

  xE32 = _mm256_max_epu8(_mm512_castsi512_si256(xEv), _mm512_extracti64x4_epi64(xEv, 1));
  xE16 = _mm_max_epu8(_mm256_castsi256_si128(xE32), _mm256_extracti128_si256(xE32, 1));
  return esl_sse_hmax_epu8(xE16);
}

#endif /*p7_ENABLE_AVX512*/
