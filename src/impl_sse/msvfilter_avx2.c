/* The MSV filter's DP on AVX2: 32 byte lanes instead of 16.
 *
 * This is the full MSV loop of msvfilter.c (the fallback after the SSV
 * filter can't decide) with __m128i replaced by __m256i, on the same
 * scores (om->rbw: om->rbv's values striped for 32 lanes, see simd.c).
 *
 * Why the result equals the SSE one: each cell M(i,k) is the same
 * function of M(i-1,k-1), xB(i-1) and the score of residue i at node k,
 * whatever lane holds it. Padding nodes (cost 255) are 0, the smallest
 * value, so they never raise a row's max nor trigger its overflow test.
 * Like the SSE code, each row's max goes to every lane before the
 * overflow test and the special states (E, J, B), which are then the
 * same saturated byte operations on the same value. simd_utest checks
 * this against the SSE code.
 *
 * The row's max is the critical path from one row to the next (row i
 * needs xB(i-1)), so it stays in vector registers: a version that took
 * it to a scalar and broadcast xB back was slower than SSE.
 *
 * Compiled with -mavx2 (P7_AVX2_CFLAGS); p7_MSVFilter() calls it only
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

/* Function:  p7_MSVFilter_avx2()
 * Synopsis:  AVX2 version of p7_MSVFilter()'s full MSV DP.
 *
 * Purpose:   The DP part of <p7_MSVFilter()> after <p7_SSVFilter()>
 *            couldn't decide, for sequence <dsq> of length <L> against
 *            the AVX2 copy of <om>'s MSV scores. Uses <ox>'s AVX2 row.
 *            Returns <eslOK> and the final J state score <*ret_xJ>, in
 *            the same offset byte units as the SSE version; or
 *            <eslERANGE> if a row overflows.
 *
 * Throws:    <eslEMEM> if <ox>'s AVX2 row can't be allocated.
 */
int
p7_MSVFilter_avx2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, uint8_t *ret_xJ)
{
  __m256i  mpv;                   /* previous row values                                       */
  __m256i  xEv;                   /* E state: keeps max for Mk->E as we go                     */
  __m256i  xBv;                   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m256i  sv;                    /* temp storage of 1 curr row value in progress              */
  __m256i  biasv;                 /* emission bias in a vector                                 */
  __m256i *dp;                    /* the one DP row, dp[0..Q-1]                                */
  const __m256i *rsc;             /* will point at om->rbw[x] for residue x[i]                 */
  __m256i  xJv;                   /* vector for states score                                   */
  __m256i  tjbmv;                 /* vector for cost of moving from either J or N through B to an M state */
  __m256i  tecv;                  /* vector for E->C  cost                                     */
  __m256i  basev;                 /* offset for scores                                         */
  __m256i  ceilingv;              /* saturated simd value used to test for overflow            */
  __m256i  rot2v, rot1v;          /* vpshufb masks: swap the words of each dword, the bytes of each word */
  int      i;                     /* counter over sequence positions 1..L                      */
  int      q;                     /* counter over vectors 0..nq-1                              */
  int      Q = p7O_NQV(om->M, 32);/* segment length: # of vectors                              */

  if ((dp = p7_omx_SIMDRow(ox, sizeof(__m256i) * Q)) == NULL) return eslEMEM;

  /* Initialization, as msvfilter.c. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base. */
  biasv = _mm256_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q] = _mm256_setzero_si256();

  ceilingv = _mm256_cmpeq_epi8(biasv, biasv);
  basev    = _mm256_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm256_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm256_set1_epi8((int8_t) om->tec_b);
  xJv      = _mm256_subs_epu8(biasv, biasv);
  xBv      = _mm256_subs_epu8(basev, tjbmv);
  rot2v    = _mm256_setr_epi8(2,3,0,1, 6,7,4,5, 10,11,8,9, 14,15,12,13, 2,3,0,1, 6,7,4,5, 10,11,8,9, 14,15,12,13);
  rot1v    = _mm256_setr_epi8(1,0,3,2, 5,4,7,6, 9,8,11,10, 13,12,15,14, 1,0,3,2, 5,4,7,6, 9,8,11,10, 13,12,15,14);

  for (i = 1; i <= L; i++)
    {
      rsc = (const __m256i *) om->rbw[dsq[i]];
      xEv = _mm256_setzero_si256();

      /* Shift the previous row's last vector up by one lane, across the two
       * 128-bit halves; zero (-infinity) comes in at lane 0. */
      mpv = _mm256_alignr_epi8(dp[Q-1], _mm256_permute2x128_si256(dp[Q-1], dp[Q-1], 0x08), 15);
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv   = _mm256_max_epu8(mpv, xBv);
	  sv   = _mm256_adds_epu8(sv, biasv);
	  sv   = _mm256_subs_epu8(sv, *rsc);   rsc++;
	  xEv  = _mm256_max_epu8(xEv, sv);

	  mpv   = dp[q];         /* Load {MDI}(i-1,q) into mpv */
	  dp[q] = sv;            /* Do delayed store of M(i,q) now that memory is usable */
	}

      /* The row's max, in every lane: the two halves, then dwords,
       * then words and bytes within dwords. */
      xEv = _mm256_max_epu8(xEv, _mm256_permute2x128_si256(xEv, xEv, 0x01));
      xEv = _mm256_max_epu8(xEv, _mm256_shuffle_epi32(xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      xEv = _mm256_max_epu8(xEv, _mm256_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1)));
      xEv = _mm256_max_epu8(xEv, _mm256_shuffle_epi8(xEv, rot2v));
      xEv = _mm256_max_epu8(xEv, _mm256_shuffle_epi8(xEv, rot1v));

      /* Overflow: msvfilter.c tests whether some lane +sat bias is 255,
       * before taking the max; on the max, it's the same test. */
      if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_adds_epu8(xEv, biasv), ceilingv)) != 0) return eslERANGE;

      xEv = _mm256_subs_epu8(xEv, tecv);
      xJv = _mm256_max_epu8(xJv,xEv);

      xBv = _mm256_max_epu8(basev, xJv);
      xBv = _mm256_subs_epu8(xBv, tjbmv);
    }

  *ret_xJ = (uint8_t) _mm_extract_epi16(_mm256_castsi256_si128(xJv), 0);
  return eslOK;
}

#endif /*p7_ENABLE_AVX2*/
