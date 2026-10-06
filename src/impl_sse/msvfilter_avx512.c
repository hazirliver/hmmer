/* The MSV filter's DP on AVX-512BW: 64 byte lanes.
 *
 * msvfilter_avx2.c with __m512i, on om->rbw striped for 64 lanes
 * (p7_oprofile_SetSIMD(om, p7_SIMD_AVX512)); see there and msvfilter.c.
 * The argument for exactness is the AVX2 one.
 *
 * Compiled with -mavx512f -mavx512bw (P7_AVX512_CFLAGS); p7_MSVFilter()
 * calls it only for profiles carrying AVX-512 copies, which
 * p7_oprofile_SetSIMD() makes only on CPUs with AVX-512BW.
 */
#include <p7_config.h>

#ifdef p7_ENABLE_AVX512

#include <immintrin.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_sse.h"

/* Function:  p7_MSVFilter_avx512()
 * Synopsis:  AVX-512BW version of p7_MSVFilter()'s full MSV DP.
 *
 * Purpose:   As p7_MSVFilter_avx2(), for a profile with AVX-512 copies
 *            (<om->simd == p7_SIMD_AVX512>).
 */
int
p7_MSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, uint8_t *ret_xJ)
{
  __m512i  mpv;                   /* previous row values                                       */
  __m512i  xEv;                   /* E state: keeps max for Mk->E as we go                     */
  __m512i  xBv;                   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m512i  sv;                    /* temp storage of 1 curr row value in progress              */
  __m512i  biasv;                 /* emission bias in a vector                                 */
  __m512i *dp;                    /* the one DP row, dp[0..Q-1]                                */
  const __m512i *rsc;             /* will point at om->rbw[x] for residue x[i]                 */
  __m512i  xJv;                   /* vector for states score                                   */
  __m512i  tjbmv;                 /* vector for cost of moving from either J or N through B to an M state */
  __m512i  tecv;                  /* vector for E->C  cost                                     */
  __m512i  basev;                 /* offset for scores                                         */
  __m512i  ceilingv;              /* saturated simd value used to test for overflow            */
  __m512i  rot2v, rot1v;          /* vpshufb masks: swap the words of each dword, the bytes of each word */
  int      i;                     /* counter over sequence positions 1..L                      */
  int      q;                     /* counter over vectors 0..nq-1                              */
  int      Q = p7O_NQV(om->M, 64);/* segment length: # of vectors                              */

  if ((dp = p7_omx_SIMDRow(ox, sizeof(__m512i) * Q)) == NULL) return eslEMEM;

  /* Initialization, as msvfilter.c. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base. */
  biasv = _mm512_set1_epi8((int8_t) om->bias_b);
  for (q = 0; q < Q; q++) dp[q] = _mm512_setzero_si512();

  ceilingv = _mm512_set1_epi8(-1);
  basev    = _mm512_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm512_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm512_set1_epi8((int8_t) om->tec_b);
  xJv      = _mm512_subs_epu8(biasv, biasv);
  xBv      = _mm512_subs_epu8(basev, tjbmv);
  rot2v    = _mm512_broadcast_i32x4(_mm_setr_epi8(2,3,0,1, 6,7,4,5, 10,11,8,9, 14,15,12,13));
  rot1v    = _mm512_broadcast_i32x4(_mm_setr_epi8(1,0,3,2, 5,4,7,6, 9,8,11,10, 13,12,15,14));

  for (i = 1; i <= L; i++)
    {
      rsc = (const __m512i *) om->rbw[dsq[i]];
      xEv = _mm512_setzero_si512();

      /* Shift the previous row's last vector up by one lane, across the
       * 128-bit lanes; zero (-infinity) comes in at lane 0. */
      mpv = _mm512_alignr_epi8(dp[Q-1], _mm512_alignr_epi64(dp[Q-1], _mm512_setzero_si512(), 6), 15);
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv   = _mm512_max_epu8(mpv, xBv);
	  sv   = _mm512_adds_epu8(sv, biasv);
	  sv   = _mm512_subs_epu8(sv, *rsc);   rsc++;
	  xEv  = _mm512_max_epu8(xEv, sv);

	  mpv   = dp[q];         /* Load {MDI}(i-1,q) into mpv */
	  dp[q] = sv;            /* Do delayed store of M(i,q) now that memory is usable */
	}

      /* The row's max, in every lane: 256-bit halves, 128-bit lanes,
       * then dwords, words and bytes within each lane. */
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_i64x2(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_i64x2(xEv, xEv, _MM_SHUFFLE(2, 3, 0, 1)));
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_epi32(xEv, (_MM_PERM_ENUM) _MM_SHUFFLE(1, 0, 3, 2)));
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_epi32(xEv, (_MM_PERM_ENUM) _MM_SHUFFLE(2, 3, 0, 1)));
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_epi8(xEv, rot2v));
      xEv = _mm512_max_epu8(xEv, _mm512_shuffle_epi8(xEv, rot1v));

      /* Overflow: some lane +sat bias is 255 (on the max, the same test) */
      if (_mm512_cmpeq_epi8_mask(_mm512_adds_epu8(xEv, biasv), ceilingv) != 0) return eslERANGE;

      xEv = _mm512_subs_epu8(xEv, tecv);
      xJv = _mm512_max_epu8(xJv,xEv);

      xBv = _mm512_max_epu8(basev, xJv);
      xBv = _mm512_subs_epu8(xBv, tjbmv);
    }

  *ret_xJ = (uint8_t) _mm_extract_epi16(_mm512_castsi512_si128(xJv), 0);
  return eslOK;
}

#endif /*p7_ENABLE_AVX512*/
