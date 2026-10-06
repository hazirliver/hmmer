/* The Viterbi filter's DP on AVX2: 16 word lanes instead of 8.
 *
 * This is p7_ViterbiFilter() of vitfilter.c with __m128i replaced by
 * __m256i, on the same scores (om->rww, om->tww: om->rwv's and om->twv's
 * values striped for 16 lanes, see simd.c). The special states use the
 * same C statements on the same int16_t variables.
 *
 * Why the result equals the SSE one: M and I cells are the same function
 * of the previous row's cells whatever lane holds them. A row's D cells
 * end up either exact (the lazy-F loop runs to its fixed point; its early
 * exit is exact) or, when the DD bound proves D->D paths can't matter,
 * the M->D values alone; which of the two depends on the row's largest
 * M->D value, the same in both versions. Either way the next row's M
 * cells are the same (that is the bound's guarantee). Padding nodes are
 * -32768 in M, I and their M->D values, so they change neither the row
 * maxima nor any real cell. simd_utest checks this against the SSE code.
 *
 * Compiled with -mavx2 (P7_AVX2_CFLAGS); p7_ViterbiFilter() calls it only
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

/* Shift the 16 words of <v> up by one lane; lane 0 gets 0. (The SSE
 * version's _mm_slli_si128(v, 2), across the two 128-bit halves.) */
#define LSHIFT_WORD(v) _mm256_alignr_epi8((v), _mm256_permute2x128_si256((v), (v), 0x08), 14)

static inline int16_t
hmax_epi16_avx2(__m256i a)
{
  return esl_sse_hmax_epi16(_mm_max_epi16(_mm256_castsi256_si128(a), _mm256_extracti128_si256(a, 1)));
}

static inline int
any_gt_epi16_avx2(__m256i a, __m256i b)
{
  return (_mm256_movemask_epi8(_mm256_cmpgt_epi16(a, b)) != 0);
}

/* Function:  p7_ViterbiFilter_avx2()
 * Synopsis:  AVX2 version of p7_ViterbiFilter()'s DP.
 *
 * Purpose:   The DP of <p7_ViterbiFilter()> for sequence <dsq> of length
 *            <L> against the AVX2 copy of <om>'s Viterbi filter scores,
 *            using <ox>'s AVX2 row. Returns <eslOK> and the final C state
 *            score <*ret_xC> (scaled words, as the SSE version's xC); or
 *            <eslERANGE> if the score overflows. <p7_ViterbiFilter()>
 *            checks the profile's mode and turns xC into a score.
 *
 * Throws:    <eslEMEM> if <ox>'s AVX2 row can't be allocated.
 */
int
p7_ViterbiFilter_avx2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, int16_t *ret_xC)
{
  __m256i  mpv, dpv, ipv;          /* previous row values                                       */
  __m256i  sv;                     /* temp storage of 1 curr row value in progress              */
  __m256i  dcv;                    /* delayed storage of D(i,q+1)                               */
  __m256i  xEv;                    /* E state: keeps max for Mk->E as we go                     */
  __m256i  xBv;                    /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m256i  Dmaxv;                  /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;     /* special states' scores                                    */
  int16_t  Dmax;                   /* maximum D cell score on row                               */
  int      i;                      /* counter over sequence positions 1..L                      */
  int      q;                      /* counter over vectors 0..nq-1                              */
  int      Q = p7O_NQV(om->M, 16); /* segment length: # of vectors                              */
  __m256i *dp;                     /* the one DP row, for the {MDI}MXo(q) macros               */
  const __m256i *rsc;              /* will point at om->rww[x] for residue x[i]                 */
  const __m256i *tsc;              /* will point into (and step thru) om->tww                   */
  __m256i  negInfv;

  if ((dp = p7_omx_SIMDRow(ox, sizeof(__m256i) * p7X_NSCELLS * Q)) == NULL) return eslEMEM;

  /* -infinity is -32768; negInfv is -32768 in lane 0 and 0 elsewhere, for an OR */
  negInfv = _mm256_setr_epi16(-32768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  /* Initialization. In unsigned arithmetic, -infinity is -32768 */
  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = _mm256_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

  for (i = 1; i <= L; i++)
    {
      rsc   = (const __m256i *) om->rww[dsq[i]];
      tsc   = (const __m256i *) om->tww;
      dcv   = _mm256_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm256_set1_epi16(-32768);
      Dmaxv = _mm256_set1_epi16(-32768);
      xBv   = _mm256_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes); replace the 0 shifted in with -32768. */
      mpv = MMXo(Q-1);  mpv = LSHIFT_WORD(mpv);  mpv = _mm256_or_si256(mpv, negInfv);
      dpv = DMXo(Q-1);  dpv = LSHIFT_WORD(dpv);  dpv = _mm256_or_si256(dpv, negInfv);
      ipv = IMXo(Q-1);  ipv = LSHIFT_WORD(ipv);  ipv = _mm256_or_si256(ipv, negInfv);

      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv   =                       _mm256_adds_epi16(xBv, *tsc);  tsc++;
	  sv   = _mm256_max_epi16 (sv, _mm256_adds_epi16(mpv, *tsc)); tsc++;
	  sv   = _mm256_max_epi16 (sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
	  sv   = _mm256_max_epi16 (sv, _mm256_adds_epi16(dpv, *tsc)); tsc++;
	  sv   = _mm256_adds_epi16(sv, *rsc);                         rsc++;
	  xEv  = _mm256_max_epi16(xEv, sv);

	  /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
	   * {MDI}MX(q) is then the current, not the prev row
	   */
	  mpv = MMXo(q);
	  dpv = DMXo(q);
	  ipv = IMXo(q);

	  /* Do the delayed stores of {MD}(i,q) now that memory is usable */
	  MMXo(q) = sv;
	  DMXo(q) = dcv;

	  /* Calculate the next D(i,q+1) partially: M->D only;
	   * delay storage, holding it in dcv
	   */
	  dcv   = _mm256_adds_epi16(sv, *tsc);  tsc++;
	  Dmaxv = _mm256_max_epi16(dcv, Dmaxv);

	  /* Calculate and store I(i,q) */
	  sv      =                       _mm256_adds_epi16(mpv, *tsc);  tsc++;
	  IMXo(q) = _mm256_max_epi16 (sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
	}

      /* Now the "special" states, which start from Mk->E (->C, ->J->B): as vitfilter.c */
      xE = hmax_epi16_avx2(xEv);
      if (xE >= 32767) return eslERANGE;	/* immediately detect overflow */
      xN = xN + om->xw[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* The "lazy F" loop, as vitfilter.c (see there for the DD bound) */
      Dmax = hmax_epi16_avx2(Dmaxv);
      if (Dmax + om->ddbound_w > xB)
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv = LSHIFT_WORD(dcv);
	  dcv = _mm256_or_si256(dcv, negInfv);
	  tsc = (const __m256i *) om->tww + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++)
	    {
	      DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
	      dcv     = _mm256_adds_epi16(DMXo(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score.
	   */
	  do {
	    dcv = LSHIFT_WORD(dcv);
	    dcv = _mm256_or_si256(dcv, negInfv);
	    tsc = (const __m256i *) om->tww + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++)
	      {
		if (! any_gt_epi16_avx2(dcv, DMXo(q))) break;
		DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
		dcv     = _mm256_adds_epi16(DMXo(q), *tsc);   tsc++;
	      }
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
	  dcv = LSHIFT_WORD(dcv);
	  DMXo(0) = _mm256_or_si256(dcv, negInfv);
	}
    } /* end loop over sequence residues 1..L */

  *ret_xC = xC;
  return eslOK;
}

#endif /*p7_ENABLE_AVX2*/
