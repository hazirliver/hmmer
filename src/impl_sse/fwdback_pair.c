/* The Forward and Backward parsers on two target sequences at once.
 *
 * The parsers (fwdback.c, do_full = FALSE) are latency-bound: most of a
 * row's time goes to the serial D->D passes, one dependent add or
 * multiply per vector of the row. Two targets have independent chains,
 * so here one 256-bit AVX2 register holds both: target A in its low
 * 128-bit half, target B in its high half. Each half runs exactly the
 * SSE code's operation sequence for its own target, so each target's
 * results are bitwise those of p7_ForwardParser()/p7_BackwardParser():
 *
 *  - element-wise ops (_mm256_add/mul_ps) act on each half as the SSE
 *    op does on its vector; no FMA (compiled with -mavx2 only), and the
 *    operands and nesting are the SSE code's;
 *  - the SSE code's vector shifts and horizontal sums use in-lane AVX2
 *    shifts and shuffles (_mm256_slli/srli_si256 and _mm256_shuffle_ps
 *    work on each 128-bit half separately), so the halves never mix;
 *  - the special states (E N J B C), rescaling and the stores to each
 *    target's own P7_OMX are the SSE code's scalar statements, per target,
 *    with that target's length-dependent transitions <xf>; a row is
 *    rescaled by blending the scaled values into the rescaled target's
 *    half only;
 *  - the Forward parser's early end of the D->D passes (M >= 100): a
 *    target whose own pass would have stopped gets no more D->D: its half
 *    of the carried <dcv> is zeroed, so its D values stay unchanged
 *    (x + 0 = x for the non-negative values here) while the other
 *    target's passes go on;
 *  - targets of different lengths: the shorter one's half goes on with
 *    zero emissions and is never read again (Backward: its half is
 *    terminated as soon as it has done row 1).
 *
 * Only the special-state rows (xmx), M, L, has_own_scales and totscale of
 * the targets' P7_OMX are written, which is all that the Backward parser,
 * the Forward filter test and domain decoding read from them; the DP row
 * lives in a buffer of the first target's P7_OMX (p7_omx_SIMDRow()).
 *
 * Compiled with -mavx2 (P7_AVX2_CFLAGS). Callers check
 * p7_fb_PairAvailable() first.
 */
#include <p7_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "hmmer.h"
#include "impl_sse.h"

#ifndef p7FWDBACK_PAIR_TESTDRIVE

/* Function:  p7_fb_PairAvailable()
 * Synopsis:  TRUE if the paired parsers can run on this host.
 */
int
p7_fb_PairAvailable(void)
{
#ifdef p7_ENABLE_AVX2
  return (p7_simd_Level() >= p7_SIMD_AVX2);
#else
  return FALSE;
#endif
}

#if defined(p7_ENABLE_AVX2) && defined(__AVX2__)

#include <math.h>
#include <immintrin.h>

/* both targets' {MDI}(q) in one row of __m256 */
#define PM(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_M])
#define PD(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_D])
#define PI(dp,q) ((dp)[(q) * p7X_NSCELLS + p7X_I])

/* a's vector in the low half, b's in the high one */
static inline __m256 two(const __m128 *a, const __m128 *b) { return _mm256_insertf128_ps(_mm256_castps128_ps256(*a), *b, 1); }
/* the same vector (a transition score) in both halves */
static inline __m256 both(const __m128 *t)                 { return _mm256_broadcast_ps(t); }
static inline __m256 halves(float a, float b)              { return _mm256_setr_m128(_mm_set1_ps(a), _mm_set1_ps(b)); }
/* per half, esl_sse_rightshiftz_float(): [x0 x1 x2 x3] -> [0 x0 x1 x2] */
static inline __m256 rshiftz(__m256 a)                     { return _mm256_castsi256_ps(_mm256_slli_si256(_mm256_castps_si256(a), 4)); }
/* per half, _mm_move_ss(a, 0) then _mm_shuffle_ps(_MM_SHUFFLE(0,3,2,1)): [x0 x1 x2 x3] -> [x1 x2 x3 0] */
static inline __m256 lshiftz(__m256 a)                     { return _mm256_castsi256_ps(_mm256_srli_si256(_mm256_castps_si256(a), 4)); }
/* all-ones lanes in the halves whose flag is set */
static inline __m256 lanes(int a, int b)                   { return _mm256_castsi256_ps(_mm256_setr_epi32(-a,-a,-a,-a, -b,-b,-b,-b)); }

/* per half, the SSE code's horizontal sum: element 0 of each half */
static inline void
hsum(__m256 v, float *ret_a, float *ret_b)
{
  v = _mm256_add_ps(v, _mm256_shuffle_ps(v, v, _MM_SHUFFLE(0, 3, 2, 1)));
  v = _mm256_add_ps(v, _mm256_shuffle_ps(v, v, _MM_SHUFFLE(1, 0, 3, 2)));
  *ret_a = _mm_cvtss_f32(_mm256_castps256_ps128(v));
  *ret_b = _mm_cvtss_f32(_mm256_extractf128_ps(v, 1));
}

/* The pair's DP row (in place, as the parsers' one row) and Q zero vectors
 * (emissions for a target that's past its end), from <ox>'s SIMD row buffer.
 */
static __m256 *
pair_row(P7_OMX *ox, int Q, __m128 **ret_zero)
{
  __m256 *dp;
  int     q;

  if ((dp = (__m256 *) p7_omx_SIMDRow(ox, sizeof(__m256) * p7X_NSCELLS * Q + sizeof(__m128) * Q)) == NULL) return NULL;
  *ret_zero = (__m128 *) (dp + p7X_NSCELLS * Q);
  for (q = 0; q < Q; q++) (*ret_zero)[q] = _mm_setzero_ps();
  return dp;
}


/* Function:  p7_ForwardParserPair()
 * Synopsis:  p7_ForwardParser() on two targets at once.
 *
 * Purpose:   The Forward parser for target A (<dsqA>, <LA>) into <oxA>, and
 *            for target B into <oxB>, with profile <om>: the same values,
 *            bit for bit, as p7_ForwardParser() on each. <xfA> and <xfB>
 *            are <om->xf> as configured for each target's length
 *            (p7_oprofile_ReconfigLength()); everything else in <om> must
 *            be the same for both. Lengths are >= 1. <oxA> and <oxB> must
 *            hold <LA>, <LB> rows of special states (p7_omx_GrowTo()).
 *            Scores (nats) go to <*opt_scA>, <*opt_scB>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslERANGE> on a score overflow, underflow or NaN, as
 *            p7_ForwardParser() does (A is checked first); <eslEMEM> if
 *            the row buffer can't be allocated.
 */
int
p7_ForwardParserPair(const P7_OPROFILE *om,
		     const ESL_DSQ *dsqA, int LA, const float (*xfA)[p7O_NXTRANS], P7_OMX *oxA, float *opt_scA,
		     const ESL_DSQ *dsqB, int LB, const float (*xfB)[p7O_NXTRANS], P7_OMX *oxB, float *opt_scB)
{
  const ESL_DSQ *dsq[2] = { dsqA, dsqB };
  const float  (*xf[2])[p7O_NXTRANS] = { xfA, xfB };
  P7_OMX        *ox[2]  = { oxA, oxB };
  float         *sc[2]  = { opt_scA, opt_scB };
  int            L[2]   = { LA, LB };
  int            Lmax   = ESL_MAX(LA, LB);
  float          xN[2], xE[2], xB[2], xC[2], xJ[2], e[2];
  float          scl[2];
  int            act[2], live[2], resc[2];
  __m256         mpv, dpv, ipv, sv, dcv, xEv, xBv, cv, zerov;
  __m256        *dp;
  __m128        *zrow;
  const __m128  *rpa, *rpb;
  const __m128  *tp;
  int            Q = p7O_NQF(om->M);
  int            i, q, j, h, m;

  if ((dp = pair_row(oxA, Q, &zrow)) == NULL) ESL_EXCEPTION(eslEMEM, "allocation failed");
  zerov = _mm256_setzero_ps();
  for (q = 0; q < Q; q++) PM(dp,q) = PI(dp,q) = PD(dp,q) = zerov;

  for (h = 0; h < 2; h++)
    {
      ox[h]->M  = om->M;
      ox[h]->L  = L[h];
      ox[h]->has_own_scales = TRUE; 	/* all forward matrices control their own scalefactors */
      xE[h] = ox[h]->xmx[p7X_E] = 0.;
      xN[h] = ox[h]->xmx[p7X_N] = 1.;
      xJ[h] = ox[h]->xmx[p7X_J] = 0.;
      xB[h] = ox[h]->xmx[p7X_B] = xf[h][p7O_N][p7O_MOVE];
      xC[h] = ox[h]->xmx[p7X_C] = 0.;
      ox[h]->xmx[p7X_SCALE] = 1.0;
      ox[h]->totscale       = 0.0;
    }

  for (i = 1; i <= Lmax; i++)
    {
      act[0] = (i <= L[0]);
      act[1] = (i <= L[1]);
      rpa    = act[0] ? om->rfv[dsq[0][i]] : zrow;
      rpb    = act[1] ? om->rfv[dsq[1][i]] : zrow;
      tp     = om->tfv;
      dcv    = zerov;
      xEv    = zerov;
      xBv    = halves(act[0] ? xB[0] : 0.0f, act[1] ? xB[1] : 0.0f);

      mpv = rshiftz(PM(dp,Q-1));
      dpv = rshiftz(PD(dp,Q-1));
      ipv = rshiftz(PI(dp,Q-1));

      for (q = 0; q < Q; q++)
	{
	  /* new M(i,q), held in sv */
	  sv   =                   _mm256_mul_ps(xBv, both(tp));  tp++;
	  sv   = _mm256_add_ps(sv, _mm256_mul_ps(mpv, both(tp))); tp++;
	  sv   = _mm256_add_ps(sv, _mm256_mul_ps(ipv, both(tp))); tp++;
	  sv   = _mm256_add_ps(sv, _mm256_mul_ps(dpv, both(tp))); tp++;
	  sv   = _mm256_mul_ps(sv, two(rpa + q, rpb + q));
	  xEv  = _mm256_add_ps(xEv, sv);

	  /* {MDI}(i-1,q); then the delayed stores of {MD}(i,q) */
	  mpv = PM(dp,q);
	  dpv = PD(dp,q);
	  ipv = PI(dp,q);
	  PM(dp,q) = sv;
	  PD(dp,q) = dcv;

	  /* M->D part of D(i,q+1), held in dcv; I(i,q) */
	  dcv      = _mm256_mul_ps(sv, both(tp)); tp++;
	  sv       =                   _mm256_mul_ps(mpv, both(tp));  tp++;
	  PI(dp,q) = _mm256_add_ps(sv, _mm256_mul_ps(ipv, both(tp))); tp++;
	}

      /* D->D paths: one full pass first */
      dcv      = rshiftz(dcv);
      PD(dp,0) = zerov;
      tp       = om->tfv + 7*Q;
      for (q = 0; q < Q; q++)
	{
	  PD(dp,q) = _mm256_add_ps(dcv, PD(dp,q));
	  dcv      = _mm256_mul_ps(PD(dp,q), both(tp)); tp++;
	}

      if (om->M < 100)
	{			/* fully serialized: three more passes for both */
	  for (j = 1; j < 4; j++)
	    {
	      dcv = rshiftz(dcv);
	      tp  = om->tfv + 7*Q;
	      for (q = 0; q < Q; q++)
		{
		  PD(dp,q) = _mm256_add_ps(dcv, PD(dp,q));
		  dcv      = _mm256_mul_ps(dcv, both(tp)); tp++;
		}
	    }
	}
      else
	{			/* up to three more passes; a target stops after a pass that changed none of its D values */
	  live[0] = act[0];
	  live[1] = act[1];
	  for (j = 1; j < 4 && (live[0] || live[1]); j++)
	    {
	      dcv = _mm256_and_ps(rshiftz(dcv), lanes(live[0], live[1]));
	      tp  = om->tfv + 7*Q;
	      cv  = zerov;
	      for (q = 0; q < Q; q++)
		{
		  sv       = _mm256_add_ps(dcv, PD(dp,q));
		  cv       = _mm256_or_ps(cv, _mm256_cmp_ps(sv, PD(dp,q), _CMP_GT_OS));
		  PD(dp,q) = sv;
		  dcv      = _mm256_mul_ps(dcv, both(tp)); tp++;
		}
	      m = _mm256_movemask_ps(cv);
	      if (! (m & 0x0f)) live[0] = FALSE;
	      if (! (m & 0xf0)) live[1] = FALSE;
	    }
	}

      /* D's into xEv; then the special states */
      for (q = 0; q < Q; q++) xEv = _mm256_add_ps(PD(dp,q), xEv);
      hsum(xEv, &e[0], &e[1]);

      for (h = 0; h < 2; h++)
	{
	  resc[h] = FALSE;
	  if (! act[h]) continue;
	  xE[h] = e[h];
	  xN[h] =  xN[h] * xf[h][p7O_N][p7O_LOOP];
	  xC[h] = (xC[h] * xf[h][p7O_C][p7O_LOOP]) +  (xE[h] * xf[h][p7O_E][p7O_MOVE]);
	  xJ[h] = (xJ[h] * xf[h][p7O_J][p7O_LOOP]) +  (xE[h] * xf[h][p7O_E][p7O_LOOP]);
	  xB[h] = (xJ[h] * xf[h][p7O_J][p7O_MOVE]) +  (xN[h] * xf[h][p7O_N][p7O_MOVE]);

	  /* Sparse rescaling, as p7_ForwardParser() */
	  if (xE[h] > 1.0e4)
	    {
	      xN[h]  = xN[h] / xE[h];
	      xC[h]  = xC[h] / xE[h];
	      xJ[h]  = xJ[h] / xE[h];
	      xB[h]  = xB[h] / xE[h];
	      scl[h] = 1.0 / xE[h];
	      resc[h] = TRUE;
	      ox[h]->xmx[i*p7X_NXCELLS+p7X_SCALE] = xE[h];
	      ox[h]->totscale += log(xE[h]);
	      xE[h] = 1.0;
	    }
	  else ox[h]->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0;

	  ox[h]->xmx[i*p7X_NXCELLS+p7X_E] = xE[h];
	  ox[h]->xmx[i*p7X_NXCELLS+p7X_N] = xN[h];
	  ox[h]->xmx[i*p7X_NXCELLS+p7X_J] = xJ[h];
	  ox[h]->xmx[i*p7X_NXCELLS+p7X_B] = xB[h];
	  ox[h]->xmx[i*p7X_NXCELLS+p7X_C] = xC[h];
	}
      if (resc[0] || resc[1])
	{
	  __m256 sv8  = halves(resc[0] ? scl[0] : 1.0f, resc[1] ? scl[1] : 1.0f);
	  __m256 mask = lanes(resc[0], resc[1]);
	  for (q = 0; q < Q; q++)
	    {
	      PM(dp,q) = _mm256_blendv_ps(PM(dp,q), _mm256_mul_ps(PM(dp,q), sv8), mask);
	      PD(dp,q) = _mm256_blendv_ps(PD(dp,q), _mm256_mul_ps(PD(dp,q), sv8), mask);
	      PI(dp,q) = _mm256_blendv_ps(PI(dp,q), _mm256_mul_ps(PI(dp,q), sv8), mask);
	    }
	}
    }

  /* C->T, and the scores back in log space (nats) */
  for (h = 0; h < 2; h++)
    {
      if       (isnan(xC[h]))           ESL_EXCEPTION(eslERANGE, "forward score is NaN");
      else if  (L[h]>0 && xC[h] == 0.0) ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");     /* if L==0, xC *should* be 0.0; J5/118 */
      else if  (isinf(xC[h]) == 1)      ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
    }
  for (h = 0; h < 2; h++)
    if (sc[h] != NULL) *sc[h] = ox[h]->totscale + log(xC[h] * xf[h][p7O_C][p7O_MOVE]);
  return eslOK;
}


/* bwd_terminate()
 * Backward, row 0 for target <h> (whose half of <dp> holds its row 1): the
 * SSE code's termination.
 */
static int
bwd_terminate(const P7_OPROFILE *om, __m256 *dp, const __m128 *zrow, int h, const ESL_DSQ *dsq, const float (*xf)[p7O_NXTRANS],
	      P7_OMX *bck, float xN, float *opt_sc)
{
  const __m128 *tp  = om->tfv;
  const __m128 *rp  = om->rfv[dsq[1]];
  __m256        xBv = _mm256_setzero_ps();
  __m256        mpv;
  float         b[2], xB;
  int           Q   = p7O_NQF(om->M);
  int           q;

  for (q = 0; q < Q; q++)
    {
      mpv = _mm256_mul_ps(PM(dp,q), (h == 0) ? two(rp + q, zrow + q) : two(zrow + q, rp + q));
      mpv = _mm256_mul_ps(mpv,      both(tp));  tp += 7;
      xBv = _mm256_add_ps(xBv,      mpv);
    }
  hsum(xBv, &b[0], &b[1]);
  xB = b[h];

  xN = (xB * xf[p7O_N][p7O_MOVE]) + (xN * xf[p7O_N][p7O_LOOP]);

  bck->xmx[p7X_B]     = xB;
  bck->xmx[p7X_C]     = 0.0;
  bck->xmx[p7X_J]     = 0.0;
  bck->xmx[p7X_N]     = xN;
  bck->xmx[p7X_E]     = 0.0;
  bck->xmx[p7X_SCALE] = 1.0;

  if       (isnan(xN))               ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if  (bck->L>0 && xN == 0.0)   ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");    /* if L==0, xN *should* be 0.0 [J5/118]*/
  else if  (isinf(xN) == 1)          ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);
  return eslOK;
}


/* Function:  p7_BackwardParserPair()
 * Synopsis:  p7_BackwardParser() on two targets at once.
 *
 * Purpose:   The Backward parser for target A, given its Forward parser
 *            matrix <fwdA>, into <bckA>; and for target B likewise: the same
 *            values, bit for bit, as p7_BackwardParser() on each. <xfA>,
 *            <xfB> as for p7_ForwardParserPair(). Lengths >= 1; <bckA>,
 *            <bckB> must hold <LA>, <LB> rows of special states. Scores
 *            (nats) go to <*opt_scA>, <*opt_scB>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslERANGE> on a score overflow, underflow or NaN, as
 *            p7_BackwardParser() does; <eslEMEM> if the row buffer can't
 *            be allocated.
 */
int
p7_BackwardParserPair(const P7_OPROFILE *om,
		      const ESL_DSQ *dsqA, int LA, const float (*xfA)[p7O_NXTRANS], const P7_OMX *fwdA, P7_OMX *bckA, float *opt_scA,
		      const ESL_DSQ *dsqB, int LB, const float (*xfB)[p7O_NXTRANS], const P7_OMX *fwdB, P7_OMX *bckB, float *opt_scB)
{
  const ESL_DSQ *dsq[2] = { dsqA, dsqB };
  const float  (*xf[2])[p7O_NXTRANS] = { xfA, xfB };
  const P7_OMX  *fwd[2] = { fwdA, fwdB };
  P7_OMX        *bck[2] = { bckA, bckB };
  float         *sc[2]  = { opt_scA, opt_scB };
  int            L[2]   = { LA, LB };
  int            Lmax   = ESL_MAX(LA, LB);
  float          xN[2], xE[2], xB[2], xC[2], xJ[2], b[2];
  float          scl[2];
  int            ii[2], act[2], resc[2];
  __m256         mpv, ipv, dpv, mcv, dcv, tmmv, timv, tdmv, xBv, xEv, zerov;
  __m256        *dp;
  __m128        *zrow;
  const __m128  *rpa, *rpb;
  const __m128  *tp;
  int            Q = p7O_NQF(om->M);
  int            s, q, j, h;
  int            status;

  if ((dp = pair_row(bckA, Q, &zrow)) == NULL) ESL_EXCEPTION(eslEMEM, "allocation failed");
  zerov = _mm256_setzero_ps();
  dcv   = zerov;

  /* initialize each target's row L */
  for (h = 0; h < 2; h++)
    {
      bck[h]->M = om->M;
      bck[h]->L = L[h];
      bck[h]->has_own_scales = FALSE;	/* backwards scale factors are *usually* given by <fwd> */
      xJ[h]  = 0.0;
      xB[h]  = 0.0;
      xN[h]  = 0.0;
      xC[h]  = xf[h][p7O_C][p7O_MOVE];      /* C<-T */
      xE[h]  = xC[h] * xf[h][p7O_E][p7O_MOVE]; /* E<-C, no tail */
    }
  xEv = halves(xE[0], xE[1]);
  for (q = 0; q < Q; q++) PM(dp,q) = PD(dp,q) = xEv;
  for (q = 0; q < Q; q++) PI(dp,q) = zerov;

  /* row L's DD paths: 1) first segment includes xE, from DMO(q) */
  tp  = om->tfv + 8*Q - 1;
  dpv = lshiftz(PD(dp,Q-1));
  for (q = Q-1; q >= 0; q--)
    {
      dcv      = _mm256_mul_ps(dpv, both(tp));      tp--;
      PD(dp,q) = _mm256_add_ps(PD(dp,q), dcv);
      dpv      = PD(dp,q);
    }
  /* 2) three more passes, only extending DD component */
  for (j = 1; j < 4; j++)
    {
      tp  = om->tfv + 8*Q - 1;
      dcv = lshiftz(dcv);
      for (q = Q-1; q >= 0; q--)
	{
	  dcv      = _mm256_mul_ps(dcv, both(tp)); tp--;
	  PD(dp,q) = _mm256_add_ps(PD(dp,q), dcv);
	}
    }
  /* MD init */
  tp  = om->tfv + 7*Q - 3;
  dcv = lshiftz(PD(dp,0));
  for (q = Q-1; q >= 0; q--)
    {
      PM(dp,q) = _mm256_add_ps(PM(dp,q), _mm256_mul_ps(dcv, both(tp))); tp -= 7;
      dcv      = PD(dp,q);
    }

  /* rescaling with the fwd matrix's scale factors; stores */
  for (h = 0; h < 2; h++)
    {
      float fs = fwd[h]->xmx[L[h]*p7X_NXCELLS+p7X_SCALE];
      resc[h] = FALSE;
      if (fs > 1.0)
	{
	  xE[h] = xE[h] / fs;
	  xN[h] = xN[h] / fs;
	  xC[h] = xC[h] / fs;
	  xJ[h] = xJ[h] / fs;
	  xB[h] = xB[h] / fs;
	  scl[h]  = 1.0 / fs;
	  resc[h] = TRUE;
	}
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_SCALE] = fs;
      bck[h]->totscale                         = log(bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_SCALE]);
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_E] = xE[h];
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_N] = xN[h];
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_J] = xJ[h];
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_B] = xB[h];
      bck[h]->xmx[L[h]*p7X_NXCELLS+p7X_C] = xC[h];
    }
  if (resc[0] || resc[1])
    {
      __m256 sv8  = halves(resc[0] ? scl[0] : 1.0f, resc[1] ? scl[1] : 1.0f);
      __m256 mask = lanes(resc[0], resc[1]);
      for (q = 0; q < Q; q++)
	{
	  PM(dp,q) = _mm256_blendv_ps(PM(dp,q), _mm256_mul_ps(PM(dp,q), sv8), mask);
	  PD(dp,q) = _mm256_blendv_ps(PD(dp,q), _mm256_mul_ps(PD(dp,q), sv8), mask);
	  PI(dp,q) = _mm256_blendv_ps(PI(dp,q), _mm256_mul_ps(PI(dp,q), sv8), mask);
	}
    }
  for (h = 0; h < 2; h++)  /* a target of length 1 has no more rows */
    if (L[h] == 1 && (status = bwd_terminate(om, dp, zrow, h, dsq[h], xf[h], bck[h], xN[h], sc[h])) != eslOK) return status;

  /* main recursion: step s is row i = L - s of each target that has one */
  for (s = 1; s < Lmax; s++)
    {
      for (h = 0; h < 2; h++) { ii[h] = L[h] - s; act[h] = (ii[h] >= 1); }

      /* phase 1. B(i) collected. Old row destroyed, new row contains
       *    complete I(i,k), partial {MD}(i,k) w/ no {MD}->{DE} paths yet.
       */
      rpa = (act[0] ? om->rfv[dsq[0][ii[0]+1]] : zrow) + Q-1;
      rpb = (act[1] ? om->rfv[dsq[1][ii[1]+1]] : zrow) + Q-1;
      tp  = om->tfv + 7*Q - 1;

      tmmv = lshiftz(both(om->tfv + 1));
      timv = lshiftz(both(om->tfv + 2));
      tdmv = lshiftz(both(om->tfv + 3));

      mpv = _mm256_mul_ps(PM(dp,0), two(rpa - (Q-1), rpb - (Q-1))); /* precalc M(i+1,k+1) * e(M_k+1, x_{i+1}) */
      mpv = lshiftz(mpv);

      xBv = zerov;
      for (q = Q-1; q >= 0; q--)
	{
	  ipv      = PI(dp,q);
	  PI(dp,q) = _mm256_add_ps(_mm256_mul_ps(ipv, both(tp)), _mm256_mul_ps(mpv, timv));   tp--;
	  PD(dp,q) =                                             _mm256_mul_ps(mpv, tdmv);
	  mcv      = _mm256_add_ps(_mm256_mul_ps(ipv, both(tp)), _mm256_mul_ps(mpv, tmmv));   tp-= 2;

	  mpv      = _mm256_mul_ps(PM(dp,q), two(rpa, rpb));  rpa--; rpb--;
	  PM(dp,q) = mcv;

	  tdmv = both(tp);   tp--;
	  timv = both(tp);   tp--;
	  tmmv = both(tp);   tp--;

	  xBv = _mm256_add_ps(xBv, _mm256_mul_ps(mpv, both(tp))); tp--;
	}

      /* phase 2: the specials */
      hsum(xBv, &b[0], &b[1]);
      for (h = 0; h < 2; h++)
	{
	  if (! act[h]) continue;
	  xB[h] = b[h];
	  xC[h] =  xC[h] * xf[h][p7O_C][p7O_LOOP];
	  xJ[h] = (xB[h] * xf[h][p7O_J][p7O_MOVE]) + (xJ[h] * xf[h][p7O_J][p7O_LOOP]); /* must come after xB */
	  xN[h] = (xB[h] * xf[h][p7O_N][p7O_MOVE]) + (xN[h] * xf[h][p7O_N][p7O_LOOP]); /* must come after xB */
	  xE[h] = (xC[h] * xf[h][p7O_E][p7O_MOVE]) + (xJ[h] * xf[h][p7O_E][p7O_LOOP]); /* must come after xJ, xC */
	}
      xEv = halves(act[0] ? xE[0] : 0.0f, act[1] ? xE[1] : 0.0f);

      /* phase 3: {MD}->E paths and one step of the D->D paths */
      tp  = om->tfv + 8*Q - 1;
      dpv = _mm256_add_ps(PD(dp,0), xEv);
      dpv = lshiftz(dpv);
      for (q = Q-1; q >= 0; q--)
	{
	  dcv      = _mm256_mul_ps(dpv, both(tp)); tp--;
	  PD(dp,q) = _mm256_add_ps(PD(dp,q), _mm256_add_ps(dcv, xEv));
	  dpv      = PD(dp,q);
	  PM(dp,q) = _mm256_add_ps(PM(dp,q), xEv);
	}

      /* phase 4: finish extending the DD paths, fully serialized */
      for (j = 1; j < 4; j++)
	{
	  dcv = lshiftz(dcv);
	  tp  = om->tfv + 8*Q - 1;
	  for (q = Q-1; q >= 0; q--)
	    {
	      dcv      = _mm256_mul_ps(dcv, both(tp)); tp--;
	      PD(dp,q) = _mm256_add_ps(PD(dp,q), dcv);
	    }
	}

      /* phase 5: add M->D paths */
      dcv = lshiftz(PD(dp,0));
      tp  = om->tfv + 7*Q - 3;
      for (q = Q-1; q >= 0; q--)
	{
	  PM(dp,q) = _mm256_add_ps(PM(dp,q), _mm256_mul_ps(dcv, both(tp))); tp -= 7;
	  dcv      = PD(dp,q);
	}

      /* Sparse rescaling, as p7_BackwardParser(); stores */
      for (h = 0; h < 2; h++)
	{
	  int i = ii[h];
	  resc[h] = FALSE;
	  if (! act[h]) continue;
	  if (xB[h] > 1.0e16) bck[h]->has_own_scales = TRUE;
	  if      (bck[h]->has_own_scales)  bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB[h] > 1.0e4) ? xB[h] : 1.0;
	  else                              bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE] = fwd[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];

	  if (bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0)
	    {
	      xE[h] /= bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      xN[h] /= bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      xJ[h] /= bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      xB[h] /= bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      xC[h] /= bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      scl[h]  = 1.0 / bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE];
	      resc[h] = TRUE;
	      bck[h]->totscale += log(bck[h]->xmx[i*p7X_NXCELLS+p7X_SCALE]);
	    }

	  bck[h]->xmx[i*p7X_NXCELLS+p7X_E] = xE[h];
	  bck[h]->xmx[i*p7X_NXCELLS+p7X_N] = xN[h];
	  bck[h]->xmx[i*p7X_NXCELLS+p7X_J] = xJ[h];
	  bck[h]->xmx[i*p7X_NXCELLS+p7X_B] = xB[h];
	  bck[h]->xmx[i*p7X_NXCELLS+p7X_C] = xC[h];
	}
      if (resc[0] || resc[1])
	{
	  __m256 sv8  = halves(resc[0] ? scl[0] : 1.0f, resc[1] ? scl[1] : 1.0f);
	  __m256 mask = lanes(resc[0], resc[1]);
	  for (q = 0; q < Q; q++)
	    {
	      PM(dp,q) = _mm256_blendv_ps(PM(dp,q), _mm256_mul_ps(PM(dp,q), sv8), mask);
	      PD(dp,q) = _mm256_blendv_ps(PD(dp,q), _mm256_mul_ps(PD(dp,q), sv8), mask);
	      PI(dp,q) = _mm256_blendv_ps(PI(dp,q), _mm256_mul_ps(PI(dp,q), sv8), mask);
	    }
	}

      /* a target that has just done its row 1 terminates now, before the other's next row overwrites its half */
      for (h = 0; h < 2; h++)
	if (ii[h] == 1 && (status = bwd_terminate(om, dp, zrow, h, dsq[h], xf[h], bck[h], xN[h], sc[h])) != eslOK) return status;
    }
  return eslOK;
}

#endif /* p7_ENABLE_AVX2 && __AVX2__ */
#endif /* ! p7FWDBACK_PAIR_TESTDRIVE */



/*****************************************************************
 * Unit tests: the pair parsers against p7_Forward/BackwardParser().
 *****************************************************************/
#ifdef p7FWDBACK_PAIR_TESTDRIVE
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"

/* every value the pair parsers write, compared bit for bit */
static void
compare_omx(const P7_OMX *a, const P7_OMX *b, int L, const char *what, int M, int La, int Lb)
{
  int i;
  if (a->M != b->M || a->L != b->L || a->has_own_scales != b->has_own_scales || memcmp(&a->totscale, &b->totscale, sizeof(float)) != 0)
    esl_fatal("fwdback_pair utest: %s: header differs (M=%d, L=%d|%d)", what, M, La, Lb);
  for (i = 0; i <= L; i++)
    if (memcmp(a->xmx + i*p7X_NXCELLS, b->xmx + i*p7X_NXCELLS, sizeof(float) * p7X_NXCELLS) != 0)
      esl_fatal("fwdback_pair utest: %s: row %d differs (M=%d, L=%d|%d)", what, i, M, La, Lb);
}

/* A target for <om>: half are random (iid, background frequencies), half
 * are emitted by the model (high scores, so rows get rescaled), of length
 * about <L>.
 */
static void
sample_target(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, int L, ESL_SQ *sq)
{
  esl_sq_Reuse(sq);
  if (esl_rnd_Roll(r, 2) == 0)
    {
      esl_sq_GrowTo(sq, L);
      esl_rsq_xfIID(r, bg->f, hmm->abc->K, L, sq->dsq);
      sq->n = L;
    }
  else
    {
      do { p7_ProfileEmit(r, hmm, gm, bg, sq, NULL); } while (sq->n == 0 || sq->n > 3*L + 10);
    }
}

/* <ddheavy>: models with strong M->D and D->D transitions (long deletion
 * chains), so the Forward parser's D->D passes run long and the two
 * targets of a pair often stop after different passes.
 */
static int nresc_fwd = 0;   /* coverage: forward rows rescaled, targets whose backward used its own scales, */
static int nown_bck  = 0;   /*           pairs of unequal length                                          */
static int nunequal  = 0;

static void
utest_pairs(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, P7_BG *bg, int M, int Lbase, int N, int ddheavy)
{
  P7_HMM      *hmm  = NULL;
  P7_PROFILE  *gm   = NULL;
  P7_OPROFILE *om   = NULL;
  ESL_SQ      *sqa  = esl_sq_CreateDigital(abc);
  ESL_SQ      *sqb  = esl_sq_CreateDigital(abc);
  P7_OMX      *fa   = p7_omx_Create(M, 0, 100), *fb  = p7_omx_Create(M, 0, 100);   /* singles */
  P7_OMX      *ba   = p7_omx_Create(M, 0, 100), *bb  = p7_omx_Create(M, 0, 100);
  P7_OMX      *pfa  = p7_omx_Create(M, 0, 100), *pfb = p7_omx_Create(M, 0, 100);   /* pairs */
  P7_OMX      *pba  = p7_omx_Create(M, 0, 100), *pbb = p7_omx_Create(M, 0, 100);
  float        xfa[p7O_NXSTATES][p7O_NXTRANS], xfb[p7O_NXSTATES][p7O_NXTRANS];
  float        sa, sb, psa, psb, bsa, bsb, pbsa, pbsb;
  int          n, La, Lb;

  if      (ddheavy)         p7_hmm_SampleUniform(r, M, abc, 0.05, 0.5, 0.3 + 0.4 * esl_random(r), 0.9 + 0.09 * esl_random(r), &hmm);
  else if (esl_rnd_Roll(r, 2)) p7_hmm_SampleUniform(r, M, abc, 0.05, 0.5, 0.05, 0.2, &hmm);  /* Pfam-like transitions */
  else                    p7_hmm_Sample(r, M, abc, &hmm);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, Lbase, p7_LOCAL);
  p7_oprofile_Convert(gm, om);

  for (n = 0; n < N; n++)
    {
      sample_target(r, hmm, gm, bg, 1 + esl_rnd_Roll(r, Lbase), sqa);
      if (esl_rnd_Roll(r, 4) == 0) esl_sq_Copy(sqa, sqb);   /* equal lengths too */
      else sample_target(r, hmm, gm, bg, 1 + esl_rnd_Roll(r, Lbase), sqb);
      La = sqa->n; Lb = sqb->n;

      /* singles, each with <om> configured for its own length */
      p7_oprofile_ReconfigLength(om, La);  memcpy(xfa, om->xf, sizeof(xfa));
      p7_omx_GrowTo(fa, M, 0, La); p7_omx_GrowTo(ba, M, 0, La);
      p7_ForwardParser (sqa->dsq, La, om, fa, &sa);
      p7_BackwardParser(sqa->dsq, La, om, fa, ba, &bsa);
      p7_oprofile_ReconfigLength(om, Lb);  memcpy(xfb, om->xf, sizeof(xfb));
      p7_omx_GrowTo(fb, M, 0, Lb); p7_omx_GrowTo(bb, M, 0, Lb);
      p7_ForwardParser (sqb->dsq, Lb, om, fb, &sb);
      p7_BackwardParser(sqb->dsq, Lb, om, fb, bb, &bsb);

      /* the pair, with <om> configured for neither on purpose (only xfa, xfb matter) */
      p7_oprofile_ReconfigLength(om, 7);
      p7_omx_GrowTo(pfa, M, 0, La); p7_omx_GrowTo(pba, M, 0, La);
      p7_omx_GrowTo(pfb, M, 0, Lb); p7_omx_GrowTo(pbb, M, 0, Lb);
      if (p7_ForwardParserPair (om, sqa->dsq, La, (const float (*)[p7O_NXTRANS]) xfa, pfa, &psa, sqb->dsq, Lb, (const float (*)[p7O_NXTRANS]) xfb, pfb, &psb) != eslOK)
	esl_fatal("fwdback_pair utest: forward pair failed");
      if (p7_BackwardParserPair(om, sqa->dsq, La, (const float (*)[p7O_NXTRANS]) xfa, pfa, pba, &pbsa, sqb->dsq, Lb, (const float (*)[p7O_NXTRANS]) xfb, pfb, pbb, &pbsb) != eslOK)
	esl_fatal("fwdback_pair utest: backward pair failed");

      if (memcmp(&sa, &psa, sizeof(float)) || memcmp(&sb, &psb, sizeof(float)))
	esl_fatal("fwdback_pair utest: forward scores differ: %a|%a vs %a|%a (M=%d, L=%d|%d)", sa, sb, psa, psb, M, La, Lb);
      if (memcmp(&bsa, &pbsa, sizeof(float)) || memcmp(&bsb, &pbsb, sizeof(float)))
	esl_fatal("fwdback_pair utest: backward scores differ: %a|%a vs %a|%a (M=%d, L=%d|%d)", bsa, bsb, pbsa, pbsb, M, La, Lb);
      compare_omx(fa, pfa, La, "forward A",  M, La, Lb);
      compare_omx(fb, pfb, Lb, "forward B",  M, La, Lb);
      compare_omx(ba, pba, La, "backward A", M, La, Lb);
      compare_omx(bb, pbb, Lb, "backward B", M, La, Lb);

      /* Backward again, given forward matrices without scale factors: the
       * values outgrow them and the Backward parser switches to its own
       * (the has_own_scales path).
       */
      if (esl_rnd_Roll(r, 4) == 0)
	{
	  int i;
	  for (i = 0; i <= La; i++) fa->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0;
	  for (i = 0; i <= Lb; i++) fb->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0;
	  p7_oprofile_ReconfigLength(om, La);  p7_BackwardParser(sqa->dsq, La, om, fa, ba, &bsa);
	  p7_oprofile_ReconfigLength(om, Lb);  p7_BackwardParser(sqb->dsq, Lb, om, fb, bb, &bsb);
	  if (p7_BackwardParserPair(om, sqa->dsq, La, (const float (*)[p7O_NXTRANS]) xfa, fa, pba, &pbsa, sqb->dsq, Lb, (const float (*)[p7O_NXTRANS]) xfb, fb, pbb, &pbsb) != eslOK)
	    esl_fatal("fwdback_pair utest: backward pair failed");
	  if (memcmp(&bsa, &pbsa, sizeof(float)) || memcmp(&bsb, &pbsb, sizeof(float)))
	    esl_fatal("fwdback_pair utest: unscaled backward scores differ (M=%d, L=%d|%d)", M, La, Lb);
	  compare_omx(ba, pba, La, "unscaled backward A", M, La, Lb);
	  compare_omx(bb, pbb, Lb, "unscaled backward B", M, La, Lb);
	  nown_bck += ba->has_own_scales + bb->has_own_scales;
	  continue;
	}

      { int i; for (i = 1; i <= La; i++) nresc_fwd += (fa->xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0);
	       for (i = 1; i <= Lb; i++) nresc_fwd += (fb->xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0); }
      nown_bck += ba->has_own_scales + bb->has_own_scales;
      nunequal += (La != Lb);
    }

  p7_omx_Destroy(fa);  p7_omx_Destroy(fb);  p7_omx_Destroy(ba);  p7_omx_Destroy(bb);
  p7_omx_Destroy(pfa); p7_omx_Destroy(pfb); p7_omx_Destroy(pba); p7_omx_Destroy(pbb);
  esl_sq_Destroy(sqa); esl_sq_Destroy(sqb);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-N",        eslARG_INT,     "40", NULL, NULL,  NULL,  NULL, NULL, "number of target pairs per model",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the paired Forward/Backward parsers";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N    = esl_opt_GetInteger(go, "-N");
  int             Ms[] = { 1, 2, 3, 4, 5, 8, 13, 31, 64, 99, 100, 101, 150, 200, 257, 400 };
  int             atype, k, rep;

  impl_Init();
  if (! p7_fb_PairAvailable()) { printf("no AVX2 on this host: nothing to test\n"); return 0; }

  for (atype = 0; atype < 2; atype++)
    {
      ESL_ALPHABET *abc = esl_alphabet_Create(atype == 0 ? eslAMINO : eslDNA);
      P7_BG        *bg  = p7_bg_Create(abc);
      for (k = 0; k < (int) (sizeof(Ms) / sizeof(int)); k++)
	for (rep = 0; rep < 3; rep++)
	  utest_pairs(r, abc, bg, Ms[k], (rep == 2 ? 600 : 150), N, FALSE);
      for (k = 0; k < (int) (sizeof(Ms) / sizeof(int)); k++)
	if (Ms[k] >= 100)
	  for (rep = 0; rep < 4; rep++)
	    utest_pairs(r, abc, bg, Ms[k], 150, N, TRUE);
      p7_bg_Destroy(bg);
      esl_alphabet_Destroy(abc);
    }
  if (nresc_fwd == 0 || nown_bck == 0 || nunequal == 0) esl_fatal("fwdback_pair utest: rescaling, backward's own scales or unequal lengths not reached");
  printf("ok (forward rows rescaled: %d; backward with own scales: %d; unequal-length pairs: %d)\n", nresc_fwd, nown_bck, nunequal);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /* p7FWDBACK_PAIR_TESTDRIVE */
