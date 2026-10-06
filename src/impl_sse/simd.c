/* Run-time choice of the vector width for the integer filters (SSV, MSV,
 * Viterbi filter).
 *
 * The SSE kernels always exist. Builds on x86 may also carry AVX2
 * versions of the same kernels (configure: p7_ENABLE_AVX2), compiled
 * with -mavx2 in their own files; they run only on CPUs that have the
 * instructions. They compute the same integer DP on the same scores,
 * striped for more lanes, so every choice gives identical results.
 *
 * The environment variable HMMER_SIMD=sse|avx2 forces one path. It
 * exists for testing and benchmarking one binary on each path
 * (dev/check.sh --simd); asking for a path that this build or CPU
 * doesn't have is an error.
 *
 * The AVX2 kernels read their own copies of the profile's filter
 * scores: the same values, striped for 32 lanes (p7_oprofile_SetSIMD()).
 *
 * Contents:
 *   1. Choosing the level.
 *   2. AVX2 copies of a profile's filter scores.
 *   3. Unit tests.
 *   4. Test driver.
 */
#include <p7_config.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"

#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. Choosing the level.
 *****************************************************************/

static int simd_level = -1;     /* p7_SIMD_*; -1 until first use. Accessed atomically. */

/* simd_cpu_max()
 * The widest level that both this build and this CPU support.
 * (__builtin_cpu_supports() also checks that the OS saves the
 * wide registers.)
 */
static int
simd_cpu_max(void)
{
  int level = p7_SIMD_SSE;
#ifdef p7_ENABLE_AVX2
  if (__builtin_cpu_supports("avx2")) level = p7_SIMD_AVX2;
#endif
  return level;
}

/* simd_parse()
 * Level asked for by <s> (the value of HMMER_SIMD, or NULL) on a host
 * whose widest level is <max>. Unset or empty means <max>.
 * Returns <eslOK> and the level in <*ret_level>; <eslEINVAL> for an
 * unknown name; <eslENORESULT> for a level above <max>.
 */
static int
simd_parse(const char *s, int max, int *ret_level)
{
  int level;

  if      (s == NULL || *s == '\0') level = max;
  else if (strcmp(s, "sse")  == 0)  level = p7_SIMD_SSE;
  else if (strcmp(s, "avx2") == 0)  level = p7_SIMD_AVX2;
  else return eslEINVAL;

  if (level > max) return eslENORESULT;
  *ret_level = level;
  return eslOK;
}


/* Function:  p7_simd_Level()
 * Synopsis:  The vector width the integer filters use in this process.
 *
 * Purpose:   Returns <p7_SIMD_SSE> or <p7_SIMD_AVX2>: the widest level
 *            this build and CPU support, or the one HMMER_SIMD asks
 *            for. Decided on the first call, then fixed. Safe to call
 *            from several threads.
 *
 * Throws:    Fails (<p7_Fail()>) if HMMER_SIMD names an unknown level
 *            or one this build or CPU doesn't have.
 */
int
p7_simd_Level(void)
{
  int   level = __atomic_load_n(&simd_level, __ATOMIC_RELAXED);
  char *s;
  int   max;
  int   status;

  if (level >= 0) return level;

  s      = getenv("HMMER_SIMD");
  max    = simd_cpu_max();
  status = simd_parse(s, max, &level);
  if (status == eslEINVAL)    p7_Fail("HMMER_SIMD=%s: expected sse or avx2", s);
  if (status == eslENORESULT) p7_Fail("HMMER_SIMD=%s: not available (this build and CPU support up to %s)", s, p7_simd_Name(max));

  __atomic_store_n(&simd_level, level, __ATOMIC_RELAXED); /* threads racing here all store the same value */
  return level;
}

/* Function:  p7_simd_Name()
 * Synopsis:  Name of a level, as HMMER_SIMD spells it.
 */
const char *
p7_simd_Name(int level)
{
  switch (level) {
  case p7_SIMD_SSE:  return "sse";
  case p7_SIMD_AVX2: return "avx2";
  default:           return "unknown";
  }
}
/*---------------- end, choosing the level ----------------------*/



/*****************************************************************
 * 2. AVX2 copies of a profile's filter scores.
 *****************************************************************/

/* stripe_index()
 * Offsets of nodes k=1..M (index k-1) in a row striped for <Q> vectors
 * of <n> lanes: vector (k-1)%Q, lane (k-1)/Q. Nested loops instead of
 * wrapping counters: the wraps were a branch miss every few nodes.
 */
static void
stripe_index(int M, int Q, int n, int *pos)
{
  int q, z, k;

  for (z = 0, k = 0; z < n && k < M; z++)
    for (q = 0; q < Q && k < M; q++, k++)
      pos[k] = q * n + z;
}

/* stripe_u8()
 * Copy the <M> node values of one SSE-striped byte row <src> (node
 * offsets <spos>) into <dst> (node offsets <dpos>, <nb> bytes). Padding
 * nodes get <pad>: what the SSE conversion gives its own padding.
 */
static void
stripe_u8(const uint8_t *src, const int *spos, uint8_t *dst, const int *dpos, int M, int nb, uint8_t pad)
{
  int k;

  memset(dst, pad, nb);
  for (k = 0; k < M; k++) dst[dpos[k]] = src[spos[k]];
}


/* Function:  p7_oprofile_SetSIMD()
 * Synopsis:  Make (or drop) the AVX2 copies of a profile's filter scores.
 *
 * Purpose:   For <level> <p7_SIMD_AVX2>, fill the AVX2 copies of
 *            <om>'s MSV and SSV scores from the SSE ones and set <om->simd>,
 *            so the filters run their AVX2 versions on <om>. For
 *            <p7_SIMD_SSE>, only reset <om->simd>: the filters use
 *            the SSE versions.
 *
 *            <p7_oprofile_Convert()> calls this with <p7_simd_Level()>.
 *            Anything that changes the SSE scores of a profile with
 *            copies has to call it again.
 *
 *            The copies hold the same values and pad nodes beyond M
 *            like the SSE conversion does, so the AVX2 kernels compute
 *            the same DP.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> for a level
 *            this build doesn't have.
 */
int
p7_oprofile_SetSIMD(P7_OPROFILE *om, int level)
{
  int     Kp = om->abc->Kp;
  int     Qs = p7O_NQB(om->M);       /* # of SSE byte vectors   */
  int     Qd = p7O_NQV(om->M, 32);   /* # of AVX2 byte vectors  */
  int    *spos, *dpos;               /* node k's offset in an SSE row, an AVX2 row */
  __m128i tmp, tmp2;
  int     x, q, j;
  int     status;

  if (level == p7_SIMD_SSE) { om->simd = p7_SIMD_SSE; return eslOK; }
#ifdef p7_ENABLE_AVX2
  if (level != p7_SIMD_AVX2) ESL_EXCEPTION(eslEINVAL, "no such SIMD level");

  /* Allocated once, for om->allocM (Convert() refuses bigger models):
   * row pointers, then 64-byte aligned rows of Q + extra vectors.
   */
  if (om->simd_mem == NULL)
    {
      size_t ptrsz = sizeof(uint8_t *) * Kp;
      size_t idxsz = sizeof(int) * 2 * om->allocM;
      size_t rbsz  = (size_t) p7O_NQV(om->allocM, 32) * 32;
      size_t sbsz  = (size_t) (p7O_NQV(om->allocM, 32) + p7O_EXTRA_SB_AVX2) * 32;
      uint8_t *data;

      om->simd_nalloc = 2 * ptrsz + idxsz + 63 + (rbsz + sbsz) * Kp;
      ESL_ALLOC(om->simd_mem, om->simd_nalloc);
      om->rbw = (uint8_t **) om->simd_mem;
      om->sbw = om->rbw + Kp;
      data    = (uint8_t *) (((uintptr_t) om->simd_mem + 2 * ptrsz + idxsz + 63) & ~((uintptr_t) 63));
      for (x = 0; x < Kp; x++) om->rbw[x] = data + x * rbsz;
      data   += Kp * rbsz;
      for (x = 0; x < Kp; x++) om->sbw[x] = data + x * sbsz;
    }
  spos = (int *) ((uint8_t *) om->simd_mem + 2 * sizeof(uint8_t *) * Kp);  /* scratch: node offsets */
  dpos = spos + om->allocM;

  /* MSV costs. SSV scores from them as sf_conversion() makes sbv from rbv (bytewise:
   * ((127 + bias) -sat rbv) ^ 127), and the extra vectors past Q that
   * ssvfilter_avx2.c reads: copies of q % Q, as in sf_conversion().
   */
  stripe_index(om->M, Qs, 16, spos);
  stripe_index(om->M, Qd, 32, dpos);
  tmp  = _mm_set1_epi8((int8_t) (om->bias_b + 127));
  tmp2 = _mm_set1_epi8(127);
  for (x = 0; x < Kp; x++)
    {
      stripe_u8((const uint8_t *) om->rbv[x], spos, om->rbw[x], dpos, om->M, Qd * 32, 255);
      for (j = 0; j < 2 * Qd; j++)
	((__m128i *) om->sbw[x])[j] = _mm_xor_si128(_mm_subs_epu8(tmp, ((__m128i *) om->rbw[x])[j]), tmp2);
      for (q = Qd; q < Qd + p7O_EXTRA_SB_AVX2; q++)
	memcpy(om->sbw[x] + q * 32, om->sbw[x] + (q % Qd) * 32, 32);
    }

  om->simd = level;
  return eslOK;

 ERROR:
  return status;
#else
  ESL_EXCEPTION(eslEINVAL, "no AVX2 in this build");
#endif
}


/* Function:  p7_omx_SIMDRow()
 * Synopsis:  The AVX2 filters' one-row DP storage in a <P7_OMX>.
 *
 * Purpose:   Returns a 64-byte aligned row of at least <nbytes> bytes
 *            in <ox> (contents undefined), growing it if needed; or
 *            NULL if that allocation fails.
 */
void *
p7_omx_SIMDRow(P7_OMX *ox, size_t nbytes)
{
  if (nbytes + 63 > ox->simd_nalloc)
    {
      free(ox->simd_mem);
      ox->simd_nalloc = 0;
      if ((ox->simd_mem = malloc(nbytes + 63)) == NULL) return NULL;
      ox->simd_nalloc = nbytes + 63;
    }
  return (void *) (((uintptr_t) ox->simd_mem + 63) & ~((uintptr_t) 63));
}
/*-------------- end, AVX2 copies of profile scores -------------*/



/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7SIMD_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

extern uint8_t get_xE(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om); /* ssvfilter.c: the SSE DP */

/* utest_parse()
 * HMMER_SIMD values on hosts of every width.
 */
static void
utest_parse(void)
{
  char msg[] = "simd parse unit test failed";
  int  max, level;

  for (max = p7_SIMD_SSE; max <= p7_SIMD_AVX2; max++)
    {
      level = -1; if (simd_parse(NULL, max, &level) != eslOK || level != max)          esl_fatal(msg);
      level = -1; if (simd_parse("",   max, &level) != eslOK || level != max)          esl_fatal(msg);
      level = -1; if (simd_parse("sse", max, &level) != eslOK || level != p7_SIMD_SSE) esl_fatal(msg);
      if (max >= p7_SIMD_AVX2) { level = -1; if (simd_parse("avx2", max, &level) != eslOK || level != p7_SIMD_AVX2) esl_fatal(msg); }
      else                     {             if (simd_parse("avx2", max, &level) != eslENORESULT)                  esl_fatal(msg); }
      if (simd_parse("SSE",    max, &level) != eslEINVAL) esl_fatal(msg);
      if (simd_parse("avx",    max, &level) != eslEINVAL) esl_fatal(msg);
      if (simd_parse("sse ",   max, &level) != eslEINVAL) esl_fatal(msg);
      if (simd_parse("avx512", max, &level) != eslEINVAL) esl_fatal(msg);
    }
}

/* utest_level()
 * The level in use is one the host supports, and the one HMMER_SIMD
 * asks for if it is set.
 */
static void
utest_level(void)
{
  char *s     = getenv("HMMER_SIMD");
  int   level = p7_simd_Level();
  int   want;

  if (level < p7_SIMD_SSE || level > simd_cpu_max())                   esl_fatal("simd level unit test failed: level %d", level);
  if (simd_parse(s, simd_cpu_max(), &want) != eslOK || level != want)  esl_fatal("simd level unit test failed: HMMER_SIMD not honored");
  if (p7_simd_Level() != level)                                        esl_fatal("simd level unit test failed: level changed");
}

/* Test sequences for comparing a wide kernel with the SSE one, cycling
 * through six kinds by <n>: iid from the background; uniform over all
 * residue codes (degenerate ones too); emitted by the model; emitted
 * with random flanks; the model's consensus with random flanks (top
 * scores: overflows); two consensus copies apart (the J state: the SSV
 * filter hands over to the full MSV). Lengths 0..2*Lcfg for the random
 * kinds. Grows <*dsq> as needed; returns the length.
 */
static int
sample_target(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg,
	      int Lcfg, int n, ESL_SQ *sq, ESL_DSQ **dsq, int *nalloc)
{
  const ESL_ALPHABET *abc = gm->abc;
  int L, f1 = 0, f2 = 0, g = 0, j;

  if (n % 6 >= 2)
    {
      esl_sq_Reuse(sq);
      if (n % 6 >= 4) { if (p7_emit_SimpleConsensus(hmm, sq)          != eslOK) esl_fatal("consensus failed"); }
      else            { if (p7_ProfileEmit(r, hmm, gm, bg, sq, NULL) != eslOK) esl_fatal("emit failed");      }
      if (n % 6 >= 3) { f1 = esl_rnd_Roll(r, Lcfg+1); f2 = esl_rnd_Roll(r, Lcfg+1); }
      if (n % 6 == 5) { g  = esl_rnd_Roll(r, Lcfg+1); f2 += g + sq->n; }  /* second copy goes into the right flank */
      L = f1 + sq->n + f2;
    }
  else L = esl_rnd_Roll(r, 2*Lcfg+1);

  if (L + 2 > *nalloc) { *nalloc = L + 2; *dsq = realloc(*dsq, sizeof(ESL_DSQ) * *nalloc); if (! *dsq) esl_fatal("malloc failed"); }

  switch (n % 6) {
  case 0:  esl_rsq_xfIID(r, bg->f, abc->K, L, *dsq);                                        break;
  case 1:  for (j = 1; j <= L; j++) (*dsq)[j] = esl_rnd_Roll(r, abc->Kp);                   break;
  default: esl_rsq_xfIID(r, bg->f, abc->K, f1, *dsq);                 /* flanks first: each writes */
           esl_rsq_xfIID(r, bg->f, abc->K, f2, *dsq + f1 + sq->n);     /* sentinels at its ends     */
           memcpy(*dsq + 1 + f1, sq->dsq + 1, sq->n);
           if (n % 6 == 5) memcpy(*dsq + 1 + f1 + sq->n + g, sq->dsq + 1, sq->n);            break;
  }
  (*dsq)[0] = (*dsq)[L+1] = eslDSQ_SENTINEL;
  return L;
}

/* utest_filters()
 * The AVX2 filters against the SSE ones on <N> targets for a random
 * model of length <M> configured for length <Lcfg>: same status, same
 * score bits. For SSV, also the same DP maximum xE unless both
 * overflowed. Counts SSV outcomes in <nout>: [0] eslOK, [1] eslERANGE,
 * [2] eslENORESULT from the J state test, [3] eslENORESULT from a
 * possible overflow, [4] eslENORESULT before any DP (score system can't
 * be trusted); and full MSV outcomes (after an SSV eslENORESULT): [5]
 * eslOK, [6] eslERANGE.
 */
static void
utest_filters(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, P7_BG *bg, int M, int Lcfg, int N, int nout[7])
{
  char         msg[] = "simd SSV unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_SQ      *sq  = esl_sq_CreateDigital(abc);
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  ESL_DSQ     *dsq = NULL;
  int          nalloc = 0;
  int          n, L;
  int          st0, st1, ssv_st;
  float        sc0, sc1;
  uint8_t      xE0, xE1;
  int          ovfl;

  if (p7_oprofile_Sample(r, abc, bg, M, Lcfg, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  if (p7_oprofile_SetSIMD(om, p7_SIMD_AVX2)                    != eslOK) esl_fatal(msg);

  for (n = 0; n < N; n++)
    {
      L = sample_target(r, hmm, gm, bg, Lcfg, n, sq, &dsq, &nalloc);
      p7_oprofile_ReconfigLength(om, L);  /* as the pipeline does; changes tjb_b */

      sc0 = sc1 = -42.0f;
      om->simd = p7_SIMD_SSE;  st0 = p7_SSVFilter(dsq, L, om, &sc0);  xE0 = get_xE(dsq, L, om);
      om->simd = p7_SIMD_AVX2; st1 = p7_SSVFilter(dsq, L, om, &sc1);  xE1 = p7_SSVFilter_xE_avx2(dsq, L, om);

      ovfl = (xE0 >= 255 - om->bias_b);
      if (st0 != st1 || memcmp(&sc0, &sc1, sizeof(float)) != 0)
	esl_fatal("%s: M=%d L=%d: SSE status %d score %g, AVX2 status %d score %g", msg, M, L, st0, sc0, st1, sc1);
      if (ovfl ? (xE1 < 255 - om->bias_b) : (xE1 != xE0))
	esl_fatal("%s: M=%d L=%d: xE %d (SSE) vs %d (AVX2), bias %d", msg, M, L, xE0, xE1, om->bias_b);

      if      (st0 == eslOK)      nout[0]++;
      else if (st0 == eslERANGE)  nout[1]++;
      else if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) nout[4]++;
      else if (ovfl)              nout[3]++;
      else                        nout[2]++;
      ssv_st = st0;

      /* MSV: the SSV filter, then the full DP if it couldn't decide */
      sc0 = sc1 = -42.0f;
      om->simd = p7_SIMD_SSE;  st0 = p7_MSVFilter(dsq, L, om, ox, &sc0);
      om->simd = p7_SIMD_AVX2; st1 = p7_MSVFilter(dsq, L, om, ox, &sc1);
      if (st0 != st1 || memcmp(&sc0, &sc1, sizeof(float)) != 0)
	esl_fatal("simd MSV unit test failed: M=%d L=%d: SSE status %d score %g, AVX2 status %d score %g", M, L, st0, sc0, st1, sc1);
      if (ssv_st == eslENORESULT) nout[st0 == eslOK ? 5 : 6]++;
    }

  free(dsq);
  p7_omx_Destroy(ox);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* utest_msv_boundary()
 * The full MSV DP at its overflow boundary. Targets are a long random
 * flank (so that the SSV filter can't trust its overflow test and hands
 * over, see ssvfilter.c) followed by the consensus's prefix 1..j, for
 * every j: the score climbs residue by residue, so the last row's max
 * passes through each value up to the overflow, often landing exactly
 * on it (max + bias == 255, the case an off-by-one in the overflow test
 * gets wrong). Counts in <nout>: [0] targets that reached the full DP,
 * [1] of those that overflowed.
 */
static void
utest_msv_boundary(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, P7_BG *bg, int M, int flank, int nout[2])
{
  char         msg[] = "simd MSV boundary unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  ESL_SQ      *sq  = esl_sq_CreateDigital(abc);
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (flank + M + 2));
  float        sc0, sc1;
  int          j, st0, st1;
  int          novfl = 0;

  if (p7_oprofile_Sample(r, abc, bg, M, flank, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  if (p7_oprofile_SetSIMD(om, p7_SIMD_AVX2)                     != eslOK) esl_fatal(msg);
  if (p7_emit_SimpleConsensus(hmm, sq)                          != eslOK) esl_fatal(msg);
  esl_rsq_xfIID(r, bg->f, abc->K, flank, dsq);
  for (j = 1; j <= sq->n; j++)
    {
      dsq[flank + j]     = sq->dsq[j];
      dsq[flank + j + 1] = eslDSQ_SENTINEL;
      p7_oprofile_ReconfigLength(om, flank + j);
      if (p7_SSVFilter(dsq, flank + j, om, &sc0) != eslENORESULT) continue;

      sc0 = sc1 = -42.0f;
      om->simd = p7_SIMD_SSE;  st0 = p7_MSVFilter(dsq, flank + j, om, ox, &sc0);
      om->simd = p7_SIMD_AVX2; st1 = p7_MSVFilter(dsq, flank + j, om, ox, &sc1);
      if (st0 != st1 || memcmp(&sc0, &sc1, sizeof(float)) != 0)
	esl_fatal("%s: M=%d L=%d: SSE status %d score %g, AVX2 status %d score %g", msg, M, flank + j, st0, sc0, st1, sc1);
      nout[0]++;
      if (st0 == eslERANGE && ++novfl > 5) break;  /* well past the boundary */
    }
  nout[1] += novfl;

  free(dsq);
  esl_sq_Destroy(sq);
  p7_omx_Destroy(ox);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* utest_copies()
 * A profile's AVX2 copies survive p7_oprofile_Copy(), and a Clone()
 * shares them; destroying either doesn't touch the other's.
 */
static void
utest_copies(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, P7_BG *bg)
{
  char         msg[] = "simd copies unit test failed";
  P7_OPROFILE *om, *om2, *om3;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * 302);
  float        sc1, sc2, sc3;
  int          i;

  if (p7_oprofile_Sample(r, abc, bg, 77, 300, NULL, NULL, &om) != eslOK) esl_fatal(msg);
  if (p7_oprofile_SetSIMD(om, p7_SIMD_AVX2) != eslOK)                    esl_fatal(msg);
  if ((om2 = p7_oprofile_Copy(om))  == NULL || om2->simd != p7_SIMD_AVX2 || om2->simd_mem == om->simd_mem) esl_fatal(msg);
  if ((om3 = p7_oprofile_Clone(om)) == NULL || om3->simd != p7_SIMD_AVX2 || om3->simd_mem != om->simd_mem) esl_fatal(msg);
  for (i = 0; i < 20; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, 300, dsq);
      sc1 = sc2 = sc3 = 0.0f;
      if (p7_SSVFilter(dsq, 300, om, &sc1) != p7_SSVFilter(dsq, 300, om2, &sc2) || sc1 != sc2) esl_fatal(msg);
      if (p7_SSVFilter(dsq, 300, om, &sc1) != p7_SSVFilter(dsq, 300, om3, &sc3) || sc1 != sc3) esl_fatal(msg);
    }
  p7_oprofile_Destroy(om3);
  p7_oprofile_Destroy(om);
  if (p7_SSVFilter(dsq, 300, om2, &sc2) != p7_SSVFilter(dsq, 300, om2, &sc2)) esl_fatal(msg); /* om2 still usable (ASan) */
  p7_oprofile_Destroy(om2);
  free(dsq);
}
#endif /*p7SIMD_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/



/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7SIMD_TESTDRIVE
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-N",        eslARG_INT,      "8", NULL, "n>0", NULL,  NULL, NULL, "number of targets per model",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the run-time choice of SIMD width";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N   = esl_opt_GetInteger(go, "-N");
  int             be_verbose = esl_opt_GetBoolean(go, "-v");
  int             atype[2]   = { eslAMINO, eslDNA };
  int             bigM[]     = { 150, 200, 255, 256, 257, 300, 448, 449, 500, 1000, 2295, 3000 };
  int             Lcfg[]     = { 10, 100, 400, 3000, 100000 };
  ESL_ALPHABET   *abc;
  P7_BG          *bg;
  int             nout[7];
  int             a, M, i, j;

  if (be_verbose) printf("SIMD level: %s (widest here: %s)\n", p7_simd_Name(p7_simd_Level()), p7_simd_Name(simd_cpu_max()));

  utest_parse();
  utest_level();

  if (simd_cpu_max() < p7_SIMD_AVX2)
    { if (be_verbose) printf("no AVX2 here: kernel comparisons skipped\n"); }
  else
    for (a = 0; a < 2; a++)
      {
	if ((abc = esl_alphabet_Create(atype[a])) == NULL) esl_fatal("failed to create alphabet");
	if ((bg  = p7_bg_Create(abc))             == NULL) esl_fatal("failed to create null model");
	esl_vec_ISet(nout, 7, 0);

	utest_copies(r, abc, bg);
	for (M = 1; M <= 130; M++)                                 utest_filters(r, abc, bg, M,       Lcfg[M % 3],            N, nout);
	for (i = 0; i < (int) (sizeof(bigM)/sizeof(int)); i++)
	  for (j = 0; j < (int) (sizeof(Lcfg)/sizeof(int)); j++)            utest_filters(r, abc, bg, bigM[i], Lcfg[j], ESL_MAX(6, N/(j+1)), nout);

	if (be_verbose) printf("%s SSV: %d ok, %d overflow, %d J state, %d possible overflow, %d untrusted; full MSV: %d ok, %d overflow\n",
			       esl_abc_DecodeType(atype[a]), nout[0], nout[1], nout[2], nout[3], nout[4], nout[5], nout[6]);
	if (nout[0] == 0 || nout[1] == 0 || nout[2] == 0 || nout[5] == 0 || nout[6] == 0) esl_fatal("simd unit test didn't reach every outcome");

	esl_vec_ISet(nout, 2, 0);
	for (i = 0; i < 12; i++) utest_msv_boundary(r, abc, bg, 100 + 25 * i, 10000, nout);
	if (be_verbose) printf("%s full MSV at the overflow boundary: %d targets, %d overflowed\n", esl_abc_DecodeType(atype[a]), nout[0], nout[1]);
	if (nout[0] == 0 || nout[1] == 0 || nout[1] == nout[0]) esl_fatal("simd MSV boundary unit test didn't reach both outcomes");

	p7_bg_Destroy(bg);
	esl_alphabet_Destroy(abc);
      }

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7SIMD_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/
