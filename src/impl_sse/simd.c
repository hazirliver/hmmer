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
 * Contents:
 *   1. Choosing the level.
 *   2. Unit tests.
 *   3. Test driver.
 */
#include <p7_config.h>

#include <stdlib.h>
#include <string.h>

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
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7SIMD_TESTDRIVE

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
#endif /*p7SIMD_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/



/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef p7SIMD_TESTDRIVE
#include <p7_config.h>

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the run-time choice of SIMD width";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  if (esl_opt_GetBoolean(go, "-v"))
    printf("SIMD level: %s (widest here: %s)\n", p7_simd_Name(p7_simd_Level()), p7_simd_Name(simd_cpu_max()));

  utest_parse();
  utest_level();

  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7SIMD_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/
