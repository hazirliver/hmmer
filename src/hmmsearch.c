/* hmmsearch: search profile HMM(s) against a sequence database.
 *
 * To do:
 *  - in MPI mode, add a check to make sure ncpus >= 2. If 1, then we
 *    only have a master, no workers. See Infernal commit r3972 on the
 *    same point; and same note in hmmscan.c's to do list.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __GLIBC__
#include <malloc.h>     /* mallopt(): see worker_Trim() */
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HMMER_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif 

#include <unistd.h>       /* sysconf(): default --dbcache_max */
#if defined(__linux__)
#include <sys/mman.h>     /* madvise(MADV_HUGEPAGE) for the --dbcache arena */
#endif

#ifdef HMMER_THREADS
#include <unistd.h>
#include <pthread.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif 

#include "hmmer.h"

#ifdef HMMER_THREADS
/* Worker threads are created once and search one pass over the target
 * db per query. Between passes they wait at this gate. To start a pass,
 * the master hands its number to as many workers as there will be blocks
 * to search (known from the previous pass), and wakes only those; at
 * the end it sets <quit>.
 */
/* Query-parallel mode: with a small target db wholly cached, one worker
 * searches all targets for a query, while other workers search other
 * queries; there is no per-query barrier. The master prints the results
 * in query order. A job holds one query and its results. The mode is
 * decided by the db size alone, so all queries after the first (which
 * fills the cache) use it or none do: switching modes would make each
 * target-parallel query wait for the slowest query-parallel one.
 */
typedef struct {
  int               nquery;      /* query number, for the output                         */
  P7_HMM           *hmm;
  P7_PROFILE       *gm;
  P7_OPROFILE      *om;          /* used (and length-reconfigured) by the job's worker    */
  P7_TOPHITS       *th;
  P7_PIPELINE       stats;       /* copy of the worker's pipeline after the search; its
                                    pointers are NULL: the output only reads its scalars */
  ESL_STOPWATCH    *w;
  int               status;      /* eslOK, or eslEINVAL from p7_pli_NewModel()            */
  char              errbuf[eslERRBUFSIZE];
  int               done;        /* TRUE: searched, waiting for the master to print it    */
} QP_JOB;

typedef struct {
  pthread_mutex_t   mutex;
  pthread_cond_t   *go;          /* [0..nworkers-1] master -> worker: new pass, or quit */
  int              *wpass;       /* [0..nworkers-1] pass each worker is to run          */
  pthread_cond_t    done;        /* workers -> master: the pass is finished             */
  int               nworkers;
  int               pass;        /* number of the current pass                          */
  int               nrunning;    /* workers still busy in the current pass              */
  int               nblocks;     /* nonempty blocks in the last pass; 0 = unknown yet   */
  int               nseqs;       /* targets in the last pass; 0 = unknown yet           */
  int               split_free;  /* TRUE: how targets are split among workers can't change
                                    which hits a worker keeps (see split_is_free())     */
  int               quit;        /* TRUE: workers exit                                  */
  QP_JOB          **wjob;        /* [0..nworkers-1] query-parallel job given to a worker, or NULL */
  int              *idle;        /* stack of workers without a job                      */
  int               nidle;
  const ESL_GETOPTS *go_opts;    /* for p7_pipeline_Reinit() in query-parallel jobs     */
  struct dbcache_s *dbc;         /* the target db cache (--dbcache), or NULL            */

  /* Streaming (after the first pass, see stream_Start()): workers stay in the
   * work queue and search blocks of two batches of queries in flight. Each
   * pool block is tagged with the batch (slot set 0 or 1) it carries.
   */
  int               stream;      /* TRUE: woken workers run stream_search()             */
  ESL_SQ_BLOCK    **pool;        /* [0..npool-1] the work queue's blocks                */
  int              *tag;         /* [0..npool-1] slot set of the batch a block carries  */
  int               npool;
  int               bnq[2];      /* per slot set: queries in its batch                  */
  int64_t           bid[2];      /*   ... the batch's number                            */
  int               bsent[2];    /*   ... blocks handed out (final once bfin is set)    */
  int               bdone[2];    /*   ... blocks searched                               */
  int               bfin[2];     /*   ... TRUE: all its blocks are handed out           */
} WORKER_GATE;

#define QP_MAXRES 2000000  /* query-parallel mode when the cached target db has up to this many residues */

/* the query-parallel jobs in flight, in query order: job[k % njobs] for head <= k < tail */
typedef struct {
  QP_JOB           *job;
  int               njobs;
  uint64_t          head;        /* next job to print */
  uint64_t          tail;        /* next job to fill  */
} QP_QUEUE;
#endif

#ifdef HMMER_THREADS
/* While the workers search, a few helper threads read the next queries
 * and build their profiles, exactly as the master would have. Reading
 * is serial and in file order (under the mutex); building the profiles
 * runs in parallel. Query number n goes to slot n % nslots, and the
 * master takes them in order.
 */
typedef struct {
  P7_HMM           *hmm;
  P7_PROFILE       *gm;          /* profiles, if hstatus is eslOK               */
  P7_OPROFILE      *om;
  int               hstatus;     /* return status of p7_hmmfile_Read()          */
  int               ready;       /* TRUE: filled in, waiting for the master     */
} PREP_SLOT;

typedef struct {
  pthread_t        *thread;      /* [0..nthreads-1] helper threads              */
  int               nthreads;
  pthread_mutex_t   mutex;
  pthread_cond_t    work;        /* master -> helpers: a slot came free, or quit */
  pthread_cond_t    ready;       /* helpers -> master: a slot is filled         */
  P7_HMMFILE       *hfp;
  ESL_ALPHABET     *abc;
  PREP_SLOT        *slot;        /* [0..nslots-1]                               */
  int               nslots;
  uint64_t          nread;       /* queries read so far                         */
  uint64_t          nused;       /* queries handed to the master                */
  int               done;        /* TRUE: a read returned something but eslOK   */
  int               quit;        /* TRUE: helpers exit                          */
} QUERY_PREP;
#endif

/* One query's state in one worker. A pass over the targets can search
 * several queries (a batch, in slots wq[qbase..qbase+nq-1]): each has its
 * own profile clone (its length is reconfigured per target), null model
 * (the bias filter HMM is set per query), hit list and pipeline. The
 * slots' pipelines are struct copies of <tmpl>, the freshly configured
 * state of the worker's own pipeline <work>: each has its own thresholds
 * and accounting, and they all share <work>'s DP matrices, domain
 * definition workspace and RNG, which only ever hold one target's work
 * (the RNG is reseeded for every region). So a worker needs one such
 * workspace however many queries are in flight, and setting up a slot
 * never touches it. There are two sets of MAXQ slots, for two batches.
 */
/* A target that passed the filters, waiting for a second one to run the
 * Forward/Backward parsers with (pair_feed()): what its later stages need
 * from the time it was filtered.
 */
typedef struct {
  const ESL_SQ     *sq;
  float             nullsc;
  float             filtersc;
  double            Z;           /* the pipeline's running Z after this target's p7_pli_NewSeq() */
  float             xf[p7O_NXSTATES][p7O_NXTRANS];  /* om->xf configured for this target's length */
} PAIR_WAIT;

typedef struct {
  P7_BG            *bg;	         /* null model                              */
  P7_PIPELINE      *pli;         /* pipeline: thresholds and accounting     */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile (a clone)       */
  PAIR_WAIT         wait;        /* paired parsers: a target waiting for a partner, if <nwait> */
  int               nwait;
} WORKER_QUERY;

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
  WORKER_GATE      *gate;
#endif 
  P7_PIPELINE      *work;        /* owns the DP matrices, domain definition workspace and RNG */
  P7_PIPELINE       tmpl;        /* <work> as configured by p7_pipeline_Create(); slots start as copies */
  WORKER_QUERY     *wq;          /* [0..2*MAXQ-1] query slots; pli NULL until first used  */
  int               qbase;       /* a one-batch pass searches slots wq[qbase..qbase+nq-1] */
  int               nq;
  int               reuse_seqs;  /* TRUE: reset each target after use. FALSE: targets are cached, read-only */
  int               pair;        /* TRUE: run the parsers on two targets at once (p7_fb_PairAvailable()) */
  P7_OMX           *pf[2];       /* the pair's Forward and Backward parser matrices, or NULL until used  */
  P7_OMX           *pb[2];
  ESL_SQ          **sqpool;      /* --cpu 0, paired parsers: targets read one at a time, kept while a query holds them */
  int               nsqpool;
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#if defined (HMMER_THREADS) && defined (HMMER_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <f>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "report domains <= this E-value threshold in output",           4 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "report domains >= this score cutoff in output",                4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "consider domains >= this score threshold as significant",      5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                             7 },

/* Other options */
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>: no autodetection",  12 },
  { "--dbcache",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "keep <seqdb> in memory: read it once for all queries",        12 },
  { "--dbcache_max",eslARG_INT,     NULL, NULL, "n>0",   NULL,"--dbcache",NULL,         "memory limit for --dbcache, in MB [default: half of RAM]",   12 },

#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT, p7_NCPU,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,      "number of parallel CPU workers to use for multithreads",      12 },
#endif
#ifdef HMMER_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                              12 },
#endif

  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   * Doesn't work with MPI
   */
  { "--restrictdb_stkey", eslARG_STRING, "0",  NULL, NULL,    NULL,  NULL,  NULL,       "Search starts at the sequence with name <s> (not with MPI)",     99 },
  { "--restrictdb_n",eslARG_INT,        "-1",  NULL, NULL,    NULL,  NULL,  NULL,       "Search <j> target sequences (starting at --restrictdb_stkey)",   99 },
  { "--ssifile",    eslARG_STRING,       NULL, NULL, NULL,    NULL,  NULL,  NULL,       "restrictdb_x values require ssi file. Override default to <s>",  99 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <hmmfile> <seqdb>";
static char banner[] = "search profile(s) against a sequence database";

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *hmmfile;           /* query HMM file                                  */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};

#define BLOCK_SIZE 1000
#define BLOCK_MIN  20     /* smallest block a small, cached target db is split into */
#define MAXQ       16     /* queries searched per pass over the targets              */
#define SLICES     4      /* cached passes: blocks per worker...                     */
#define SLICES_FREE 16    /*   ... or this many, when splitting can't change results  */
#ifndef DP_KEEP_CELLS
#define DP_KEEP_CELLS (1 << 18) /* a worker keeps its domain matrices up to this many cells (3 MB each) */
#endif
#ifndef MMAP_FLOOR
#define MMAP_FLOOR  (1 << 20)   /* glibc: blocks from 1 MB up get their own mapping, freed to the system */
#endif

/* DBCACHE: the target database kept in memory (--dbcache).
 *
 * The first query's pass reads the (possibly --restrictdb'ed) database
 * into blocks that are kept. Later queries are served from those blocks
 * and parse nothing, so a non-rewindable target (.gz, stdin) also works
 * with several queries.
 *
 * The first pass reads exactly as without a cache, into the usual reused
 * block, and appends a compact copy of each block to the cache: one
 * ESL_SQ array per block, with all strings and digital sequences packed
 * into large arena chunks instead of Easel's per-field growable buffers.
 * Cached sequences are never written again.
 *
 * Threaded mode keeps the work queue's protocol unchanged. The same pool
 * of blocks circulates between reader and workers. In a cached pass each
 * pool block acts as a carrier: while it's in flight, its <list> points at
 * a cached block's sequences. carrier_list[] holds each pool block's own
 * list, which is restored at the end of every pass.
 */
#define DBCACHE_CHUNK (64 * 1024 * 1024)  /* arena chunk size (bytes); a multiple of 2MB      */
#define DBCACHE_ALIGN (2 * 1024 * 1024)   /* chunk alignment, so chunks can use huge pages   */

typedef struct dbcache_s {
  ESL_SQ_BLOCK  **blk;          /* cached blocks, in file order                         */
  int             nblk;
  int             nalloc;
  int             complete;     /* TRUE once a pass has read the whole (restricted) db  */
  int             next;         /* in a cached pass: next block to hand out             */
  int             pass_eof;     /* the current pass has reached its end                 */
  int64_t         nseq;         /* number of sequences cached                           */
  int64_t         nres;         /* number of residues cached                            */
  uint64_t        nbytes;       /* bytes of the arena in use                            */
  uint64_t        maxbytes;     /* limit (--dbcache_max); above it the cache is dropped */
  int             dropped;      /* TRUE: over the limit; targets are read per query     */

  char          **chunk;        /* arena: strings and digital sequences of cached seqs  */
  int             nchunk;
  int             nchunkalloc;
  size_t          used;         /* bytes used in the last chunk                         */
  size_t          chunksize;    /* size of the last chunk                               */

  ESL_SQ_BLOCK   *scratch;      /* serial mode: the reused block the first pass reads into */
  ESL_SQ_BLOCK  **carrier;      /* threaded mode: the work queue's pool blocks          */
  ESL_SQ        **carrier_list; /*   ... and their own sequence lists                   */
  int             ncarrier;
  const ESL_ALPHABET *abc;
} DBCACHE;

static DBCACHE *dbcache_Create(const ESL_ALPHABET *abc, int ncarrier, uint64_t maxbytes);
static uint64_t dbcache_DefaultMax(void);
static void     dbcache_StartPass(DBCACHE *dbc);
static void     dbcache_EndPass(DBCACHE *dbc);
static int      dbcache_NextBlock(DBCACHE *dbc, ESL_SQFILE *dbfp, int *n_targetseqs, ESL_SQ_BLOCK *rdblk, ESL_SQ_BLOCK **ret_blk);
static void     dbcache_Destroy(DBCACHE *dbc);

/* the output streams for a search's results */
typedef struct {
  FILE *ofp;         /* main output                         */
  FILE *afp;         /* -A alignments, or NULL              */
  FILE *tblfp;       /* --tblout, or NULL                   */
  FILE *domtblfp;    /* --domtblout, or NULL                */
  FILE *pfamtblfp;   /* --pfamtblout, or NULL               */
} OUTPUTS;

static int  output_query_header(FILE *ofp, const P7_HMM *hmm);
static int  output_results(const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw, const P7_HMM *hmm, int nquery,
			   P7_TOPHITS *th, P7_PIPELINE *pli, ESL_STOPWATCH *w);
#ifdef HMMER_THREADS
static void qp_Create    (QP_QUEUE *qp, int njobs);
static int  qp_Submit    (QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw,
			  P7_HMM *hmm, int nquery, P7_PROFILE **next_gm, P7_OPROFILE **next_om, const P7_BG *bg);
static int  qp_OutputNext(QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw);
static int  qp_Drain     (QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw);
static void qp_Destroy   (QP_QUEUE *qp);
#endif
static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc);

#ifdef HMMER_THREADS
static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, WORKER_GATE *gate, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc, int nq);
static void pipeline_thread(void *arg);
static void search_pass(WORKER_INFO *info);
static int  split_is_free(const ESL_GETOPTS *go);
static int  worker_SetQuery(WORKER_INFO *info, int k, const P7_OPROFILE *om, const ESL_GETOPTS *go, const ESL_ALPHABET *abc);
static void worker_Fresh(WORKER_INFO *info);
static void worker_Trim(WORKER_INFO *info, int M);
static void worker_Search(WORKER_INFO *info, int base, int nq, ESL_SQ *dbsq, int n);
static void pair_feed (WORKER_INFO *info, WORKER_QUERY *q, const ESL_SQ *sq);
static void pair_flush(WORKER_INFO *info, WORKER_QUERY *q);
static int  output_batch(WORKER_INFO *info, int infocnt, int base, int nb, P7_HMM **bhmm, P7_PROFILE **bgm, P7_OPROFILE **bom, int *bnum, ESL_STOPWATCH **w,
                         const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw);
#ifdef HMMER_THREADS
static int  pass_blocksize(const WORKER_GATE *gate, const DBCACHE *dbc, int nq);
static void stream_Start  (ESL_WORK_QUEUE *queue, WORKER_GATE *gate, void **ret_block);
static int  stream_HandOut(ESL_WORK_QUEUE *queue, WORKER_GATE *gate, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc, int set, int nq, int64_t bid, void **io_block);
static void stream_Wait   (WORKER_GATE *gate, int set);
static void stream_Stop   (ESL_WORK_QUEUE *queue, WORKER_GATE *gate, DBCACHE *dbc, void **io_block);
static void stream_search (WORKER_INFO *info);
#endif
static void search_query(WORKER_INFO *info, WORKER_GATE *gate, QP_JOB *job);
static void query_prep_Start(QUERY_PREP *prep, P7_HMMFILE *hfp, ESL_ALPHABET *abc, int nthreads);
static void query_prep_Next (QUERY_PREP *prep, P7_HMM **ret_hmm, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om, int *ret_hstatus);
static void query_prep_Stop (QUERY_PREP *prep);
static void query_prep_thread(void *arg);
#endif 

#ifdef HMMER_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif 


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions controlling model-specific thresholding:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0 && strcmp(*ret_seqfile, "-") == 0) 
    { if (puts("Either <hmmfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere most common options are:")                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")           && fprintf(ofp, "# output directed to file:         %s\n",             esl_opt_GetString(go, "-o"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-A")           && fprintf(ofp, "# MSA of all hits saved to file:   %s\n",             esl_opt_GetString(go, "-A"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")     && fprintf(ofp, "# per-seq hits tabular output:     %s\n",             esl_opt_GetString(go, "--tblout"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domtblout")  && fprintf(ofp, "# per-dom hits tabular output:     %s\n",             esl_opt_GetString(go, "--domtblout"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pfamtblout") && fprintf(ofp, "# pfam-style tabular hit output:   %s\n",             esl_opt_GetString(go, "--pfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")        && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")      && fprintf(ofp, "# show alignments in output:       no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")    && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")      && fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")           && fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "-E"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")           && fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(go, "-T"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domE")       && fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--domE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domT")       && fprintf(ofp, "# domain reporting threshold:      score >= %g\n",    esl_opt_GetReal(go, "--domT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")       && fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "--incE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")       && fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal(go, "--incT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomE")    && fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--incdomE"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomT")    && fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",    esl_opt_GetReal(go, "--incdomT"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")     && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_nc")     && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_tc")     && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")        && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")         && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F1"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")         && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F2"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")         && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F3"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")     && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") && fprintf(ofp, "# Restrict db to start at seq key: %s\n",            esl_opt_GetString(go, "--restrictdb_stkey"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_n")     && fprintf(ofp, "# Restrict db to # target seqs:    %d\n",            esl_opt_GetInteger(go, "--restrictdb_n")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ssifile")          && fprintf(ofp, "# Override ssi file to:            %s\n",            esl_opt_GetString(go, "--ssifile"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# sequence search space set to:    %.0f\n",           esl_opt_GetReal(go, "-Z"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")       && fprintf(ofp, "# domain search space set to:      %.0f\n",           esl_opt_GetReal(go, "--domZ"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                               fprintf(ofp, "# random number seed set to:       %d\n",             esl_opt_GetInteger(go, "--seed"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--tformat")    && fprintf(ofp, "# targ <seqfile> format asserted:  %s\n",             esl_opt_GetString(go, "--tformat"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dbcache")    && fprintf(ofp, "# target db cached in memory:      yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dbcache_max")&& fprintf(ofp, "# target db cache limit:           %d MB\n",         esl_opt_GetInteger(go, "--dbcache_max")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")        && fprintf(ofp, "# number of worker threads:        %d\n",             esl_opt_GetInteger(go, "--cpu"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
#ifdef HMMER_MPI
  if (esl_opt_IsUsed(go, "--mpi")        && fprintf(ofp, "# MPI:                             on\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;	
  struct cfg_s     cfg;        
  int              status   = eslOK;

  impl_Init();                  /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

/* is the range restricted? */
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") )
    if ((cfg.firstseq_key = esl_opt_GetString(go, "--restrictdb_stkey")) == NULL)  p7_Fail("Failure capturing --restrictdb_stkey\n");

  if (esl_opt_IsUsed(go, "--restrictdb_n") )
    cfg.n_targetseq = esl_opt_GetInteger(go, "--restrictdb_n");

  if ( cfg.n_targetseq != -1 && cfg.n_targetseq < 1 )
    p7_Fail("--restrictdb_n must be >= 1\n");


  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HMMER_MPI
  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HMMER_MPI*/
    {
      status = serial_master(go, &cfg);
    }

  esl_getopts_Destroy(go);

  return status;
}


/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().  We also use the
 * ESL_EXCEPTION and ERROR: mechanisms, but only because we know we're
 * using a fatal exception handler.
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam tabular output (--pfamtblout)    */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  P7_BG           *bg       = NULL;              /* null model, for building profiles here          */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w[2][MAXQ];                    /* per slot set and query of a batch: elapsed time */
  P7_HMM          *bhmm[2][MAXQ];                 /* the queries of a batch, and their profiles      */
  P7_PROFILE      *bgm[2][MAXQ];
  P7_OPROFILE     *bom[2][MAXQ];
  int              bnum[2][MAXQ];                 /* their query numbers                             */
  int              nbs[2]   = { 0, 0 };           /* queries in each slot set's batch                */
  int              nb;                            /* queries in this batch                           */
  int              maxb;                          /* most queries this batch may take                */
  int              set      = 0;                  /* slot set of this batch                          */
  int              pending  = -1;                 /* streaming: slot set of the batch handed out
                                                     before this one and not printed yet, or -1     */
  int64_t          nbatch   = 0;                  /* streaming: batches handed out                   */
  int              streaming = FALSE;
  void            *sblock   = NULL;               /* streaming: the empty block the master holds     */
  int              textw    = 0;
  OUTPUTS          outs;                          /* the output streams, for output_results()       */
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i, k;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
  DBCACHE         *dbc      = NULL;              /* in-memory target db (--dbcache), or NULL        */
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
  WORKER_GATE      gate;
  QUERY_PREP       prep;
  int              prep_started = FALSE;
  QP_QUEUE         qp;
#endif
  P7_PROFILE      *next_gm  = NULL;              /* profiles of the next query, if prefetched       */
  P7_OPROFILE     *next_om  = NULL;
  P7_HMM          *hmm_err  = NULL;              /* a query the pipeline can't use (p7_pli_NewModel()) */
  char             errbuf[eslERRBUFSIZE];

  for (k = 0; k < 2*MAXQ; k++) w[k/MAXQ][k%MAXQ] = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  


  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
    if (esl_opt_IsUsed(go, "--ssifile"))
      esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
    else
      esl_sqfile_OpenSSI(dbfp, NULL);
  }



  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout")); }
  if (esl_opt_IsOn(go, "--pfamtblout")){ if ((pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)  esl_fatal("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout")); }
  outs.ofp       = ofp;
  outs.afp       = afp;
  outs.tblfp     = tblfp;
  outs.domtblfp  = domtblfp;
  outs.pfamtblfp = pfamtblfp;

#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN( esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
  if (ncpus > 0)
    {
#ifdef __GLIBC__
      if (MMAP_FLOOR > 0)
	{ /* trimmed matrices go back to the system: see worker_Trim() */
	  mallopt(M_MMAP_THRESHOLD, MMAP_FLOOR);
	  mallopt(M_TRIM_THRESHOLD, 2 * MMAP_FLOOR);
	}
#endif
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
      ESL_ALLOC(gate.go,    sizeof(pthread_cond_t) * ncpus);
      ESL_ALLOC(gate.wpass, sizeof(int)            * ncpus);
      ESL_ALLOC(gate.wjob,  sizeof(QP_JOB *)       * ncpus);
      ESL_ALLOC(gate.idle,  sizeof(int)            * ncpus);
      if (pthread_mutex_init(&gate.mutex, NULL) != 0 || pthread_cond_init(&gate.done, NULL) != 0)
        esl_fatal("Failed to initialize the worker thread gate");
      for (i = 0; i < ncpus; i++)
        {
          if (pthread_cond_init(&gate.go[i], NULL) != 0) esl_fatal("Failed to initialize the worker thread gate");
          gate.wpass[i] = 0;
          gate.wjob[i]  = NULL;
          gate.idle[i]  = ncpus - 1 - i;   /* worker 0 on top */
        }
      gate.nidle    = ncpus;
      gate.go_opts  = go;
      gate.dbc      = NULL;                /* set once the cache exists */
      gate.nworkers = ncpus;
      gate.pass     = 0;
      gate.nrunning = 0;
      gate.nblocks  = 0;
      gate.nseqs    = 0;
      gate.split_free = split_is_free(go);
      gate.quit     = FALSE;
      gate.stream   = FALSE;
      gate.npool    = 0;
      ESL_ALLOC(gate.pool, sizeof(ESL_SQ_BLOCK *) * ncpus * 2);
      ESL_ALLOC(gate.tag,  sizeof(int)            * ncpus * 2);
      qp_Create(&qp, 2 * ncpus);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * infocnt);

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      esl_sqfile_SetDigital(dbfp, abc); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks
      bg = p7_bg_Create(abc);           /* p7_ProfileConfig() only reads its background frequencies */

      if (esl_opt_GetBoolean(go, "--dbcache"))
	{
	  uint64_t maxbytes = esl_opt_IsOn(go, "--dbcache_max") ? (uint64_t) esl_opt_GetInteger(go, "--dbcache_max") << 20 : dbcache_DefaultMax();
	  if ((dbc = dbcache_Create(abc, ncpus * 2, maxbytes)) == NULL) p7_Fail("Failed to allocate the target sequence cache");
#ifdef HMMER_THREADS
	  if (ncpus > 0) gate.dbc = dbc;
#endif
	}

      for (i = 0; i < infocnt; ++i)
	{
	  ESL_ALLOC(info[i].wq, sizeof(WORKER_QUERY) * 2 * MAXQ);
	  for (k = 0; k < 2*MAXQ; k++) info[i].wq[k].pli = NULL;   /* made as they're needed (worker_SetQuery()) */
	  info[i].work       = NULL;
	  info[i].qbase      = 0;
	  info[i].nq         = 0;
	  info[i].pair       = p7_fb_PairAvailable();
	  info[i].pf[0] = info[i].pf[1] = info[i].pb[0] = info[i].pb[1] = NULL;
	  info[i].sqpool     = NULL;
	  info[i].nsqpool    = 0;
	  info[i].reuse_seqs = (dbc == NULL);
#ifdef HMMER_THREADS
	  info[i].queue = queue;
	  info[i].gate  = &gate;
#endif
	}

#ifdef HMMER_THREADS
      for (i = 0; i < ncpus * 2; ++i)
	{
	  block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
	  if (block == NULL) 	      esl_fatal("Failed to allocate sequence block");

 	  status = esl_workqueue_Init(queue, block);
	  if (status != eslOK)	      esl_fatal("Failed to add block to work queue");
	  gate.pool[gate.npool]  = block;
	  gate.tag[gate.npool++] = -1;

	  if (dbc) { /* pool blocks may carry cached sequences; remember their own lists */
	    dbc->carrier[dbc->ncarrier]      = block;
	    dbc->carrier_list[dbc->ncarrier] = block->list;
	    dbc->ncarrier++;
	  }
	}

      if (ncpus > 0)
	{ /* helpers that read and prepare the next queries during the searches */
	  query_prep_Start(&prep, hfp, abc, (ncpus + 1) / 2);
	  prep_started = TRUE;
	}
#endif
    }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      nquery++;

#ifdef HMMER_THREADS
      /* With a small target db wholly cached, each query is searched by
       * one worker, while other workers search the next queries; results
       * are printed in query order. Otherwise the targets of each query are
       * split across all workers (below).
       */
      if (ncpus > 0 && dbc && dbc->complete && dbc->nres <= QP_MAXRES)
	{
	  if ((status = qp_Submit(&qp, &gate, &outs, go, abc, textw, hmm, nquery, &next_gm, &next_om, bg)) != eslOK) return status;
	  query_prep_Next(&prep, &hmm, &next_gm, &next_om, &hstatus);
	  continue;
	}
      if (ncpus > 0 && (status = qp_Drain(&qp, &gate, &outs, go, abc, textw)) != eslOK) return status;
#endif

      /* seqfile may need to be rewound (multiquery mode), unless every target is cached */
      if (nquery > 1 && ! (dbc && dbc->complete))
      {
        if (! esl_sqfile_IsRewindable(dbfp))
        {
          if (dbc && dbc->dropped)
            esl_fatal("Target sequence file %s isn't rewindable and doesn't fit in the --dbcache limit (%" PRIu64 " MB); can't search it with multiple queries", cfg->dbfile, dbc->maxbytes >> 20);
          esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
        }

        if (! esl_opt_IsUsed(go, "--restrictdb_stkey") )
          esl_sqfile_Position(dbfp, 0); //only re-set current position to 0 if we're not planning to set it in a moment
      }

      if ( cfg->firstseq_key != NULL && ! (dbc && dbc->complete)) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
        sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
        if (sstatus != eslOK)
          p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
      }

      /* One pass over the targets searches a batch of queries: this one
       * and, after the first pass, up to MAXQ-1 more that follow it in the
       * file (the first pass is alone: it learns the number of targets,
       * fills the --dbcache, and a target file that can't be rewound fails
       * at the second query as before). Every block of targets is then read (or
       * taken from the --dbcache) once for the whole batch, and each
       * worker searches it with each query. Per query, the targets and the
       * pipeline calls are those of a one-query pass, so the results don't
       * change; the queries' outputs are printed in order after the pass.
       * A query the pipeline can't use (p7_pli_NewModel() fails, e.g.
       * --cut_ga without GA lines) ends the batch; it's reported after
       * the queries before it, as a one-query pass would.
       */
      maxb    = (nquery > 1) ? MAXQ : 1;
      nb      = 0;
      hmm_err = NULL;
      set     = (ncpus > 0 && nquery > 1) ? (int) ((nbatch + 1) & 1) : 0;   /* streaming batches alternate slot sets */
      for (;;)
	{
	  /* Convert to an optimized model, unless the prefetch thread already did */
	  if (next_om)
	    {
	      gm = next_gm;  next_gm = NULL;
	      om = next_om;  next_om = NULL;
	    }
	  else
	    {
	      gm = p7_profile_Create (hmm->M, abc);
	      om = p7_oprofile_Create(hmm->M, abc);
	      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
	      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */
	    }
	  esl_stopwatch_Start(w[set][nb]);

	  /* the query's pipeline, hit list and profile clone in each worker */
	  for (i = 0; i < infocnt; ++i)
	    if ((status = worker_SetQuery(&info[i], set*MAXQ + nb, om, go, abc)) == eslEINVAL) break;
	  if (status == eslEINVAL)
	    {
	      for (k = 0; k <= i && k < infocnt; k++) p7_oprofile_Destroy(info[k].wq[set*MAXQ + nb].om);
	      hmm_err = hmm;
	      p7_oprofile_Destroy(om);
	      p7_profile_Destroy(gm);
	      break;
	    }
	  bhmm[set][nb] = hmm;
	  bgm[set][nb]  = gm;
	  bom[set][nb]  = om;
	  bnum[set][nb] = nquery;
	  nb++;
	  if (nb == maxb) break;
#ifdef HMMER_THREADS
	  if (ncpus > 0) query_prep_Next(&prep, &hmm, &next_gm, &next_om, &hstatus);
	  else           hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
#else
	  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
#endif
	  if (hstatus != eslOK) break;
	  nquery++;
	}

      nbs[set] = nb;
      if (nb > 0)
	{
	  /* (while streaming, the workers' fields stay as stream_Start() found them) */
	  if (! streaming)
	    for (i = 0; i < infocnt; ++i)
	      {
		info[i].qbase = set * MAXQ;
		info[i].nq    = nb;
		/* workers reset the reused blocks they read into, but never cached targets */
		info[i].reuse_seqs = (dbc == NULL || ! dbc->complete);
#ifdef HMMER_THREADS
		if (ncpus > 0 && bnum[set][0] == 1) esl_threads_AddThread(threadObj, &info[i]); /* workers persist across queries */
#endif
	      }

#ifdef HMMER_THREADS
	  if (ncpus > 0 && bnum[set][0] == 1) esl_threads_WaitForStart(threadObj);
	  if (ncpus > 0 && bnum[set][0] > 1)
	    { /* after the first pass: hand this batch out, then print the one before it (see stream_Start()) */
	      if (! streaming) { stream_Start(queue, &gate, &sblock); streaming = TRUE; }
	      sstatus = stream_HandOut(queue, &gate, dbfp, cfg->n_targetseq, dbc, set, nb, ++nbatch, &sblock);
	      if (pending >= 0)
		{
		  stream_Wait(&gate, pending);
		  if ((status = output_batch(info, infocnt, pending*MAXQ, nbs[pending], bhmm[pending], bgm[pending], bom[pending], bnum[pending], w[pending],
					     &outs, go, abc, textw)) != eslOK) return status;
		}
	      pending = set;
	    }
	  else if (ncpus > 0) sstatus = thread_loop(threadObj, queue, &gate, dbfp, cfg->n_targetseq, dbc, nb);
	  else                sstatus = serial_loop(info, dbfp, cfg->n_targetseq, dbc);
#else
	  sstatus = serial_loop(info, dbfp, cfg->n_targetseq, dbc);
#endif
	  /* a read error: as a one-query pass would, the failing pass's query
	   * header comes after the queries before it, then the error
	   */
	  if (sstatus != eslEOF && (status = output_query_header(ofp, bhmm[set][0])) != eslOK) return status;
	  switch(sstatus)
	  {
	  case eslEFORMAT:
	    esl_fatal("Parse failed (sequence file %s):\n%s\n",
		dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
	    break;
	  case eslEOF:
	    /* do nothing */
	    break;
	  default:
	    esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
	  }
	}

      /* a pass that ran here (not streamed): print its batch now */
      if (pending != set && (status = output_batch(info, infocnt, set*MAXQ, nb, bhmm[set], bgm[set], bom[set], bnum[set], w[set], &outs, go, abc, textw)) != eslOK) return status;

      if (hmm_err)
	{
#ifdef HMMER_THREADS
	  if (pending >= 0)
	    {
	      stream_Wait(&gate, pending);
	      if ((status = output_batch(info, infocnt, pending*MAXQ, nbs[pending], bhmm[pending], bgm[pending], bom[pending], bnum[pending], w[pending],
					 &outs, go, abc, textw)) != eslOK) return status;
	    }
#endif
	  if ((status = output_query_header(ofp, hmm_err)) != eslOK) return status;
	  p7_Fail(info->wq[set*MAXQ + nb].pli->errbuf);
	}

      if (nb == maxb)   /* else the batch ended at a read error or EOF, in <hstatus> */
	{
#ifdef HMMER_THREADS
	  if (ncpus > 0) query_prep_Next(&prep, &hmm, &next_gm, &next_om, &hstatus);
	  else           hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
#else
	  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
#endif
	}
    } /* end outer loop over query HMMs */

#ifdef HMMER_THREADS
  if (pending >= 0)
    {
      stream_Wait(&gate, pending);
      if ((status = output_batch(info, infocnt, pending*MAXQ, nbs[pending], bhmm[pending], bgm[pending], bom[pending], bnum[pending], w[pending],
				 &outs, go, abc, textw)) != eslOK) return status;
      pending = -1;
    }
  if (streaming) stream_Stop(queue, &gate, dbc, &sblock);
  if (ncpus > 0 && (status = qp_Drain(&qp, &gate, &outs, go, abc, textw)) != eslOK) return status;   /* before any read error */
#endif

  switch(hstatus) {
  case eslEOD:       p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
  case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);      break;
  case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
  case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
  default:           p7_Fail("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
  }


  /* Terminate outputs... any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (pfamtblfp) p7_tophits_TabularTail(pfamtblfp,"hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit
   */
  for (i = 0; i < infocnt; ++i)
    {
      for (k = 0; k < 2*MAXQ; k++)
	if (info[i].wq[k].pli)
	  {
	    free(info[i].wq[k].pli);   /* a struct copy; the workspace is <work>'s */
	    p7_tophits_Destroy(info[i].wq[k].th);
	    p7_bg_Destroy(info[i].wq[k].bg);
	  }
      p7_pipeline_Destroy(info[i].work);
      for (k = 0; k < 2; k++) { p7_omx_Destroy(info[i].pf[k]); p7_omx_Destroy(info[i].pb[k]); }
      for (k = 0; k < info[i].nsqpool; k++) esl_sq_Destroy(info[i].sqpool[k]);
      free(info[i].sqpool);
      free(info[i].wq);
    }

#ifdef HMMER_THREADS
  if (prep_started) query_prep_Stop(&prep);

  if (ncpus > 0)
    {
      pthread_mutex_lock(&gate.mutex);   /* release the parked workers, and wait for them to exit */
      gate.quit = TRUE;
      for (i = 0; i < gate.nworkers; i++) pthread_cond_signal(&gate.go[i]);
      pthread_mutex_unlock(&gate.mutex);
      esl_threads_WaitForFinish(threadObj);
      for (i = 0; i < gate.nworkers; i++) pthread_cond_destroy(&gate.go[i]);
      pthread_cond_destroy(&gate.done);
      pthread_mutex_destroy(&gate.mutex);
      free(gate.go);
      free(gate.wpass);
      free(gate.wjob);
      free(gate.idle);
      free(gate.pool);
      free(gate.tag);
      qp_Destroy(&qp);

      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
	esl_sq_DestroyBlock(block);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  dbcache_Destroy(dbc);
  p7_bg_Destroy(bg);
  free(info);
  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  for (k = 0; k < 2*MAXQ; k++) esl_stopwatch_Destroy(w[k/MAXQ][k%MAXQ]);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

#ifdef HMMER_MPI

/* Define common tags used by the MPI master/slave processes */
#define HMMER_ERROR_TAG          1
#define HMMER_HMM_TAG            2
#define HMMER_SEQUENCE_TAG       3
#define HMMER_BLOCK_TAG          4
#define HMMER_PIPELINE_TAG       5
#define HMMER_TOPHITS_TAG        6
#define HMMER_HIT_TAG            7
#define HMMER_TERMINATING_TAG    8
#define HMMER_READY_TAG          9

/* mpi_failure()
 * Generate an error message.  If the clients rank is not 0, a
 * message is created with the error message and sent to the
 * master process for handling.
 */
static void
mpi_failure(char *format, ...)
{
  va_list  argp;
  int      status = eslFAIL;
  int      len;
  int      rank;
  char     str[512];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* format the error mesg */
  va_start(argp, format);
  len = vsnprintf(str, sizeof(str), format, argp);
  va_end(argp);

  /* make sure the error string is terminated */
  str[sizeof(str)-1] = '\0';

  /* if the caller is the master, print the results and abort */
  if (rank == 0)
    {
      if (fprintf(stderr, "\nError: ") < 0) exit(eslEWRITE);
      if (fprintf(stderr, "%s", str)   < 0) exit(eslEWRITE);
      if (fprintf(stderr, "\n")        < 0) exit(eslEWRITE);
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD, status);
      exit(1);
    }
  else
    {
      MPI_Send(str, len, MPI_CHAR, 0, HMMER_ERROR_TAG, MPI_COMM_WORLD);
      pause();
    }
}

#define MAX_BLOCK_SIZE (512*1024)

typedef struct {
  uint64_t  offset;
  uint64_t  length;
  uint64_t  count;
} SEQ_BLOCK;

typedef struct {
  int        complete;
  int        size;
  int        current;
  int        last;
  SEQ_BLOCK *blocks;
} BLOCK_LIST;

/* this routine parses the database keeping track of the blocks
 * offset within the file, number of sequences and the length
 * of the block.  These blocks are passed as work units to the
 * MPI workers.  If multiple hmm's are in the query file, the
 * blocks are reused without parsing the database a second time.
 */
int next_block(ESL_SQFILE *sqfp, ESL_SQ *sq, BLOCK_LIST *list, SEQ_BLOCK *block, int n_targetseqs)
{
  int      status   = eslOK;

  /* if the list has been calculated, use it instead of parsing the database */
  if (list->complete)
    {
      if (list->current == list->last)
      {
        block->offset = 0;
        block->length = 0;
        block->count  = 0;

        status = eslEOF;
      }
      else
      {
        int inx = list->current++;

        block->offset = list->blocks[inx].offset;
        block->length = list->blocks[inx].length;
        block->count  = list->blocks[inx].count;

        status = eslOK;
      }

      return status;
    }

  block->offset = 0;
  block->length = 0;
  block->count = 0;

  esl_sq_Reuse(sq);
  if (n_targetseqs == 0) status = eslEOF; //this is to handle the end-case of a restrictdb scenario, where no more targets are required, and we want to mark the list as complete
  while (block->length < MAX_BLOCK_SIZE && (n_targetseqs <0 || block->count < n_targetseqs) && (status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      if (block->count == 0) block->offset = sq->roff;
      block->length = sq->eoff - block->offset + 1;
      block->count++;
      esl_sq_Reuse(sq);
    }

  if (block->count > 0)
    if (status == eslEOF || block->count == n_targetseqs)
      status = eslOK;
  if (status == eslEOF) list->complete = 1;

  /* add the block to the list of known blocks */
  if (status == eslOK)
    {
      int inx;

      if (list->last >= list->size)
	{
	  void *tmp;
	  list->size += 500;
	  ESL_RALLOC(list->blocks, tmp, sizeof(SEQ_BLOCK) * list->size);
	}

      inx = list->last++;
      list->blocks[inx].offset = block->offset;
      list->blocks[inx].length = block->length;
      list->blocks[inx].count  = block->count;
    }

  return status;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of hmmbuild.
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam-style tabular output  (--pfamtblout) */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  BLOCK_LIST      *list     = NULL;
  SEQ_BLOCK        block;

  int              i;
  int              size;
  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  int              n_targets;

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
      if (esl_opt_IsUsed(go, "--ssifile"))
        esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
      else
        esl_sqfile_OpenSSI(dbfp, NULL);
  }


  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o") && (ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o"));

  if (esl_opt_IsOn(go, "-A") && (afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) 
    mpi_failure("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));

  if (esl_opt_IsOn(go, "--tblout") && (tblfp = fopen(esl_opt_GetString(go, "--tblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout"));

  if (esl_opt_IsOn(go, "--domtblout") && (domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout"));

  if (esl_opt_IsOn(go, "--pfamtblout") && (pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)
    mpi_failure("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout"));

  ESL_ALLOC(list, sizeof(BLOCK_LIST));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
    }
  

  if ( cfg->firstseq_key != NULL ) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
    sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
    if (sstatus != eslOK)
      p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
  }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;
      int              seq_cnt = 0;
      nquery++;
      esl_stopwatch_Start(w);

      n_targets = cfg->n_targetseq;

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)   list->current = 0;

      if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
      if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, hmm->M, 100, FALSE, p7_SEARCH_SEQS);
      p7_pli_NewModel(pli, om, bg);

      /* Main loop: */
      while ((n_targets==-1 || seq_cnt<=n_targets) && (sstatus = next_block(dbfp, dbsq, list, &block, n_targets-seq_cnt)) == eslOK )
      {
        seq_cnt += block.count;

        if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0)
          mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

        MPI_Get_count(&mpistatus, MPI_PACKED, &size);
        if (mpi_buf == NULL || size > mpi_size) {
          void *tmp;
          ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
          mpi_size = size;
        }

        dest = mpistatus.MPI_SOURCE;
        MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

        if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
          mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
        if (mpistatus.MPI_TAG != HMMER_READY_TAG)
          mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);

        MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
      }

      if (n_targets!=-1 && seq_cnt==n_targets)
        sstatus = eslEOF;

      switch(sstatus)
      {
      case eslEFORMAT:
        mpi_failure("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOF:
        break;
      default:
        mpi_failure("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      block.offset = 0;
      block.length = 0;
      block.count  = 0;

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  P7_PIPELINE     *mpi_pli   = NULL;
	  P7_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  p7_tophits_Merge(th, mpi_th);
	  p7_pipeline_Merge(pli, mpi_pli);

	  p7_pipeline_Destroy(mpi_pli);
	  p7_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      p7_tophits_SortBySortkey(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, th, pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, th, pli, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, hmm->name, hmm->acc, th, pli);

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (p7_tophits_Alignment(th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
	  {
	    esl_msa_SetName     (msa, hmm->name, -1);
	    esl_msa_SetAccession(msa, hmm->acc,  -1);
	    esl_msa_SetDesc     (msa, hmm->desc, -1);
	    esl_msa_FormatAuthor(msa, "hmmsearch (HMMER %s)", HMMER_VERSION);

	    if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
	  } 
	else { if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
	  
	esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus) {
  case eslEOD:       mpi_failure("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
  case eslEFORMAT:   mpi_failure("bad file format in HMM file %s",             cfg->hmmfile);      break;
  case eslEINCOMPAT: mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
  case eslEOF:       /* EOF is good, that's what we expect here */                                 break;
  default:           mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
  }

  /* monitor all the workers to make sure they have ended */
  for (i = 1; i < cfg->nproc; ++i)
    {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

      MPI_Get_count(&mpistatus, MPI_PACKED, &size);
      if (mpi_buf == NULL || size > mpi_size) {
	void *tmp;
	ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	mpi_size = size; 
      }

      dest = mpistatus.MPI_SOURCE;
      MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != HMMER_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

  /* Terminate outputs... any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,     "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp,  "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (pfamtblfp)p7_tophits_TabularTail(pfamtblfp, "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)     { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit
   */
  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);

  return eslOK;

 ERROR:
  return eslEMEM;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures           */
  int              mpi_size = 0;                 /* size of the allocated buffer                    */

  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
      esl_sqfile_SetDigital(dbfp, abc);
    }
  
  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;

      SEQ_BLOCK        block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      p7_pli_NewModel(pli, om, bg);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  status = esl_sqfile_Position(dbfp, block.offset);
	  if (status != eslOK) mpi_failure("Cannot position sequence database to %ld\n", block.offset);

	  while (count > 0 && (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	    {
	      length = dbsq->eoff - block.offset + 1;

	      p7_pli_NewSeq(pli, dbsq);
	      p7_bg_SetLength(bg, dbsq->n);
	      p7_oprofile_ReconfigLength(om, dbsq->n);
      
	      p7_Pipeline(pli, om, bg, dbsq, NULL, th);

	      esl_sq_Reuse(dbsq);
	      p7_pipeline_Reuse(pli);

	      --count;
	    }

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (count > 0)              mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n",  block.count,  block.count - count, block.offset);
	  if (block.length != length) mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length,              block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus)
    {
    case eslEOF:
      /* do nothing */
      break;
    case eslEFORMAT:
      mpi_failure("bad file format in HMM file %s", cfg->hmmfile);
      break;
    case eslEINCOMPAT:
      mpi_failure("HMM file %s contains different alphabets", cfg->hmmfile);
      break;
    default:
      mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
    }

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, HMMER_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  return eslOK;
}
#endif /*HMMER_MPI*/

/* dbcache_*: the --dbcache in-memory target database (see DBCACHE). */
static DBCACHE *
dbcache_Create(const ESL_ALPHABET *abc, int ncarrier, uint64_t maxbytes)
{
  DBCACHE *dbc = NULL;
  int      status;

  ESL_ALLOC(dbc, sizeof(DBCACHE));
  dbc->blk          = NULL;
  dbc->chunk        = NULL;
  dbc->scratch      = NULL;
  dbc->carrier      = NULL;
  dbc->carrier_list = NULL;
  dbc->nalloc       = 64;
  dbc->nblk         = 0;
  dbc->complete     = FALSE;
  dbc->next         = 0;
  dbc->pass_eof     = FALSE;
  dbc->nseq         = 0;
  dbc->nres         = 0;
  dbc->nbytes       = 0;
  dbc->maxbytes     = maxbytes;
  dbc->dropped      = FALSE;
  dbc->nchunk       = 0;
  dbc->nchunkalloc  = 16;
  dbc->used         = 0;
  dbc->chunksize    = 0;
  dbc->ncarrier     = 0;
  dbc->abc          = abc;
  ESL_ALLOC(dbc->blk,          sizeof(ESL_SQ_BLOCK *) * dbc->nalloc);
  ESL_ALLOC(dbc->chunk,        sizeof(char *)         * dbc->nchunkalloc);
  ESL_ALLOC(dbc->carrier,      sizeof(ESL_SQ_BLOCK *) * ESL_MAX(1, ncarrier));
  ESL_ALLOC(dbc->carrier_list, sizeof(ESL_SQ *)       * ESL_MAX(1, ncarrier));
  if (ncarrier == 0 && (dbc->scratch = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc)) == NULL) goto ERROR;
  return dbc;

 ERROR:
  dbcache_Destroy(dbc);
  return NULL;
}

/* dbcache_DefaultMax()
 * Default --dbcache limit: half of the memory available to this process,
 * i.e. physical RAM or, if lower, a cgroup (container) memory limit.
 */
static uint64_t
dbcache_DefaultMax(void)
{
  uint64_t mem = 0;
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  long     npages = sysconf(_SC_PHYS_PAGES);
  long     psize  = sysconf(_SC_PAGESIZE);
  if (npages > 0 && psize > 0) mem = (uint64_t) npages * (uint64_t) psize;
#endif
#if defined(__linux__)
  {
    const char *cgfile[] = { "/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes" };
    unsigned long long lim;
    FILE *fp;
    int   i;
    for (i = 0; i < 2; i++)
      if ((fp = fopen(cgfile[i], "r")) != NULL)
        {
          if (fscanf(fp, "%llu", &lim) == 1 && lim > 0 && (mem == 0 || lim < mem)) mem = lim; /* "max" = no limit */
          fclose(fp);
        }
  }
#endif
  return (mem > 0) ? mem / 2 : (uint64_t) 4 << 30;   /* unknown: 4 GB */
}

/* Over the limit: free what's cached; from now on, targets are read per query. */
static void
dbcache_drop(DBCACHE *dbc)
{
  int i;
  for (i = 0; i < dbc->nchunk; i++) free(dbc->chunk[i]);
  dbc->nchunk  = 0;
  dbc->nblk    = 0;
  dbc->nseq    = 0;
  dbc->nres    = 0;
  dbc->used    = dbc->chunksize = 0;
  dbc->nbytes  = 0;
  dbc->dropped = TRUE;
  fprintf(stderr, "Warning: the target database needs more than the --dbcache limit (%" PRIu64 " MB);\n"
                  "         it is read again for each query instead. Raise the limit with --dbcache_max <MB>.\n",
          dbc->maxbytes >> 20);
}

static void
dbcache_StartPass(DBCACHE *dbc)
{
  dbc->next     = 0;
  dbc->pass_eof = FALSE;
}

/* After a threaded pass: give the pool blocks their own lists back. */
static void
dbcache_EndPass(DBCACHE *dbc)
{
  int i;
  for (i = 0; i < dbc->ncarrier; i++)
    {
      dbc->carrier[i]->list  = dbc->carrier_list[i];
      dbc->carrier[i]->count = 0;
    }
}

/* <n> bytes from the arena, aligned to <align> (a power of 2); NULL on allocation failure. */
static void *
dbcache_alloc(DBCACHE *dbc, size_t n, size_t align)
{
  size_t off = (dbc->used + align - 1) & ~(align - 1);
  void  *p;
  int    status;

  if (dbc->nchunk == 0 || off + n > dbc->chunksize)
    {
      if (dbc->nchunk == dbc->nchunkalloc)
        {
          ESL_RALLOC(dbc->chunk, p, sizeof(char *) * dbc->nchunkalloc * 2);
          dbc->nchunkalloc *= 2;
        }
      dbc->chunksize = ESL_MAX(DBCACHE_CHUNK, (n + DBCACHE_ALIGN - 1) & ~((size_t) DBCACHE_ALIGN - 1));
      if (posix_memalign(&p, DBCACHE_ALIGN, dbc->chunksize) != 0) { status = eslEMEM; goto ERROR; }
#if defined(__linux__) && defined(MADV_HUGEPAGE)
      madvise(p, dbc->chunksize, MADV_HUGEPAGE);   /* fewer page faults; no effect on contents */
#endif
      dbc->chunk[dbc->nchunk++] = p;
      off = 0;
    }
  dbc->nbytes += (off + n) - dbc->used;
  dbc->used    = off + n;
  return dbc->chunk[dbc->nchunk-1] + off;

 ERROR:
  return NULL;
}

static char *
dbcache_strdup(DBCACHE *dbc, const char *s, int *ret_alloc)
{
  size_t n = strlen(s) + 1;
  char  *p = dbcache_alloc(dbc, n, 1);
  if (p) { memcpy(p, s, n); *ret_alloc = (int) n; }
  return p;
}

/* Compact copy of digital sequence <src> into <dst>: every field has the
 * same value; strings and residue arrays live in the arena, exactly sized.
 */
static int
dbcache_copyseq(DBCACHE *dbc, const ESL_SQ *src, ESL_SQ *dst)
{
  int x;

  *dst = *src;   /* every scalar field (n, L, start, end, C, W, tax_id, idx, offsets, abc...) */
  dst->salloc = src->n + 2;   /* digital [0..n+1], sentinels included; ss and xr likewise */
  if ((dst->name   = dbcache_strdup(dbc, src->name,   &dst->nalloc))   == NULL) return eslEMEM;
  if ((dst->acc    = dbcache_strdup(dbc, src->acc,    &dst->aalloc))   == NULL) return eslEMEM;
  if ((dst->desc   = dbcache_strdup(dbc, src->desc,   &dst->dalloc))   == NULL) return eslEMEM;
  if ((dst->source = dbcache_strdup(dbc, src->source, &dst->srcalloc)) == NULL) return eslEMEM;
  if ((dst->dsq    = dbcache_alloc(dbc, dst->salloc, 1))               == NULL) return eslEMEM;
  memcpy(dst->dsq, src->dsq, dst->salloc);
  if (src->ss)
    {
      if ((dst->ss = dbcache_alloc(dbc, dst->salloc, 1)) == NULL) return eslEMEM;
      memcpy(dst->ss, src->ss, dst->salloc);
    }
  if (src->nxr > 0)
    {
      int tagalloc;
      if ((dst->xr_tag = dbcache_alloc(dbc, sizeof(char *) * src->nxr, sizeof(char *))) == NULL) return eslEMEM;
      if ((dst->xr     = dbcache_alloc(dbc, sizeof(char *) * src->nxr, sizeof(char *))) == NULL) return eslEMEM;
      for (x = 0; x < src->nxr; x++)
        {
          if ((dst->xr_tag[x] = dbcache_strdup(dbc, src->xr_tag[x], &tagalloc)) == NULL) return eslEMEM;
          if ((dst->xr[x]     = dbcache_alloc(dbc, dst->salloc, 1))            == NULL) return eslEMEM;
          memcpy(dst->xr[x], src->xr[x], dst->salloc);
        }
    }
  return eslOK;
}

/* Append a compact copy of block <src> to the cache. */
static int
dbcache_addblock(DBCACHE *dbc, const ESL_SQ_BLOCK *src)
{
  ESL_SQ_BLOCK *b = NULL;
  void         *p;
  int           i;
  int           status;

  if (dbc->nblk == dbc->nalloc)
    {
      ESL_RALLOC(dbc->blk, p, sizeof(ESL_SQ_BLOCK *) * dbc->nalloc * 2);
      dbc->nalloc *= 2;
    }
  /* the block and its ESL_SQ array live in the arena too */
  if ((b       = dbcache_alloc(dbc, sizeof(ESL_SQ_BLOCK), 16))                      == NULL) return eslEMEM;
  if ((b->list = dbcache_alloc(dbc, sizeof(ESL_SQ) * ESL_MAX(1, src->count), 16)) == NULL) return eslEMEM;
  b->count        = b->listSize = src->count;
  b->complete     = src->complete;
  b->first_seqidx = src->first_seqidx;

  for (i = 0; i < src->count; i++)
    {
      if ((status = dbcache_copyseq(dbc, src->list + i, b->list + i)) != eslOK) return status;
      dbc->nres += src->list[i].n;
    }
  dbc->blk[dbc->nblk++] = b;
  dbc->nseq += src->count;
  return eslOK;

 ERROR:
  return status;
}

/* dbcache_NextBlock()
 * Next block of targets for this pass, in file order. While the cache is
 * incomplete (first pass) it reads into <rdblk> exactly as without a
 * cache (esl_sqio_ReadBlock(), honoring <*n_targetseqs> the same way),
 * keeps a compact copy, and returns <rdblk> itself. In later passes it
 * returns cached blocks. Returns eslOK, eslEOF at the end of the pass, or
 * a read/allocation error.
 */
static int
dbcache_NextBlock(DBCACHE *dbc, ESL_SQFILE *dbfp, int *n_targetseqs, ESL_SQ_BLOCK *rdblk, ESL_SQ_BLOCK **ret_blk)
{
  int status;

  *ret_blk = NULL;
  if (dbc->pass_eof) return eslEOF;

  if (dbc->complete)
    {
      if (dbc->next < dbc->nblk) { *ret_blk = dbc->blk[dbc->next++]; return eslOK; }
      dbc->pass_eof = TRUE;
      return eslEOF;
    }

  if (*n_targetseqs == 0) { rdblk->count = 0; dbc->pass_eof = TRUE; dbc->complete = ! dbc->dropped; return eslEOF; }

  status = esl_sqio_ReadBlock(dbfp, rdblk, -1, *n_targetseqs, /*max_init_window=*/FALSE, FALSE);
  if (status == eslEOF) { dbc->pass_eof = TRUE; dbc->complete = ! dbc->dropped; }
  if (status != eslOK)  return status;

  *n_targetseqs -= rdblk->count;
  if (! dbc->dropped)
    {
      if ((status = dbcache_addblock(dbc, rdblk)) != eslOK) return status;
      if (dbc->nbytes > dbc->maxbytes) dbcache_drop(dbc);
    }
  *ret_blk = rdblk;
  return eslOK;
}

static void
dbcache_Destroy(DBCACHE *dbc)
{
  int i;

  if (dbc == NULL) return;
  for (i = 0; i < dbc->nchunk; i++) free(dbc->chunk[i]);   /* everything cached lives in the arena */
  if (dbc->scratch) esl_sq_DestroyBlock(dbc->scratch);
  free(dbc->blk);
  free(dbc->chunk);
  free(dbc->carrier);
  free(dbc->carrier_list);
  free(dbc);
}


/* output_query_header()
 * The lines that open a query's results.
 */
static int
output_query_header(FILE *ofp, const P7_HMM *hmm)
{
  if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
  if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
  return eslOK;
}

/* output_results()
 * Print one query's results: per-target and per-domain hit lists, the
 * tabular outputs, the pipeline statistics and the -A alignment. Stops
 * the query's stopwatch <w> before the statistics.
 *
 * Throws:    <eslEWRITE> on a write failure.
 */
static int
output_results(const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw, const P7_HMM *hmm, int nquery,
	       P7_TOPHITS *th, P7_PIPELINE *pli, ESL_STOPWATCH *w)
{
  FILE *ofp = outs->ofp;
  FILE *afp = outs->afp;

  p7_tophits_SortBySortkey(th);
  p7_tophits_Threshold(th, pli);
  p7_tophits_Targets(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  p7_tophits_Domains(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (outs->tblfp)     p7_tophits_TabularTargets(outs->tblfp,    hmm->name, hmm->acc, th, pli, (nquery == 1));
  if (outs->domtblfp)  p7_tophits_TabularDomains(outs->domtblfp, hmm->name, hmm->acc, th, pli, (nquery == 1));
  if (outs->pfamtblfp) p7_tophits_TabularXfam(outs->pfamtblfp, hmm->name, hmm->acc, th, pli);

  esl_stopwatch_Stop(w);
  p7_pli_Statistics(ofp, pli, w);
  if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  /* Output the results in an MSA (-A option) */
  if (afp) {
    ESL_MSA *msa = NULL;

    if (p7_tophits_Alignment(th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
      {
	esl_msa_SetName     (msa, hmm->name, -1);
	esl_msa_SetAccession(msa, hmm->acc,  -1);
	esl_msa_SetDesc     (msa, hmm->desc, -1);
	esl_msa_FormatAuthor(msa, "hmmsearch (HMMER %s)", HMMER_VERSION);

	if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);

	if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      }
    else { if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

    esl_msa_Destroy(msa);
  }
  return eslOK;
}

/* serial_sq()
 * --cpu 0 with paired parsers: an ESL_SQ to read the next target into,
 * reset, that no query slot holds as its waiting target (the pool grows
 * to at most one more than the number of queries in a pass).
 */
static ESL_SQ *
serial_sq(WORKER_INFO *info)
{
  void *p;
  int   i, k, held;
  int   status;

  for (i = 0; i < info->nsqpool; i++)
    {
      for (held = FALSE, k = info->qbase; k < info->qbase + info->nq && ! held; k++)
	held = (info->wq[k].nwait && info->wq[k].wait.sq == info->sqpool[i]);
      if (! held) { esl_sq_Reuse(info->sqpool[i]); return info->sqpool[i]; }
    }
  ESL_RALLOC(info->sqpool, p, sizeof(ESL_SQ *) * (info->nsqpool + 1));
  if ((info->sqpool[info->nsqpool] = esl_sq_CreateDigital(info->wq[info->qbase].om->abc)) == NULL) goto ERROR;
  return info->sqpool[info->nsqpool++];

 ERROR:
  esl_fatal("Failed to allocate a target sequence");
  return NULL;
}

static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc)
{
  int      sstatus;
  ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */
  int seq_cnt = 0;
  int i, k;

  worker_Fresh(info);
  if (dbc)
    { /* --dbcache: the same targets in the same order, one block at a time */
      ESL_SQ_BLOCK *blk;
      int           reading = ! dbc->complete;  /* first pass: <blk> is the reused scratch block */

      dbcache_StartPass(dbc);
      while ((sstatus = dbcache_NextBlock(dbc, dbfp, &n_targetseqs, dbc->scratch, &blk)) == eslOK)
        {
          if (info->pair && ! reading)   /* cached targets stay put: pairs can span blocks */
            { for (k = info->qbase; k < info->qbase + info->nq; k++) for (i = 0; i < blk->count; i++) pair_feed(info, &(info->wq[k]), blk->list + i); }
          else worker_Search(info, info->qbase, info->nq, blk->list, blk->count);
          if (reading) for (i = 0; i < blk->count; i++) esl_sq_Reuse(blk->list + i);
        }
      if (info->pair && ! reading) for (k = info->qbase; k < info->qbase + info->nq; k++) pair_flush(info, &(info->wq[k]));
      return sstatus;
    }

  if (! info->pair)
    {
      dbsq = esl_sq_CreateDigital(info->wq[info->qbase].om->abc);

      /* Main loop: */
      while ( (n_targetseqs==-1 || seq_cnt<n_targetseqs) &&  (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
      {
          worker_Search(info, info->qbase, info->nq, dbsq, 1);
          seq_cnt++;
          esl_sq_Reuse(dbsq);
      }

      if (n_targetseqs!=-1 && seq_cnt==n_targetseqs)
        sstatus = eslEOF;

      esl_sq_Destroy(dbsq);
      return sstatus;
    }

  /* With paired parsers: targets are read one at a time as above, but a
   * target a query holds (waiting for a partner) keeps its ESL_SQ; the next
   * one is read into another from the pool. Every query finishes its
   * waiting target at the end of the pass.
   */
  dbsq = serial_sq(info);
  while ( (n_targetseqs==-1 || seq_cnt<n_targetseqs) &&  (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
    {
      for (k = info->qbase; k < info->qbase + info->nq; k++) pair_feed(info, &(info->wq[k]), dbsq);
      seq_cnt++;
      dbsq = serial_sq(info);
    }
  if (n_targetseqs!=-1 && seq_cnt==n_targetseqs)
    sstatus = eslEOF;
  for (k = info->qbase; k < info->qbase + info->nq; k++) pair_flush(info, &(info->wq[k]));
  return sstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, WORKER_GATE *gate, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc, int nq)
{
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  int  nactive;                 /* workers taking part in this pass                     */
  int  nblocks  = 0;            /* nonempty blocks handed out                           */
  int  nseqs    = 0;            /* targets handed out                                   */
  int  bs       = 0;            /* max targets per block; 0 = as many as a block holds  */
  int  inplace;
  int  off      = 0;            /* --dbcache: next target in <cur> to hand out          */
  int  w;
  ESL_SQ_BLOCK *block;
  ESL_SQ_BLOCK *cblk;
  ESL_SQ_BLOCK *cur     = NULL; /* --dbcache: cached block being handed out in slices   */
  void         *newBlock;

  esl_workqueue_Reset(queue);   /* safe: the workers are parked at the gate, not in the queue */
  if (dbc) dbcache_StartPass(dbc);

  /* Every pass hands out the same targets, so after the first one we
   * know how many there are. When they come from the --dbcache, split a
   * small db into SLICES blocks per worker (but at least BLOCK_MIN targets
   * each), so they all get a share; a big one keeps full blocks. Finer
   * slices even out the end of a pass better, but they change how many
   * targets each worker has seen when it decides to keep a hit, so they're
   * used only when that can't matter (split_is_free()). Without
   * the cache, the reader sets the pace of a one-query pass and small
   * blocks don't help; a pass with several queries reads blocks of the
   * same size as the slices, when splitting can't change results.
   * Then wake only as many workers as there will be blocks: the rest
   * would only collect an end-of-pass block. Which worker searches which
   * target doesn't change the results.
   */
  bs = pass_blocksize(gate, dbc, nq);
  if      (bs > 0)          nactive = ESL_MIN((gate->nseqs + bs - 1) / bs, gate->nworkers);
  else if (gate->nblocks)   nactive = ESL_MIN(gate->nblocks, gate->nworkers);
  else                      nactive = gate->nworkers;
  pthread_mutex_lock(&gate->mutex);
  gate->nrunning = nactive;
  gate->pass++;
  for (w = 0; w < nactive; w++) { gate->wpass[w] = gate->pass; pthread_cond_signal(&gate->go[w]); }
  pthread_mutex_unlock(&gate->mutex);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  /* Main loop: */
  while (sstatus == eslOK )
    {
      block = (ESL_SQ_BLOCK *) newBlock;

      if (dbc)
      { /* --dbcache: first pass reads into the pool block; later, it carries cached sequences, up to <bs> at a time */
        inplace = FALSE;
        while (cur == NULL || off >= cur->count)
          {
            sstatus = dbcache_NextBlock(dbc, dbfp, &n_targetseqs, block, &cblk);
            if (sstatus != eslOK) break;
            if (cblk == block) { inplace = TRUE; break; }
            cur = cblk;
            off = 0;
          }
        if (sstatus != eslOK) block->count = 0;
        else if (! inplace)
          {
            block->list  = cur->list + off;
            block->count = (bs > 0) ? ESL_MIN(bs, cur->count - off) : cur->count - off;
            off += block->count;
          }
      } else if (n_targetseqs == 0)
      {
        block->count = 0;
        sstatus = eslEOF;
      } else {
        sstatus = esl_sqio_ReadBlock(dbfp, block, -1, (bs > 0 && (n_targetseqs < 0 || n_targetseqs > bs)) ? bs : n_targetseqs, /*max_init_window=*/FALSE, FALSE);
        n_targetseqs -= block->count;
      }

      if (sstatus == eslEOF)
      {
        if (eofCount < nactive) sstatus = eslOK;
        ++eofCount;
      }

      if (sstatus == eslOK)
      {
        if (block->count > 0) { nblocks++; nseqs += block->count; }
        status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
        if (status != eslOK) esl_fatal("Work queue reader failed");
      }
    }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the workers to finish this pass */
      pthread_mutex_lock(&gate->mutex);
      while (gate->nrunning > 0) pthread_cond_wait(&gate->done, &gate->mutex);
      pthread_mutex_unlock(&gate->mutex);
      esl_workqueue_Complete(queue);
      if (dbc) dbcache_EndPass(dbc);
      gate->nblocks = nblocks;
      gate->nseqs   = nseqs;
    }

  return sstatus;
}

/* pass_blocksize()
 * Most targets per block in a pass with <nq> queries (0: as many as a
 * block holds), from the number of targets seen in the first pass; see
 * the comment in thread_loop().
 */
static int
pass_blocksize(const WORKER_GATE *gate, const DBCACHE *dbc, int nq)
{
  int bs = 0;
  int nper;

  if (gate->nseqs > 0 && ((dbc && dbc->complete) || (nq > 1 && gate->split_free)))
    {
      nper = (gate->split_free ? SLICES_FREE : SLICES) * gate->nworkers;
      bs   = (gate->nseqs + nper - 1) / nper;
      bs   = ESL_MAX(bs, BLOCK_MIN);
      if (bs >= BLOCK_SIZE) bs = 0;
    }
  return bs;
}

/* Streaming: after the first pass, consecutive batches of queries
 * overlap. The master hands out all blocks of batch b (read, or sliced
 * from the --dbcache), each tagged with b's slot set, and only then waits
 * for batch b-1 to finish and prints it; workers that run out of b-1's
 * blocks go on with b's instead of waiting at a barrier. Two batches are
 * in flight at most, in slot sets 0 and 1. Which worker searches which
 * block, and when, doesn't change any result (PLAN.md §2), and each
 * query's results are printed in order, as before.
 *
 * stream_Start() puts all workers into stream_search(); stream_HandOut()
 * hands out one batch; stream_Wait() waits for one; stream_Stop() sends
 * each worker an empty block to leave the stream and waits for them.
 */
static void
stream_Start(ESL_WORK_QUEUE *queue, WORKER_GATE *gate, void **ret_block)
{
  int w;

  esl_workqueue_Reset(queue);   /* safe: the workers are parked at the gate */
  pthread_mutex_lock(&gate->mutex);
  gate->stream   = TRUE;
  gate->nrunning = gate->nworkers;
  gate->pass++;
  gate->bfin[0]  = gate->bfin[1]  = TRUE;
  gate->bsent[0] = gate->bsent[1] = gate->bdone[0] = gate->bdone[1] = 0;
  for (w = 0; w < gate->nworkers; w++) { gate->wpass[w] = gate->pass; pthread_cond_signal(&gate->go[w]); }
  pthread_mutex_unlock(&gate->mutex);
  if (esl_workqueue_ReaderUpdate(queue, NULL, ret_block) != eslOK) esl_fatal("Work queue reader failed");
}

/* the slot set a pool block carries */
static int
stream_tag(const WORKER_GATE *gate, const ESL_SQ_BLOCK *block)
{
  int i;
  for (i = 0; i < gate->npool; i++) if (gate->pool[i] == block) return gate->tag[i];
  esl_fatal("stream: block not in the pool");
  return -1;
}

static void
stream_settag(WORKER_GATE *gate, const ESL_SQ_BLOCK *block, int set)
{
  int i;
  for (i = 0; i < gate->npool; i++) if (gate->pool[i] == block) { gate->tag[i] = set; return; }
  esl_fatal("stream: block not in the pool");
}

/* stream_HandOut()
 * Hand out every block of batch <bid> (<nq> queries in slot set <set>):
 * read from <dbfp> (rewound by the caller) or sliced from the complete
 * --dbcache, as thread_loop() does for a pass. <*io_block> is the empty
 * block the master holds; on return, it holds the next one. Returns
 * eslEOF at the end of the targets, or a read error.
 */
static int
stream_HandOut(ESL_WORK_QUEUE *queue, WORKER_GATE *gate, ESL_SQFILE *dbfp, int n_targetseqs, DBCACHE *dbc, int set, int nq, int64_t bid, void **io_block)
{
  ESL_SQ_BLOCK *block;
  ESL_SQ_BLOCK *cblk;
  ESL_SQ_BLOCK *cur     = NULL;   /* --dbcache: cached block being handed out in slices */
  int           off     = 0;
  int           bs      = pass_blocksize(gate, dbc, nq);
  int           nsent   = 0;
  int           sstatus = eslOK;
  int           inplace;

  /* the set's previous batch is printed already: its slots are free */
  pthread_mutex_lock(&gate->mutex);
  gate->bnq[set]   = nq;
  gate->bid[set]   = bid;
  gate->bsent[set] = gate->bdone[set] = 0;
  gate->bfin[set]  = FALSE;
  pthread_mutex_unlock(&gate->mutex);
  if (dbc) dbcache_StartPass(dbc);

  for (;;)
    {
      block = (ESL_SQ_BLOCK *) *io_block;
      if (dbc)
	{
	  inplace = FALSE;
	  while (cur == NULL || off >= cur->count)
	    {
	      sstatus = dbcache_NextBlock(dbc, dbfp, &n_targetseqs, block, &cblk);
	      if (sstatus != eslOK) break;
	      if (cblk == block) { inplace = TRUE; break; }
	      cur = cblk;
	      off = 0;
	    }
	  if (sstatus != eslOK) break;
	  if (! inplace)
	    {
	      block->list  = cur->list + off;
	      block->count = (bs > 0) ? ESL_MIN(bs, cur->count - off) : cur->count - off;
	      off += block->count;
	    }
	}
      else
	{
	  if (n_targetseqs == 0) { sstatus = eslEOF; break; }
	  sstatus = esl_sqio_ReadBlock(dbfp, block, -1, (bs > 0 && (n_targetseqs < 0 || n_targetseqs > bs)) ? bs : n_targetseqs, /*max_init_window=*/FALSE, FALSE);
	  if (sstatus != eslOK) break;
	  n_targetseqs -= block->count;
	}
      if (block->count == 0) continue;

      stream_settag(gate, block, set);
      nsent++;
      if (esl_workqueue_ReaderUpdate(queue, block, io_block) != eslOK) esl_fatal("Work queue reader failed");
    }

  pthread_mutex_lock(&gate->mutex);
  gate->bsent[set] = nsent;
  gate->bfin[set]  = TRUE;
  pthread_mutex_unlock(&gate->mutex);
  return sstatus;
}

/* stream_Wait()
 * Wait until every block of the batch in slot set <set> is searched.
 */
static void
stream_Wait(WORKER_GATE *gate, int set)
{
  pthread_mutex_lock(&gate->mutex);
  while (! (gate->bfin[set] && gate->bdone[set] == gate->bsent[set])) pthread_cond_wait(&gate->done, &gate->mutex);
  pthread_mutex_unlock(&gate->mutex);
}

/* stream_Stop()
 * After the last batch: one empty block per worker takes it out of the
 * stream; wait until all are parked at the gate again.
 */
static void
stream_Stop(ESL_WORK_QUEUE *queue, WORKER_GATE *gate, DBCACHE *dbc, void **io_block)
{
  ESL_SQ_BLOCK *block;
  int           w;

  for (w = 0; w < gate->nworkers; w++)
    {
      block = (ESL_SQ_BLOCK *) *io_block;
      block->count = 0;
      if (esl_workqueue_ReaderUpdate(queue, block, (w < gate->nworkers - 1) ? io_block : NULL) != eslOK) esl_fatal("Work queue reader failed");
    }
  *io_block = NULL;
  pthread_mutex_lock(&gate->mutex);
  while (gate->nrunning > 0) pthread_cond_wait(&gate->done, &gate->mutex);
  gate->stream = FALSE;
  pthread_mutex_unlock(&gate->mutex);
  esl_workqueue_Complete(queue);
  if (dbc) dbcache_EndPass(dbc);
}

/* stream_search()
 * A worker's part of the stream: search each block with the queries of
 * the batch it's tagged with, until an empty block arrives.
 */
static void
stream_search(WORKER_INFO *info)
{
  WORKER_GATE  *gate = info->gate;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  int64_t       last = -1;
  int           set, i;

  if (esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock) != eslOK) esl_fatal("Work queue worker failed");
  block = (ESL_SQ_BLOCK *) newBlock;
  while (block->count > 0)
    {
      set = stream_tag(gate, block);
      if (gate->bid[set] != last) { worker_Fresh(info); last = gate->bid[set]; }

      worker_Search(info, set * MAXQ, gate->bnq[set], block->list, block->count);
      if (info->reuse_seqs)   /* not with --dbcache: cached targets are read-only */
	for (i = 0; i < block->count; ++i) esl_sq_Reuse(block->list + i);

      pthread_mutex_lock(&gate->mutex);
      if (++gate->bdone[set] == gate->bsent[set] && gate->bfin[set]) pthread_cond_signal(&gate->done);
      pthread_mutex_unlock(&gate->mutex);

      if (esl_workqueue_WorkerUpdate(info->queue, block, &newBlock) != eslOK) esl_fatal("Work queue worker failed");
      block = (ESL_SQ_BLOCK *) newBlock;
    }
  if (esl_workqueue_WorkerUpdate(info->queue, block, NULL) != eslOK) esl_fatal("Work queue worker failed");
}

static void 
pipeline_thread(void *arg)
{
  int            workeridx;
  int            pass = 0;
  int            stream;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  WORKER_GATE   *gate;
  QP_JOB        *job;

  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);
  gate = info->gate;

  /* a share of one pass over the targets per query, or a whole query
   * in query-parallel mode, until the master says quit
   */
  for (;;)
    {
      pthread_mutex_lock(&gate->mutex);
      while (gate->wpass[workeridx] == pass && gate->wjob[workeridx] == NULL && ! gate->quit) pthread_cond_wait(&gate->go[workeridx], &gate->mutex);
      if (gate->quit) { pthread_mutex_unlock(&gate->mutex); break; }

      if ((job = gate->wjob[workeridx]) != NULL)
	{
	  pthread_mutex_unlock(&gate->mutex);

	  search_query(info, gate, job);

	  pthread_mutex_lock(&gate->mutex);
	  job->done              = TRUE;
	  gate->wjob[workeridx]  = NULL;
	  gate->idle[gate->nidle++] = workeridx;
	  pthread_cond_signal(&gate->done);
	  pthread_mutex_unlock(&gate->mutex);
	  continue;
	}

      pass   = gate->wpass[workeridx];
      stream = gate->stream;
      pthread_mutex_unlock(&gate->mutex);

      if (stream) stream_search(info);
      else        search_pass(info);

      pthread_mutex_lock(&gate->mutex);
      if (--gate->nrunning == 0) pthread_cond_signal(&gate->done);
      pthread_mutex_unlock(&gate->mutex);
    }

  esl_threads_Finished(obj, workeridx);
  return;
}

/* worker_SetQuery()
 * Prepare query slot <k> of a worker for query profile <om>: a hit
 * list, the pipeline (a copy of the worker's freshly configured one,
 * sharing its workspace), a clone of <om>, and p7_pli_NewModel(). The
 * worker's own pipeline is created on first use. A slot's setup touches
 * nothing the worker may be using for other slots, so the master can set
 * up one batch while the workers search another. Returns eslOK, or
 * eslEINVAL if the pipeline can't use this model (e.g. --cut_ga without
 * GA lines); the message is in the slot's pli->errbuf.
 */
static int
worker_SetQuery(WORKER_INFO *info, int k, const P7_OPROFILE *om, const ESL_GETOPTS *go, const ESL_ALPHABET *abc)
{
  WORKER_QUERY *q = &(info->wq[k]);
  int           status;

  if (info->work == NULL)
    {
      if ((info->work = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS)) == NULL) goto ERROR; /* L_hint = 100 is just a dummy for now */
      info->tmpl = *(info->work);
    }
  if (q->pli == NULL)
    {
      if ((q->th = p7_tophits_Create())  == NULL) goto ERROR;
      if ((q->bg = p7_bg_Create(abc))    == NULL) goto ERROR;
      ESL_ALLOC(q->pli, sizeof(P7_PIPELINE));
    }
  else p7_tophits_Reuse(q->th);

  *(q->pli) = info->tmpl;
  q->om    = p7_oprofile_Clone(om);
  q->nwait = 0;
  return p7_pli_NewModel(q->pli, q->om, q->bg);

 ERROR:
  esl_fatal("Failed to allocate a query's pipeline");
  return eslEMEM;
}

/* worker_Fresh()
 * Called by the worker itself when it starts on a new batch: lay out the
 * DP matrices' rows for the next model length again, as a new pipeline
 * would (p7_pipeline_Reinit() does the same), instead of keeping the
 * widest model's row stride. Only speed depends on it.
 */
static void
worker_Fresh(WORKER_INFO *info)
{
  if (info->work == NULL) return;
  info->work->oxf->validR = info->work->oxb->validR = info->work->fwd->validR = info->work->bck->validR = 0;
}

/* worker_Trim()
 * After a target's domain stage: if a long region made the worker's
 * domain definition matrices (fwd, bck: a row per residue of the region,
 * a column per model node) bigger than DP_KEEP_CELLS, give the memory
 * back. Upstream does that by making a new pipeline for every query;
 * reused matrices would otherwise keep, in every worker, the largest
 * region times the longest model it has met (16 threads, 250 human
 * proteins of 696-13,477 residues against Pfam-A: 740 MB peak, upstream
 * 192 MB). The trim is in place: the slots' pipelines are struct copies
 * that hold these P7_OMX pointers, so a fresh matrix's contents go into
 * the old object and the old buffers are freed with the fresh one.
 * glibc would keep freed blocks in the worker's arena once its mmap
 * threshold has grown, so threaded runs fix the threshold at MMAP_FLOOR
 * (mallopt() in serial_master()): big matrices then have their own
 * mappings and are unmapped when trimmed. Fixing it also stops glibc's
 * own adjustment of the trim threshold, which would stay at 128 KB and
 * make arenas give back and refault memory all the time (Pfam-A x 5
 * sequences: 10x the madvise() calls, 5% slower); it is set as glibc's
 * adjustment would set it, to twice the mmap threshold. Only memory and
 * speed depend on any of this.
 */
static void
worker_Trim(WORKER_INFO *info, int M)
{
  P7_OMX *ox[2], *nx, tmp;
  int     h;

  ox[0] = info->work->fwd;
  ox[1] = info->work->bck;
  for (h = 0; h < 2; h++)
    if (ox[h]->ncells > DP_KEEP_CELLS)
      {
	if ((nx = p7_omx_Create(M, 100, 100)) == NULL) esl_fatal("Failed to allocate a DP matrix");
	tmp = *(ox[h]); *(ox[h]) = *nx; *nx = tmp;
	p7_omx_Destroy(nx);
      }
}

/* worker_Search()
 * Search targets <dbsq[0..n-1]> with queries wq[base..base+nq-1], in
 * order: per query, exactly the calls a one-query search makes.
 */
static void
worker_Search(WORKER_INFO *info, int base, int nq, ESL_SQ *dbsq, int n)
{
  WORKER_QUERY *q;
  int           k, i;

  for (k = base; k < base + nq; k++)
    {
      q = &(info->wq[k]);
      if (info->pair)
	{
	  for (i = 0; i < n; i++) pair_feed(info, q, dbsq + i);
	  pair_flush(info, q);    /* the block's targets don't outlive it */
	  continue;
	}
      for (i = 0; i < n; i++)
	{
	  p7_pli_NewSeq(q->pli, dbsq + i);
	  p7_bg_SetLength(q->bg, dbsq[i].n);
	  p7_oprofile_ReconfigLength(q->om, dbsq[i].n);

	  p7_Pipeline(q->pli, q->om, q->bg, dbsq + i, NULL, q->th);

	  p7_pipeline_Reuse(q->pli);
	  worker_Trim(info, q->om->M);
	}
    }
}

/* pair_finish()
 * The rest of p7_Pipeline() for a target whose parsers have run (Forward
 * score <fwdsc>, matrices <fwd>, <bck>): <om> and <bg> configured for its
 * length again and the pipeline's Z as it was for it, then the domain
 * stage and the hit.
 */
static void
pair_finish(WORKER_INFO *info, WORKER_QUERY *q, const PAIR_WAIT *t, float fwdsc, P7_OMX *fwd, P7_OMX *bck)
{
  double Z = q->pli->Z;

  p7_oprofile_ReconfigLength(q->om, t->sq->n);
  p7_bg_SetLength(q->bg, t->sq->n);
  q->pli->Z = t->Z;
  p7_pli_Domains(q->pli, q->om, q->bg, t->sq, NULL, t->nullsc, fwdsc, fwd, bck, q->th);
  q->pli->Z = Z;
  p7_pipeline_Reuse(q->pli);
  worker_Trim(info, q->om->M);
}

/* pair_one()
 * A target without a partner: the single-target parsers, as p7_Pipeline().
 */
static void
pair_one(WORKER_INFO *info, WORKER_QUERY *q, const PAIR_WAIT *t)
{
  float fwdsc;

  p7_oprofile_ReconfigLength(q->om, t->sq->n);
  p7_omx_GrowTo(info->pf[0], q->om->M, 0, t->sq->n);
  p7_ForwardParser(t->sq->dsq, t->sq->n, q->om, info->pf[0], &fwdsc);
  if (! p7_pli_ForwardPass(q->pli, q->om, t->filtersc, fwdsc)) return;
  p7_omx_GrowTo(info->pb[0], q->om->M, 0, t->sq->n);
  p7_BackwardParser(t->sq->dsq, t->sq->n, q->om, info->pf[0], info->pb[0], NULL);
  pair_finish(info, q, t, fwdsc, info->pf[0], info->pb[0]);
}

/* pair_two()
 * Two targets, <a> before <b> in the db: both parsers on the pair (each
 * target's results are bitwise those of the single-target parsers, see
 * fwdback_pair.c), then each target's domain stage in db order.
 */
static void
pair_two(WORKER_INFO *info, WORKER_QUERY *q, const PAIR_WAIT *a, const PAIR_WAIT *b)
{
  const P7_OPROFILE *om = q->om;
  int                La = a->sq->n, Lb = b->sq->n;
  float              sa, sb;
  int                pa, pb;

  p7_omx_GrowTo(info->pf[0], om->M, 0, La);
  p7_omx_GrowTo(info->pf[1], om->M, 0, Lb);
  p7_ForwardParserPair(om, a->sq->dsq, La, (const float (*)[p7O_NXTRANS]) a->xf, info->pf[0], &sa,
		           b->sq->dsq, Lb, (const float (*)[p7O_NXTRANS]) b->xf, info->pf[1], &sb);
  pa = p7_pli_ForwardPass(q->pli, om, a->filtersc, sa);
  pb = p7_pli_ForwardPass(q->pli, om, b->filtersc, sb);
  if (! pa && ! pb) return;

  p7_omx_GrowTo(info->pb[0], om->M, 0, La);
  p7_omx_GrowTo(info->pb[1], om->M, 0, Lb);
  if (pa && pb)
    p7_BackwardParserPair(om, a->sq->dsq, La, (const float (*)[p7O_NXTRANS]) a->xf, info->pf[0], info->pb[0], NULL,
			      b->sq->dsq, Lb, (const float (*)[p7O_NXTRANS]) b->xf, info->pf[1], info->pb[1], NULL);
#ifdef p7_CROSSCHECK
  { /* shadow-run the single-target parsers; any bit difference is fatal */
    const PAIR_WAIT *t[2] = { a, b };
    float            ps[2] = { sa, sb };
    int              ok[2] = { pa, pb };
    int              h, i;
    for (h = 0; h < 2; h++)
      {
	P7_OMX *xf  = p7_omx_Create(om->M, 0, t[h]->sq->n);
	P7_OMX *xb  = p7_omx_Create(om->M, 0, t[h]->sq->n);
	float   s1;
	p7_oprofile_ReconfigLength(q->om, t[h]->sq->n);
	p7_ForwardParser(t[h]->sq->dsq, t[h]->sq->n, om, xf, &s1);
	if (memcmp(&s1, &ps[h], sizeof(float)) != 0) esl_fatal("p7_CROSSCHECK: paired Forward score %a != %a (L=%d)", ps[h], s1, t[h]->sq->n);
	for (i = 0; i <= t[h]->sq->n; i++)
	  if (memcmp(xf->xmx + i*p7X_NXCELLS, info->pf[h]->xmx + i*p7X_NXCELLS, sizeof(float) * p7X_NXCELLS) != 0)
	    esl_fatal("p7_CROSSCHECK: paired Forward row %d differs (L=%d)", i, t[h]->sq->n);
	if (ok[0] && ok[1])
	  {
	    p7_BackwardParser(t[h]->sq->dsq, t[h]->sq->n, om, xf, xb, NULL);
	    if (xb->has_own_scales != info->pb[h]->has_own_scales || memcmp(&xb->totscale, &info->pb[h]->totscale, sizeof(float)) != 0)
	      esl_fatal("p7_CROSSCHECK: paired Backward scales differ (L=%d)", t[h]->sq->n);
	    for (i = 0; i <= t[h]->sq->n; i++)
	      if (memcmp(xb->xmx + i*p7X_NXCELLS, info->pb[h]->xmx + i*p7X_NXCELLS, sizeof(float) * p7X_NXCELLS) != 0)
		esl_fatal("p7_CROSSCHECK: paired Backward row %d differs (L=%d)", i, t[h]->sq->n);
	  }
	p7_omx_Destroy(xf);
	p7_omx_Destroy(xb);
      }
  }
#endif
  if (pa)
    {
      if (! pb) { p7_oprofile_ReconfigLength(q->om, La); p7_BackwardParser(a->sq->dsq, La, q->om, info->pf[0], info->pb[0], NULL); }
      pair_finish(info, q, a, sa, info->pf[0], info->pb[0]);
    }
  if (pb)
    {
      if (! pa) { p7_oprofile_ReconfigLength(q->om, Lb); p7_BackwardParser(b->sq->dsq, Lb, q->om, info->pf[1], info->pb[1], NULL); }
      pair_finish(info, q, b, sb, info->pf[1], info->pb[1]);
    }
}

/* pair_feed()
 * Search target <sq> with query slot <q>, running the Forward/Backward
 * parsers on two targets at once: a target that passes the filters waits
 * (in q->wait) for the next one that does; the caller keeps a waiting
 * target's <sq> valid until pair_flush(). Per target, the calls are those
 * of p7_Pipeline(), in its order: p7_pli_Filter(), the parsers,
 * p7_pli_ForwardPass(), p7_pli_Domains() (fwdback_pair.c: the paired
 * parsers give each target's single-target results bit for bit). Hits are
 * added in db order. A waiting target's later stages get its own length
 * configuration and running Z back (pair_finish()).
 */
static void
pair_feed(WORKER_INFO *info, WORKER_QUERY *q, const ESL_SQ *sq)
{
  PAIR_WAIT  b;
  PAIR_WAIT *t = q->nwait ? &b : &(q->wait);
  int        h, pass;

  if (info->pf[0] == NULL)
    for (h = 0; h < 2; h++)
      {
	if ((info->pf[h] = p7_omx_Create(q->om->M, 0, 100)) == NULL) esl_fatal("Failed to allocate parser matrices");
	if ((info->pb[h] = p7_omx_Create(q->om->M, 0, 100)) == NULL) esl_fatal("Failed to allocate parser matrices");
      }

  p7_pli_NewSeq(q->pli, sq);
  p7_bg_SetLength(q->bg, sq->n);
  p7_oprofile_ReconfigLength(q->om, sq->n);

  p7_pli_Filter(q->pli, q->om, q->bg, sq, &(t->nullsc), &(t->filtersc), &pass);
  p7_pipeline_Reuse(q->pli);
  if (! pass) return;

  t->sq = sq;
  t->Z  = q->pli->Z;
  memcpy(t->xf, q->om->xf, sizeof(t->xf));
  if (! q->nwait) { q->nwait = 1; return; }
  pair_two(info, q, &(q->wait), &b);
  q->nwait = 0;
}

/* pair_flush()
 * Finish slot <q>'s waiting target, if any, without a partner.
 */
static void
pair_flush(WORKER_INFO *info, WORKER_QUERY *q)
{
  if (q->nwait) pair_one(info, q, &(q->wait));
  q->nwait = 0;
}

/* output_batch()
 * After a batch's pass: for each of its queries (slots base..base+nb-1),
 * merge the workers' hit lists and accounting into worker 0's, print the
 * results, and free the query.
 */
static int
output_batch(WORKER_INFO *info, int infocnt, int base, int nb, P7_HMM **bhmm, P7_PROFILE **bgm, P7_OPROFILE **bom, int *bnum, ESL_STOPWATCH **w,
	     const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw)
{
  int i, k;
  int status;

  for (k = 0; k < nb; k++)
    {
      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
	{
	  p7_tophits_Merge (info[0].wq[base+k].th,  info[i].wq[base+k].th);
	  p7_pipeline_Merge(info[0].wq[base+k].pli, info[i].wq[base+k].pli);
	  p7_oprofile_Destroy(info[i].wq[base+k].om);
	}

      /* Print the results.  */
      if ((status = output_query_header(outs->ofp, bhmm[k])) != eslOK) return status;
      if ((status = output_results(outs, go, abc, textw, bhmm[k], bnum[k], info->wq[base+k].th, info->wq[base+k].pli, w[k])) != eslOK) return status;

      p7_oprofile_Destroy(info->wq[base+k].om);
      p7_oprofile_Destroy(bom[k]);
      p7_profile_Destroy(bgm[k]);
      p7_hmm_Destroy(bhmm[k]);
    }
  return eslOK;
}

/* split_is_free()
 * TRUE if splitting the targets among workers differently can't change
 * any output. A worker keeps a hit for the final thresholding when
 * p7_pli_TargetReportable() says so, and with E-value thresholds and no
 * -Z that uses the worker's own running target count as Z: which
 * borderline hits are kept, and so the widths of name columns, then
 * depend on the split (upstream behaviour; see dev/profiling/PLAN.md §2).
 * With -Z or bit score thresholds nothing depends on it.
 */
static int
split_is_free(const ESL_GETOPTS *go)
{
  return (esl_opt_IsOn(go, "-Z") || esl_opt_IsOn(go, "-T") ||
          esl_opt_GetBoolean(go, "--cut_ga") || esl_opt_GetBoolean(go, "--cut_nc") || esl_opt_GetBoolean(go, "--cut_tc"));
}

/* search_query()
 * Query-parallel mode: search all cached targets, in file order, for
 * one query, with this worker's pipeline and null model, as the serial
 * search does.
 */
static void
search_query(WORKER_INFO *info, WORKER_GATE *gate, QP_JOB *job)
{
  P7_PIPELINE  *pli = info->wq[0].pli;
  P7_BG        *bg  = info->wq[0].bg;
  WORKER_QUERY  jq;              /* paired parsers: this job's query as a slot */
  ESL_SQ_BLOCK *blk;
  ESL_SQ       *dbsq;
  int           b, i;

  *pli = info->tmpl;
  worker_Fresh(info);
  jq.bg = bg; jq.pli = pli; jq.th = job->th; jq.om = job->om; jq.nwait = 0;
  job->status = p7_pli_NewModel(pli, job->om, bg);
  if (job->status == eslEINVAL) { strcpy(job->errbuf, pli->errbuf); return; }

  for (b = 0; b < gate->dbc->nblk; b++)
    {
      blk = gate->dbc->blk[b];
      if (info->pair)
	{ /* this job's query as a slot: its profile and hit list, the worker's pipeline and null model;
	     cached targets stay put, so pairs can span blocks */
	  for (i = 0; i < blk->count; i++) pair_feed(info, &jq, blk->list + i);
	  continue;
	}
      for (i = 0; i < blk->count; i++)
	{
	  dbsq = blk->list + i;
	  p7_pli_NewSeq(pli, dbsq);
	  p7_bg_SetLength(bg, dbsq->n);
	  p7_oprofile_ReconfigLength(job->om, dbsq->n);

	  p7_Pipeline(pli, job->om, bg, dbsq, NULL, job->th);

	  p7_pipeline_Reuse(pli);   /* cached targets are read-only: no esl_sq_Reuse() */
	  worker_Trim(info, job->om->M);
	}
    }
  if (info->pair) pair_flush(info, &jq);

  job->stats     = *pli;
  job->stats.oxf = job->stats.oxb = job->stats.fwd = job->stats.bck = NULL;
  job->stats.r   = NULL;
  job->stats.ddef = NULL;
  job->stats.hfp = NULL;
}

/* search_pass()
 * A worker's part of one query's search: process blocks from the work
 * queue until the reader's empty end-of-pass block arrives.
 */
static void
search_pass(WORKER_INFO *info)
{
  int i;
  int status;

  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;

  worker_Fresh(info);
  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (ESL_SQ_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      worker_Search(info, info->qbase, info->nq, block->list, block->count);
      if (info->reuse_seqs)   /* not with --dbcache: cached targets are read-only */
	for (i = 0; i < block->count; ++i) esl_sq_Reuse(block->list + i);

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (ESL_SQ_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");
}

/* qp_Create()
 * An empty queue of <njobs> query-parallel jobs.
 */
static void
qp_Create(QP_QUEUE *qp, int njobs)
{
  int i;
  int status;

  qp->njobs = njobs;
  qp->head  = qp->tail = 0;
  ESL_ALLOC(qp->job, sizeof(QP_JOB) * njobs);
  for (i = 0; i < njobs; i++)
    {
      qp->job[i].hmm = NULL;
      qp->job[i].gm  = NULL;
      qp->job[i].om  = NULL;
      if ((qp->job[i].th = p7_tophits_Create())    == NULL) goto ERROR;
      if ((qp->job[i].w  = esl_stopwatch_Create()) == NULL) goto ERROR;
    }
  return;

 ERROR:
  esl_fatal("Failed to allocate the query-parallel jobs");
}

/* qp_Submit()
 * Hand query <hmm> (number <nquery>) to an idle worker, with its profiles
 * from the prefetch helpers (<*next_gm>, <*next_om>; taken over), or built
 * here. While waiting for a free job or an idle worker, print the jobs
 * that are done, in order.
 */
static int
qp_Submit(QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw,
	  P7_HMM *hmm, int nquery, P7_PROFILE **next_gm, P7_OPROFILE **next_om, const P7_BG *bg)
{
  QP_JOB *job;
  int     w;
  int     ready;
  int     status;

  if (qp->tail - qp->head == (uint64_t) qp->njobs && (status = qp_OutputNext(qp, gate, outs, go, abc, textw)) != eslOK) return status;

  job = &(qp->job[qp->tail % qp->njobs]);
  job->nquery = nquery;
  job->hmm    = hmm;
  if (*next_om)
    {
      job->gm  = *next_gm;  *next_gm = NULL;
      job->om  = *next_om;  *next_om = NULL;
    }
  else
    {
      job->gm = p7_profile_Create (hmm->M, abc);
      job->om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, job->gm, 100, p7_LOCAL); /* as in serial_master() */
      p7_oprofile_Convert(job->gm, job->om);
    }
  p7_tophits_Reuse(job->th);
  job->status = eslOK;
  job->done   = FALSE;
  esl_stopwatch_Start(job->w);

  pthread_mutex_lock(&gate->mutex);
  while (gate->nidle == 0)
    {
      if (qp->head < qp->tail && qp->job[qp->head % qp->njobs].done)
	{
	  pthread_mutex_unlock(&gate->mutex);
	  if ((status = qp_OutputNext(qp, gate, outs, go, abc, textw)) != eslOK) return status;
	  pthread_mutex_lock(&gate->mutex);
	}
      else pthread_cond_wait(&gate->done, &gate->mutex);
    }
  w = gate->idle[--gate->nidle];
  gate->wjob[w] = job;
  pthread_cond_signal(&gate->go[w]);
  pthread_mutex_unlock(&gate->mutex);
  qp->tail++;

  /* print what's done already, without waiting */
  for (;;)
    {
      pthread_mutex_lock(&gate->mutex);
      ready = (qp->head < qp->tail && qp->job[qp->head % qp->njobs].done);
      pthread_mutex_unlock(&gate->mutex);
      if (! ready) break;
      if ((status = qp_OutputNext(qp, gate, outs, go, abc, textw)) != eslOK) return status;
    }
  return eslOK;
}

/* qp_OutputNext()
 * Wait for the oldest job in flight, print its results as the
 * target-parallel path would, and free its query.
 */
static int
qp_OutputNext(QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw)
{
  QP_JOB *job = &(qp->job[qp->head % qp->njobs]);
  FILE   *ofp = outs->ofp;
  int     status;

  pthread_mutex_lock(&gate->mutex);
  while (! job->done) pthread_cond_wait(&gate->done, &gate->mutex);
  pthread_mutex_unlock(&gate->mutex);

  if (fprintf(ofp, "Query:       %s  [M=%d]\n", job->hmm->name, job->hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (job->hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", job->hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
  if (job->hmm->desc) { if (fprintf(ofp, "Description: %s\n", job->hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
  if (job->status == eslEINVAL) p7_Fail(job->errbuf);

  if ((status = output_results(outs, go, abc, textw, job->hmm, job->nquery, job->th, &(job->stats), job->w)) != eslOK) return status;

  p7_oprofile_Destroy(job->om);  job->om  = NULL;
  p7_profile_Destroy (job->gm);  job->gm  = NULL;
  p7_hmm_Destroy     (job->hmm); job->hmm = NULL;
  qp->head++;
  return eslOK;
}

/* qp_Drain()
 * Print all jobs in flight, in order.
 */
static int
qp_Drain(QP_QUEUE *qp, WORKER_GATE *gate, const OUTPUTS *outs, const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int textw)
{
  int status;

  while (qp->head < qp->tail)
    if ((status = qp_OutputNext(qp, gate, outs, go, abc, textw)) != eslOK) return status;
  return eslOK;
}

static void
qp_Destroy(QP_QUEUE *qp)
{
  int i;

  for (i = 0; i < qp->njobs; i++)
    {
      p7_tophits_Destroy   (qp->job[i].th);
      esl_stopwatch_Destroy(qp->job[i].w);
      p7_oprofile_Destroy  (qp->job[i].om);
      p7_profile_Destroy   (qp->job[i].gm);
      p7_hmm_Destroy       (qp->job[i].hmm);
    }
  free(qp->job);
}

/* query_prep_Start()
 * Start <nthreads> helpers that read queries from <hfp> after the ones
 * the master has read, and prepare them, up to <2*nthreads> ahead.
 */
static void
query_prep_Start(QUERY_PREP *prep, P7_HMMFILE *hfp, ESL_ALPHABET *abc, int nthreads)
{
  int i;
  int status;

  prep->hfp      = hfp;
  prep->abc      = abc;
  prep->nthreads = nthreads;
  prep->nslots   = 2 * nthreads;
  prep->nread    = 0;
  prep->nused    = 0;
  prep->done     = FALSE;
  prep->quit     = FALSE;
  ESL_ALLOC(prep->thread, sizeof(pthread_t) * prep->nthreads);
  ESL_ALLOC(prep->slot,   sizeof(PREP_SLOT) * prep->nslots);
  for (i = 0; i < prep->nslots; i++) prep->slot[i].ready = FALSE;

  if (pthread_mutex_init(&prep->mutex, NULL) != 0 || pthread_cond_init(&prep->work, NULL) != 0 || pthread_cond_init(&prep->ready, NULL) != 0)
    esl_fatal("Failed to initialize the query prefetch threads");
  for (i = 0; i < prep->nthreads; i++)
    if (pthread_create(&prep->thread[i], NULL, (void *(*)(void *)) query_prep_thread, prep) != 0)
      esl_fatal("Failed to start the query prefetch threads");
  return;

 ERROR:
  esl_fatal("Failed to allocate the query prefetch threads");
}

/* query_prep_Next()
 * Wait for the next query in file order; hand over its HMM, profiles
 * (NULL unless <*ret_hstatus> is eslOK) and the p7_hmmfile_Read() status.
 */
static void
query_prep_Next(QUERY_PREP *prep, P7_HMM **ret_hmm, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om, int *ret_hstatus)
{
  PREP_SLOT *sl;

  pthread_mutex_lock(&prep->mutex);
  sl = &(prep->slot[prep->nused % prep->nslots]);
  while (! sl->ready) pthread_cond_wait(&prep->ready, &prep->mutex);
  *ret_hmm     = sl->hmm;
  *ret_gm      = sl->gm;
  *ret_om      = sl->om;
  *ret_hstatus = sl->hstatus;
  sl->ready    = FALSE;
  prep->nused++;
  pthread_cond_signal(&prep->work);  /* a slot came free */
  pthread_mutex_unlock(&prep->mutex);
}

/* query_prep_Stop()
 * Stop the helpers, and free the queries they prepared that weren't used.
 */
static void
query_prep_Stop(QUERY_PREP *prep)
{
  int i;

  pthread_mutex_lock(&prep->mutex);
  prep->quit = TRUE;
  pthread_cond_broadcast(&prep->work);
  pthread_mutex_unlock(&prep->mutex);
  for (i = 0; i < prep->nthreads; i++) pthread_join(prep->thread[i], NULL);

  for (i = 0; i < prep->nslots; i++)
    if (prep->slot[i].ready)
      {
        p7_hmm_Destroy     (prep->slot[i].hmm);
        p7_profile_Destroy (prep->slot[i].gm);
        p7_oprofile_Destroy(prep->slot[i].om);
      }
  pthread_cond_destroy(&prep->work);
  pthread_cond_destroy(&prep->ready);
  pthread_mutex_destroy(&prep->mutex);
  free(prep->slot);
  free(prep->thread);
}

/* query_prep_thread()
 * A query prefetch helper: read the next query (serially, in file order)
 * and build its profiles with the same calls the master makes.
 */
static void
query_prep_thread(void *arg)
{
  QUERY_PREP  *prep = (QUERY_PREP *) arg;
  P7_BG       *bg;              /* own null model: p7_ProfileConfig() only reads bg->f */
  PREP_SLOT   *sl;
  P7_HMM      *hmm;
  P7_PROFILE  *gm;
  P7_OPROFILE *om;
  int          hstatus;

  impl_Init();                  /* same floating point modes as the master */
  if ((bg = p7_bg_Create(prep->abc)) == NULL) esl_fatal("Failed to allocate a null model");

  pthread_mutex_lock(&prep->mutex);
  for (;;)
    {
      while (! prep->quit && (prep->done || prep->nread - prep->nused >= (uint64_t) prep->nslots))
        pthread_cond_wait(&prep->work, &prep->mutex);
      if (prep->quit) break;

      sl  = &(prep->slot[prep->nread % prep->nslots]);
      prep->nread++;
      hmm = NULL;
      hstatus = p7_hmmfile_Read(prep->hfp, &(prep->abc), &hmm);   /* under the mutex: in file order */
      if (hstatus != eslOK) prep->done = TRUE;
      pthread_mutex_unlock(&prep->mutex);

      gm = NULL;
      om = NULL;
      if (hstatus == eslOK)
	{
	  gm = p7_profile_Create (hmm->M, prep->abc);
	  om = p7_oprofile_Create(hmm->M, prep->abc);
	  p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL); /* as in serial_master() */
	  p7_oprofile_Convert(gm, om);
	}

      pthread_mutex_lock(&prep->mutex);
      sl->hmm     = hmm;
      sl->gm      = gm;
      sl->om      = om;
      sl->hstatus = hstatus;
      sl->ready   = TRUE;
      pthread_cond_broadcast(&prep->ready);
    }
  pthread_mutex_unlock(&prep->mutex);
  p7_bg_Destroy(bg);
}
#endif   /* HMMER_THREADS */
 


