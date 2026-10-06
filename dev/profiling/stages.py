#!/usr/bin/env python3
"""Attribute perf samples of hmmsearch to pipeline stages (innermost known frame wins).
usage: perf script -i X.data -F comm,pid,tid,ip,sym | stages.py [--threads]
  --threads splits main (reader/output) thread from worker threads.
Needs a build with -g -fno-omit-frame-pointer (see dev/profiling/README.md)."""
import sys, re, collections

# innermost-first stage mapping: symbol -> stage
INNER = [
    (r'^(calc_band_(avx2_|avx512_)?\d+|p7_SSVFilter\w*|p7_MSVFilter\w*|get_xE)$', 'S1 MSV/SSV filter'),
    (r'^(esl_hmm_Forward|p7_bg_FilterScore)$',                         'S2 bias filter (bg HMM fwd)'),
    (r'^p7_bg_NullOne$',                                               'S0 null1 score'),
    (r'^p7_ViterbiFilter\w*$',                                         'S3 Viterbi filter'),
    (r'^p7_ForwardParser$',                                            'S4 Forward parser (filter)'),
    (r'^p7_BackwardParser$',                                           'S5 Backward parser'),
    (r'^p7_DomainDecoding$',                                           'D1 domain decoding (posterior B/E)'),
    (r'^(p7_StochasticTrace|p7_trace_Index|p7_spensemble_\w+|link_spsamples|esl_cluster_\w+|region_trace_ensemble|p7_Null2_ByTrace)$', 'D2 stochastic trace clustering'),
    (r'^p7_Forward$',                                                  'D3 envelope Forward'),
    (r'^p7_Backward$',                                                 'D3 envelope Backward'),
    (r'^p7_Decoding$',                                                 'D4 envelope posterior decoding'),
    (r'^p7_Null2_ByExpectation$',                                      'D5 null2'),
    (r'^(p7_OptimalAccuracy)$',                                        'D6 optimal-accuracy DP'),
    (r'^(p7_OATrace)$',                                                'D6 OA traceback'),
    (r'^p7_alidisplay_\w+$',                                           'D7 alidisplay'),
    (r'^(p7_oprofile_Reconfig\w+|p7_bg_SetLength|p7_pli_NewSeq|p7_omx_GrowTo|p7_domaindef_GrowTo|p7_pipeline_Reuse|p7_domaindef_Reuse|esl_sq_Reuse)$', 'P per-seq setup/reconfig'),
    (r'^(p7_tophits_Targets|p7_tophits_Domains|p7_alidisplay_Print|p7_tophits_Tabular\w+|p7_tophits_Alignment|esl_msafile_Write|p7_pli_Statistics)$', 'O output'),
    (r'^(p7_tophits_Sort\w*|p7_tophits_Merge|p7_tophits_Threshold|p7_tophits_CreateNextHit)$', 'O hit list sort/merge'),
    (r'^(p7_hmmfile_Read|p7_ProfileConfig|p7_oprofile_Convert|p7_oprofile_SetSIMD|p7_oprofile_Clone|p7_pipeline_Create|p7_pli_NewModel|p7_oprofile_Create|p7_profile_Create)$', 'Q per-query model setup'),
    (r'^(esl_sqio_ReadBlock|esl_sqio_Read|header_fasta|header_run\w*|seeaddbuf|seebuf|addbuf|skipbuf|loadbuf|loadmem|end_fasta|read_fasta|nextchar|esl_sqfile_Position|sqascii_\w+|esl_abc_Digitize\w*|header_embl|end_embl|check_buffers|esl_sq_Grow\w*|esl_sq_SetName|esl_sq_SetDesc|esl_sq_SetAccession|esl_sq_FetchFromMSA|esl_buffer_\w+|esl_sq_GrowTo)$', 'R seq DB read/parse'),
    (r'^(esl_workqueue_\w+|esl_threads_\w+)$',                         'T thread sync'),
]
INNER = [(re.compile(a), b) for a, b in INNER]
def classify(frames):
    # frames: leaf first
    indd = any(f in ('p7_domaindef_ByPosteriorHeuristics','rescore_isolated_domain','region_trace_ensemble') for f in frames)
    for f in frames:
        if f == 'forward_engine':  return 'D3 envelope Forward' if indd else 'S4 Forward parser (filter)'
        if f == 'backward_engine': return 'D3 envelope Backward' if indd else 'S5 Backward parser'
        if f in ('esl_vec_FNorm','p7_trace_Reverse','p7_trace_Append','esl_rnd_FChoose','esl_random','p7_trace_Grow') and indd: return 'D2 stochastic trace clustering'
        if f in ('p7_omx_Reuse','esl_vec_FSet','esl_vec_FSum') and not indd: return 'P per-seq setup/reconfig'
        for rx, st in INNER:
            if rx.match(f):
                return st
    if indd: return 'D9 domaindef other'
    if 'p7_Pipeline' in frames: return 'P pipeline glue (gumbel/exp/etc)'
    if 'thread_loop' in frames or 'serial_loop' in frames: return 'R seq DB read/parse'
    if 'pipeline_thread' in frames: return 'P worker loop other'
    return None

def main():
    by_thread = '--threads' in sys.argv
    tot = collections.Counter(); unk = collections.Counter(); thr = collections.Counter()
    n = 0
    frames = []; head = None
    def flush():
        nonlocal n
        if head is None: return
        tid, pid, comm = head
        st = classify(frames) if comm == 'hmmsearch' else f'R {comm} (child process, e.g. gzip -dc)'
        leaf = frames[0] if frames else '?'
        if st is None:
            # libc leaf: try to find a known caller deeper (already done) -> unknown
            st = 'X other/libc: ' + leaf
            unk[leaf] += 1
        role = ('main' if tid == pid else 'worker') if comm == 'hmmsearch' else 'child'
        key = (role, st) if by_thread else st
        tot[key] += 1; thr[role] += 1; n += 1
    for line in sys.stdin:
        if not line.strip():
            flush(); head = None; frames = []; continue
        if not line.startswith((' ', '\t')):
            m = re.match(r'\s*(\S+)\s+(\d+)/(\d+)', line)
            if m:
                head = (int(m.group(3)), int(m.group(2)), m.group(1))
            continue
        p = line.split()
        if len(p) >= 2:
            sym = p[1].split('+')[0]
            frames.append(sym)
    flush()
    print(f'samples: {n}   ' + '  '.join(f'{k}={v} ({100*v/n:.1f}%)' for k, v in thr.items()))
    for k, v in sorted(tot.items(), key=lambda kv: (-kv[1])):
        if 100*v/n < 0.05: continue
        lab = f'{k[0]:6s} {k[1]}' if by_thread else k
        print(f'  {100*v/n:6.2f}%  {lab}')

if __name__ == '__main__':
    main()
