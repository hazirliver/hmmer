#!/usr/bin/env bash
# Derive the benchmark datasets from the pinned raw downloads.
#   prepare.sh <rawdir> <outdir> <easel-miniapps-dir>
# Deterministic: uses only awk/sort/gzip -n and Easel miniapps with fixed seeds,
# never HMMER programs, so the data cannot drift when HMMER code changes.
# Output checksums are pinned in dev/bench/data.sha256. Changing this file
# changes the data version (new derived-<hash> dir) -> re-pin + re-hash goldens.
set -euo pipefail
RAW=$1; O=$2; E=$3
export LC_ALL=C
trap 'echo "prepare.sh: failed at line $LINENO: $BASH_COMMAND" >&2' ERR
cd "$O"
log() { printf '  %-28s %s\n' "$1" "$(date +%T)"; }

# stride <k> <offset> : keep every k-th record (records end with a line "//")
stride_rec() { awk -v k="$1" -v o="${2:-0}" '{buf = buf $0 "\n"} /^\/\/$/ {if (n % k == o) printf "%s", buf; n++; buf=""}'; }
# stride over FASTA records
stride_fa()  { awk -v k="$1" -v o="${2:-0}" '/^>/{n++} (n-1) % k == o'; }
# keep FASTA records whose sequence length satisfies lo <= L <= hi
len_fa()     { awk -v lo="$1" -v hi="$2" '/^>/{if (h!="" && L>=lo && L<=hi) printf "%s%s", h, s; h=$0"\n"; s=""; L=0; next} {s=s $0 "\n"; L+=length($0)} END{if (h!="" && L>=lo && L<=hi) printf "%s%s", h, s}'; }
# fetch FASTA records by exact name, in the order given
names_fa()   { awk -v want="$1" 'BEGIN{n=split(want,w," "); for(i=1;i<=n;i++) ord[w[i]]=i} /^>/{name=substr($1,2); keep=(name in ord); if(keep) cur=name} keep{rec[cur]=rec[cur] $0 "\n"} END{for(i=1;i<=n;i++) printf "%s", rec[w[i]]}'; }
# HMM/Stockholm records by NAME / #=GF ID, in the order given
names_rec()  { awk -v want="$1" -v key="$2" 'BEGIN{n=split(want,w," "); for(i=1;i<=n;i++) ord[w[i]]=i} {buf=buf $0 "\n"} $0 ~ key {nm=$NF} /^\/\/$/ {if (nm in ord) rec[nm]=buf; buf=""; nm=""} END{for(i=1;i<=n;i++) printf "%s", rec[w[i]]}'; }
# HMM records with lo <= LENG <= hi, every k-th of those, max m
leng_rec()   { awk -v lo="$1" -v hi="$2" -v k="$3" -v m="$4" '{buf=buf $0 "\n"} /^LENG/{L=$2} /^\/\/$/ {if (L>=lo && L<=hi) {if (n%k==0 && c<m) {printf "%s", buf; c++}; n++}; buf=""}'; }

# ---------------------------------------------------------------- protein sequences
log "sprot (UniProt text)"
zcat "$RAW/sprot_dat.gz" > sprot.dat
log "sprot.fa"
"$E/esl-reformat" --informat embl fasta sprot.dat > sprot.fa
log "sprot subsets"
stride_rec 11  < sprot.dat > sprot_50k.dat; "$E/esl-reformat" --informat embl fasta sprot_50k.dat > sprot_50k.fa
gzip -n -6 -c sprot.fa > sprot.fa.gz
stride_rec 57  < sprot.dat > sprot_10k.dat; "$E/esl-reformat" --informat embl fasta sprot_10k.dat > sprot_10k.fa
stride_rec 571 < sprot.dat > sprot_1k.dat;  "$E/esl-reformat" --informat embl fasta sprot_1k.dat  > sprot_1k.fa
gzip -n -6 -c sprot_10k.fa > sprot_10k.fa.gz
gzip -n -6 -c sprot_50k.fa > sprot_50k.fa.gz
rm sprot_1k.dat
len_fa 3000 1000000 < sprot.fa > sprot_long.fa                 # ~very long proteins (titin, mucins...)
len_fa 1 60 < sprot.fa | stride_fa 2 > sprot_short.fa         # short peptides
log "human proteome (sprot OX 9606)"
awk '{buf=buf $0 "\n"} /^OX   NCBI_TaxID=9606[; ]/{h=1} /^\/\/$/ {if (h) printf "%s", buf; buf=""; h=0}' sprot.dat \
  | "$E/esl-reformat" --informat embl fasta - > human.fa
rm sprot.dat
stride_fa 2000 7 < human.fa > q_human_10.fa
stride_fa 200 3  < human.fa > q_human_100.fa
stride_fa 20 1   < human.fa > q_human_1000.fa
# hand-picked queries: tiny, typical, kinase (many hits), multidomain, huge (titin)
names_fa "INS_HUMAN HBA_HUMAN P53_HUMAN EGFR_HUMAN BRCA1_HUMAN" < sprot.fa > q_pick5.fa
names_fa "HBA_HUMAN" < sprot.fa > q_hba.fa
names_fa "EGFR_HUMAN" < sprot.fa > q_egfr.fa
names_fa "TITIN_HUMAN" < sprot.fa > q_titin.fa
log "decoys"
"$E/esl-shuffle" --seed 1 -d sprot_10k.fa | sed 's/^>/>dshuf_/' > sprot_10k_dishuffled.fa
"$E/esl-shuffle" --seed 2 -G --amino -N 100000 -L 350 > random_100k.fa

# ---------------------------------------------------------------- Pfam profiles
log "Pfam-A.hmm"
zcat "$RAW/pfam_hmm.gz" > Pfam-A.hmm
stride_rec 22  < Pfam-A.hmm > pfam_1k.hmm
stride_rec 220 < Pfam-A.hmm > pfam_100.hmm
stride_rec 2198 < Pfam-A.hmm > pfam_10.hmm
leng_rec 1   60   100 10 < Pfam-A.hmm > pfam_len_short.hmm
leng_rec 100 200  300 10 < Pfam-A.hmm > pfam_len_med.hmm
leng_rec 400 700  30  10 < Pfam-A.hmm > pfam_len_long.hmm
leng_rec 1000 9999 5  10 < Pfam-A.hmm > pfam_len_xlong.hmm
cat pfam_len_short.hmm pfam_len_med.hmm pfam_len_long.hmm pfam_len_xlong.hmm > pfam_len_mix.hmm
# families with very many hits: stress domain definition, alignment and output
names_rec "Pkinase 7tm_1 ABC_tran WD40 zf-C2H2 RRM_1 Ank_2 TPR_1 Response_reg HATPase_c" '^NAME ' < Pfam-A.hmm > pfam_hitters.hmm
# single models for kernel microbenchmarks
leng_rec 45 55 1 1     < Pfam-A.hmm > kern_M50.hmm
leng_rec 195 205 1 1   < Pfam-A.hmm > kern_M200.hmm
leng_rec 395 405 1 1   < Pfam-A.hmm > kern_M400.hmm
leng_rec 990 1010 1 1  < Pfam-A.hmm > kern_M1000.hmm

# ---------------------------------------------------------------- Pfam seed alignments
log "Pfam-A.seed"
zcat "$RAW/pfam_seed.gz" > Pfam-A.seed.sto
stride_rec 22  < Pfam-A.seed.sto > seed_1k.sto
stride_rec 220 < Pfam-A.seed.sto > seed_100.sto
# largest seeds by #=GF SQ
big=$(awk '/^#=GF ID/{id=$3} /^#=GF SQ/{print $3, id}' Pfam-A.seed.sto | sort -k1,1nr -k2,2 | awk 'NR<=20 {printf "%s ", $2}')
names_rec "$big" '^#=GF ID ' < Pfam-A.seed.sto > seed_big20.sto
# one large alignment in every MSA format HMMER reads
names_rec "LRR_1" '^#=GF ID ' < Pfam-A.seed.sto > msa_LRR_1.sto
for f in afa a2m psiblast clustal phylip selex; do "$E/esl-reformat" -o msa_LRR_1.$f $f msa_LRR_1.sto; done
# unaligned seqs of the biggest seed family, for hmmalign
names_rec "GFO_IDH_MocA_C3" '^#=GF ID ' < Pfam-A.seed.sto > msa_GFO.sto
"$E/esl-reformat" fasta msa_GFO.sto > seqs_GFO.fa

# ---------------------------------------------------------------- Rfam (RNA alignments)
log "Rfam seed"
zcat "$RAW/rfam_seed.gz" > Rfam.seed.sto
stride_rec 41 < Rfam.seed.sto > rfam_100.sto
names_rec "5S_rRNA tRNA SSU_rRNA_bacteria 6S RNaseP_bact_a tmRNA" '^#=GF ID ' < Rfam.seed.sto > rfam_pick.sto
rm Rfam.seed.sto

# ---------------------------------------------------------------- DNA
log "genomes"
zcat "$RAW/hg38_chr21.gz" > chr21.fa
zcat "$RAW/ecoli_k12.gz" | sed 's/^>\(\S*\).*/>\1/' > ecoli.fa
"$E/esl-sfetch" --index chr21.fa > /dev/null
"$E/esl-sfetch" -n chr21_30M-35M -c 30000001..35000000 chr21.fa chr21 > chr21_5M.fa
"$E/esl-sfetch" -n chr21_30M-31M -c 30000001..31000000 chr21.fa chr21 > chr21_1M.fa
for i in 0 1 2 3 4 5 6 7 8 9; do
  s=$((15000001 + i * 2500000)); "$E/esl-sfetch" -n chr21_${s}_100k -c $s..$((s + 99999)) chr21.fa chr21
done > q_chr21_10x100k.fa
"$E/esl-sfetch" -n chr21_query_6k -c 33000001..33006000 chr21.fa chr21 > q_chr21_6k.fa
"$E/esl-sfetch" -n chr21_query_2k -c 25000001..25002000 chr21.fa chr21 > q_chr21_2k.fa       # few hits
"$E/esl-sfetch" -n chr21_query_2k_rep -c 20000001..20002000 chr21.fa chr21 > q_chr21_2k_rep.fa  # repeat-rich: many hits
rm chr21.fa.ssi
gzip -n -6 -c chr21_5M.fa > chr21_5M.fa.gz
"$E/esl-shuffle" --seed 3 -G --dna -N 10 -L 1000000 | sed 's/^>/>rnd_/' > random_dna_10M.fa
log "Dfam (human curated families)"
zcat "$RAW/dfam_curated.gz" | awk '{buf=buf $0 "\n"} /^TH +TaxId:9606;/{h=1} /^\/\/$/ {if (h) printf "%s", buf; buf=""; h=0}' > dfam_human.hmm
stride_rec 14 < dfam_human.hmm > dfam_100.hmm
names_rec "AluY AluSx AluJb MIR THE1B MER5A MLT1A0 L1HS_3end L1PA2_3end LTR12C SVA_A Tigger1 Charlie1 L1HS_5end HERVK" '^NAME ' < dfam_human.hmm > dfam_pick15.hmm
names_rec "AluY" '^NAME ' < dfam_human.hmm > dfam_AluY.hmm
names_rec "L1HS_5end" '^NAME ' < dfam_human.hmm > dfam_L1HS_5end.hmm
log "done"
