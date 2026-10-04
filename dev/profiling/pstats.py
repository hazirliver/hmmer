#!/usr/bin/env python3
# aggregate per-query "Internal pipeline statistics" blocks in hmmsearch output
import sys, re
q=0; nodes=0; tseq=0; res=0; msv=bias=vit=fwd=0; cpu=0.0
for l in open(sys.argv[1]):
    if l.startswith('Query model(s):'):
        q+=1; nodes+=int(re.search(r'\((\d+) nodes', l).group(1))
    elif l.startswith('Target sequences:'):
        m=re.search(r':\s+(\d+)\s+\((\d+) residues',l); tseq+=int(m.group(1)); res+=int(m.group(2))
    elif l.startswith('Passed MSV'): msv+=int(l.split()[3])
    elif l.startswith('Passed bias'): bias+=int(l.split()[3])
    elif l.startswith('Passed Vit'): vit+=int(l.split()[3])
    elif l.startswith('Passed Fwd'): fwd+=int(l.split()[3])
cells = None
print(f"queries={q} sumM={nodes} targets/query={tseq//max(q,1)} res/query={res//max(q,1)} "
      f"pass MSV={msv/tseq:.4f} bias={bias/tseq:.4f} Vit={vit/tseq:.5f} Fwd={fwd/tseq:.6f} "
      f"(n: {msv} {bias} {vit} {fwd})  total_cells={nodes*(res//max(q,1))/1e9:.1f}G")
