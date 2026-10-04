#!/usr/bin/env bash
# cells.sh <hmmsearch args...>   (needs root for bpftrace uprobes; uses sudo)
#   Counts calls and summed target length L handed to each pipeline stage.
#   For a single-model run, DP cells per stage = M * sum(L). Combine with the stage
#   times from profile_stage.sh to get in-situ throughput (cells/s).
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
H=${HMMSEARCH:-$ROOT/build-prof/src/hmmsearch}
BT=$(mktemp); trap 'rm -f "$BT"' EXIT
cat > "$BT" <<BTEOF
uprobe:$H:p7_MSVFilter       { @n["1 MSV"]=count();     @L["1 MSV"]=sum(arg1); }
uprobe:$H:p7_bg_FilterScore  { @n["2 bias"]=count();    @L["2 bias"]=sum(arg2); }
uprobe:$H:p7_ViterbiFilter   { @n["3 Vit"]=count();     @L["3 Vit"]=sum(arg1); }
uprobe:$H:p7_ForwardParser   { @n["4 FwdP"]=count();    @L["4 FwdP"]=sum(arg1); }
uprobe:$H:p7_BackwardParser  { @n["5 BckP"]=count();    @L["5 BckP"]=sum(arg1); }
uprobe:$H:p7_Forward         { @n["D3 Fwd"]=count();    @L["D3 Fwd"]=sum(arg1); }
uprobe:$H:p7_Backward        { @n["D3 Bck"]=count();    @L["D3 Bck"]=sum(arg1); }
uprobe:$H:p7_StochasticTrace { @n["D2 stotrace"]=count(); @L["D2 stotrace"]=sum(arg2); }
uprobe:$H:p7_OptimalAccuracy { @n["D6 OA"]=count(); }
uprobe:$H:p7_domaindef_ByPosteriorHeuristics { @n["D0 domaindef"]=count(); }
BTEOF
sudo bpftrace "$BT" -c "$H --cpu 0 -o /dev/null $*" 2>&1 | grep -v '^Attach\|^$'
