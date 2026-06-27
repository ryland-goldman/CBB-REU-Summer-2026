#!/bin/bash
# Live dashboard for the LOCAL CNSGA optimizer (sim/optimize.py via local_run.sh).
#
#   watch -n 5 bash local_monitor.sh
#
# All local (no ssh). Shows run state, eval progress + ETA, best-so-far, per-stage in-flight evals,
# and host load / RAM / disk. Safe to run while the optimizer is running.

BASE=/Users/rylandgoldman/Documents/Coding/Cornell2
PY=~/miniforge3/envs/CBB/bin/python
RUNS=/tmp/linac_runs
LOG="$BASE/logs/opt/local_run.log"
cd "$BASE" 2>/dev/null || exit 1

TOTAL=$(awk '/^[[:space:]]*max_evaluations:/{print $2; exit}' config/xopt.yaml)
POP=$(awk '/^[[:space:]]*population_size:/{print $2; exit}' config/xopt.yaml)
MAXW=$(awk '/^[[:space:]]*max_workers:/{print $2; exit}' config/xopt.yaml)
NCPU=$(sysctl -n hw.ncpu)
hms(){ printf "%dh%02dm" $(( ${1:-0}/3600 )) $(( (${1:-0}%3600)/60 )); }

pid=$(pgrep -f 'sim/optimize.py' | head -1)
main=$(pgrep -f 'sim/main.py' | wc -l | tr -d ' ')

# elapsed from the optimize.py process start (NOT the log birth time -- `>` restarts reuse the inode,
# so the log's birth time persists across relaunches and would overstate elapsed)
el=0
if [ -n "$pid" ]; then
  st=$(ps -o lstart= -p "$pid" 2>/dev/null)
  s=$(date -j -f "%a %b %e %T %Y" "$st" +%s 2>/dev/null)
  [ -n "$s" ] && el=$(( $(date +%s) - s ))
fi

echo "===== Local CNSGA optimizer =====   $(date '+%a %H:%M:%S')"
if [ -n "$pid" ]; then
  echo "controller : pid $pid  RUNNING   elapsed $(hms $el)"
else
  echo "controller : NOT RUNNING   (restart: nohup bash local_run.sh > logs/opt/local_run.log 2>&1 &)"
fi
echo "config     : pop $POP x up to $((TOTAL/POP)) gens = $TOTAL evals   |   workers $main / $MAXW   |   cores $NCPU"
echo

# ---- eval progress / ETA / best-so-far from the checkpoint ----
EL=$el TOTAL=${TOTAL:-0} "$PY" - <<'PY'
import os, sys
ck=os.path.expanduser("~/Documents/Coding/Cornell2/logs/opt/data.csv")
el=int(os.environ.get("EL","0")); total=int(os.environ.get("TOTAL","0"))
def hms(s): s=int(max(s,0)); return "%dh%02dm"%(s//3600,(s%3600)//60)
if not os.path.exists(ck):
    print("evals      : 0 / %d   (no checkpoint yet -- first generation still running)"%total); sys.exit()
try:
    import numpy as np, pandas as pd
    d=pd.read_csv(ck); done=len(d)
    q=pd.to_numeric(d["q_out_C"],errors="coerce") if "q_out_C" in d else pd.Series([],dtype=float)
    fin=int(np.isfinite(q).sum()) if len(q) else 0
    rate=done/el if el>0 else 0; eta=(total-done)/rate if rate>0 else 0
    print("evals      : %d / %d  (%.0f%%)   valid %d   %.1f/hr   ETA ~%s"
          %(done,total,100*done/total if total else 0,fin,rate*3600,hms(eta)))
    if fin:
        sub=d[np.isfinite(q)]
        bq=sub.loc[pd.to_numeric(sub["q_out_C"],errors="coerce").idxmax()]
        be=sub.loc[pd.to_numeric(sub["eps_n"],errors="coerce").idxmin()]
        print("best       : max q_out_C=%.3e C (eps_n=%.2e ke=%.0f)   min eps_n=%.3e (q=%.2e ke=%.0f)"
              %(float(bq["q_out_C"]),float(bq["eps_n"]),float(bq["ke_out_mev"]),
                float(be["eps_n"]),float(be["q_out_C"]),float(be["ke_out_mev"])))
except Exception as e:
    print("evals      : checkpoint parse error:", type(e).__name__, e)
PY
echo

# ---- host load / RAM / disk / sandbox usage ----
load=$(sysctl -n vm.loadavg | awk '{print $2}')
mem=$(top -l 1 -n 0 2>/dev/null | awk -F'PhysMem: ' '/PhysMem/{print $2}')
echo "system     : load $load (of $NCPU cores)   |   mem ${mem:-n/a}"
echo "disk       : $(df -h "$BASE" | awk 'NR==2{print $4" free, "$5" used"}')   |   sandboxes $(du -sh "$RUNS" 2>/dev/null | cut -f1 || echo 0B)"
echo

# ---- in-flight evals by stage (.py-suffixed patterns are mutually exclusive) ----
echo "in-flight evals by stage:"
any=0
for s in injector linac1-4 autophase.py autophase_impact.py converter linac5-8; do
  c=$(pgrep -f "sim/$s" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$c" -gt 0 ]; then printf "   %-18s %d\n" "$(echo "$s" | sed 's/\.py$//')" "$c"; any=1; fi
done
[ "$any" = 0 ] && echo "   (workers between stages / idle)"
[ "$main" -gt 0 ] && printf "   %-18s %d\n" "evals in flight" "$main"
