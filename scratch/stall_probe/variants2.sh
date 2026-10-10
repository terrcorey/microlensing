cd "$(dirname "$0")"; export PYTHONWARNINGS=ignore
run() { printf "%-52s " "$*"; timeout 15 /sharedscratch/gyc1/microlensing/.venv/bin/python -W ignore one.py "$@" 2>/dev/null | tail -1; [ ${PIPESTATUS[0]} = 124 ] && echo "HANG >15s"; }
q=5.623413e-4
for s in 2 3 4 5 6 7 8 9 10; do x=$(python3 -c "print(-$s*$q/(1+$q) + 5.479e-5)"); run $s $q $x 3.94096e-05 9.682674e-05 1e-3 1e-2; done
echo "-- s=10, RelTol scan"; x=-0.00557321
for r in 2e-3 3e-3 5e-3; do run 10 $q $x 3.94096e-05 9.682674e-05 $r 1e-2; done
echo "-- s=10 farther from the primary (same y)"; for dx in 2e-4 5e-4 1e-3 3e-3; do x=$(python3 -c "print(-10*$q/(1+$q) + $dx)"); run 10 $q $x 3.94096e-05 9.682674e-05 1e-3 1e-2; done
