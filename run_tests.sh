#!/bin/sh
# Runs every check suite in the repo. Needs Python 3.10+ and numpy.
# Usage: ./run_tests.sh            (uses python3 on PATH)
#        PYTHON=.venv/bin/python ./run_tests.sh
set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
status=0
for t in model/test_applicability.py model/test_strata.py \
         spec/test_ensemble.py prototype/eval-harness/test_harness.py; do
    echo "== $t"
    "$PY" "$t" || status=1
    echo
done
[ "$status" -eq 0 ] && echo "All suites passed." || echo "At least one suite failed."
exit "$status"
