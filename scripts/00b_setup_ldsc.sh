#!/usr/bin/env bash
# Install ldsc WITHOUT conda, into the currently active python env (e.g. your (ds) env).
# Uses the Python-3 port belowlab/ldsc (UNVERIFIED by me: run the check at the end before trusting results).
set -euo pipefail
D="${1:-$HOME/ldsc}"
[ -d "$D" ] || git clone https://github.com/belowlab/ldsc.git "$D"
pip install --default-timeout=60 numpy scipy pandas bitarray pybedtools || pip install --default-timeout=60 numpy scipy pandas bitarray
python "$D/ldsc.py" -h | head -5
python "$D/munge_sumstats.py" -h | head -3
echo "export LDSC_PY=\"python $D/ldsc.py\""
echo "export MUNGE_PY=\"python $D/munge_sumstats.py\""
