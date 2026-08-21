#!/usr/bin/env bash
# Wait for GMrepo to come back, verify the profile payload actually parses,
# then backfill genus-level abundances.
#
# Staged on purpose: the response shape of getFullTaxonomicProfileByRunID has
# never been observed, so a small batch has to prove itself before ~11,870
# requests are sent. The raw payload is saved for after-the-fact review.
set -uo pipefail
cd /Users/lucyshaa/Desktop/CML
source .venv/bin/activate

# Credentials are read from .env (which is gitignored) so this script can be
# committed without carrying a live database password.
set -a; . ./.env; set +a
MYSQL="mysql -h${GUTDB_HOST} -P${GUTDB_PORT} -u${GUTDB_USER} -p${GUTDB_PASSWORD} ${GUTDB_NAME}"

RAW=/tmp/gmrepo_raw_profile.json
MAX_ATTEMPTS=160     # 160 x 90s ~= 4 hours
PROBE_N=20

echo "=== stage 1: wait for GMrepo ==="
UP=0
for i in $(seq 1 $MAX_ATTEMPTS); do
  code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 20 \
    -X POST "https://gmrepo.humangut.info/api/get_all_phenotypes" \
    -H "Content-Type: application/json" -d '{}')
  if [ "$code" = "200" ]; then echo "$(date +%H:%M:%S) GMrepo UP after $i attempts"; UP=1; break; fi
  [ $((i % 10)) -eq 0 ] && echo "$(date +%H:%M:%S) attempt $i -> HTTP $code"
  sleep 90
done
if [ "$UP" -ne 1 ]; then echo "GMREPO_NEVER_CAME_UP"; exit 1; fi

echo ""
echo "=== stage 2: capture a raw profile payload for review ==="
RUN=$($MYSQL -N -e \
  "SELECT s.run_accession FROM samples s
   LEFT JOIN sample_taxon_abundances a ON a.sample_id=s.id
   WHERE a.sample_id IS NULL AND s.sex IS NOT NULL AND s.age_years IS NOT NULL
   GROUP BY s.id, s.run_accession ORDER BY s.id LIMIT 1;" 2>/dev/null)
echo "probe run: $RUN"
curl -s --max-time 60 -X POST "https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID" \
  -H "Content-Type: application/json" -d "{\"run_id\":\"$RUN\"}" -o "$RAW"
echo "saved raw payload -> $RAW ($(wc -c < "$RAW") bytes)"
python3 -c "
import json
try:
    d=json.load(open('$RAW'))
    print('top-level type:', type(d).__name__)
    print('top-level keys:', sorted(d)[:15] if isinstance(d,dict) else 'list')
except Exception as e:
    print('payload is not JSON:', str(e)[:120])
"

echo ""
echo "=== stage 3: probe load ($PROBE_N runs, genus) ==="
BEFORE=$($MYSQL -N -e \
  "SELECT COUNT(*) FROM sample_taxon_abundances;" 2>/dev/null)
python populate_database.py sync-gmrepo-abundances \
  --rank genus --only-complete-demographics --limit $PROBE_N --workers 3 \
  --output data/gmrepo_run_abundances_probe.csv
PROBE_RC=$?
AFTER=$($MYSQL -N -e \
  "SELECT COUNT(*) FROM sample_taxon_abundances;" 2>/dev/null)
echo "probe exit=$PROBE_RC  abundance rows: $BEFORE -> $AFTER"

if [ "$PROBE_RC" -ne 0 ] || [ "$AFTER" -le "$BEFORE" ]; then
  echo "PROBE_FAILED - not scaling up. Inspect $RAW."
  exit 2
fi

# sanity gate: every newly profiled sample must sum to ~1.0
BAD=$($MYSQL -N -e \
  "SELECT COUNT(*) FROM (
     SELECT sample_id, SUM(relative_abundance) s
     FROM sample_taxon_abundances GROUP BY sample_id
     HAVING ABS(s-1.0) > 0.01
   ) t;" 2>/dev/null)
echo "samples whose abundances do NOT sum to 1.0: $BAD"
if [ "$BAD" -ne 0 ]; then echo "SANITY_GATE_FAILED - not scaling up."; exit 3; fi

echo ""
echo "=== stage 4: full genus backfill for the demographically-complete cohort ==="
python populate_database.py sync-gmrepo-abundances \
  --rank genus --only-complete-demographics --workers 4 \
  --output data/gmrepo_run_abundances.csv
echo "stage 4 exit=$?"

echo ""
echo "=== final state ==="
$MYSQL -e "
SELECT COUNT(*) total_samples, SUM(has_demographics) with_demo,
       SUM(has_features) with_features, SUM(ml_ready) ml_ready
FROM v_ml_ready_samples;" 2>/dev/null
