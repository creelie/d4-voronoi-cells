#!/bin/bash
# The sum-of-squares half of the four point constraint.
#
# The blocks are grouped by size: the seven large ones go one per process,
# because a single one of them fills four gigabytes, and the small ones share a
# process, because folding a partial result into the running total costs more
# than computing one of them.  Each group is merged into the running total and
# its file deleted.
#
#   bash run_sos4.sh /path/to/proofs/4_24 ps.txt [workdir]
#
# ps.txt is the output of psker; workdir (default: a directory work4/ beside
# this script) holds the running total and the list of finished groups, so an
# interrupted run resumes where it stopped.  The zonal half of the constraint
# is  python3 verify45.py DATA ps.txt 4 zonal WORK/z4.pkl , and
# python3 combine4.py WORK/s4_total.pkl WORK/z4.pkl  adds the two halves and
# reports whether the constraint vanishes.
set -u
cd "$(dirname "$0")"
DATA=${1:?usage: run_sos4.sh DATA ps.txt [workdir]}
PS=${2:?usage: run_sos4.sh DATA ps.txt [workdir]}
WORK=${3:-work4}
mkdir -p "$WORK"
TOTAL=$WORK/s4_total.pkl
DONE=$WORK/s4_done.txt
PART=$WORK/s4_part.pkl
touch "$DONE"

run_group () {          # $1 = "lo,hi"  $2 = tag
  if grep -qx "g$2" "$DONE"; then echo "skip group $2"; return; fi
  SOS_BLOCKS=$1 python3 verify45.py "$DATA" "$PS" 4 sos "$PART" 2>&1 \
      | grep -E '^    block|Error|error'
  if [ ! -s "$PART" ]; then
    echo "GROUP $2 FAILED"
    return
  fi
  python3 merge4.py "$TOTAL" "$PART"
  rm -f "$PART"
  echo "g$2" >> "$DONE"
}

# The groups must cover the fifty sum-of-squares blocks 0..49 exactly once.
SOS_GROUPS="2,3:b2 1,2:b1 3,4:b3 4,5:b4 19,20:b19 20,21:b20 21,22:b21
        0,1:s0 5,6:s5 6,13:s6 13,19:m13 22,28:m22 28,45:s28 45,50:m45"
python3 - $SOS_GROUPS <<'EOF' || exit 1
import sys
seen = []
for g in sys.argv[1:]:
    lo, hi = map(int, g.split(':')[0].split(','))
    seen += range(lo, hi)
if sorted(seen) != list(range(50)):
    sys.exit('run_sos4.sh: the groups do not cover blocks 0..49 exactly once')
EOF

run_group 2,3   b2       # 395x350, the largest
run_group 1,2   b1       # 312x280
run_group 3,4   b3       # 292x250
run_group 4,5   b4       # 195x164
run_group 19,20 b19      # 171x139
run_group 20,21 b20      # 228x183
run_group 21,22 b21      # 171x129
run_group 0,1   s0       # 103x89
run_group 5,6   s5       # 13x3
run_group 6,13  s6       # the small ones, 13x3 up to 63x22
run_group 13,19 m13      # 29x26 up to 120x89
run_group 22,28 m22
run_group 28,45 s28
run_group 45,50 m45
echo "ALL GROUPS DONE"
