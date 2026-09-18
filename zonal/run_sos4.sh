#!/bin/bash
# The sum-of-squares half of the four point constraint.
#
# The blocks are grouped by size: the seven large ones go one per process,
# because a single one of them fills four gigabytes, and the small ones share a
# process, because folding a partial result into the running total costs more
# than computing one of them.  Each group is merged into the running total and
# its file deleted.
cd /home/claude/aud
TOTAL=/tmp/claude-0/s4_total.pkl
DONE=/tmp/claude-0/s4_done.txt
touch $DONE

run_group () {          # $1 = "lo,hi"  $2 = tag
  if grep -qx "g$2" $DONE; then echo "skip group $2"; return; fi
  SOS_BLOCKS=$1 python3 verify45.py llm/LasserreSphericalCodes/proofs/4_24 \
      /tmp/claude-0/ps_all.txt 4 sos /tmp/claude-0/s4_part.pkl 2>&1 \
      | grep -E '^    block|Error|error'
  if [ ! -s /tmp/claude-0/s4_part.pkl ]; then
    echo "GROUP $2 FAILED"
    return
  fi
  python3 merge4.py $TOTAL /tmp/claude-0/s4_part.pkl
  rm -f /tmp/claude-0/s4_part.pkl
  echo "g$2" >> $DONE
}

run_group 2,3   b2       # 395x350, the largest
run_group 19,20 b19      # 171x139
run_group 20,21 b20      # 228x183
run_group 21,22 b21      # 171x129
run_group 6,13  s6       # the small ones, 13x3 up to 63x22
run_group 13,19 m13      # 29x26 up to 120x89
run_group 22,28 m22
run_group 28,45 s28
run_group 45,50 m45
echo "ALL GROUPS DONE"
python3 combine4.py $TOTAL /tmp/claude-0/z4.pkl
