#!/bin/sh
# Check the six standalone Lean files with the pinned toolchain and print
# the axiom report of each.  Needs `lean` on the PATH (elan installs the
# version named in lean-toolchain).  The two Lake projects are built
# separately:  (cd cell600 && lake build)  and  (cd certificate && lake build).
set -e
cd "$(dirname "$0")"
for f in D4Stress D4Meet D4Certificate D4InnerProducts D4RootLattices D4Closure; do
  echo "== $f.lean"
  lean "$f.lean"
done
echo "all six files check"
