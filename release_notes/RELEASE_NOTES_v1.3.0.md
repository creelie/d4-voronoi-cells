# v1.3.0

Supplementary package for *The Sphere Packing Problem in Dimension 4 and the
Twenty-Four-Cell Conjecture*, revised so that the twenty-four-contact case is
proved inside the paper.

## The manuscript (paper/)

- The classification of 24-point contact configurations, previously cited
  from the preprint of de Laat, Leijenhorst and de Muinck Keizer, is now
  proved in the text (Theorem "Twenty-four points"): the bound of the second
  level of the Lasserre hierarchy with its equality case (Lemma "The bound of
  the second level"), the positivity of Gram kernels (Lemma "Gram kernels are
  positive"), the certificate kernel written out from the zonal matrices, and
  Proposition "What the verification establishes", the four finite properties
  of the deposited data that the argument consumes.  What is taken from that
  work is the certificate itself, which is data, and the construction of the
  kernels, which is a definition; no statement of the preprint is used.
- Abstract, introduction, the account of the contact-count argument, the
  conclusion and the data availability statement are updated accordingly.
- A seventh Lean verification, lean/D4RootLattices.lean: the finite content
  of the lemma that carries the four inner products to the D_4 root system
  (the census of norm-2 Gram matrices in rank at most 4, the exact box, the
  counts 2, 6, 12, 24 and a D_4 basis in every lattice with 24 roots), by
  native_decide.  lean/run_all.sh checks the five standalone files; the log
  of the run behind this release is lean/runs/run_all_2026-09-21.log.
- Three figures are new (the two-point polynomial and the certificate's
  blocks; the twenty-four directions and their Gram matrix; the free-volume
  identity on real geometry) and seven were redrawn from their data so that no
  label touches a curve, a bar, a legend or another label.  Every figure
  script runs a collision test before it writes (paper/figures_new/figstyle.py)
  and reproduces the shipped PNG byte for byte.
- Every DOI of the bibliography was resolved and compared with the cited
  metadata (authors, title, venue, volume, issue, year, pages); none needed a
  change.

## The package

- paper/: D4.tex, figures/, figures_new/.  D4.pdf is delivered separately
  (it is 21 MB); copy it into paper/ before tagging.
- third_party/llm24-certificate/: the deposited certificate under its MIT
  licence, with fetch_certificate.py (MD5 02acd5270f7b3fa799abdeb5291706fd)
  and verify.sh.
- multi_cap/p2_zeroset_check.py: exact-arithmetic check of the zeros of the
  two-point polynomial, independent of python-flint.
- CITATION.cff and .zenodo.json carry version v1.3.0.

## Before tagging

1. In third_party/llm24-certificate/ run `python3 fetch_certificate.py` and
   then `python3 fetch_certificate.py --split`; the first writes
   LasserreSphericalCodes.zip and checks the MD5, the second cuts it into
   80 MB parts, which are what gets committed (GitHub refuses single files
   above 100 MB, and Zenodo archives the repository snapshot, not release
   assets).  `python3 fetch_certificate.py --join` restores the archive.
2. Tag v1.3.0 and publish the release; Zenodo mints a new version DOI under
   the concept record 10.5281/zenodo.22766562.  The paper cites the concept
   DOI, which resolves to the latest version, so nothing in the paper has to
   change after minting; if a version DOI is wanted in the Data availability
   statement, add it there once it exists.

## Archive

Published on Zenodo: version DOI 10.5281/zenodo.22880335, under the concept record 10.5281/zenodo.22766562. The certificate travels in the repository as third_party/llm24-certificate/LasserreSphericalCodes.zip.part-00 and .part-01; fetch_certificate.py --join restores the archive and checks the MD5.
