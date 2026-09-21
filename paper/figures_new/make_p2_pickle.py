"""
make_p2_pickle.py -- the exact coefficients of the two-point polynomial of the
certificate, read from multi_cap/llm24_out/llm24_p2.txt (written by
llm24_certificate_check.py, step 7) and stored as a dictionary
{degree: (numerator, denominator)} for fig_c_sigma2.py.
"""
import os
import pickle
import sys
from fractions import Fraction

sys.set_int_max_str_digits(2_000_000)
SRC = None
for cand in ("../../multi_cap/llm24_out/llm24_p2.txt",):
    if os.path.exists(cand):
        SRC = cand
        break
coef = {}
for line in open(SRC):
    if line.startswith("#") or not line.strip():
        continue
    i, c = line.split()
    f = Fraction(c)
    coef[int(i)] = (f.numerator, f.denominator)
assert sorted(coef) == list(range(17)), sorted(coef)
pickle.dump(coef, open("p2.pkl", "wb"))
print("p2.pkl written: degree", max(coef), "leading coefficient of", len(str(coef[16][0])), "digits")
