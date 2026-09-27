/-
D4OmegaMain.lean: the checks of D4LabelledMain again, with nothing taken from
outside Lean but the certificate: the tables of omega are computed in D4Omega
from the closed forms, the slab constant and the enclosures of the cell integrals
of the data are checked against values computed there, and so is the comparison
of A_* with the bound of the certificate of Theorem 7.73.  Settled by
native_decide.  D4OmegaRegionI.lean does the same for region I (Theorem 7.73).
-/
import D4Omega

/-- A_* > 8 - 0.0928555703 (the bound of the certificate, D4Certificate.lean); s(D) > 8 - A_*;
fr(1, 1/2) < 1/22; the slab constant of the data at most 4000 r kappa; and the 256 cell
integrals of the data enclosing the values computed here. -/
theorem omega_constants : labelledConstantsOk = true := by native_decide

/-- The two tables: increasing grids of step 2^-16 up to 1/2 and up to a_D, a positive
bound on |omega''|, and d1lo <= d1hi throughout. -/
theorem omega_tables : (tableOk leanTabI ⟨1, 1⟩ && tableOk leanTabL aDtop) = true := by native_decide

def verifySlabLean : Bool := (statSlabWith leanTabL 20000000).1 == 0
def verifyGammaLean : Bool := (statGammaWith leanTabL 20000000).1 == 0

/-- II_s with the tables of omega computed in Lean. -/
theorem region_IIs_lean : verifySlabLean = true := by native_decide

/-- II_f with the tables of omega computed in Lean. -/
theorem region_IIf_lean : verifyGammaLean = true := by native_decide

#eval (leanTabI.us.size, leanTabL.us.size, leanTabI.m2.toRat.num.toFloat / leanTabI.m2.toRat.den.toFloat,
       leanTabL.m2.toRat.num.toFloat / leanTabL.m2.toRat.den.toFloat)
#eval (statSlabWith leanTabL 20000000, statGammaWith leanTabL 20000000)
#print axioms region_IIf_lean
