/-
D4OmegaRegionI.lean: the branch and bound of Theorem 7.73 (region I of
Theorem 21.8) again, with the table of omega computed in Lean (D4Omega) in place
of the tables of D4CertData.  Settled by native_decide.
-/
import D4Omega

def verifyDomainLean : Bool := (verifyStatWith leanTabI 4000000).1 == 0

theorem domain_ok_lean_tables : verifyDomainLean = true := by native_decide

#print axioms domain_ok_lean_tables
