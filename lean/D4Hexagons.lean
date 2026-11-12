/-
D4Hexagons.lean: the finite and algebraic content of lem:hexagon-identity,
thm:hexagon-rigidity and cor:hexagon-close.

  * The 24 integral roots (+-1, +-1, 0, 0) in the numbering of the paper's
    appendix, the 96 pairs T with integral inner product 1 (inner product 1/2
    after normalisation), and the 16 hexagons a, b, b - a, -a, -b, a - b.
    Proved by `decide`: every hexagon satisfies u_{k-1} + u_{k+1} = u_k, its
    six consecutive pairs lie in T, the 96 hexagon pairs are distinct and so
    partition T, and every root lies in exactly four hexagons.
  * Over every commutative ring (`grind`): the six-cycle identities behind
    -1/2 sum e^2 <= Psi_H <= 3/2 sum e^2, and the end value of the convex
    bound in step (i) of the proof of thm:hexagon-rigidity.
  * Exact integer arithmetic: F(1/48, 17/20000) < 1/48 and F(7/100, 9/5000) <
    7/100 for the first form of F in thm:hexagon-rigidity, the two constants
    that compare it with sqrt480 s + (sqrt325/2) eta^2, the radius 2/sqrt325,
    and 2/sqrt(1 - 2s) > 2.0017 at s = 17/20000.

No Mathlib; nothing uses `native_decide`.
-/

namespace D4Hexagons

/-! ## Roots, pairs and hexagons -/

/-- The roots `eps_i e_i + eps_j e_j`, pairs `i < j` lexicographic, signs
`(+,+), (+,-), (-,+), (-,-)`; root `4k + l`. -/
def roots : List (List Int) :=
  [[1,1,0,0],[1,-1,0,0],[-1,1,0,0],[-1,-1,0,0],
   [1,0,1,0],[1,0,-1,0],[-1,0,1,0],[-1,0,-1,0],
   [1,0,0,1],[1,0,0,-1],[-1,0,0,1],[-1,0,0,-1],
   [0,1,1,0],[0,1,-1,0],[0,-1,1,0],[0,-1,-1,0],
   [0,1,0,1],[0,1,0,-1],[0,-1,0,1],[0,-1,0,-1],
   [0,0,1,1],[0,0,1,-1],[0,0,-1,1],[0,0,-1,-1]]

def dot (a b : List Int) : Int := (List.zipWith (· * ·) a b).foldl (· + ·) 0
def add (a b : List Int) : List Int := List.zipWith (· + ·) a b
def root (i : Nat) : List Int := roots.getD i []

/-- The pairs `{i, j}`, `i < j`, of roots at inner product `1`. -/
def T : List (Nat × Nat) :=
  (List.range 24).flatMap fun i =>
    ((List.range 24).filter fun j => i < j && dot (root i) (root j) == 1).map fun j => (i, j)

theorem roots_length : roots.length = 24 := by decide
theorem roots_norm : roots.all (fun r => dot r r == 2) = true := by decide
theorem T_length : T.length = 96 := by decide

/-- The sixteen hexagons, each in cyclic order `a, b, b - a, -a, -b, a - b`. -/
def hexagons : List (List Nat) :=
  [[0,4,14,3,7,13],[0,5,15,3,6,12],[0,8,18,3,11,17],[0,9,19,3,10,16],
   [1,4,12,2,7,15],[1,5,13,2,6,14],[1,8,16,2,11,19],[1,9,17,2,10,18],
   [4,8,22,7,11,21],[4,9,23,7,10,20],[5,8,20,6,11,23],[5,9,21,6,10,22],
   [12,16,22,15,19,21],[12,17,23,15,18,20],[13,16,20,14,19,23],[13,17,21,14,18,22]]

def at6 (h : List Nat) (k : Nat) : Nat := h.getD (k % 6) 0

/-- `u_{k-1} + u_{k+1} = u_k` for every `k` mod 6. -/
def cyclicRelation (h : List Nat) : Bool :=
  (List.range 6).all fun k => add (root (at6 h (k + 5))) (root (at6 h (k + 1))) == root (at6 h k)

def norm (p : Nat × Nat) : Nat × Nat := if p.1 ≤ p.2 then p else (p.2, p.1)

def hexPairs (h : List Nat) : List (Nat × Nat) :=
  (List.range 6).map fun k => norm (at6 h k, at6 h (k + 1))

def allHexPairs : List (Nat × Nat) := hexagons.flatMap hexPairs

theorem hexagons_length : hexagons.length = 16 := by decide
theorem hexagons_cyclic : hexagons.all cyclicRelation = true := by decide
theorem hexPairs_in_T : allHexPairs.all (fun p => T.contains p) = true := by decide
theorem hexPairs_count : allHexPairs.length = 96 := by decide
theorem hexPairs_distinct : allHexPairs.Nodup := by decide
theorem each_root_in_four :
    (List.range 24).all (fun i => (hexagons.filter fun h => h.contains i).length == 4) = true := by
  decide

/-! ## Identities over every commutative ring -/

section Identities
open Lean.Grind
variable {α : Type} [CommRing α]

/-- `2 Psi_H = sum x_k^2 - 2 sum x_k x_{k+1}` in one coordinate; then
`3 sum x^2 - 2 Psi_H` is a sum of squares, which is `Psi_H <= 3/2 sum e^2`. -/
theorem psi_upper (x0 x1 x2 x3 x4 x5 : α) :
    3*(x0^2+x1^2+x2^2+x3^2+x4^2+x5^2)
      - ((x0^2+x1^2+x2^2+x3^2+x4^2+x5^2) - 2*(x0*x1+x1*x2+x2*x3+x3*x4+x4*x5+x5*x0))
    = (x0+x1)^2+(x1+x2)^2+(x2+x3)^2+(x3+x4)^2+(x4+x5)^2+(x5+x0)^2 := by
  grind

/-- `2 Psi_H + sum x^2` is a sum of squares, which is `Psi_H >= -1/2 sum e^2`. -/
theorem psi_lower (x0 x1 x2 x3 x4 x5 : α) :
    ((x0^2+x1^2+x2^2+x3^2+x4^2+x5^2) - 2*(x0*x1+x1*x2+x2*x3+x3*x4+x4*x5+x5*x0))
      + (x0^2+x1^2+x2^2+x3^2+x4^2+x5^2)
    = (x0-x1)^2+(x1-x2)^2+(x2-x3)^2+(x3-x4)^2+(x4-x5)^2+(x5-x0)^2 := by
  grind

/-- Step (i): the end value of `6 s^2 - 2 s Y + Y^2` at `Y = 6 s + P`. -/
theorem endpoint (s P : α) :
    6*s^2 - 2*s*(6*s+P) + (6*s+P)^2 = 30*s^2 + 10*s*P + P^2 := by
  grind

/-- Step (i): with `y_p = s - x_p` and `Y = sum y_p`, `sum x_p^2 = 6 s^2 - 2 s Y + sum y_p^2`. -/
theorem squares_shift (s y0 y1 y2 y3 y4 y5 : α) :
    (s-y0)^2+(s-y1)^2+(s-y2)^2+(s-y3)^2+(s-y4)^2+(s-y5)^2
      = 6*s^2 - 2*s*(y0+y1+y2+y3+y4+y5) + (y0^2+y1^2+y2^2+y3^2+y4^2+y5^2) := by
  grind

end Identities

/-! ## Exact arithmetic

For `eta = p/q` and `s = u/v` put `D = 4 q^4 v^2`,
`A D = 4 (480 u^2 q^4 + 60 u v p^2 q^2 + 36 p^4 v^2)` (`A = 480 s^2 + 60 s eta^2 + 36 eta^4`),
`B D = 4 p^2 q^2 v^2 - p^4 v^2` (`B = eta^2 - eta^4/4`) and `E D = 4 p^4 v^2` (`E = eta^4`).
`F(eta, s) < eta` is `(sqrt A + 3 eta^2)^2 < B`, that is `sqrt A < sqrt B - 3 eta^2`;
it holds when `B > 9 E`, `R = B + 9 E - A > 0` and `36 E B < R^2`. -/

def hexCheck (p q u v : Int) : Bool :=
  let A := 4*(480*u^2*q^4 + 60*u*v*p^2*q^2 + 36*p^4*v^2)
  let B := 4*p^2*q^2*v^2 - p^4*v^2
  let E := 4*p^4*v^2
  let R := B + 9*E - A
  decide (B > 9*E) && decide (R > 0) && decide (36*E*B < R*R)

/-- `F(1/48, 17/20000) < 1/48`: `d(W) /= 1/48` for `s <= 17/20000`. -/
theorem eta_1_48 : hexCheck 1 48 17 20000 = true := by decide
/-- `F(7/100, 9/5000) < 7/100`: `d(W) /= 7/100` for `s <= 9/5000`. -/
theorem eta_7_100 : hexCheck 7 100 9 5000 = true := by decide
/-- The constant `17/20000` is close to sharp for this check. -/
theorem eta_1_48_sharp : hexCheck 1 48 18 20000 = false := by decide

/-- `60 <= 12 sqrt480` and `18 <= sqrt325`, `9^2 + 1/4 = 325/4`: the comparison of
`F` with `sqrt480 s + (sqrt325/2) eta^2`. -/
theorem compare_constants : 60^2 ≤ 12^2*480 ∧ 18^2 ≤ 325 ∧ 4*81 + 1 = 325 := by decide

/-- `0.110940 < 2/sqrt325 < 0.110941`. -/
theorem radius_325 : 325*110940^2 < 4*10^12 ∧ 4*10^12 < 325*110941^2 := by decide

/-- `2/sqrt(1 - 2 s) > 2.0017` at `s = 17/20000`: `4 > 2.0017^2 (1 - 17/10000)`. -/
theorem radius_close : 4*10^8*10^4 > 20017^2*(10^4 - 17) := by decide

end D4Hexagons

#print axioms D4Hexagons.hexPairs_distinct
#print axioms D4Hexagons.psi_upper
#print axioms D4Hexagons.eta_1_48
