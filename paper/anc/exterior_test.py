"""Exact exterior-cut test in four dimensions, using the standard library."""

from fractions import Fraction as Q
from itertools import combinations, product


ZERO = (Q(0),) * 4
AXES = tuple(tuple(Q(i == j) for i in range(4)) for j in range(4))


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def shift(u, v, amount):
    return tuple(a + amount * b for a, b in zip(u, v))


def scale(v, amount):
    return tuple(amount * x for x in v)


def sign(x):
    return (x > 0) - (x < 0)


def radical_sign(a, b, d):
    """Sign of a + b sqrt(d), for rational coefficients and d >= 0."""
    if b == 0 or d == 0:
        return sign(a)
    if a == 0 or sign(a) == sign(b):
        return sign(b)
    return sign(a) * sign(a * a - b * b * d)


def affine_space(equations):
    """A point and orthogonal normals for linear equations."""
    point = ZERO
    normals = []
    for row, rhs in equations:
        residual = rhs - dot(row, point)
        normal = project(row, normals)
        if normal == ZERO:
            if residual:
                return None
            continue
        length = dot(normal, normal)
        point = shift(point, normal, residual / length)
        normals.append((normal, length))
    return point, tuple(normals)


def project(row, normals):
    for normal, length in normals:
        coefficient = dot(row, normal)
        if coefficient:
            row = shift(row, normal, -coefficient / length)
    return row


def stratum(equations):
    """Return a point, or a sphere in an affine space, for a boundary subset."""
    planes = [(b, c) for a, b, c in equations if not a]
    spheres = [(scale(b, 1 / a), c / a) for a, b, c in equations if a]
    if spheres:
        row, rhs = spheres[0]
        planes.extend((shift(b, row, -1), c - rhs) for b, c in spheres[1:])
    affine = affine_space(planes)
    if affine is None:
        return None
    point, normals = affine
    if not spheres:
        return (point, (), Q(0)) if len(normals) == 4 else None
    centre = scale(row, Q(-1, 2))
    difference = shift(centre, point, -1)
    tangent = project(difference, normals)
    perpendicular = shift(difference, tangent, -1)
    radius2 = rhs + dot(centre, centre) - dot(perpendicular, perpendicular)
    if radius2 < 0 or (radius2 > 0 and len(normals) == 4):
        return None
    return shift(point, tangent, 1), normals if radius2 else (), radius2


class ExteriorTest:
    def __init__(self, shell):
        # Store constraints as a |p|^2 + b . p <= c.
        self.constraints = ((Q(1), ZERO, Q(2, 3)),) + tuple(
            (4 - dot(y, y), scale(y, Q(4)), Q(4)) for y in shell
        )
        self.strata = set()
        for size in range(1, 5):
            for equations in combinations(self.constraints, size):
                data = stratum(equations)
                if data is not None:
                    self.strata.add(data)

    def feasible(self, p, v, d):
        # Candidate p + v sqrt(d). Squaring uses the same radical throughout.
        square = dot(p, p) + dot(v, v) * d
        cross = 2 * dot(p, v)
        for a, b, c in self.constraints:
            if radical_sign(c - a * square - dot(b, p), -a * cross - dot(b, v), d) < 0:
                return False
        return True

    def cutting_witness(self, x):
        """Return p, v, d for a feasible cut p + v sqrt(d), or None."""
        for p, normals, radius2 in self.strata:
            if radius2 == 0:
                directions = (ZERO,)
                d = Q(0)
            else:
                for functional in (x, *AXES):
                    v = project(functional, normals)
                    if v != ZERO:
                        break
                d = radius2 / dot(v, v)
                directions = (v, scale(v, Q(-1))) if len(normals) == 3 else (v,)
            excess = dot(x, p) - 1
            for v in directions:
                if radical_sign(excess, dot(x, v), d) > 0 and self.feasible(p, v, d):
                    return p, v, d
        return None


def verify_packing(shell):
    assert all(4 <= dot(y, y) < 6 for y in shell)
    assert all(dot(d, d) >= 4 for y, z in combinations(shell, 2)
               for d in (shift(y, z, -1),))


def rational_roots():
    return tuple(scale(axis, Q(s)) for axis in AXES for s in (-1, 1)) + tuple(
        tuple(Q(s, 2) for s in signs) for signs in product((-1, 1), repeat=4)
    )


def layer_shell(n):
    equator = []
    for i, j in combinations(range(3), 2):
        for s, t in product((-1, 1), repeat=2):
            equator.append(shift(scale(AXES[i], Q(17 * s, 10)), AXES[j], Q(17 * t, 10)))
    upper = [(*map(Q, signs), Q(17, 10)) for signs in product((-1, 1), repeat=3)]
    lower = [(*map(Q, signs), Q(-17, 10)) for signs in product((-1, 1), repeat=3)]
    return tuple(equator + upper + lower[:n - 20])


def run():
    a = Q(1001, 1000)
    c = Q(99999, 100000)
    shell = tuple(scale(u, 2 * a) for u in rational_roots())
    verify_packing(shell)
    test = ExteriorTest(shell)
    # Rational coordinates: Q is cut out by |x_i| <= 1 and sum |x_i| <= 2.
    # Truncate each vertex of aQ at x_i +/- x_j = 2c. The clipped
    # vertex below is one representative of a single orbit.
    vertex = (a, 2 * c - a, 2 * (a - c), Q(0))
    assert all(2 * dot(vertex, y) <= dot(y, y) for y in shell)
    assert test.cutting_witness(vertex) is None
    assert test.cutting_witness((a, a, Q(0), Q(0))) is not None
    lower = 8 * a**4 - 192 * (a-c)**4
    assert lower > Q(8032, 1000)
    # The protected polytope fits inside B(2); apply the stability theorem.
    delta = Q(1, 2000)
    assert (1 - Q(3, 2) * delta)**4 * lower > Q(8007, 1000)
    print('24 neighbours: clipped orbit protected; untrimmed vertex cut', flush=True)
    print('24 neighbours: exact clipped volume > 8.032; perturbations of 1/2000 give > 8.007',
          flush=True)

    for n in (25, 26, 27):
        shell = layer_shell(n)
        verify_packing(shell)
        test = ExteriorTest(shell)
        near = scale(AXES[3], Q(5, 4))
        far = scale(AXES[3], Q(3, 2))
        assert all(2 * dot(far, y) <= dot(y, y) for y in shell)
        assert test.cutting_witness(near) is None
        assert test.cutting_witness(far) is not None
        # A ball of radius 6/5 survives perturbations up to 1/500.
        assert all(dot(y, y) >= (Q(12, 5) + Q(1, 500))**2 for y in shell)
        assert Q(1, 2) * Q(314, 100)**2 * Q(6, 5)**4 > Q(51, 5)
        print(f'{n} neighbours: x=(0,0,0,5/4) protected; x=(0,0,0,3/2) cut; '
              'perturbations of 1/500 give volume > 10.2',
              flush=True)


if __name__ == '__main__':
    run()
