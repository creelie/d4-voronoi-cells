"""Rational interval checks for the pair and triple cap comparisons."""

from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __add__(self, other):
        other = enclose(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-enclose(other))

    def __rsub__(self, other):
        return enclose(other) + (-self)

    def __mul__(self, other):
        other = enclose(other)
        products = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = enclose(other)
        assert other.lo > 0
        return self * Interval(1 / other.hi, 1 / other.lo)


def enclose(number):
    return number if isinstance(number, Interval) else Interval(Q(number), Q(number))


def sqrt_enclosure(value):
    assert value >= 0
    grid = 10**24
    integer = isqrt(value.numerator * grid**2 // value.denominator)
    return Interval(Q(integer, grid), Q(integer + 1, grid))


def atan_endpoint(value):
    assert 0 <= value < 1
    lower = sum((-1) ** k * value ** (2 * k + 1) / (2 * k + 1) for k in range(48))
    return Interval(lower, lower + value**97 / 97)


def atan_enclosure(value):
    value = enclose(value)
    return Interval(atan_endpoint(value.lo).lo, atan_endpoint(value.hi).hi)


PI = 16 * atan_enclosure(Q(1, 5)) - 4 * atan_enclosure(Q(1, 239))


def single_cap(radius2, h):
    radius4 = radius2**2
    perpendicular = sqrt_enclosure(radius2 - h**2)
    angle = atan_enclosure(perpendicular / h)
    return (4 * PI / 3) * (Q(3, 8) * radius4 * angle
                          - h * (5 * radius2 - 2 * h**2) * perpendicular / 8)


def cap_data(radius2, h):
    radius4 = radius2**2
    tangent = sqrt_enclosure(radius2 - h**2) / h
    angle = atan_enclosure(tangent)
    # Equal offsets and normals at angle pi/3. Integrate one half-lens
    # from pi/6 to arctan(tangent), then double by reflection.
    def primitive(arg, theta):
        return (radius4 * theta - 2 * radius2 * h**2 * arg
                + h**4 * (arg + arg * arg * arg / 3))

    pair = PI / 2 * (primitive(tangent, angle)
                     - primitive(sqrt_enclosure(Q(1, 3)), PI / 6))
    return PI * PI * radius4 / 2, single_cap(radius2, h), pair


def decimal_interval(interval, digits=8):
    scale = 10**digits
    lower = interval.lo * scale
    upper = interval.hi * scale
    low = lower.numerator // lower.denominator
    high = -(-upper.numerator // upper.denominator)
    return f'[{low // scale}.{low % scale:0{digits}d}, {high // scale}.{high % scale:0{digits}d}]'


def run():
    dilation = Q(53, 40)
    assert (dilation**2 + 8 * dilation) / 9 < Q(11, 8)
    assert Q(32 * 27 * 6**3, 105 * 4**7) / Q(314, 100) < Q(7, 200)
    assert 5 * Q(7, 200) / 3 == Q(7, 120)
    assert Q(29, 9) < Q(15, 4)
    print('PASS pair dilation, triple-to-pair ratio and five-completion constants')

    for h in (Q(1), Q(10031, 10000)):
        b0, s0, p0 = cap_data(Q(3, 2), h)
        b1, s1, p1 = cap_data(Q(25, 16), h)
        truncated = b0 - 24 * s0 + 96 * p0
        common = b1 - b0 - 24 * (s1 - s0)
        beta = dilation**3 - 1 - Q(7, 120)
        coarse = truncated + Q(77, 100) * (common + 96 * beta * p0)
        sharp = truncated + Q(77, 100) * (common + 96 * (p1 - Q(127, 120) * p0))
        if h == 1:
            assert sharp.hi < Q(79361, 10000)
        else:
            assert truncated.hi < Q(79963, 10000)
            assert coarse.lo > Q(80008, 10000)
            assert sharp.lo > Q(80291, 10000)
        print(f'h={h} T={decimal_interval(truncated)} '
              f'one-radius lower expression={decimal_interval(coarse)} '
              f'two-radius lower expression={decimal_interval(sharp)}')


if __name__ == '__main__':
    run()
