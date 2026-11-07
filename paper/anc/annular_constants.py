"""Exact rational comparisons in the annular-volume appendix."""

from fractions import Fraction as F


def atan_bounds(x):
    lower = sum((-1) ** j * x ** (2 * j + 1) / (2 * j + 1) for j in range(8))
    return lower, lower + x**17 / 17


def run():
    a, a_upper = F(12247, 10000), F(122475, 100000)
    q_lower, q_upper = F(12649, 10000), F(1265, 1000)
    slope_lower, slope_upper = F(77, 100), F(31, 40)
    assert a**2 < F(3, 2) < a_upper**2
    assert q_lower**2 < F(8, 5) < q_upper**2
    assert slope_lower**2 < F(3, 5) < slope_upper**2

    rho, b = F(5, 4), F(1, 4)
    loss_square = (F(2, 5) * F(96, 95)**2)**2 * (2 * rho * (rho - a))**3 / (a**6 * b**6)
    assert loss_square == F(22535460391546257408000000000, 439736697911706067259903971009)
    assert loss_square < F(23, 100)**2
    print('PASS exterior removal fraction < 23/100')

    assert (q_upper - a) / (a * slope_lower) < F(1, 20)
    tail_lower = 2 * F(314, 100) * (q_lower - a_upper)**4 / slope_upper**3
    assert tail_lower == F(65277247497517, 1861937500000000000)
    assert tail_lower > F(1, 28600)
    print('PASS vertex-cone slope < 1/20 and volume correction > 1/28600')

    l5, u5 = atan_bounds(F(1, 5))
    l239, u239 = atan_bounds(F(1, 239))
    pi_lower, pi_upper = F(31415926, 10000000), F(31415927, 10000000)
    assert pi_lower < 16 * l5 - 4 * u239
    assert 16 * u5 - 4 * l239 < pi_upper
    denominator = 9 * pi_lower**2 / 8 - F(33274, 10000) + F(1, 28600)
    assert 0 < denominator < 8
    assert pi_lower**2 / (2 * denominator) < F(634625, 1000000)
    print('PASS 3.1415926 < pi < 3.1415927 and density bound < 0.634625')


if __name__ == '__main__':
    run()
