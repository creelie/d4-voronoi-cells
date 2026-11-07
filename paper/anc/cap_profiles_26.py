"""Combine the recorded 26-centre certificate with single-cap upper bounds."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from itertools import accumulate

from cap_comparison import PI, single_cap


LOWER_EDGES = tuple(map(Q, ('2', '2.0161', '2.05', '2.1', '2.15', '2.2', '2.444')))
POINT = tuple(map(Q, (
    '0.10730785', '0.10730646', '-0.05347062', '-0.05363574',
    '-0.05425642', '-0.05347062', '-0.03946852',
)))
PAIR = {
    pair: Q(logged) + Q('0.0000001')
    for pair, logged in zip(('AA', 'AB', 'AF', 'BB', 'BF', 'FF'),
                            ('-.02327176', '-.01036330', '-.01503736',
                             '-.001713338', '-.004512281', '-.003681549'))
}
TRIPLE = {
    triple: Q(logged) + Q('0.0000001')
    for triple, logged in zip(
        ('AAA', 'AAB', 'AAF', 'ABB', 'ABF', 'AFF', 'BBB', 'BBF', 'BFF'),
        ('.001497925', '.0009598990', '.001322105', '.0005966748',
         '.001043879', '.001145076', '.0004094797', '.0006564131', '.001043982'),
    )
}
T = Q('8.282667610')


def compositions(total, bins):
    if bins == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, bins - 1):
            yield (first, *tail)


@lru_cache(maxsize=None)
def angular_terms(a, b, f):
    counts = {'A': a, 'B': b, 'F': f}
    assert f <= 2
    total = T / 2
    for word, coefficient in (*PAIR.items(), *TRIPLE.items()):
        multiplicity = 1
        for letter in set(word):
            multiplicity *= comb(counts[letter], word.count(letter))
        total += multiplicity * coefficient
    return total


def run():
    def round_up(value):
        grid = 10**10
        return Q(-(-value.numerator * grid // value.denominator), grid)

    caps = tuple(round_up(single_cap(Q(3, 2), edge / 2).hi) for edge in LOWER_EDGES)
    target = (Q(9, 8) * PI * PI - 8).lo
    best = (Q(0), ())
    profile_count = excluded = 0
    for n in compositions(26, 7):
        cumulative = tuple(accumulate(n))
        if n[0] > 24 or any(actual < minimum for actual, minimum in
                             zip(cumulative[1:6], (15, 20, 22, 23, 24))):
            continue
        profile_count += 1
        a, b, f = cumulative[1], cumulative[5] - cumulative[1], n[-1]
        certificate = angular_terms(a, b, f) + sum(x * m for x, m in zip(n, POINT))
        caps_only = sum(x * cap for x, cap in zip(n, caps))
        joint = min(certificate, caps_only)
        if joint < target:
            excluded += 1
        if joint > best[0]:
            best = joint, n
    assert profile_count == 17845
    assert best[0] < Q('3.288')
    print(f'PASS U < 3.288 on all {profile_count} count profiles; '
          f'largest upper bound at {best[1]}')
    print(f'{excluded} profiles clear T > 8 by these two bounds alone')
    previous_max = (15, 0, 9, 0, 0, 0, 2)
    cap_sum = sum(x * cap for x, cap in zip(previous_max, caps))
    assert cap_sum < Q('3.027')
    print(f'PASS former maximising profile {previous_max} has cap sum < 3.027')
    level_25 = Q('1.819065983') / 2 + 23 * Q('0.097085101') + 2 * Q('0.087814019')
    assert level_25 < Q('3.3182')
    print('PASS 25-neighbour level < 3.3182 if at most 23 lie within 2.0161')


if __name__ == '__main__':
    run()
