"""Small regression suite: python -m unittest test_interval_transfer."""

from math import isqrt
import unittest

from adjoining_interval_transfer import measure_transfer
from interval_accounting.interval_exact import process_prime_interval
from rolling_segmented_sieve import primes_up_to, rolling_counts


def prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


class RollingTests(unittest.TestCase):
    def test_root_sieve(self):
        for n in range(100):
            self.assertEqual(primes_up_to(n), [i for i in range(n + 1) if prime(i)])

    def test_against_existing_sieve(self):
        roots = [p for p in range(2, 110) if prime(p)]
        counts = [
            (p, q, process_prime_interval(p, q, verbose=False)[0])
            for p, q in zip(roots, roots[1:])
        ]
        for block in (1, 5, 89, 1 << 20):
            with self.subTest(block=block):
                actual = list(rolling_counts(2, 109, block_odds=block))
                self.assertEqual([(r.p1, r.p2, r.survivors) for r in actual], counts)
                for row in actual:
                    self.assertEqual(row.initial, row.p2**2 - row.p1**2 - 1)
                    self.assertEqual(row.removed + row.survivors, row.initial)

    def test_range_clipping(self):
        self.assertEqual(list(rolling_counts(8, 12, block_odds=7)), [])
        rows = list(rolling_counts(8, 14, block_odds=7))
        self.assertEqual([(r.p1, r.p2, r.survivors) for r in rows], [(11, 13, 9)])
        rows = list(rolling_counts(880, 888, block_odds=109))
        self.assertEqual([r.survivors for r in rows], [262, 524])


class TransferTests(unittest.TestCase):
    def test_all_small_triples_by_trial_division(self):
        roots = [p for p in range(2, 110) if prime(p)]
        for p, q, r in zip(roots, roots[1:], roots[2:]):
            with self.subTest(p=p, q=q, r=r):
                row = measure_transfer(p, q, r)
                old = [prime(n) for n in range(p*p + 1, q*q)]
                new = [prime(n) for n in range(q*q + 1, r*r)]
                prefix = list(zip(old, new))
                self.assertEqual(row.matched, len(prefix))
                self.assertEqual(row.old_total, sum(old))
                self.assertEqual(row.new_total, sum(new))
                self.assertEqual(row.kept, prefix.count((True, True)))
                self.assertEqual(row.lost, prefix.count((True, False)))
                self.assertEqual(row.gained, prefix.count((False, True)))
                self.assertEqual(row.old_tail, sum(old[row.matched:]))
                self.assertEqual(row.new_tail, sum(new[row.matched:]))

    def test_manuscript_example(self):
        row = measure_transfer(881, 883, 887)
        self.assertEqual(row.stable_primes, "2;3;7")
        self.assertEqual(
            (row.old_total, row.kept, row.lost, row.gained, row.new_tail, row.new_total),
            (262, 66, 196, 195, 263, 524),
        )

    def test_reject_nonconsecutive(self):
        for triple in ((2, 3, 7), (7, 11, 14), (7, 7, 11)):
            with self.subTest(triple=triple), self.assertRaises(ValueError):
                measure_transfer(*triple)


if __name__ == "__main__":
    unittest.main()
