"""Measure exact offset-by-offset transfer across adjoining prime-square intervals."""

import argparse
from bisect import bisect_left
import csv
import sys
from typing import NamedTuple, Sequence

from rolling_segmented_sieve import primes_up_to


class TransferRow(NamedTuple):
    p: int
    q: int
    r: int
    shift: int
    old_width: int
    new_width: int
    matched: int
    stable_primes: str
    old_prefix: int
    kept: int
    lost: int
    gained: int
    new_prefix: int
    net_prefix: int
    old_tail: int
    new_tail: int
    old_total: int
    new_total: int


def first_offset(start: int, prime: int) -> int:
    # Mark divisibility by this prime, except for the prime itself.
    return max(2 * prime, ((start + prime - 1) // prime) * prime) - start


def tail_count(start: int, width: int, divisors: Sequence[int]) -> int:
    if not width:
        return 0
    alive = bytearray(b"\x01") * width
    for prime in divisors:
        first = first_offset(start, prime)
        if first < width:
            alive[first::prime] = bytes(1 + (width - 1 - first) // prime)
    return alive.count(1)


def _measure(p: int, q: int, r: int, divisors: Sequence[int]) -> TransferRow:
    a, b, c = p * p, q * q, r * r
    old_width, new_width = b - a - 1, c - b - 1
    matched = min(old_width, new_width)
    shift = b - a

    # Mark multiplicities are needed: removing one mark may leave another.
    # With a bytearray-sized interval, the number of distinct prime factors
    # of a position here is far below 255.
    marks = bytearray(matched)
    moves = []
    stable = []
    for prime in divisors:
        old = first_offset(a + 1, prime)
        new = first_offset(b + 1, prime)
        for j in range(old, matched, prime):
            marks[j] += 1
        if shift % prime == 0 and old == new:
            stable.append(prime)
        elif old != new and min(old, new) < matched:
            moves.append((prime, old, new))

    previous = marks.copy()
    for prime, old, new in moves:
        for j in range(old, matched, prime):
            marks[j] -= 1
        for j in range(new, matched, prime):
            marks[j] += 1
    old_prefix = previous.count(0)
    new_prefix = marks.count(0)
    lost = sum(before == 0 and after != 0 for before, after in zip(previous, marks))
    kept = old_prefix - lost
    gained = new_prefix - kept

    old_tail = tail_count(a + 1 + matched, old_width - matched, divisors)
    new_tail = tail_count(b + 1 + matched, new_width - matched, divisors)
    return TransferRow(p, q, r, shift, old_width, new_width, matched,
                       ";".join(map(str, stable)), old_prefix, kept, lost, gained,
                       new_prefix, gained - lost, old_tail, new_tail,
                       old_prefix + old_tail, new_prefix + new_tail)


def measure_transfer(p: int, q: int, r: int) -> TransferRow:
    """Compare (p², q²) and (q², r²); reject nonconsecutive prime triples."""
    if not 2 <= p < q < r:
        raise ValueError("require p < q < r, with p >= 2")
    primes = primes_up_to(r)
    i = bisect_left(primes, p)
    if primes[i:i + 3] != [p, q, r]:
        raise ValueError("p, q, r must be consecutive primes")
    return _measure(p, q, r, primes[:i + 2])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", metavar="ROOT", type=int, nargs="*",
                        help="three consecutive primes p q r")
    parser.add_argument("--range", type=int, nargs=2, metavar=("START", "STOP"),
                        help="all consecutive triples within inclusive root bounds")
    parser.add_argument("--csv", action="store_true", help="stream CSV to stdout")
    args = parser.parse_args()

    if args.range is not None:
        start, stop = args.range
        if args.roots or not 2 <= start <= stop:
            parser.error("give either three roots or --range START STOP (2 <= START <= STOP)")
        primes = primes_up_to(stop)
        i = bisect_left(primes, start)
        rows = (_measure(*primes[j:j + 3], primes[:j + 2])
                for j in range(i, len(primes) - 2))
    else:
        if len(args.roots) != 3:
            parser.error("give three consecutive roots or --range START STOP")
        try:
            rows = iter((measure_transfer(*args.roots),))
        except ValueError as exc:
            parser.error(str(exc))

    if args.csv:
        writer = csv.writer(sys.stdout, lineterminator="\n")
        writer.writerow(TransferRow._fields)
        writer.writerows(rows)
    else:
        for row in rows:
            print(f"({row.p}^2, {row.q}^2) -> ({row.q}^2, {row.r}^2)")
            print(f"  Widths: {row.old_width:,} -> {row.new_width:,}; "
                  f"matched: {row.matched:,}; stable primes: {row.stable_primes or 'none'}")
            print(f"  Matched: {row.old_prefix:,} - {row.lost:,} lost "
                  f"+ {row.gained:,} gained = {row.new_prefix:,} "
                  f"({row.kept:,} kept; net {row.net_prefix:+,})")
            print(f"  Tail primes: {row.old_tail:,} -> {row.new_tail:,}; "
                  f"total primes: {row.old_total:,} -> {row.new_total:,}")


if __name__ == "__main__":
    main()
