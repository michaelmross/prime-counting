"""Exact prime counts in consecutive prime-square intervals, with a rolling sieve."""

import argparse
import csv
from math import isqrt
import sys
from typing import Iterator, NamedTuple


class IntervalCount(NamedTuple):
    p1: int
    p2: int
    gap: int
    initial: int
    removed: int
    survivors: int


def primes_up_to(limit: int) -> list[int]:
    """Sieve the roots, using one byte for each odd candidate."""
    if limit < 2:
        return []
    size = (limit + 1) // 2
    flags = bytearray(b"\x01") * size
    flags[0] = 0  # Index i represents 2*i + 1.
    for p in range(3, isqrt(limit) + 1, 2):
        if flags[p // 2]:
            j = p * p // 2
            flags[j::p] = bytes(1 + (size - 1 - j) // p)
    return [2] + [2 * i + 1 for i, live in enumerate(flags) if live]


def rolling_counts(
    start: int, stop: int, *, block_odds: int = 1 << 22
) -> Iterator[IntervalCount]:
    """Yield rows for p1 >= start and p2 <= stop, with inclusive root bounds."""
    if not 2 <= start <= stop:
        raise ValueError("require 2 <= start <= stop")
    if block_odds <= 0:
        raise ValueError("block_odds must be positive")

    primes = primes_up_to(stop)
    selected = [p for p in primes if p >= start]
    if len(selected) < 2:
        return
    edges = [p * p // 2 for p in selected]
    left, right = edges[0], edges[-1]

    divisors = [p for p in primes if 2 < p < selected[-1]]
    squares = [p * p // 2 for p in divisors]
    next_hits = [s + max(0, (left - s + p - 1) // p) * p
                 for s, p in zip(squares, divisors)]
    active = cell = survivors = 0
    for lo in range(left, right, block_odds):
        hi = min(lo + block_odds, right)
        size = hi - lo
        flags = bytearray(b"\x01") * size
        while active < len(squares) and squares[active] < hi:
            active += 1
        for i in range(active):
            j = next_hits[i] - lo
            if j < size:
                p = divisors[i]
                hits = 1 + (size - 1 - j) // p
                flags[j::p] = bytes(hits)
                next_hits[i] += hits * p

        # Squares are composite, so half-open [p1², p2²) gives the open count.
        cursor = lo
        while cursor < hi:
            boundary = edges[cell + 1]
            end = min(hi, boundary)
            survivors += flags.count(1, cursor - lo, end - lo)
            cursor = end
            if end == boundary:
                p1, p2 = selected[cell], selected[cell + 1]
                initial = p2 * p2 - p1 * p1 - 1
                yield IntervalCount(p1, p2, p2 - p1, initial,
                                    initial - survivors, survivors)
                cell += 1
                survivors = 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", type=int, help="inclusive lower root bound")
    parser.add_argument("stop", type=int, help="inclusive upper root bound")
    parser.add_argument("--block-odds", type=int, default=1 << 22,
                        help="maximum odd positions per block (default: 4194304)")
    parser.add_argument("--csv", action="store_true", help="stream CSV to stdout")
    args = parser.parse_args()
    if not 2 <= args.start <= args.stop or args.block_odds <= 0:
        parser.error("require 2 <= start <= stop and a positive --block-odds")
    rows = rolling_counts(args.start, args.stop, block_odds=args.block_odds)
    if args.csv:
        writer = csv.writer(sys.stdout, lineterminator="\n")
        writer.writerow(IntervalCount._fields)
        writer.writerows(rows)
    else:
        fmt = "{:>9} {:>9} {:>5} {:>15} {:>15} {:>15}"
        print(fmt.format("p1", "p2", "Gap", "Initial", "Removed", "Survivors"))
        for row in rows:
            print(fmt.format(*row))


if __name__ == "__main__":
    main()
