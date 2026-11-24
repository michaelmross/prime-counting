"""Measure the smallest core that certifies four primes (requires NumPy).

Example: python core_diagnostic.py 375019 --start 375017
This scans interval positions; it is a diagnostic, not a floor-only proof.
"""

import argparse
from bisect import bisect_left
import csv
from math import comb, isqrt, log
import sys
from time import perf_counter

import numpy as np


def prime_list(n: int) -> list[int]:
    flags = np.ones(n + 1, dtype=bool)
    flags[:2] = False
    for p in range(2, isqrt(n) + 1):
        if flags[p]:
            flags[p*p::p] = False
    return np.flatnonzero(flags).tolist()


def diagnose(p: int, q: int, primes: list[int], block_size: int, reference: int):
    started = perf_counter()
    lo, hi = p*p + 1, q*q - 1
    ceiling = isqrt(isqrt(hi))
    losses = [0] * (ceiling + 1)
    labels = [r for r in primes if r <= q]
    survivors = 0

    for a in range(lo, hi + 1, block_size):
        length = min(block_size, hi - a + 1)
        marks = np.zeros(length, dtype=np.uint8)
        least = np.zeros(length, dtype=np.uint32)
        for r in labels:
            first = max((-a) % r, 2*r-a)
            if first >= length:
                continue
            marks[first::r] += 1
            # Four distinct factors imply that the least is <= hi**(1/4).
            if r <= ceiling:
                view = least[first::r]
                view[view == 0] = r
        survivors += length - int(np.count_nonzero(marks))
        heavy = marks >= 4
        if not np.any(heavy):
            continue
        heavy_marks = marks[heavy]
        heavy_least = least[heavy]
        for k in range(4, int(heavy_marks.max()) + 1):
            counts = np.bincount(heavy_least[heavy_marks == k])
            credit = comb(k-1, 3)
            for factor in np.flatnonzero(counts):
                losses[int(factor)] += int(counts[factor]) * credit

    # For any fixed core, eligible rows have least factor above its cutoff.
    reference_loss = sum(losses[reference+1:])
    reference_raw = survivors - reference_loss
    cutoff = None
    at_cutoff = None
    prior_core = None
    prior_raw = None
    if survivors >= 4:
        loss = sum(losses[3:])
        cutoff = 2
        while survivors-loss < 4:
            cutoff += 1
            loss -= losses[cutoff]
        at_cutoff = survivors-loss
        if cutoff > 2:
            prior_core = primes[bisect_left(primes, cutoff) - 1]
            prior_raw = survivors - sum(losses[prior_core+1:])

    return (
        p, q, hi-lo+1, survivors, reference, reference_loss, reference_raw,
        max(0, reference_raw), cutoff,
        None if cutoff is None else f"{log(cutoff)/log(p):.8f}",
        at_cutoff, "CERTIFIED" if cutoff is not None else "FEWER_THAN_FOUR",
        f"{perf_counter()-started:.3f}",
        prior_core, prior_raw,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, help="Inclusive upper root bound")
    parser.add_argument("--start", type=int, default=3, help="Inclusive lower root bound")
    parser.add_argument("--step", type=int, default=1, help="Test every Nth eligible pair")
    parser.add_argument("--all-pairs", action="store_true", help="Include non-twin pairs")
    parser.add_argument("--reference-core", type=int, default=7, help="Also report the bound at this core")
    parser.add_argument("--block-size", type=int, default=8_000_000,
                        help="Positions per block; smaller values reduce memory")
    args = parser.parse_args()
    if not 2 <= args.start <= args.n:
        parser.error("require 2 <= --start <= n")
    if args.step < 1 or args.block_size < 1 or args.reference_core < 2:
        parser.error("require positive --step and --block-size, and --reference-core >= 2")

    primes = prime_list(args.n)
    pairs = [(p,q) for p,q in zip(primes, primes[1:])
             if p >= args.start and (args.all_pairs or q-p == 2)]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow((
        "p1", "p2", "Initial", "Survivors", "Reference Core", "Reference Loss",
        "Reference Raw LB", "Reference Prime LB", "Minimum Core", "Core Exponent",
        "Prime LB at Minimum", "Status", "Seconds",
        "Core Below Minimum", "Raw LB Below Minimum",
    ))
    sys.stdout.flush()
    for p,q in pairs[::args.step]:
        writer.writerow(diagnose(p,q,primes,args.block_size,args.reference_core))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
