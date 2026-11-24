"""Certify a lower bound on primes between two squares using floor counts.

Example: python overlap_bound.py 101 103

The roots need not be prime. A lower bound below four is inconclusive.
"""

import argparse
from functools import lru_cache
from math import isqrt

CORE_LIMITS = (2, 3, 5, 7, 11, 13)


def primes_up_to(limit: int) -> list[int]:
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start::p] = b"\x00" * ((limit - start) // p + 1)
    return [p for p in range(2, limit + 1) if flags[p]]


@lru_cache(maxsize=None)
def wheel_prefix(core: tuple[int, ...]) -> tuple[int, ...]:
    """Count integers coprime to the core within one repeating period."""
    period = 1
    for p in core:
        period *= p
    prefix = [0]
    for n in range(1, period + 1):
        prefix.append(prefix[-1] + all(n % p for p in core))
    return tuple(prefix)


def overlap_bound(p1: int, p2: int, core_limit: int = 7) -> dict[str, int]:
    """Count hits and certify a bound without marking interval positions."""
    if not 1 <= p1 < p2:
        raise ValueError("Require 1 <= p1 < p2.")
    if core_limit not in CORE_LIMITS:
        raise ValueError(f"core_limit must be one of {CORE_LIMITS}.")
    a, b = p1**2, p2**2 - 1
    h = b - a

    def multiples(d: int) -> int:
        return b // d - a // d

    # Match the exact sieve's list, including the upper root if it is prime.
    primes = primes_up_to(p2)
    core = tuple(p for p in primes if p <= core_limit)
    tail = [p for p in primes if p > core_limit]
    prefix = wheel_prefix(core)
    period = len(prefix) - 1

    def coprime_through(n: int) -> int:
        periods, remainder = divmod(n, period)
        return periods * prefix[-1] + prefix[remainder]

    def clean_multiples(d: int) -> int:
        # d is coprime to the core, so only its cofactor needs screening.
        return coprime_through(b // d) - coprime_through(a // d)

    all_hits = sum(multiples(p) - (a < p <= b) for p in primes)
    core_pool = clean_multiples(1) + sum(a < p <= b for p in core)
    tail_singles = sum(clean_multiples(p) - (a < p <= b) for p in tail)
    # Each core-covered position accounts for one first hit. All its other
    # hits (including hits by tail primes) are repeats counted exactly here.
    core_repeats = all_hits - tail_singles - (h - core_pool)

    tail_pairs = tail_triples = 0
    size = len(tail)
    for i, r in enumerate(tail):
        for j in range(i + 1, size):
            rs = r * tail[j]
            if rs > b:
                break
            hits = clean_multiples(rs)
            if not hits:
                continue  # An extending triple cannot have a clean hit either.
            tail_pairs += hits
            for k in range(j + 1, size):
                rst = rs * tail[k]
                if rst > b:
                    break
                tail_triples += clean_multiples(rst)

    # On a still-eligible position with k marks, pairs minus triples is at
    # most max(0, k-1). The repeats on the other positions are exact.
    lower_overlap = max(0, core_repeats + tail_pairs - tail_triples)
    return {
        "initial": h,
        "all_hits": all_hits,
        "core_repeats": core_repeats,
        "tail_pairs": tail_pairs,
        "tail_triples": tail_triples,
        "overlap_lower_bound": lower_overlap,
        "four_prime_threshold": all_hits - h + 4,
        "survivor_lower_bound": max(0, h - all_hits + lower_overlap),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p1", type=int, help="Lower square root")
    parser.add_argument("p2", type=int, help="Upper square root")
    parser.add_argument("--core", type=int, choices=CORE_LIMITS, default=7,
                        help="Presieve through this prime (default: 7; 2 reproduces the old bound)")
    args = parser.parse_args()
    try:
        result = overlap_bound(args.p1, args.p2, core_limit=args.core)
    except ValueError as exc:
        parser.error(str(exc))
    for label, key in (
        ("Initial positions", "initial"),
        ("All hits", "all_hits"),
        ("Repeat hits on positions removed by the presieve", "core_repeats"),
        ("Pair hits after presieving", "tail_pairs"),
        ("Triple hits after presieving", "tail_triples"),
        ("Minimum certified overlap", "overlap_lower_bound"),
        ("Overlap needed for four primes", "four_prime_threshold"),
        ("Minimum certified prime count", "survivor_lower_bound"),
    ):
        print(f"{label}: {result[key]:,}")
    if result["survivor_lower_bound"] >= 4:
        print("At least four primes certified.")
    else:
        print("The lower bound is inconclusive for four primes.")
