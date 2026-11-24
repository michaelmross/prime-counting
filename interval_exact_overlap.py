"""Exact sieve and overlap accounting for the open interval (p1**2, p2**2).

Run: python interval_exact.py 881 883
Use --ascending to assign each composite to its least prime factor.
Use --quiet to print only the final totals.
"""

from math import isqrt


def sieve_of_eratosthenes(max_n: int) -> list[int]:
    """Return primes <= max_n in descending order."""
    if max_n < 2:
        return []
    flags = bytearray(b"\x01") * (max_n + 1)
    flags[0:2] = b"\x00\x00"
    for p in range(2, isqrt(max_n) + 1):
        if flags[p]:
            for multiple in range(p * p, max_n + 1, p):
                flags[multiple] = 0
    return [p for p in range(max_n, 1, -1) if flags[p]]


def process_prime_interval(
    p1: int = 881,
    p2: int = 883,
    *,
    descending: bool = True,
    verbose: bool = True,
):
    """Return (exact prime count, per-prime removal and overlap history).

    total_hits counts all composite multiples visited at the step.
    overlap counts the hits already removed by an earlier sieve prime.
    removed counts the newly removed positions.
    """
    if not isinstance(p1, int) or not isinstance(p2, int):
        raise TypeError("The two square roots must be integers.")
    if not 1 <= p1 < p2:
        raise ValueError("Require 1 <= p1 < p2.")

    lower, upper = p1 * p1 + 1, p2 * p2 - 1
    total = upper - lower + 1
    alive = bytearray(b"\x01") * total
    remaining = total
    history = []
    primes = sieve_of_eratosthenes(p2)
    if not descending:
        primes.reverse()

    if verbose:
        print(f"Integers: {lower} through {upper}; initial count: {total}")
        print("Prime |     Before |   All hits |    Overlap | Newly removed |  Remaining")

    for p in primes:
        before = remaining
        # Include all in-range composite multiples, but never erase p itself.
        first = max(2 * p, ((lower + p - 1) // p) * p)
        multiples = range(first, upper + 1, p)
        total_hits = len(multiples)
        removed = 0
        for n in multiples:
            index = n - lower
            if alive[index]:
                alive[index] = 0
                removed += 1
        overlap = total_hits - removed
        remaining -= removed
        history.append({
            "prime": p,
            "prev_interval": before,
            "total_hits": total_hits,
            "overlap": overlap,
            "removed": removed,
            "new_interval": remaining,
        })
        if verbose:
            print(f"{p:5d} | {before:10d} | {total_hits:10d} | {overlap:10d} | "
                  f"{removed:13d} | {remaining:10d}")

    return remaining, history


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p1", type=int, nargs="?", default=881)
    parser.add_argument("p2", type=int, nargs="?", default=883)
    parser.add_argument("--ascending", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    count, history = process_prime_interval(
        args.p1, args.p2,
        descending=not args.ascending,
        verbose=not args.quiet,
    )
    initial = args.p2**2 - args.p1**2 - 1
    print(f"Initial integers: {initial}")
    print(f"Total sieve hits: {sum(row['total_hits'] for row in history)}")
    print(f"Overlap hits (repeat hits): {sum(row['overlap'] for row in history)}")
    print(f"Composites removed: {initial - count}")
    print(f"Prime survivors: {count}")
