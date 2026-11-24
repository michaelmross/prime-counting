"""Exact step-by-step accounting for the open interval (p1**2, p2**2).

Run: python interval_exact.py 881 883
Use --ascending to count composites by their least prime factor.
Use --quiet to print only the final totals.

The interval flags use one byte per integer. All arithmetic is integral.
The supplied square roots need not themselves be prime for this sieve to work.
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
    """Return (exact prime count, per-prime removal history).

Each history entry counts only multiples still present at that step.
Descending order preserves the original program's order, including p2
when p2 is prime. Ascending order gives least-prime-factor frequencies.
"""
    if not isinstance(p1, int) or not isinstance(p2, int):
        raise TypeError("The two square roots must be integers.")
    if not 2 <= p1 < p2:
        raise ValueError("Require 2 <= p1 < p2.")

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
        print("Prime | Before | Newly removed | Remaining")

    for p in primes:
        before = remaining
        # Start at the first multiple in the interval, but never erase p
        # itself. Unlike p*p, 2*p also preserves descending-step accounting.
        first = max(2 * p, ((lower + p - 1) // p) * p)
        removed = 0
        for n in range(first, upper + 1, p):
            index = n - lower
            if alive[index]:
                alive[index] = 0
                removed += 1
        remaining -= removed
        history.append({
            "prime": p,
            "prev_interval": before,
            "removed": removed,
            "new_interval": remaining,
        })
        if verbose:
            print(f"{p:5d} | {before:6d} | {removed:13d} | {remaining:9d}")

    # Every composite <= upper has a prime divisor <= isqrt(upper) < p2.
    # Every prime is retained, so the remaining flags count exactly primes.
    return remaining, history


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p1", type=int, nargs="?", default=881)
    parser.add_argument("p2", type=int, nargs="?", default=883)
    parser.add_argument("--ascending", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    count, _ = process_prime_interval(
        args.p1, args.p2,
        descending=not args.ascending,
        verbose=not args.quiet,
    )
    initial = args.p2**2 - args.p1**2 - 1
    print(f"Initial integers: {initial}")
    print(f"Composites removed: {initial - count}")
    print(f"Prime survivors: {count}")
