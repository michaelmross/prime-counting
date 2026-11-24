import argparse
from math import isqrt

from interval_exact_overlap import process_prime_interval
from overlap_bound import CORE_LIMITS, overlap_bound


def get_consecutive_prime_pairs(max_n: int) -> list[tuple[int, int]]:
    """Return consecutive prime pairs (p1, p2) with p2 <= max_n."""
    if max_n < 3:
        return []

    is_prime = [True] * (max_n + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, isqrt(max_n) + 1):
        if is_prime[p]:
            for multiple in range(p * p, max_n + 1, p):
                is_prime[multiple] = False

    primes = [p for p in range(2, max_n + 1) if is_prime[p]]
    return list(zip(primes, primes[1:]))


def run_intervals(max_n: int, ascending: bool = False, core_limit: int = 7, *, start: int = 2):
    """Print and return exact sieve totals and overlap certificates."""
    pairs = [(p1, p2) for p1, p2 in get_consecutive_prime_pairs(max_n) if p1 >= start]
    if not pairs:
        print(f"No consecutive prime pairs found in the root range [{start}, {max_n}].")
        return

    results = []
    header = (
        f"{'p1':>6} | {'p2':>6} | {'Gap':>5} | {'Initial':>10} | "
        f"{'Removed':>10} | {'Survivors':>10} | {'All hits':>10} | {'Overlap':>10} | "
        f"{'Core repeats':>12} | {'Tail pairs':>10} | {'Tail triples':>12} | "
        f"{'Overlap LB':>10} | {'Overlap for 4':>13} | {'Prime LB':>10} | Status"
    )
    print(header)
    print("-" * len(header))

    for p1, p2 in pairs:
        count, history = process_prime_interval(
            p1, p2, descending=not ascending, verbose=False
        )
        initial = p2**2 - p1**2 - 1
        removed = initial - count
        gap = p2 - p1
        total_hits = sum(step["total_hits"] for step in history)
        overlap = sum(step["overlap"] for step in history)
        bound = overlap_bound(p1, p2, core_limit=core_limit)
        certified = bound["survivor_lower_bound"] >= 4
        status = "CERTIFIED" if certified else "INCONCLUSIVE"

        results.append(
            {
                **bound,
                "p1": p1,
                "p2": p2,
                "gap": gap,
                "initial": initial,
                "removed": removed,
                "survivors": count,
                "total_hits": total_hits,
                "overlap": overlap,
                "certified_four": certified,
                "status": status,
            }
        )
        print(
            f"{p1:6d} | {p2:6d} | {gap:5d} | {initial:10d} | {removed:10d} | "
            f"{count:10d} | {total_hits:10d} | {overlap:10d} | "
            f"{bound['core_repeats']:12d} | {bound['tail_pairs']:10d} | "
            f"{bound['tail_triples']:12d} | {bound['overlap_lower_bound']:10d} | "
            f"{bound['four_prime_threshold']:13d} | {bound['survivor_lower_bound']:10d} | {status}"
        )

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Report exact sieve totals and four-prime certificates for p2 <= n."
    )
    parser.add_argument("n", type=int, help="Upper bound for p2 (e.g. 103)")
    parser.add_argument("--start", type=int,
                        help="Inclusive lower bound for p1 (default: 2)")
    parser.add_argument(
        "--ascending", action="store_true", help="Use ascending prime removal order"
    )
    parser.add_argument("--core", type=int, choices=CORE_LIMITS, default=7,
                        help="Presieve through this prime (default: 7)")
    args = parser.parse_args()
    if args.start is not None and not 2 <= args.start <= args.n:
        parser.error("--start must be between 2 and n")
    run_intervals(args.n, ascending=args.ascending, core_limit=args.core,
                  start=2 if args.start is None else args.start)
