import argparse
import logging
from interval_exact import process_prime_interval


def get_consecutive_prime_pairs(start: int, end: int) -> list[tuple[int, int]]:
    """Return consecutive prime pairs (p1, p2) where p1 >= start and p2 <= end."""
    if end < 3:
        return []

    # Sieve of Eratosthenes up to `end`
    is_prime = [True] * (end + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(end**0.5) + 1):
        if is_prime[p]:
            for multiple in range(p * p, end + 1, p):
                is_prime[multiple] = False

    primes = [p for p in range(2, end + 1) if is_prime[p]]

    # Filter pairs matching start and end bounds
    pairs = []
    for i in range(len(primes) - 1):
        p1, p2 = primes[i], primes[i + 1]
        if p1 >= start and p2 <= end:
            pairs.append((p1, p2))

    return pairs


def setup_logging(log_file: str | None = None):
    """Configure logging output to console and an optional log file."""
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w", encoding="utf-8"))

    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=handlers,
    )


def run_intervals(start: int, end: int, ascending: bool = False):
    """Run process_prime_interval on all consecutive prime pairs in range [start, end]."""
    pairs = get_consecutive_prime_pairs(start, end)

    if not pairs:
        logging.info(f"No consecutive prime pairs found in range [{start}, {end}].")
        return

    header = f"{'p1':>6} | {'p2':>6} | {'Gap':>5} | {'Initial':>10} | {'Removed':>10} | {'Survivors':>10}"
    divider = "-" * 62

    logging.info(f"Processing prime pairs with p1 >= {start} and p2 <= {end}\n")
    logging.info(header)
    logging.info(divider)

    results = []
    for p1, p2 in pairs:
        count, _ = process_prime_interval(
            p1, p2, descending=not ascending, verbose=False
        )
        initial = p2**2 - p1**2 - 1
        removed = initial - count
        gap = p2 - p1

        results.append({
            "p1": p1,
            "p2": p2,
            "gap": gap,
            "initial": initial,
            "removed": removed,
            "survivors": count,
        })

        logging.info(
            f"{p1:6d} | {p2:6d} | {gap:5d} | {initial:10d} | {removed:10d} | {count:10d}"
        )

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Feed consecutive prime pairs in [start, end] into interval_exact in quiet mode."
    )
    parser.add_argument("start", type=int, help="Lower bound for p1")
    parser.add_argument("end", type=int, help="Upper bound for p2")
    parser.add_argument(
        "-l", "--log", type=str, default=None, help="Log file path (e.g., results.log)"
    )
    parser.add_argument(
        "--ascending",
        action="store_true",
        help="Use ascending prime removal order",
    )
    args = parser.parse_args()

    setup_logging(args.log)
    run_intervals(args.start, args.end, ascending=args.ascending)