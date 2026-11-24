# Prime Populations Between Prime Squares: Companion Code

This snapshot accompanies *An Algorithmic Account of Prime Populations Between Prime Squares* (Ross, M. M., 2026). It contains the algorithms and selected verification inputs, not the complete log of the computational sweep.

| Path | What it does |
| --- | --- |
| `interval_accounting/interval_exact.py` | Counts primes in the open interval and records removals by sieve prime. `--ascending` gives the least-prime-factor bins. |
| `interval_accounting/run_consecutive_intervals.py` | Runs the exact sieve on consecutive prime pairs in a specified range of roots. |
| `interval_overlap/interval_exact_overlap.py` | Adds all divisor hits and repeat hits to the exact sieve ledger. |
| `interval_overlap/overlap_bound.py` | Computes a floor-count lower bound with a presieved core. An inconclusive lower bound does not refute the four-prime claim. |
| `interval_overlap/run_interval_overlaps.py` | Reports exact and certified lower counts for a range of consecutive prime roots. |
| `interval_overlap/core_diagnostic.py` | Diagnoses how the overlap certificate depends on the core cutoff; requires NumPy. |
| `research_cycle_1032/verify.cjs` | Reproduces the middle-semiprime strip checks using archived JSON inputs from `research_cycle_1029` and `research_cycle_1031`. |
| `verify_final_interval.py` | Checks the final pair of roots and independently counts primes between their squares. |
| `rolling_segmented_sieve.py` | Streams exact prime counts at square boundaries while carrying each sieve prime's next multiple from block to block. |
| `adjoining_interval_transfer.py` | Updates the common sieve at aligned offsets and reports survivor gains, losses, and unmatched tails. |
| `test_interval_transfer.py` | Regression tests for the rolling and adjoining-transfer scripts. |

## Quick checks

Run from this directory with Python 3.10+ and Node.js 18+:

```text
python interval_overlap/interval_exact_overlap.py 881 883 --quiet
python interval_overlap/overlap_bound.py 881 883 --core 7
node research_cycle_1032/verify.cjs
```

The first command gives 3,527 positions, 7,696 total divisor hits, 4,431 repeated hits, 3,265 composites, and 262 primes. The second certifies a lower bound of 254 primes.

To test a small range:

```text
python interval_accounting/run_consecutive_intervals.py 800 900
```

## Adjoining-interval scripts

Both scripts use only the Python standard library. Keep them in the same folder:

```text
python rolling_segmented_sieve.py 881 887
python adjoining_interval_transfer.py 881 883 887
python adjoining_interval_transfer.py --range 9767 9859 --csv
python -m unittest test_interval_transfer
```

The rolling script takes inclusive root bounds and prints one row per complete consecutive-prime pair in range. Add `--csv` for machine-readable output; `--block-odds` controls the number of odd positions in each sieve block, independently of the output interval widths. The default is 4,194,304 odd positions. This is the script for large sweeps.

The transfer script accepts three consecutive prime roots or `--range START STOP` to compare each adjoining pair of cells in range. It matches prefixes of equal length and prints both unmatched tails, so it handles shrinking as well as growing intervals. On a matched prefix, `new_prefix = old_prefix - lost + gained`; each full count adds its own tail. For `881 883 887`, the first prefix has 262 primes: 66 offsets are kept, 196 are lost, 195 are gained, and the destination gets 263 more primes in its tail, for a total of 524. This is an exact diagnostic, not a certificate derived from the previous prime count alone. It keeps mark counts for each matched offset and takes more space on wide intervals.

## Reported sweep endpoint

The study counted primes in the open intervals between consecutive prime squares through the pair $p=1{,}224{,}709$, $q=1{,}224{,}739$. The last interval is $(1{,}499{,}912{,}134{,}681,\;1{,}499{,}985{,}618{,}121)$ and was independently checked to contain 2,618,304 primes. Run `python verify_final_interval.py` to repeat that check. This package does not include the full sweep output or assert that it has been independently rerun. Running the range driver over the full range is a substantial computation.

## Citation and license

Code: MIT. Note and data: CC-BY 4.0.
