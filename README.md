# Prime Populations Between Prime Squares: Companion Code

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24488933-blue.svg)](https://doi.org/10.5281/zenodo.24488933)

This snapshot accompanies [*An Algorithmic Account of Prime Populations Between Prime Squares*](https://doi.org/10.5281/zenodo.24488671) (Ross, M. M., 2026). It contains the algorithms, selected verification inputs, and dataset.

| Path | What it does |
| --- | --- |
| `interval_exact.py` | Counts primes in the open interval and records removals by sieve prime. `--ascending` gives the least-prime-factor bins. |
| `run_consecutive_intervals.py` | Runs the exact sieve on consecutive prime pairs in a specified range of roots. |
| `interval_exact_overlap.py` | Adds all divisor hits and repeat hits to the exact sieve ledger. |
| `overlap_bound.py` | Computes a floor-count lower bound with a presieved core. An inconclusive lower bound does not refute the four-prime claim. |
| `run_interval_overlaps.py` | Reports exact and certified lower counts for a range of consecutive prime roots. |
| `core_diagnostic.py` | Diagnoses how the overlap certificate depends on the core cutoff; requires NumPy. |
| `verify.cjs` | Reproduces the middle-semiprime strip checks using archived JSON inputs. |
| `verify_final_interval.py` | Independently checks the earlier $1.5\times10^{12}$ sweep endpoint. |
| `rolling_segmented_sieve.py` | Streams exact prime counts at square boundaries while carrying each sieve prime's next multiple from block to block. |
| `adjoining_interval_transfer.py` | Updates the common sieve at aligned offsets and reports survivor gains, losses, and unmatched tails. |
| `test_interval_transfer.py` | Regression tests for the rolling and adjoining-transfer scripts. |
| `data/4x10e12_sweep.csv` | Sweep data through upper root $2{,}000{,}029$.|

## Quick checks

Run from this directory with Python 3.10 or later:

```text
python interval_exact.py 881 883 --quiet
python interval_exact_overlap.py 881 883 --quiet
python overlap_bound.py 881 883 --core 7
```

In the first command, the **Newly removed** column is the per-prime removal histogram: each composite is counted under the first sieve prime that marks it. The default order is descending; `--ascending` shows least-prime-factor bins. The overlap ledger also counts later hits: it gives 3,527 positions, 7,696 total divisor hits, 4,431 repeated hits, 3,265 composites, and 262 primes. The floor-count command certifies a lower bound of 254 primes without finding their positions.

The optional `node verify.cjs` check (Node.js 18 or later) requires two archived JSON inputs. The legacy verifier imports `../research_cycle_1029/results.json` and `../research_cycle_1031/results.json`. Without those files at the expected paths, the command does not run from the flat checkout alone.

To test a small range:

```text
python run_consecutive_intervals.py 800 900
```

## Adjoining-interval scripts

Both scripts use only the Python standard library. Keep them in the same folder:

```text
python rolling_segmented_sieve.py 881 887
python adjoining_interval_transfer.py 881 883 887
python adjoining_interval_transfer.py --range 9767 9859 --csv
```

For `python -m unittest test_interval_transfer`, the flat-layout test must import `process_prime_interval` from `interval_exact` rather than `interval_accounting.interval_exact`.

The rolling script takes inclusive root bounds and prints one row per complete consecutive-prime pair in range. Add `--csv` for machine-readable output; `--block-odds` controls the number of odd positions in each sieve block, independently of the output interval widths. The default is 4,194,304 odd positions. This is the script for large sweeps.

The transfer script accepts three consecutive prime roots or `--range START STOP` to compare each adjoining pair of cells in range. It matches prefixes of equal length and prints both unmatched tails, so it handles shrinking as well as growing intervals. On a matched prefix, `new_prefix = old_prefix - lost + gained`; each full count adds its own tail. For `881 883 887`, the first prefix has 262 primes: 66 offsets are kept, 196 are lost, 195 are gained, and the destination gets 263 more primes in its tail, for a total of 524. This is an exact diagnostic, not a certificate derived from the previous prime count alone. It keeps mark counts for each matched offset and takes more space on wide intervals.


## Citation and license

Code: MIT. Note and data: CC-BY 4.0.
