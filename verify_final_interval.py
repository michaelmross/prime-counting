"""Independently check the reported sweep's last interval, not the whole sweep."""

from math import isqrt


def verify_final_interval() -> None:
    p, q = 1_224_709, 1_224_739
    a, b = p * p, q * q
    cap = isqrt(1_500_000_000_000)

    roots = bytearray(b"\x01") * (cap + 1)
    roots[:2] = b"\x00\x00"
    for r in range(2, isqrt(cap) + 1):
        if roots[r]:
            roots[r * r :: r] = bytes((cap - r * r) // r + 1)
    assert [r for r in range(p, cap + 1) if roots[r]] == [p, q]

    # Only odd candidates can be prime here. Slot j represents a + 2 + 2*j.
    lo = a + 2
    candidates = bytearray(b"\x01") * ((b - a) // 2 - 1)
    for r in range(3, q, 2):
        if not roots[r]:
            continue
        first = ((lo + r - 1) // r) * r
        if first % 2 == 0:
            first += r
        j = (first - lo) // 2
        if j < len(candidates):
            candidates[j::r] = bytes((len(candidates) - 1 - j) // r + 1)

    count = candidates.count(1)
    assert count == 2_618_304
    print(f"({a:,}, {b:,}): {count:,} primes")


if __name__ == "__main__":
    verify_final_interval()
