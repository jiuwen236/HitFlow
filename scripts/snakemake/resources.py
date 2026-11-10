"""Memory estimates for selected tasks, based on measured workflow runs."""

import math
import re


def memory_mb(value):
    """Parse binary K/M/G/T sizes; bare numbers are MiB."""
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)(?:i?B)?\s*", str(value), re.I)
    if not match:
        raise ValueError(f"Invalid memory size: {value!r}")
    amount, unit = match.groups()
    result = float(amount) * {"": 1, "K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 ** 2}[unit.upper()]
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"Memory must be positive: {value!r}")
    return math.ceil(result)


def retry_memory_mb(base_mb, attempt=1, cap_mb=128 * 1024):
    """Grow by 1.5 per attempt, without lowering the initial estimate."""
    # 1.5x per retry and the 128 GiB default cap are empirical choices;
    # set the cap via cluster.retry_max_mem_mb to match the node memory.
    base_mb, cap_mb = memory_mb(base_mb), memory_mb(cap_mb)
    if not isinstance(attempt, int) or attempt < 1:
        raise ValueError("attempt must be a positive integer")
    if base_mb > cap_mb:
        raise ValueError(
            f"Initial memory {base_mb} MiB exceeds cluster.retry_max_mem_mb={cap_mb}; "
            "increase the limit or reduce the tool's buffers/concurrency."
        )
    result = base_mb
    for _ in range(attempt - 1):
        result = min(cap_mb, math.ceil(result * 1.5))
        if result == cap_mb:
            break
    return result


def samtools_sort_memory_mb(sort_mem, threads, headroom_mb=2048):
    """Budget per-thread buffers plus the main thread and process overhead."""
    # samtools -m is per thread; +1 covers the main thread. Headroom is
    # empirical: 2 GiB for sort alone; the minimap2 | sort pipe passes 24 GiB.
    # That ONT pipe (3G x 16 threads) peaked at 65,561 MiB, so the old 64 GiB
    # request failed; with 24 GiB headroom this gives 75 GiB.
    return memory_mb(sort_mem) * (int(threads) + 1) + memory_mb(headroom_mb)


def sambamba_sort_memory_mb(sort_mem, threads, floor_mem):
    """Sambamba -m is a total budget; allow for piped view and compression."""
    # Unlike samtools, -m is not multiplied by threads. The 2 GiB + 256 MiB
    # per thread margin is empirical.
    return max(memory_mb(floor_mem), memory_mb(sort_mem) + 2048 + int(threads) * 256)


def deepvariant_memory_mb(floor_mem, shards):
    """Reserve 8 GiB for inference/overhead and 2 GiB per concurrent shard."""
    # Empirical: a 16-shard HiFi run peaked at 33,496 MiB; this gives 40 GiB.
    return max(memory_mb(floor_mem), (8 + 2 * int(shards)) * 1024)
