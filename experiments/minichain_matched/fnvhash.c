/* FNV-1a over the little-endian bytes of int64 values: the same hash the instrumented
   Minichain (vendor/minichain-dp/graphUtils.cpp, namespace recomb_dp) writes to its DP dump. */
#include <stdint.h>

uint64_t fnv_int64(const int64_t *a, long n, uint64_t h)
{
    for (long i = 0; i < n; ++i) {
        uint64_t u = (uint64_t)a[i];
        for (int b = 0; b < 8; ++b) {
            h ^= (u >> (8 * b)) & 0xff;
            h *= 1099511628211ULL;
        }
    }
    return h;
}
