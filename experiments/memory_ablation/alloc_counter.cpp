// Heap accounting for the memory measurement (experiments/memory_ablation, 2026-10-20).
// Linked into a separate copy of research/current/new/gap_native.cpp; the production library is unchanged.
// Every replaceable global operator new/delete of this image is counted by malloc_size, so the counter
// measures the native engine's own heap (std::vector, caches, records), not Python's.
#include <atomic>
#include <cstdlib>
#include <new>
#include <malloc/malloc.h>

namespace {
std::atomic<long long> current{0}, peak{0}, calls{0};
inline void* counted(void* p) {
    if (!p) throw std::bad_alloc();
    long long now = current.fetch_add((long long)malloc_size(p)) + (long long)malloc_size(p);
    long long seen = peak.load();
    while (now > seen && !peak.compare_exchange_weak(seen, now)) {}
    calls.fetch_add(1);
    return p;
}
inline void release(void* p) noexcept {
    if (p) { current.fetch_sub((long long)malloc_size(p)); std::free(p); }
}
void* aligned(std::size_t n, std::size_t a) {
    void* p = nullptr;
    if (posix_memalign(&p, a < sizeof(void*) ? sizeof(void*) : a, n ? n : 1)) return nullptr;
    return p;
}
}  // namespace

void* operator new(std::size_t n) { return counted(std::malloc(n ? n : 1)); }
void* operator new[](std::size_t n) { return counted(std::malloc(n ? n : 1)); }
void* operator new(std::size_t n, std::align_val_t a) { return counted(aligned(n, std::size_t(a))); }
void* operator new[](std::size_t n, std::align_val_t a) { return counted(aligned(n, std::size_t(a))); }
void operator delete(void* p) noexcept { release(p); }
void operator delete[](void* p) noexcept { release(p); }
void operator delete(void* p, std::size_t) noexcept { release(p); }
void operator delete[](void* p, std::size_t) noexcept { release(p); }
void operator delete(void* p, std::align_val_t) noexcept { release(p); }
void operator delete[](void* p, std::align_val_t) noexcept { release(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { release(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { release(p); }

extern "C" {
void mem_reset_peak() { peak.store(current.load()); }
long long mem_current() { return current.load(); }
long long mem_peak() { return peak.load(); }
long long mem_calls() { return calls.load(); }
}
