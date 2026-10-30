// Run by test_gap_budget.py. Include the actual implementation to exercise interruption
// inside rehashing without exposing a test-only mutation API to Python callers.
#include "gap_native.cpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {
    PairCache original;
    original.reset(512);  // 1024 slots; all these keys collide under the existing hash
    for (unsigned u = 1; u <= 512; ++u) {
        U64 key = (U64(u) << 32) | 1000;
        size_t p = original.slot(key);
        original.keys[p] = key; original.lo[p] = u; original.hi[p] = u + 17;
        original.plane[p] = int(u - 1); ++original.count;
        original.pool.push_back(u);
    }
    auto validate = [](const PairCache& c) {
        assert(c.count == 512 && c.pool.size() == 512);
        for (unsigned u = 1; u <= 512; ++u) {
            U64 key = (U64(u) << 32) | 1000;
            size_t p = c.slot(key);
            assert(c.keys[p] == key && c.lo[p] == u && c.hi[p] == u + 17);
            assert(c.pool[c.plane[p]] == u);
        }
    };
    WorkBudget<true> full(1e9);
    PairCache successful = original;
    successful.grow(full);
    validate(successful);
    assert(successful.cap == 2048);
    int aborted = 0;
    for (double cap : {0.0, 3071.0, 3072.0, full.work / 2, full.work - 1, full.work}) {
        PairCache c = original;
        WorkBudget<true> limited(cap);
        bool stopped = false;
        try { c.grow(limited); } catch (const BudgetExceeded&) { stopped = true; }
        assert(limited.work <= cap);
        validate(c);
        if (stopped) {
            ++aborted;
            assert(c.cap == original.cap && c.keys == original.keys);
            assert(c.lo == original.lo && c.hi == original.hi && c.plane == original.plane);
            assert(c.pool == original.pool);
            assert(limited.rehashes == 1 && limited.rehashes_done == 0);
            // Reuse after interruption: the old cache remains usable and can still grow.
            WorkBudget<true> retry(1e9);
            c.grow(retry);
            validate(c);
        } else assert(c.cap == 2048 && limited.rehashes_done == 1);
    }
    assert(aborted == 5);
    // A hit at the end of a collision chain must itself be interruptible.
    WorkBudget<true> lookup(8);
    try { original.slot((U64(512) << 32) | 1000, lookup); assert(false); }
    catch (const BudgetExceeded&) { assert(lookup.work == 8 && lookup.cache_probes == 8); }
    validate(original);
    // Exercise complete native calls under sanitizers too: branching geometry, multiple
    // mask widths, starts/ends, local/remote/same-anchor pointers, and warm reuse.
    int native_calls = 0;
    for (int H : {2, 65, 129}) {
        int W = 1; while (64 * W < H) W *= 2;
        std::vector<I> lengths{1, 2, 5, 2}, phi(4 * H), dist(4 * H);
        std::vector<U64> masks(4 * W);
        std::vector<int> last(4 * H);
        for (int h = 0; h < H; ++h) {
            int middle = h % 2 ? 2 : 1;
            for (int v : {0, middle, 3}) masks[v * W + h / 64] |= U64(1) << (h % 64);
            phi[middle * H + h] = 1; phi[3 * H + h] = 1 + lengths[middle];
            last[h] = 0; last[H + h] = middle == 1 ? 1 : 0;
            last[2 * H + h] = middle == 2 ? 2 : 0; last[3 * H + h] = 3;
            dist[(middle == 1 ? 2 : 1) * H + h] = 1;
        }
        double st[40]{};
        void* graph = gap_graph_new(4, H, W, lengths.data(), masks.data(), phi.data(), last.data(), dist.data(), st);
        assert(graph);
        constexpr int n = 120;
        std::mt19937 rng(12345);
        std::vector<I> anchors(6 * n);
        std::vector<int> offsets(n + 1);
        for (int i = 0; i < n; ++i) {
            int v = rng() % 4; I x = 1 + rng() % lengths[v], y = x + rng() % (lengths[v] - x + 1);
            I c = 1 + rng() % 100, d = c + rng() % 6;
            I vals[]{v, x, y, c, d, I(rng() % 100)};
            std::copy(vals, vals + 6, &anchors[6 * i]);
            int deg = v == 0 || v == 3 ? H : (v == 1 ? (H + 1) / 2 : H / 2);
            offsets[i + 1] = offsets[i] + deg;
        }
        int K = offsets[n];
        std::vector<I> C(K), back(K), B(n), R(n), fc(K), fb(K), fB(n), fR(n), cc(K), cb(K), cB(n), cR(n);
        for (I gamma : {I(0), I(3), I(1000)}) {
            assert(gap_ph_solve(graph,n,anchors.data(),gamma,offsets.data(),fc.data(),fb.data(),fB.data(),fR.data(),st,128)==0);
            assert(gap_class_packed(graph,n,anchors.data(),gamma,offsets.data(),cc.data(),cb.data(),cB.data(),cR.data(),st,0)==0);
            assert(cc == fc && cB == fB && cR == fR);
            for (int warm : {0, 1}) for (double factor : {0.0, .1, 1.0, 10.0, 1000.0}) {
                assert(gap_auto_packed(graph,n,anchors.data(),gamma,offsets.data(),C.data(),back.data(),B.data(),R.data(),st,warm,factor,128)==0);
                assert(C == fc && B == fB && R == fR);
                assert(back == (st[7] ? fb : cb));
                assert(st[36] <= st[38]);
                ++native_calls;
            }
        }
        gap_graph_free(graph);
    }
    std::cout << "{\"rehash_abort_boundaries\":" << aborted
              << ",\"successful_rehash_work\":" << full.work
              << ",\"successful_rehash_probes\":" << full.cache_probes
              << ",\"interrupted_lookup\":true,\"native_auto_calls\":" << native_calls << "}\n";
}
