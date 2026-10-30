// Exact gap-sensitive haplotype-aware chaining: native kernels for the contract of
// docs/MATHEMATICAL_MODEL.md ("Gap-sensitive haplotype contract", transcribed from
// Algorithm 2 of the 2024 preprint). C ABI consumed by gap_native.py. This file is new;
// native_chaining.* and class_chaining.* are unchanged.
//
// Two engines share one prepared graph and one per-read setup (validated anchors,
// global END ranks sorted by (vertex, y, id), thr(i), query events with STARTs before
// ENDs at equal coordinates):
//
//  * class: the origin-class engine of gap_chaining.py with the certified fast paths
//    of gap_fastpath.py, "full" plan. Records (K, origin, source, mask[W]) live in one Fenwick tree over END
//    ranks. Merge pairs are decided by (i) a common origin vertex, (ii) per-vertex offset
//    bounds, (iii) the exact range of route lengths between the two origin vertices,
//    cached per vertex pair; only undecided pairs use one bit-sliced threshold on cached
//    route planes. Floors, switches and best states use the same ranges; remote switches
//    are answered once per last2reach group; prefix records that cannot beat the floor
//    are dropped before merging. Decisions, and hence records, equal the specification's.
//  * per_haplotype: query-major per-haplotype Fenwick trees over the same global END
//    ranks, node-major lanes (one row of H lanes per node). Each lane packs the key and
//    the publishing anchor as key*2^24 + (2^24-1-anchor), so one int64 max keeps an
//    attaining witness. The local query is at thr(i) and the update at er(i), over all H
//    lanes (vectorized) when the anchor has many labels and over its own labels otherwise.
//    Each missing label g is queried in its own lane at the last END on last(v_i, g); labels
//    sharing that vertex are walked together, reading node rows contiguously.
//
//  * auto (added 2026-10-06): the class engine under a deterministic work budget, a factor
//    times the per-haplotype engine's lane operations on the read; if exceeded, the read is
//    re-solved by the per-haplotype engine. Work counting is compiled only into this mode
//    (ClassEngine<W, true>), so the class modes run unchanged code. See gap_auto_packed.
//
// Integers only. All arithmetic is int64; the adapter rejects inputs whose explicit bound
// exceeds 2^60 (rechecked here with overflow-checked sums), and the per-haplotype lanes
// additionally need keys < 2^39 and fewer than 2^24 anchors. Nothing is rounded.
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using I = int64_t;
using U64 = uint64_t;
constexpr I NEG = -(I(1) << 61);
constexpr I LIMIT = I(1) << 60;
constexpr int SH = 24;                       // per-haplotype lane: key << SH | (JM - anchor)
constexpr I JM = (I(1) << SH) - 1;
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
void require(bool ok, const char* why) { if (!ok) throw std::invalid_argument(why); }
thread_local std::string last_error;

I checked_add(I a, I b) {
    I r;
    if (__builtin_add_overflow(a, b, &r) || r >= LIMIT || r <= -LIMIT) throw std::invalid_argument("value bound exceeds 2^60");
    return r;
}

// ---------------------------------------------------------------- graph index
struct Graph {
    int V = 0, H = 0, W = 0, b = 0;
    I S = 0;
    std::vector<I> len;           // V
    std::vector<U64> vmask;       // V*W, bit h: haplotype h visits v
    std::vector<I> phi;           // V*H, characters before v on h (0 off-path)
    // When S < 2^29: int32 offsets in rows of 64W lanes, with sentinels P = 2^30-1 off-path
    // (+P in phi_plus, -P in phi_minus), so masked range scans are plain reductions.
    std::vector<int32_t> phi_plus, phi_minus;
    std::vector<I> vlo, vhi;      // min/max phi over A_v
    std::vector<U64> planes;      // (V*b + t)*W: bit t of phi_h(v), h in A_v
    std::vector<int> loff, lab;   // sorted label lists per vertex
    std::vector<int> wrank;       // V*W: popcount of the words before q
    std::vector<int> goff, gw;    // class remote groups per vertex: last2reach vertex,
    std::vector<I> gd;            //   D(w, v),
    std::vector<U64> gmask;       //   and the missing labels with that last vertex (W words)
    std::vector<int> roff, rg, rw;// per-haplotype remote entries per vertex: label, last vertex,
    std::vector<I> rcost;         //   phi_g(w) + D(w, v)
    void* warm = nullptr;         // optional graph-persistent pair cache (PairCache*)
    int degree(int v) const { return loff[v + 1] - loff[v]; }
    const U64* mask(int v) const { return &vmask[size_t(v) * W]; }
    int rank(int v, int h) const {  // position of label h in v's sorted label list
        int q = h >> 6;
        return wrank[size_t(v) * W + q] + std::popcount(vmask[size_t(v) * W + q] & ((U64(1) << (h & 63)) - 1));
    }
};

Graph* build_graph(int V, int H, int W, const I* lengths, const U64* vmask, const I* phi,
                   const int* last, const I* dist, double* stats) {
    auto t = Clock::now();
    require(V > 0 && H > 0 && (W == 1 || W == 2 || W == 4 || W == 8 || W == 16) && H <= 64 * W, "invalid graph dimensions");
    auto g = new Graph();
    try {
        g->V = V; g->H = H; g->W = W;
        g->len.assign(lengths, lengths + V);
        g->vmask.assign(vmask, vmask + size_t(V) * W);
        I S = 0;
        for (int v = 0; v < V; ++v) { require(lengths[v] > 0, "vertex lengths must be positive"); S = checked_add(S, lengths[v]); }
        require(S < (I(1) << 40), "total graph length exceeds 2^40");
        g->S = S;
        g->b = std::max(1, int(std::bit_width(uint64_t(S))));
        int b = g->b;
        g->phi.assign(size_t(V) * H, 0);
        g->vlo.assign(V, 0); g->vhi.assign(V, 0);
        g->planes.assign(size_t(V) * b * W, 0);
        g->loff.assign(V + 1, 0); g->wrank.assign(size_t(V) * W, 0);
        for (int v = 0; v < V; ++v) {
            const U64* row = g->mask(v);
            int r = 0; bool first = true;
            for (int q = 0; q < W; ++q) {
                g->wrank[size_t(v) * W + q] = r;
                for (U64 bits = row[q]; bits; bits &= bits - 1) {
                    int h = q * 64 + std::countr_zero(bits);
                    require(h < H, "mask bit beyond H");
                    I off = phi[size_t(v) * H + h];
                    require(off >= 0 && off + lengths[v] <= S, "offset out of range");
                    g->phi[size_t(v) * H + h] = off;
                    if (first || off < g->vlo[v]) g->vlo[v] = off;
                    if (first || off > g->vhi[v]) g->vhi[v] = off;
                    first = false;
                    for (int tt = 0; tt < b; ++tt)
                        if ((off >> tt) & 1) g->planes[(size_t(v) * b + tt) * W + q] |= U64(1) << (h & 63);
                    g->lab.push_back(h); ++r;
                }
            }
            g->loff[v + 1] = int(g->lab.size());
        }
        if (S < (I(1) << 29)) {
            constexpr int32_t BIG = (int32_t(1) << 30) - 1;
            g->phi_plus.assign(size_t(V) * 64 * W, BIG); g->phi_minus.assign(size_t(V) * 64 * W, -BIG);
            for (int v = 0; v < V; ++v)
                for (int h = 0; h < H; ++h)
                    if ((g->vmask[size_t(v) * W + (h >> 6)] >> (h & 63)) & 1) {
                        g->phi_plus[size_t(v) * 64 * W + h] = g->phi_minus[size_t(v) * 64 * W + h] = int32_t(g->phi[size_t(v) * H + h]);
                    }
        }
        // Remote structures. For g not on v with w = last(v, g) >= 0, dist[v][g] is the
        // shortest distance D(w, v); it must agree within each w group.
        g->goff.assign(V + 1, 0); g->roff.assign(V + 1, 0);
        std::vector<int> ws; std::vector<I> ds; std::vector<U64> ms;
        struct RemoteEntry { int w, g; I cost; };
        std::vector<RemoteEntry> rent;
        for (int v = 0; v < V; ++v) {
            const U64* row = g->mask(v);
            ws.clear(); ds.clear(); ms.clear();
            for (int h = 0; h < H; ++h) {
                if ((row[h >> 6] >> (h & 63)) & 1) continue;
                int w = last[size_t(v) * H + h];
                if (w < 0) continue;
                require(w < v && ((g->mask(w)[h >> 6] >> (h & 63)) & 1), "last2reach vertex must precede v on the path");
                I d = dist[size_t(v) * H + h];
                require(d > 0 && d <= S, "last2reach distance out of range");
                size_t k = 0;
                while (k < ws.size() && ws[k] != w) ++k;
                if (k == ws.size()) { ws.push_back(w); ds.push_back(d); ms.resize(ms.size() + W, 0); }
                require(ds[k] == d, "last2reach distance must depend on the last vertex only");
                ms[k * W + (h >> 6)] |= U64(1) << (h & 63);
                rent.push_back({w, h, checked_add(g->phi[size_t(w) * H + h], d)});
            }
            // Per-haplotype remote entries grouped by last vertex (then label), so the baseline
            // walks each group's prefix nodes once and reads node rows contiguously.
            std::sort(rent.begin(), rent.end(), [](const auto& p, const auto& q) { return std::tie(p.w, p.g) < std::tie(q.w, q.g); });
            for (const auto& r : rent) { g->rg.push_back(r.g); g->rw.push_back(r.w); g->rcost.push_back(r.cost); }
            rent.clear();
            std::vector<int> order(ws.size());
            for (size_t k = 0; k < ws.size(); ++k) order[k] = int(k);
            std::sort(order.begin(), order.end(), [&](int a, int c) { return ws[a] < ws[c]; });
            for (int k : order) {
                g->gw.push_back(ws[k]); g->gd.push_back(ds[k]);
                g->gmask.insert(g->gmask.end(), ms.begin() + size_t(k) * W, ms.begin() + size_t(k + 1) * W);
            }
            g->goff[v + 1] = int(g->gw.size());
            g->roff[v + 1] = int(g->rg.size());
        }
        stats[0] = elapsed(t); stats[1] = double(S); stats[2] = b; stats[3] = double(g->gw.size()); stats[4] = double(g->rg.size());
        return g;
    } catch (...) { delete g; throw; }
}

// ---------------------------------------------------------------- per-read setup
struct Anchor { int v; I x, y, c, d, w; };
struct Read {
    int n = 0;
    std::vector<Anchor> A;
    std::vector<int> events;      // 2i: START(i), 2i+1: END(i), in query order
    std::vector<int> er, thr;     // 1-based END rank; number of ENDs before (v_i, x_i)
    std::vector<int> endv;        // vertex of the END with rank r+1
    I maxkey = 0;                 // bound on every published key C + y + d + phi
    int tau(int w) const {        // ENDs on vertices <= w
        return int(std::upper_bound(endv.begin(), endv.end(), w) - endv.begin());
    }
};

// LSD radix sort of distinct packed keys using their low `bits` bits.
void radix_sort(std::vector<U64>& a, unsigned bits) {
    if (a.size() < 256) { std::sort(a.begin(), a.end()); return; }
    constexpr unsigned D = 11; constexpr size_t NB = size_t(1) << D;
    std::vector<U64> tmp(a.size());
    std::vector<uint32_t> cnt(NB);
    for (unsigned shift = 0; shift < bits; shift += D) {
        std::fill(cnt.begin(), cnt.end(), 0u);
        for (U64 x : a) ++cnt[(x >> shift) & (NB - 1)];
        uint32_t sum = 0;
        for (auto& c : cnt) { uint32_t t = c; c = sum; sum += t; }
        for (U64 x : a) tmp[cnt[(x >> shift) & (NB - 1)]++] = x;
        a.swap(tmp);
    }
}

void sort_by_keys(std::vector<int>& items, int n, const std::vector<std::tuple<I, I, I>>& key) {
    // key[k] = (major, minor, phase) for item k; ties by item id. Packs into one word when it fits.
    I m0 = 0, m1 = 0, m2 = 0;
    for (auto& [a, c, d] : key) { m0 = std::max(m0, a); m1 = std::max(m1, c); m2 = std::max(m2, d); }
    unsigned b0 = std::bit_width(uint64_t(m0)), b1 = std::bit_width(uint64_t(m1)), b2 = std::bit_width(uint64_t(m2));
    unsigned bi = std::bit_width(unsigned(std::max(1, int(items.size()))));
    (void)n;
    if (b0 + b1 + b2 + bi <= 64) {
        std::vector<U64> packed(items.size());
        for (size_t k = 0; k < items.size(); ++k) {
            auto& [a, c, d] = key[k];
            packed[k] = ((((U64(a) << b1) | U64(c)) << b2 | U64(d)) << bi) | U64(k);
        }
        radix_sort(packed, b0 + b1 + b2 + bi);
        std::vector<int> out(items.size());
        U64 low = bi == 64 ? ~U64(0) : ((U64(1) << bi) - 1);
        for (size_t k = 0; k < items.size(); ++k) out[k] = items[packed[k] & low];
        items.swap(out);
        return;
    }
    std::vector<int> idx(items.size());
    for (size_t k = 0; k < idx.size(); ++k) idx[k] = int(k);
    std::sort(idx.begin(), idx.end(), [&](int p, int q) { return std::tie(std::get<0>(key[p]), std::get<1>(key[p]), std::get<2>(key[p]), p) <
                                                                 std::tie(std::get<0>(key[q]), std::get<1>(key[q]), std::get<2>(key[q]), q); });
    std::vector<int> out(items.size());
    for (size_t k = 0; k < idx.size(); ++k) out[k] = items[idx[k]];
    items.swap(out);
}

void prepare_read(Read& R, const Graph& G, int n, const I* input, I gamma, double* stats) {
    auto t = Clock::now();
    require(n >= 0 && n < (1 << 27), "anchor count out of range");
    require(gamma >= 0 && gamma < LIMIT, "gamma must be a nonnegative integer below 2^60");
    R.n = n; R.A.resize(n);
    I sumw = 0, mx = 0, my = 0, mc = 0, md = 0, mw = 0;
    for (int i = 0; i < n; ++i) {
        const I* q = input + 6 * size_t(i);
        require(q[0] >= 0 && q[0] < G.V, "anchor vertex out of bounds");
        Anchor a{int(q[0]), q[1], q[2], q[3], q[4], q[5]};
        require(G.degree(a.v) > 0, "anchor vertex lies on no haplotype");
        require(1 <= a.x && a.x <= a.y && a.y <= G.len[a.v], "invalid graph interval");
        require(1 <= a.c && a.c <= a.d && a.d < LIMIT, "invalid query interval");
        require(a.w >= 0 && a.w < LIMIT, "weights must be nonnegative integers below 2^60");
        R.A[i] = a;
        sumw = checked_add(sumw, a.w);
        mx = std::max(mx, a.x); my = std::max(my, a.y); mc = std::max(mc, a.c); md = std::max(md, a.d); mw = std::max(mw, a.w);
    }
    // Every intermediate value lies within +-4*bound (see gap_native.py); bound < 2^60.
    I bound = checked_add(checked_add(checked_add(sumw, 2 * G.S), checked_add(checked_add(mx, my), checked_add(mc, md))), checked_add(checked_add(mw, gamma), 16));
    (void)bound;
    R.maxkey = checked_add(checked_add(sumw, G.S), checked_add(my, md));
    double K = 0;
    for (int i = 0; i < n; ++i) K += G.degree(R.A[i].v);
    stats[5] = K;
    // Graph events by (vertex, coordinate, phase START=0 < END=1, id): END ranks follow
    // (vertex, y, id), and thr(i) = #ENDs before START(i) = #ENDs with (v, y) < (v_i, x_i).
    std::vector<int> gev(2 * size_t(n));
    std::vector<std::tuple<I, I, I>> gk(2 * size_t(n));
    for (int i = 0; i < n; ++i) {
        gev[2 * i] = 2 * i; gk[2 * i] = {R.A[i].v, R.A[i].x, 0};
        gev[2 * i + 1] = 2 * i + 1; gk[2 * i + 1] = {R.A[i].v, R.A[i].y, 1};
    }
    sort_by_keys(gev, 2 * n, gk);
    R.er.assign(n, 0); R.thr.assign(n, 0); R.endv.assign(n, 0);
    int ends = 0;
    for (int e : gev) {
        int i = e >> 1;
        if (e & 1) { R.endv[ends] = R.A[i].v; R.er[i] = ++ends; } else R.thr[i] = ends;
    }
    // Query events by (coordinate, phase START=0 < END=1, id).
    R.events.resize(2 * size_t(n));
    std::vector<std::tuple<I, I, I>> qk(2 * size_t(n));
    for (int i = 0; i < n; ++i) {
        R.events[2 * i] = 2 * i; qk[2 * i] = {R.A[i].c, 0, 0};
        R.events[2 * i + 1] = 2 * i + 1; qk[2 * i + 1] = {R.A[i].d, 1, 0};
    }
    // sort_by_keys breaks ties by position; positions 2i/2i+1 order ties by id then phase,
    // and the phase is already part of the key, so ties are by id.
    sort_by_keys(R.events, 2 * n, qk);
    stats[1] = elapsed(t);
}

// ---------------------------------------------------------------- masks
template<int W> struct Mask {
    U64 w[W];
    static Mask zero() { Mask m; for (int k = 0; k < W; ++k) m.w[k] = 0; return m; }
    static Mask load(const U64* p) { Mask m; for (int k = 0; k < W; ++k) m.w[k] = p[k]; return m; }
    bool any() const { U64 z = 0; for (int k = 0; k < W; ++k) z |= w[k]; return z != 0; }
    Mask operator&(const Mask& o) const { Mask m; for (int k = 0; k < W; ++k) m.w[k] = w[k] & o.w[k]; return m; }
    Mask operator|(const Mask& o) const { Mask m; for (int k = 0; k < W; ++k) m.w[k] = w[k] | o.w[k]; return m; }
    Mask andnot(const Mask& o) const { Mask m; for (int k = 0; k < W; ++k) m.w[k] = w[k] & ~o.w[k]; return m; }
    Mask& operator|=(const Mask& o) { for (int k = 0; k < W; ++k) w[k] |= o.w[k]; return *this; }
    bool operator==(const Mask& o) const { U64 z = 0; for (int k = 0; k < W; ++k) z |= w[k] ^ o.w[k]; return z == 0; }
    int lowest() const { for (int k = 0; k < W; ++k) if (w[k]) return k * 64 + std::countr_zero(w[k]); return -1; }
};

// Node / candidate record: key K (or state coefficient), origin anchor and its vertex,
// source anchor (the publishing anchor of a node record).
template<int W> struct Rec { I k; int o; int s; int ov; Mask<W> m; };
enum Kind { RESTART = 0, LOCAL = 1, REMOTE = 2, SAME = 3 };
// State record at an anchor: C(i,h) = k - (phi_h(v_i) - phi_h(v_o)) for h in m. During
// START(i), [lo, hi] is the range of l_h(v_o, v_i) over A_{v_o} & A_{v_i}.
template<int W> struct St { I k; int o, kind, j, g, w, ov; I lo, hi; Mask<W> m; };

enum Counter {
    C_MERGES = 0, C_CONTESTED, C_PAIRS, C_L1, C_L2, C_L3, C_L4, C_THRESH, C_MINIMA, C_RANGES, C_LABELS,
    C_PLANES, C_CUT, C_PREFIX, C_GROUPS, C_RCANDS, C_RSLICED, C_PUBLISHED, C_OUTPUT, C_MAXNODE, C_COUNT
};

// Thrown by the class engine when its work exceeds the budget set by gap_auto_packed.
struct BudgetExceeded {};

// Exact route-length ranges and route planes per vertex pair u < w (open addressing).
struct PairCache {
    std::vector<U64> keys; std::vector<I> lo, hi; std::vector<int> plane;
    std::vector<U64> pool; size_t count = 0, cap = 0;
    void reset(size_t want) {
        cap = 1024; while (cap < 2 * want) cap <<= 1;
        keys.assign(cap, ~U64(0)); lo.assign(cap, 0); hi.assign(cap, 0); plane.assign(cap, -1); pool.clear(); count = 0;
    }
    size_t slot(U64 key) const {
        size_t p = size_t((key * 0x9E3779B97F4A7C15ULL) >> 20) & (cap - 1);
        while (keys[p] != key && keys[p] != ~U64(0)) p = (p + 1) & (cap - 1);
        return p;
    }
    void grow() {
        std::vector<U64> k2; std::vector<I> l2, h2; std::vector<int> p2;
        k2.swap(keys); l2.swap(lo); h2.swap(hi); p2.swap(plane);
        size_t old = cap; cap <<= 1;
        keys.assign(cap, ~U64(0)); lo.assign(cap, 0); hi.assign(cap, 0); plane.assign(cap, -1);
        for (size_t s = 0; s < old; ++s) if (k2[s] != ~U64(0)) { size_t p = slot(k2[s]); keys[p] = k2[s]; lo[p] = l2[s]; hi[p] = h2[s]; plane[p] = p2[s]; }
    }
};

template<typename T, typename Less> void insertion_sort(std::vector<T>& a, Less less) {
    for (size_t k = 1; k < a.size(); ++k) {
        T x = a[k]; size_t j = k;
        while (j > 0 && less(x, a[j - 1])) { a[j] = a[j - 1]; --j; }
        a[j] = x;
    }
}

template<int W, bool Budget = false> struct ClassEngine {
    const Graph& G; const Read& R; I gamma;
    int n, H, b;
    std::vector<std::vector<Rec<W>>> node;
    std::vector<Mask<W>> nodeU;
    std::vector<St<W>> pool;      // final state records of every anchor
    std::vector<int> soff, scnt;  // anchor -> pool range
    std::vector<I> B, Rv;
    PairCache own;                // cold: a fresh cache per read
    PairCache& cache;             // or the graph's persistent (warm) cache
    double cnt[C_COUNT] = {0};
    // Work in machine words (a W-word mask operation counts W; scalar steps count 1), checked
    // against `budget` once per event. The default budget is unlimited; auto mode sets it.
    double work = 0, budget = std::numeric_limits<double>::infinity();
    // Counting is compiled in only for the auto mode, so the class modes run unchanged code.
    void spend(double x) { if constexpr (Budget) work += x; }
    void check() const { if constexpr (Budget) { if (work > budget) throw BudgetExceeded{}; } }
    // scratch
    std::vector<Rec<W>> acc, filt, inc;
    std::vector<St<W>> pre, post;
    std::vector<unsigned char> flag;
    struct Cand { I ub, k; int order, s, vo; I lo, hi; Mask<W> m; };
    std::vector<Cand> cands;
    std::vector<std::pair<I, int>> best_order;
    // Direct-mapped memo of range(ov, v_i) during the current START.
    int memo_i[64], memo_ov[64]; I memo_lo[64], memo_hi[64];

    ClassEngine(const Graph& g, const Read& r, I gm, PairCache* shared)
        : G(g), R(r), gamma(gm), n(r.n), H(g.H), b(g.b), cache(shared ? *shared : own) {
        for (int k = 0; k < 64; ++k) memo_i[k] = -1;
        node.assign(n + 1, {}); nodeU.assign(n + 1, Mask<W>::zero());
        soff.assign(n, 0); scnt.assign(n, 0); B.assign(n, NEG); Rv.assign(n, NEG);
        if (!shared) cache.reset(2048);
        pool.reserve(size_t(n) * 2);
    }
    int vert(int anchor) const { return R.A[anchor].v; }

    // ---- exact pair ranges and route planes (u < w) ----
    size_t pair_slot(int u, int w) {
        U64 key = (U64(unsigned(u)) << 32) | unsigned(w);
        size_t p = cache.slot(key);
        if (cache.keys[p] == key) return p;
        if (2 * (cache.count + 1) > cache.cap) { cache.grow(); p = cache.slot(key); }
        const U64* mu = G.mask(u); const U64* mw = G.mask(w);
        I lo = std::numeric_limits<I>::max(), hi = std::numeric_limits<I>::min(); long labels = 0;
        if (!G.phi_plus.empty()) {
            // On-path differences lie in (-S, S), S < 2^29. A lane off either path gives at
            // least P - S > 2^29 to the minimum scan and at most S - P < -2^29 to the maximum
            // scan, so it is neutral; sentinel sums stay within +-(2P) < 2^31.
            const size_t row = size_t(64) * W;
            const int32_t* up = &G.phi_plus[u * row]; const int32_t* um = &G.phi_minus[u * row];
            const int32_t* wp = &G.phi_plus[w * row]; const int32_t* wm = &G.phi_minus[w * row];
            int32_t l32 = std::numeric_limits<int32_t>::max(), h32 = std::numeric_limits<int32_t>::min();
            for (int q = 0; q < W; ++q) {
                U64 bits = mu[q] & mw[q];
                if (!bits) continue;
                labels += std::popcount(bits);
                int32_t l = l32, h = h32;
                for (int k = 64 * q; k < 64 * q + 64; ++k) {
                    l = std::min(l, wp[k] - um[k]);
                    h = std::max(h, wm[k] - up[k]);
                }
                l32 = l; h32 = h;
            }
            lo = l32; hi = h32;
        } else {
            const I* pu = &G.phi[size_t(u) * H]; const I* pw = &G.phi[size_t(w) * H];
            for (int q = 0; q < W; ++q)
                for (U64 bits = mu[q] & mw[q]; bits; bits &= bits - 1) {
                    int hh = q * 64 + std::countr_zero(bits);
                    I d = pw[hh] - pu[hh];
                    lo = std::min(lo, d); hi = std::max(hi, d); ++labels;
                }
        }
        spend(G.phi_plus.empty() ? double(labels) : double(64 * W));
        require(labels > 0, "pair range requested for vertices without a common haplotype");
        cache.keys[p] = key; cache.lo[p] = lo; cache.hi[p] = hi; cache.plane[p] = -1; ++cache.count;
        cnt[C_RANGES] += 1; cnt[C_LABELS] += double(labels);
        return p;
    }
    std::pair<I, I> target_range(int ov, int i, int vi) {  // range(ov, vi) at START(i)
        int k = ov & 63;
        if (memo_i[k] == i && memo_ov[k] == ov) return {memo_lo[k], memo_hi[k]};
        auto r = range(ov, vi);
        memo_i[k] = i; memo_ov[k] = ov; memo_lo[k] = r.first; memo_hi[k] = r.second;
        return r;
    }
    std::pair<I, I> range(int u, int w) {  // range of phi_h(w) - phi_h(u) over A_u & A_w, any order
        if (u == w) return {0, 0};
        if (u < w) { size_t p = pair_slot(u, w); return {cache.lo[p], cache.hi[p]}; }
        size_t p = pair_slot(w, u); return {-cache.hi[p], -cache.lo[p]};
    }
    const U64* route_planes(int u, int w) {  // b planes of phi_h(w) - phi_h(u), u < w
        require(u < w, "route planes need u < w (invariant R+)");
        size_t p = pair_slot(u, w);
        if (cache.plane[p] < 0) {
            spend(double(b) * W);
            int base = int(cache.pool.size());
            cache.pool.resize(cache.pool.size() + size_t(b) * W);
            U64* out = &cache.pool[base];
            const U64* X = &G.planes[size_t(w) * b * W]; const U64* Y = &G.planes[size_t(u) * b * W];
            U64 borrow[W];
            for (int q = 0; q < W; ++q) borrow[q] = 0;
            for (int t = 0; t < b; ++t)
                for (int q = 0; q < W; ++q) {
                    U64 x = X[t * W + q], y = Y[t * W + q], c = borrow[q];
                    out[t * W + q] = x ^ y ^ c;
                    borrow[q] = (~x & (y | c)) | (y & c);
                }
            cache.plane[p] = base; cnt[C_PLANES] += 1;
        }
        return &cache.pool[cache.plane[p]];
    }
    // {h in M : value(h) < theta}, values in [0, S] stored in b planes.
    Mask<W> threshold(const U64* P, I theta, const Mask<W>& M) {
        if (theta <= 0) return Mask<W>::zero();
        if (theta > G.S) return M;
        cnt[C_THRESH] += 1;
        spend(double(b) * W);
        Mask<W> lt = Mask<W>::zero(), eq = M;
        for (int t = b - 1; t >= 0; --t) {
            const Mask<W> pl = Mask<W>::load(P + size_t(t) * W);
            if ((theta >> t) & 1) { lt |= eq.andnot(pl); eq = eq & pl; }
            else eq = eq.andnot(pl);
            if (!eq.any()) break;
        }
        return lt;
    }
    // Masked minimum and its lowest attaining label.
    std::pair<I, int> minimum(const U64* P, const Mask<W>& M) {
        cnt[C_MINIMA] += 1;
        spend(double(b) * W);
        Mask<W> live = M; I val = 0;
        for (int t = b - 1; t >= 0; --t) {
            const Mask<W> pl = Mask<W>::load(P + size_t(t) * W);
            Mask<W> z = live.andnot(pl);
            if (z.any()) live = z; else { val |= I(1) << t; live = live & pl; }
        }
        return {val, live.lowest()};
    }

    // ---- merge: L := pointwise max(L, R), ties to L ----
    // Right wins at h iff delta(h) = phi_h(ub) - phi_h(ua) > t = K_a - K_b.
    int decide(int ua, int ub, I t) {
        if (ua == ub) { cnt[C_L1] += 1; return t < 0 ? 1 : 0; }
        I lo = G.vlo[ub] - G.vhi[ua], hi = G.vhi[ub] - G.vlo[ua];
        if (lo > t) { cnt[C_L2] += 1; return 1; }
        if (hi <= t) { cnt[C_L2] += 1; return 0; }
        auto [rlo, rhi] = range(ua, ub);
        if (rlo > t) { cnt[C_L3] += 1; return 1; }
        if (rhi <= t) { cnt[C_L3] += 1; return 0; }
        cnt[C_L4] += 1;
        return -1;
    }
    Mask<W> sliced_right_wins(int ua, int ub, I t, const Mask<W>& inter) {
        if (ua < ub) return inter.andnot(threshold(route_planes(ua, ub), t + 1, inter));  // l > t
        return threshold(route_planes(ub, ua), -t, inter);                                 // -l > t
    }
    void merge(std::vector<Rec<W>>& L, Mask<W>& LU, const Rec<W>* Rr, int nr, const Mask<W>& RU) {
        cnt[C_MERGES] += 1;
        spend(W);
        Mask<W> ov = LU & RU;
        if (!ov.any()) { L.insert(L.end(), Rr, Rr + nr); LU |= RU; return; }
        cnt[C_CONTESTED] += 1;
        Mask<W> win = RU.andnot(LU);
        int nl = int(L.size());
        if (int(flag.size()) < nl) flag.resize(nl);
        for (int a = 0; a < nl; ++a) {
            check();
            Mask<W> ia = L[a].m & ov;
            spend(W);
            flag[a] = ia.any();
            if (!flag[a]) continue;
            int ua = L[a].ov;
            for (int c = 0; c < nr; ++c) {
                Mask<W> inter = ia & Rr[c].m;
                spend(W);
                if (!inter.any()) continue;
                cnt[C_PAIRS] += 1;
                I t = L[a].k - Rr[c].k;
                int ub = Rr[c].ov;
                int verdict = decide(ua, ub, t);
                if (verdict > 0) win |= inter;
                else if (verdict < 0) win |= sliced_right_wins(ua, ub, t, inter);
            }
        }
        spend(double(W) * double(nl + nr));
        int out = 0;
        for (int a = 0; a < nl; ++a) {
            if (flag[a]) { L[a].m = L[a].m.andnot(win); if (!L[a].m.any()) continue; }
            if (out != a) L[out] = L[a];
            ++out;
        }
        L.resize(out);
        for (int c = 0; c < nr; ++c) {
            Mask<W> keep = Rr[c].m & win;
            if (keep.any()) { L.push_back(Rr[c]); L.back().m = keep; }
        }
        LU |= RU;
    }

    // ---- remote switches, one pass per last2reach group ----
    struct RemoteBest { bool have = false; I cand = 0; int g = 0, order = 0, j = -1, w = -1; };
    RemoteBest remote(int i, I alpha) {
        RemoteBest best;
        int v = R.A[i].v;
        for (int gi = G.goff[v]; gi < G.goff[v + 1]; ++gi) {
            cnt[C_GROUPS] += 1;
            spend(W);
            check();
            int w = G.gw[gi]; I cst = alpha - G.gd[gi] - gamma;
            const Mask<W> Gm = Mask<W>::load(&G.gmask[size_t(gi) * W]);
            cands.clear();
            int order = 0;
            for (int x = R.tau(w); x > 0; x -= x & -x, ++order) {
                spend(W);
                if (!(nodeU[x] & Gm).any()) continue;
                for (const auto& rec : node[x]) {
                    Mask<W> inter = rec.m & Gm;
                    spend(W);
                    if (!inter.any()) continue;
                    auto [lo, hi] = range(rec.ov, w);
                    cands.push_back({rec.k - lo, rec.k, order, rec.s, rec.ov, lo, hi, inter});
                }
            }
            cnt[C_RCANDS] += double(cands.size());
            insertion_sort(cands, [](const Cand& p, const Cand& q) { return p.ub > q.ub; });
            spend(0.5 * double(cands.size()) * double(cands.size()) + double(cands.size()));
            for (const auto& c : cands) {
                int low = c.m.lowest();
                I ubc = c.ub + cst;
                if (best.have && std::make_tuple(ubc, -low, -c.order) < std::make_tuple(best.cand, -best.g, -best.order)) continue;
                I val; int g;
                if (c.lo == c.hi) { val = c.ub; g = low; }
                else { auto [mn, gg] = minimum(route_planes(c.vo, w), c.m); val = c.k - mn; g = gg; cnt[C_RSLICED] += 1; }
                I cand = val + cst;
                if (!best.have || std::make_tuple(cand, -g, -c.order) > std::make_tuple(best.cand, -best.g, -best.order))
                    best = {true, cand, g, c.order, c.s, w};
            }
        }
        return best;
    }

    // Keep r on {h : l_h(v_o, v_i) < theta}, using r's stored range.
    Mask<W> keep_below(const St<W>& r, int vi, I theta) {
        if (r.hi < theta) return r.m;
        if (r.lo >= theta) return Mask<W>::zero();
        return threshold(route_planes(r.ov, vi), theta, r.m);
    }

    void start(int i) {
        const Anchor& a = R.A[i];
        int vi = a.v;
        const Mask<W> Ai = Mask<W>::load(G.mask(vi));
        I alpha = a.w - a.x - a.c + 2;
        RemoteBest rb = remote(i, alpha);
        Rv[i] = rb.have ? rb.cand : NEG;
        I f = a.w; int fk = RESTART, fj = -1, fg = -1, fw = -1;
        if (rb.have && rb.cand > f) { f = rb.cand; fk = REMOTE; fj = rb.j; fg = rb.g; fw = rb.w; }
        // Local prefix query; records that cannot beat the floor are cut first.
        I bound = f - alpha;
        acc.clear(); Mask<W> accU = Mask<W>::zero();
        for (int x = R.thr[i]; x > 0; x -= x & -x) {
            spend(1);
            check();
            filt.clear(); Mask<W> fU = Mask<W>::zero();
            for (const auto& z : node[x]) {
                cnt[C_PREFIX] += 1;
                spend(W);
                if (z.k <= bound) { cnt[C_CUT] += 1; continue; }
                Mask<W> m = z.m & Ai;
                if (!m.any()) continue;
                if (z.k - target_range(z.ov, i, vi).first <= bound) { cnt[C_CUT] += 1; continue; }
                filt.push_back({z.k, z.o, z.s, z.ov, m}); fU |= m;
            }
            if (filt.empty()) continue;
            if (acc.empty()) { acc.swap(filt); accU = fU; }
            else merge(acc, accU, filt.data(), int(filt.size()), fU);
        }
        // Floor: keep a local record where N - l > f.
        pre.clear(); Mask<W> used = Mask<W>::zero();
        for (const auto& r : acc) {
            auto [lo, hi] = target_range(r.ov, i, vi);
            St<W> s{r.k + alpha, r.o, LOCAL, r.s, -1, -1, r.ov, lo, hi, r.m};
            Mask<W> keep = keep_below(s, vi, s.k - f);
            if (keep.any()) { s.m = keep; pre.push_back(s); used |= keep; }
        }
        Mask<W> rest = Ai.andnot(used);
        if (rest.any()) pre.push_back({f, i, fk, fj, fg, fw, vi, 0, 0, rest});
        spend(double(W) * double(acc.size()));
        // Best state, best-first with range bounds; ties to the smallest label.
        best_order.clear();
        spend(0.5 * double(pre.size()) * double(pre.size()) + double(W + 1) * double(pre.size()));
        for (int k = 0; k < int(pre.size()); ++k) best_order.push_back({pre[k].k - pre[k].lo, k});
        insertion_sort(best_order, [](const std::pair<I, int>& p, const std::pair<I, int>& q) { return p.first > q.first; });
        bool have = false; I Bi = 0; int bh = -1;
        for (auto [ub, k] : best_order) {
            const auto& r = pre[k];
            int low = r.m.lowest();
            if (have && (ub < Bi || (ub == Bi && low > bh))) continue;
            I val; int h;
            if (r.lo == r.hi) { val = ub; h = low; }
            else { auto [mn, hh] = minimum(route_planes(r.ov, vi), r.m); val = r.k - mn; h = hh; }
            if (!have || val > Bi || (val == Bi && h < bh)) { have = true; Bi = val; bh = h; }
        }
        B[i] = Bi;
        // Same-anchor switch: keep h where N - l >= s; the old state wins ties.
        I s = Bi - gamma;
        post.clear(); used = Mask<W>::zero();
        for (const auto& r : pre) {
            Mask<W> keep = keep_below(r, vi, r.k - s + 1);
            if (keep.any()) { post.push_back(r); post.back().m = keep; used |= keep; }
        }
        rest = Ai.andnot(used);
        if (rest.any()) post.push_back({s, i, SAME, i, bh, -1, vi, 0, 0, rest});
        // Coalesce equal (value, origin, witness); store the anchor's final states.
        spend(double(post.size()) * double(post.size()));
        soff[i] = int(pool.size());
        for (const auto& r : post) {
            bool merged = false;
            for (int k = soff[i]; k < int(pool.size()); ++k) {
                auto& q = pool[k];
                if (q.k == r.k && q.o == r.o && q.kind == r.kind && q.j == r.j && q.g == r.g && q.w == r.w) { q.m |= r.m; merged = true; break; }
            }
            if (!merged) pool.push_back(r);
        }
        scnt[i] = int(pool.size()) - soff[i];
        cnt[C_OUTPUT] += scnt[i];
    }

    void end(int i) {
        const Anchor& a = R.A[i];
        inc.clear(); Mask<W> incU = Mask<W>::zero();
        for (int k = soff[i]; k < soff[i] + scnt[i]; ++k) {
            const auto& r = pool[k];
            I K = r.k + a.y + a.d;
            bool merged = false;
            for (auto& q : inc) if (q.k == K && q.o == r.o) { q.m |= r.m; merged = true; break; }
            if (!merged) inc.push_back({K, r.o, i, r.ov, r.m});
            incU |= r.m;
        }
        cnt[C_PUBLISHED] += double(inc.size());
        spend(double(scnt[i]) * double(scnt[i]));
        for (int x = R.er[i]; x <= n; x += x & -x) {
            if (node[x].empty()) node[x].reserve(4);
            merge(node[x], nodeU[x], inc.data(), int(inc.size()), incU);
            cnt[C_MAXNODE] = std::max(cnt[C_MAXNODE], double(node[x].size()));
        }
    }

    void run() {
        for (int e : R.events) {
            check();
            if (e & 1) end(e >> 1); else start(e >> 1);
        }
    }

    // Per-state values and packed pointers (state indices; -1 restart).
    void expand(const int* offsets, I* C, I* back) const {
        for (int i = 0; i < n; ++i) {
            int vi = R.A[i].v;
            const U64* row = G.mask(vi);
            for (int k = soff[i]; k < soff[i] + scnt[i]; ++k) {
                const auto& r = pool[k];
                const I* pi = &G.phi[size_t(vi) * H]; const I* po = &G.phi[size_t(r.ov) * H];
                I fixed = -1;
                if (r.kind == REMOTE) fixed = offsets[r.j] + G.rank(R.A[r.j].v, r.g);
                else if (r.kind == SAME) fixed = offsets[i] + G.rank(vi, r.g);
                int vj = r.kind == LOCAL ? R.A[r.j].v : 0;
                for (int q = 0; q < W; ++q) {
                    int base = offsets[i] + G.wrank[size_t(vi) * W + q];
                    for (U64 bits = r.m.w[q]; bits; bits &= bits - 1) {
                        int bit = std::countr_zero(bits), h = q * 64 + bit;
                        size_t st = size_t(base + std::popcount(row[q] & ((U64(1) << bit) - 1)));
                        C[st] = r.k - (pi[h] - po[h]);
                        back[st] = r.kind == LOCAL ? I(offsets[r.j] + G.rank(vj, h)) : fixed;
                    }
                }
            }
        }
    }
};

struct Compact {
    int W = 1;
    std::vector<int> off, origin, kind, wj, wg, ww;
    std::vector<I> value;
    std::vector<U64> masks;
};

template<int W, bool Budget = false> void class_solve(const Graph& G, const Read& R, I gamma, I* B, I* Rout, double* stats,
                                 Compact* out, const int* offsets, I* C, I* back, bool warm,
                                 double budget = std::numeric_limits<double>::infinity()) {
    auto t = Clock::now();
    PairCache* shared = nullptr;
    if (warm) {
        if (!G.warm) { auto c = new PairCache(); c->reset(1 << 14); const_cast<Graph&>(G).warm = c; }
        shared = static_cast<PairCache*>(G.warm);
    }
    ClassEngine<W, Budget> E(G, R, gamma, shared);
    E.budget = budget;
    stats[2] = elapsed(t); t = Clock::now();
    try { E.run(); } catch (const BudgetExceeded&) { stats[36] = E.work; throw; }
    stats[36] = E.work;
    stats[3] = elapsed(t); t = Clock::now();
    for (int i = 0; i < R.n; ++i) { B[i] = E.B[i]; Rout[i] = E.Rv[i]; }
    if (out) {
        int n = R.n;
        out->W = W; out->off.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) out->off[i + 1] = out->off[i] + E.scnt[i];
        size_t total = size_t(out->off[n]);
        out->value.resize(total); out->origin.resize(total); out->kind.resize(total);
        out->wj.resize(total); out->wg.resize(total); out->ww.resize(total); out->masks.resize(total * W);
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < E.scnt[i]; ++k) {
                const auto& r = E.pool[E.soff[i] + k];
                size_t z = size_t(out->off[i] + k);
                out->value[z] = r.k; out->origin[z] = r.o; out->kind[z] = r.kind;
                out->wj[z] = r.j; out->wg[z] = r.g; out->ww[z] = r.w;
                for (int q = 0; q < W; ++q) out->masks[z * W + q] = r.m.w[q];
            }
    }
    if (C) E.expand(offsets, C, back);
    stats[4] = elapsed(t);
    for (int k = 0; k < C_COUNT; ++k) stats[8 + k] = E.cnt[k];
}

// ---------------------------------------------------------------- per-haplotype
// dense_cut (in 1/256 units): an anchor with deg labels uses dense all-H lanes iff
// 256*deg >= dense_cut*H, and otherwise loops over its own labels only. 0 = always dense.
void ph_solve(const Graph& G, const Read& R, I gamma, const int* offsets, I* C, I* back, I* B, I* Rout, double* stats,
              int dense_cut) {
    auto t = Clock::now();
    int n = R.n, H = G.H;
    require(n < (1 << SH), "per-haplotype lanes support fewer than 2^24 anchors");
    require(R.maxkey < (I(1) << (62 - SH)), "per-haplotype lanes need published keys below 2^38");
    require(dense_cut >= 0 && dense_cut <= 257, "dense_cut must lie in [0, 257]");
    std::vector<I> tree(size_t(n + 1) * H, -1);  // -1: empty lane
    std::vector<I> acc(H, -1), kv(H, -1), racc, keys;
    stats[2] = elapsed(t); t = Clock::now();
    double steps = 0, queries = 0, sparse_anchors = 0;
    for (int e : R.events) {
        int i = e >> 1; const Anchor& a = R.A[i]; int vi = a.v;
        int deg = G.degree(vi); const int* labs = &G.lab[G.loff[vi]];
        bool dense = 256 * I(deg) >= I(dense_cut) * H;
        if (!(e & 1)) {
            I alpha = a.w - a.x - a.c + 2;
            // Remote: one prefix walk per last2reach group, reading the group's labels per row.
            bool have = false; I rv = 0, rptr = -1; int rgbest = -1;
            for (int k0 = G.roff[vi]; k0 < G.roff[vi + 1];) {
                int w = G.rw[k0], k1 = k0;
                while (k1 < G.roff[vi + 1] && G.rw[k1] == w) ++k1;
                int m = k1 - k0;
                if (int(racc.size()) < m) racc.resize(m);
                std::fill(racc.begin(), racc.begin() + m, I(-1));
                const int* gs = &G.rg[k0];
                for (int x = R.tau(w); x > 0; x -= x & -x) {
                    const I* row = &tree[size_t(x) * H];
                    for (int k = 0; k < m; ++k) racc[k] = std::max(racc[k], row[gs[k]]);
                    steps += m;
                }
                queries += m;
                for (int k = 0; k < m; ++k) {
                    I best = racc[k];
                    if (best < 0) continue;
                    int g = gs[k];
                    I cand = (best >> SH) + alpha - G.rcost[k0 + k] - gamma;
                    if (!have || cand > rv || (cand == rv && g < rgbest)) {
                        int j = int(JM - (best & JM));
                        have = true; rv = cand; rgbest = g; rptr = offsets[j] + G.rank(R.A[j].v, g);
                    }
                }
                k0 = k1;
            }
            Rout[i] = have ? rv : NEG;
            I f = a.w; I fptr = -1;
            if (have && rv > f) { f = rv; fptr = rptr; }
            // Local: prefix maxima at thr(i) over all lanes (dense) or the anchor's labels.
            if (dense) {
                std::fill(acc.begin(), acc.end(), I(-1));
                for (int x = R.thr[i]; x > 0; x -= x & -x) {
                    const I* row = &tree[size_t(x) * H];
                    for (int h = 0; h < H; ++h) acc[h] = std::max(acc[h], row[h]);
                    steps += deg;
                }
            } else {
                sparse_anchors += 1;
                for (int k = 0; k < deg; ++k) acc[labs[k]] = -1;
                for (int x = R.thr[i]; x > 0; x -= x & -x) {
                    const I* row = &tree[size_t(x) * H];
                    for (int k = 0; k < deg; ++k) { int h = labs[k]; acc[h] = std::max(acc[h], row[h]); }
                    steps += deg;
                }
            }
            const I* phi = &G.phi[size_t(vi) * H];
            I* Ci = C + offsets[i]; I* bi = back + offsets[i];
            I Bi = NEG; int bk = -1;
            for (int k = 0; k < deg; ++k) {
                int h = labs[k]; I P = acc[h];
                I val = f, ptr = fptr;
                if (P >= 0) {
                    I local = (P >> SH) + alpha - phi[h];
                    if (local > f) { val = local; int j = int(JM - (P & JM)); ptr = offsets[j] + G.rank(R.A[j].v, h); }
                }
                Ci[k] = val; bi[k] = ptr;
                if (val > Bi) { Bi = val; bk = k; }
            }
            B[i] = Bi;
            I s = Bi - gamma, sptr = offsets[i] + bk;
            for (int k = 0; k < deg; ++k) if (s > Ci[k]) { Ci[k] = s; bi[k] = sptr; }
        } else {
            const I* phi = &G.phi[size_t(vi) * H];
            const I* Ci = C + offsets[i];
            I tag = JM - i;
            if (dense) {
                for (int k = 0; k < deg; ++k) { int h = labs[k]; kv[h] = ((Ci[k] + a.y + a.d + phi[h]) << SH) | tag; }
                for (int x = R.er[i]; x <= n; x += x & -x) {
                    I* row = &tree[size_t(x) * H];
                    for (int h = 0; h < H; ++h) row[h] = std::max(row[h], kv[h]);
                    steps += deg;
                }
                for (int k = 0; k < deg; ++k) kv[labs[k]] = -1;
            } else {
                if (int(keys.size()) < deg) keys.resize(deg);
                for (int k = 0; k < deg; ++k) keys[k] = ((Ci[k] + a.y + a.d + phi[labs[k]]) << SH) | tag;
                for (int x = R.er[i]; x <= n; x += x & -x) {
                    I* row = &tree[size_t(x) * H];
                    for (int k = 0; k < deg; ++k) { int h = labs[k]; row[h] = std::max(row[h], keys[k]); }
                    steps += deg;
                }
            }
        }
    }
    stats[3] = elapsed(t);
    stats[8 + C_COUNT] = steps; stats[9 + C_COUNT] = queries; stats[10 + C_COUNT] = sparse_anchors;
}
// Lane operations ph_solve performs on this read: Fenwick rows times the lanes read or
// written (all H lanes for dense anchors, the anchor's labels otherwise), for local queries,
// updates and the grouped remote queries. The base of the auto mode's budget.
double ph_work_estimate(const Graph& G, const Read& R, int dense_cut) {
    int n = R.n, H = G.H; double est = 0;
    for (int i = 0; i < n; ++i) {
        int vi = R.A[i].v, deg = G.degree(vi);
        double lanes = 256 * I(deg) >= I(dense_cut) * H ? double(H) : double(deg);
        for (int k0 = G.roff[vi]; k0 < G.roff[vi + 1];) {
            int w = G.rw[k0], k1 = k0;
            while (k1 < G.roff[vi + 1] && G.rw[k1] == w) ++k1;
            est += double(k1 - k0) * std::popcount(unsigned(R.tau(w)));
            k0 = k1;
        }
        est += lanes * std::popcount(unsigned(R.thr[i]));
        int up = 0;
        for (int x = R.er[i]; x <= n; x += x & -x) ++up;
        est += lanes * up;
    }
    return est;
}
}  // namespace

extern "C" {
const char* gap_error() { return last_error.c_str(); }

// stats[5]: build seconds, S, b, groups, remote entries.
void* gap_graph_new(int V, int H, int W, const I* lengths, const U64* vmask, const I* phi, const int* last,
                    const I* dist, double* stats) {
    try { return build_graph(V, H, W, lengths, vmask, phi, last, dist, stats); }
    catch (const std::exception& e) { last_error = e.what(); return nullptr; }
}
void gap_graph_free(void* g) {
    auto G = static_cast<Graph*>(g);
    if (G) delete static_cast<PairCache*>(G->warm);
    delete G;
}
// Warm pair cache: entries and route-plane words held; reset discards it.
void gap_warm_info(void* g, double* out) {
    auto c = static_cast<PairCache*>(static_cast<Graph*>(g)->warm);
    out[0] = c ? double(c->count) : 0; out[1] = c ? double(c->pool.size()) : 0;
}
void gap_warm_reset(void* g) {
    auto G = static_cast<Graph*>(g);
    delete static_cast<PairCache*>(G->warm); G->warm = nullptr;
}

// stats has 40 slots: [0] kernel, [1] read setup, [2] engine init, [3] sweep, [4] export
// (compact copy or expansion), [5] K, [6] n, [8..8+C_COUNT) class counters; the
// per-haplotype engine uses [8+C_COUNT] Fenwick label-steps, [9+C_COUNT] remote label
// queries and [10+C_COUNT] anchors that took the label-sparse path.
void class_dispatch(const Graph& G, int n, const I* anchors, I gamma, I* B, I* R, double* stats,
                    Compact* out, const int* offsets, I* C, I* back, bool warm) {
    auto total = Clock::now();
    std::fill(stats, stats + 40, 0.0);
    Read rd; prepare_read(rd, G, n, anchors, gamma, stats);
    if (out) { out->W = G.W; out->off.assign(1, 0); }
    if (n) switch (G.W) {
        case 1: class_solve<1>(G, rd, gamma, B, R, stats, out, offsets, C, back, warm); break;
        case 2: class_solve<2>(G, rd, gamma, B, R, stats, out, offsets, C, back, warm); break;
        case 4: class_solve<4>(G, rd, gamma, B, R, stats, out, offsets, C, back, warm); break;
        case 8: class_solve<8>(G, rd, gamma, B, R, stats, out, offsets, C, back, warm); break;
        case 16: class_solve<16>(G, rd, gamma, B, R, stats, out, offsets, C, back, warm); break;
        default: require(false, "unsupported mask width");
    }
    stats[6] = n; stats[0] = elapsed(total);
}
// Compact output: returns a handle (nullptr on error); B and R are written.
void* gap_class_compact(void* graph, int n, const I* anchors, I gamma, I* B, I* R, double* stats, int warm) {
    Compact* out = nullptr;
    try {
        out = new Compact();
        class_dispatch(*static_cast<Graph*>(graph), n, anchors, gamma, B, R, stats, out, nullptr, nullptr, nullptr, warm != 0);
        return out;
    } catch (const std::exception& e) { delete out; last_error = e.what(); return nullptr; }
}
// Expanded output: packed per-state values and pointers written in the same call.
int gap_class_packed(void* graph, int n, const I* anchors, I gamma, const int* offsets, I* C, I* back, I* B, I* R,
                     double* stats, int warm) {
    try {
        class_dispatch(*static_cast<Graph*>(graph), n, anchors, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0);
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); return -1; }
}
long long gap_compact_size(void* h) { return h ? (long long)static_cast<Compact*>(h)->value.size() : -1; }
void gap_compact_copy(void* h, int* off, I* value, int* origin, int* kind, int* wj, int* wg, int* ww, U64* masks) {
    auto c = static_cast<Compact*>(h);
    std::copy(c->off.begin(), c->off.end(), off);
    std::copy(c->value.begin(), c->value.end(), value); std::copy(c->origin.begin(), c->origin.end(), origin);
    std::copy(c->kind.begin(), c->kind.end(), kind); std::copy(c->wj.begin(), c->wj.end(), wj);
    std::copy(c->wg.begin(), c->wg.end(), wg); std::copy(c->ww.begin(), c->ww.end(), ww);
    std::copy(c->masks.begin(), c->masks.end(), masks);
}
void gap_compact_free(void* h) { delete static_cast<Compact*>(h); }

// Expand compact records (arrays as copied by gap_compact_copy) into packed states.
int gap_expand(void* graph, int n, const int* vertex, const int* offsets, const int* off, const I* value,
               const int* origin, const int* kind, const int* wj, const int* wg, const U64* masks, I* C, I* back) {
    try {
        const Graph& G = *static_cast<Graph*>(graph);
        int W = G.W, H = G.H;
        for (int i = 0; i < n; ++i) {
            int vi = vertex[i];
            require(vi >= 0 && vi < G.V, "vertex out of bounds");
            const U64* row = G.mask(vi);
            const I* pi = &G.phi[size_t(vi) * H];
            for (int z = off[i]; z < off[i + 1]; ++z) {
                int o = origin[z], kd = kind[z], j = wj[z], g = wg[z];
                require(o >= 0 && o < n && kd >= 0 && kd <= 3, "invalid compact record");
                require(kd == RESTART || (j >= 0 && j < n && g >= -1 && g < H), "invalid compact witness");
                require(kd != SAME || (j == i && g >= 0), "invalid same-anchor witness");
                require(kd != REMOTE || g >= 0, "invalid remote witness");
                const I* po = &G.phi[size_t(vertex[o]) * H];
                I fixed = -1;
                if (kd == REMOTE) fixed = offsets[j] + G.rank(vertex[j], g);
                else if (kd == SAME) fixed = offsets[i] + G.rank(vi, g);
                const U64* m = masks + size_t(z) * W;
                for (int q = 0; q < W; ++q) {
                    require(!(m[q] & ~row[q]), "compact mask outside the anchor's labels");
                    int base = offsets[i] + G.wrank[size_t(vi) * W + q];
                    for (U64 bits = m[q]; bits; bits &= bits - 1) {
                        int bit = std::countr_zero(bits), h = q * 64 + bit;
                        size_t st = size_t(base + std::popcount(row[q] & ((U64(1) << bit) - 1)));
                        C[st] = value[z] - (pi[h] - po[h]);
                        back[st] = kd == LOCAL ? I(offsets[j] + G.rank(vertex[j], h)) : fixed;
                    }
                }
            }
        }
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); return -1; }
}

// Auto mode, packed output: the class engine under a work budget of
// factor * (the per-haplotype engine's lane operations on this read) + 64*W*n words. If the
// budget is exceeded the read is re-solved by the per-haplotype engine, so the total work is
// at most about (1 + factor) times the per-haplotype engine's, plus one pass that computes the
// estimate. Without the per-haplotype domain (keys < 2^38, n < 2^24) the budget is unlimited.
// stats[7] = 1 if the fallback ran; [36] class work in words (at completion or abort);
// [37] the per-haplotype estimate; [38] the budget; [39] 1 if the fallback was available.
int gap_auto_packed(void* graph, int n, const I* anchors, I gamma, const int* offsets, I* C, I* back, I* B, I* R,
                    double* stats, int warm, double factor, int dense_cut) {
    try {
        auto total = Clock::now();
        std::fill(stats, stats + 40, 0.0);
        require(factor >= 0, "auto factor must be nonnegative");
        require(dense_cut >= 0 && dense_cut <= 257, "dense_cut must lie in [0, 257]");
        const Graph& G = *static_cast<Graph*>(graph);
        Read rd; prepare_read(rd, G, n, anchors, gamma, stats);
        bool ph_ok = n < (1 << SH) && rd.maxkey < (I(1) << (62 - SH));
        double est = n ? ph_work_estimate(G, rd, dense_cut) : 0.0;
        double budget = ph_ok ? factor * est + 64.0 * G.W * double(n) : std::numeric_limits<double>::infinity();
        bool fell_back = false;
        if (n) {
            try {
                switch (G.W) {
                    case 1: class_solve<1, true>(G, rd, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0, budget); break;
                    case 2: class_solve<2, true>(G, rd, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0, budget); break;
                    case 4: class_solve<4, true>(G, rd, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0, budget); break;
                    case 8: class_solve<8, true>(G, rd, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0, budget); break;
                    case 16: class_solve<16, true>(G, rd, gamma, B, R, stats, nullptr, offsets, C, back, warm != 0, budget); break;
                    default: require(false, "unsupported mask width");
                }
            } catch (const BudgetExceeded&) {
                fell_back = true;
                ph_solve(G, rd, gamma, offsets, C, back, B, R, stats, dense_cut);
            }
        }
        stats[7] = fell_back ? 1.0 : 0.0;
        stats[37] = est; stats[38] = budget; stats[39] = ph_ok ? 1.0 : 0.0;
        stats[6] = n; stats[0] = elapsed(total);
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); return -1; }
}

int gap_ph_solve(void* graph, int n, const I* anchors, I gamma, const int* offsets, I* C, I* back, I* B, I* R,
                 double* stats, int dense_cut) {
    try {
        auto total = Clock::now();
        std::fill(stats, stats + 40, 0.0);
        const Graph& G = *static_cast<Graph*>(graph);
        Read rd; prepare_read(rd, G, n, anchors, gamma, stats);
        if (n) ph_solve(G, rd, gamma, offsets, C, back, B, R, stats, dense_cut);
        stats[6] = n; stats[0] = elapsed(total);
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); return -1; }
}
}
