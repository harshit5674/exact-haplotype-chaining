// Emulation of vg Giraffe v1.77.0 recombination-aware chaining (single-record DP) and the
// exact per-haplotype DP for the objective it approximates. See SPEC_GIRAFFE.md.
//
// Build: clang++ -std=c++20 -O3 -shared -fPIC giraffe_emul.cpp -o <lib>
//
// Graph: vertices 0..V-1 already in topological order, lengths, directed edges, and one
// W-word haplotype mask per vertex. D(u,v) = shortest start-of-u to start-of-v distance
// (sum of lengths of u and every intermediate vertex on a shortest path; D(u,u) = 0).
//
// Anchors must arrive sorted by (c ascending, d descending, tie order chosen by the caller);
// the array index is the "source" number used by Giraffe's tie rule (SPEC section 5).
// Anchor tuple: v, x, y, c, d (1-based closed), weight (ignored here; item points are given).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int32_t NEG = -(1 << 29);
constexpr int64_t NOWHERE = std::numeric_limits<int64_t>::max();
thread_local std::string g_error;

struct Graph {
    int V = 0, H = 0, W = 0;
    std::vector<int64_t> len;
    std::vector<int32_t> dist;   // V*V, -1 unreachable
    std::vector<uint64_t> vmask; // V*W
};

struct Params {
    int64_t graph_lookback;   // keep iff Z <= graph_lookback (and Z > 0)
    int64_t read_lookback;    // keep iff read gap <= read_lookback
    int64_t indel_limit;      // keep iff |read gap - graph gap| <= indel_limit
    int32_t base_seed_length; // L in score_chain_gap
    int32_t rho;              // recombination penalty
    double gap_scale;
    int32_t n_bonus;          // number of heuristic variants (<= 4)
    int32_t bonus[4];         // consistency bonus of each variant
    int32_t do_exact;         // run the exact DP
    int32_t do_heuristic;     // run the heuristic DP (requires W == 1)
    int32_t prune;            // exact DP: skip sources that cannot improve any label (values unchanged)
};

// vg: score_chain_gap (chain_items.cpp 365-373), then jump0 = (int)(-g * gap_scale) (line 509).
inline int score_chain_gap(uint64_t d, uint64_t L) {
    if (d == 0) return 0;
    return 0.01 * (double)L * (double)d + 0.5 * std::log2((double)d);
}
inline int jump0_of(uint64_t indel, uint64_t L, double gap_scale) {
    int g = score_chain_gap(indel, L);
    int j = -g * gap_scale;  // int * double -> double -> int (truncation), as in vg
    return j;
}

}  // namespace

extern "C" {

const char* ge_error() { return g_error.c_str(); }

void* ge_graph_new(int V, int H, int W, const int64_t* lengths, int E, const int32_t* eu, const int32_t* ev,
                   const uint64_t* vmask) {
    try {
        auto* g = new Graph();
        g->V = V; g->H = H; g->W = W;
        g->len.assign(lengths, lengths + V);
        g->vmask.assign(vmask, vmask + (size_t)V * W);
        std::vector<std::vector<int>> out(V);
        for (int e = 0; e < E; ++e) {
            if (eu[e] < 0 || eu[e] >= V || ev[e] < 0 || ev[e] >= V || eu[e] >= ev[e]) {
                g_error = "edge out of range or not topological"; delete g; return nullptr;
            }
            out[eu[e]].push_back(ev[e]);
        }
        g->dist.assign((size_t)V * V, -1);
        std::vector<int64_t> d(V);
        for (int u = 0; u < V; ++u) {
            std::fill(d.begin(), d.end(), -1);
            d[u] = 0;
            for (int v = u; v < V; ++v) {
                if (d[v] < 0) continue;
                int64_t nd = d[v] + g->len[v];
                for (int w : out[v])
                    if (d[w] < 0 || nd < d[w]) d[w] = nd;
            }
            for (int v = 0; v < V; ++v) {
                if (d[v] > std::numeric_limits<int32_t>::max()) { g_error = "distance overflow"; delete g; return nullptr; }
                g->dist[(size_t)u * V + v] = (int32_t)d[v];
            }
        }
        return g;
    } catch (const std::exception& e) {
        g_error = e.what();
        return nullptr;
    }
}

void ge_graph_free(void* gp) { delete static_cast<Graph*>(gp); }

int32_t ge_graph_dist(void* gp, int u, int v) {
    auto* g = static_cast<Graph*>(gp);
    return g->dist[(size_t)u * g->V + v];
}

// DP state shared by the geometric solver and hand-built examples. Destinations are processed in
// order 0..n-1; process(i, ...) receives every incoming transition of i (sources ascending).
// Outputs (null when the corresponding DP is off):
//   heuristic, per variant b: hs[b*n+i], hsrc[b*n+i] (-1 = nowhere), hmask[b*n+i], hrec[b*n+i]
//   exact: C[i*H+h] (NEG if h not in M_i), Pj[i*H+h] (-1 = fresh start), Pg[i*H+h], B[i]
struct DP {
    int n, H, W;
    const int32_t* item;
    const uint64_t* amask;
    Params p;
    int nb;
    int32_t *hs, *hsrc, *hrec, *C, *Pj, *B;
    uint64_t* hmask;
    int16_t* Pg;
    std::vector<int64_t> src;
    std::vector<int32_t> evb, best, bj, bl, amax;
    double vector_updates = 0;

    DP(int n_, int H_, int W_, const int32_t* item_, const uint64_t* amask_, const Params& p_, int32_t* hs_,
       int32_t* hsrc_, uint64_t* hmask_, int32_t* hrec_, int32_t* C_, int32_t* Pj_, int16_t* Pg_, int32_t* B_)
        : n(n_), H(H_), W(W_), item(item_), amask(amask_), p(p_), hs(hs_), hsrc(hsrc_), hrec(hrec_), C(C_), Pj(Pj_),
          B(B_), hmask(hmask_), Pg(Pg_) {
        if (p.do_heuristic && W != 1) throw std::runtime_error("heuristic needs at most 64 haplotypes");
        if (p.n_bonus < 0 || p.n_bonus > 4) throw std::runtime_error("n_bonus out of range");
        nb = p.do_heuristic ? p.n_bonus : 0;
        src.assign((size_t)nb * n, NOWHERE);
        evb.assign((size_t)nb * n, 0);
        for (int b = 0; b < nb; ++b)
            for (int i = 0; i < n; ++i) {
                hs[(size_t)b * n + i] = item[i];
                hmask[(size_t)b * n + i] = amask[i];
                hrec[(size_t)b * n + i] = 0;
                evb[(size_t)b * n + i] = p.bonus[b];
            }
        best.assign(H, NEG); bj.assign(H, -1); bl.assign(H, 0);
        amax.assign(p.do_exact ? n : 0, -1);
    }

    void process(int i, size_t m, const int32_t* tsrc, const int32_t* tj0) {
        const int32_t it = item[i];
        // ---- heuristic (vg chain_items.cpp 438-582) ----
        for (int b = 0; b < nb; ++b) {
            const int cb = p.bonus[b];
            int32_t* S = hs + (size_t)b * n;
            int64_t* SR = src.data() + (size_t)b * n;
            uint64_t* PM = hmask + (size_t)b * n;
            int32_t* RC = hrec + (size_t)b * n;
            int32_t* EB = evb.data() + (size_t)b * n;
            const uint64_t Mi = amask[i];
            for (size_t t = 0; t < m; ++t) {
                const int j = tsrc[t];
                // from nowhere, re-offered on every transition (lines 462-474)
                {
                    const int64_t eval_now = (int64_t)it + cb, eval_cur = (int64_t)S[i] + EB[i];
                    const bool greater = it > S[i] || (it == S[i] && NOWHERE > SR[i]);
                    if (eval_now > eval_cur || (eval_now == eval_cur && greater)) {
                        S[i] = it; SR[i] = NOWHERE; PM[i] = Mi; RC[i] = 0; EB[i] = cb;
                    }
                }
                // candidate from j (lines 509-548)
                const uint64_t pj = PM[j];
                const bool reset = (pj & Mi) == 0;
                const int jump = tj0[t] - (reset ? p.rho : 0);
                const int32_t cscore = S[j] + jump + it;
                const uint64_t cpaths = reset ? Mi : (pj & Mi);
                const int32_t crec = RC[j] + (reset ? 1 : 0);
                int bonus = 0;
                if (cb > 0) {
                    const int pre = __builtin_popcountll(pj);
                    if (pre > 0 && !reset) bonus = (cb * __builtin_popcountll(cpaths)) / pre;
                }
                const int64_t eval_from = (int64_t)cscore + bonus, eval_best = (int64_t)S[i] + EB[i];
                const bool greater = cscore > S[i] || (cscore == S[i] && (int64_t)j > SR[i]);
                if (eval_from > eval_best || (eval_from == eval_best && greater)) {
                    S[i] = cscore; SR[i] = j; PM[i] = cpaths; RC[i] = crec; EB[i] = bonus;
                }
            }
        }
        // ---- exact per-haplotype DP (SPEC section 7) ----
        // Exact pruning: u_t = B_j + jump0 bounds every label's candidate from source j. After the
        // source with the largest u is applied, every label holds >= u_max - rho (switch term), so a
        // source with u_t <= u_max - rho cannot strictly improve any label and is skipped. Values are
        // unchanged; only which of several tied predecessors is recorded can differ.
        if (p.do_exact) {
            std::fill(best.begin(), best.end(), NEG);
            std::fill(bj.begin(), bj.end(), -1);
            std::fill(bl.begin(), bl.end(), 0);
            int32_t* bestp = best.data();
            int32_t* bjp = bj.data();
            int32_t* blp = bl.data();
            auto apply = [&](size_t t) {
                const int j = tsrc[t];
                const int32_t j0 = tj0[t];
                const int32_t bjr = B[j] - p.rho;
                const int32_t* Cj = C + (size_t)j * H;
                for (int h = 0; h < H; ++h) {
                    const int32_t vj = Cj[h];
                    const int32_t loc = vj >= bjr ? 1 : 0;
                    const int32_t val = (loc ? vj : bjr) + j0;
                    const bool upd = val > bestp[h];
                    bestp[h] = upd ? val : bestp[h];
                    bjp[h] = upd ? j : bjp[h];
                    blp[h] = upd ? loc : blp[h];
                }
                ++vector_updates;
            };
            if (!p.prune) {
                // Original loop (all sources, ascending): the engine used for the Phase 2 outputs.
                for (size_t t = 0; t < m; ++t) apply(t);
            } else if (m > 0) {
                int64_t umax = std::numeric_limits<int64_t>::min();
                size_t tmax = 0;
                for (size_t t = 0; t < m; ++t) {
                    const int64_t u = (int64_t)B[tsrc[t]] + tj0[t];
                    if (u > umax) { umax = u; tmax = t; }
                }
                apply(tmax);
                const int64_t thr = umax - p.rho;
                for (size_t t = 0; t < m; ++t)
                    if (t != tmax && (int64_t)B[tsrc[t]] + tj0[t] > thr) apply(t);
            }
            int32_t bi = NEG;
            int ai = -1;
            int32_t* Ci = C + (size_t)i * H;
            for (int h = 0; h < H; ++h) {
                const size_t k = (size_t)i * H + h;
                const bool member = (amask[(size_t)i * W + (h >> 6)] >> (h & 63)) & 1;
                if (!member) { Ci[h] = NEG; Pj[k] = -1; Pg[k] = -1; continue; }
                if (best[h] > 0) {
                    Ci[h] = it + best[h];
                    Pj[k] = bj[h];
                    Pg[k] = (int16_t)(bl[h] ? h : amax[bj[h]]);
                } else {
                    Ci[h] = it; Pj[k] = -1; Pg[k] = -1;
                }
                if (Ci[h] > bi) { bi = Ci[h]; ai = h; }
            }
            if (ai < 0) throw std::runtime_error("anchor with empty haplotype mask");
            B[i] = bi;
            amax[i] = ai;
        }
    }

    void finish() {
        for (int b = 0; b < nb; ++b)
            for (int i = 0; i < n; ++i) {
                const int64_t s = src[(size_t)b * n + i];
                hsrc[(size_t)b * n + i] = s == NOWHERE ? -1 : (int32_t)s;
            }
    }
};

// Geometric transition generation (SPEC section 3) followed by both DPs, one destination at a time.
// anchors: n*6 int64 (v, x, y, c, d, w), sorted by (c asc, d desc, caller's tie order).
// item[n]: item points. stats[10]: transitions, pairs scanned (d_j < c_i), rejected: graph lookback,
// graph gap < 0, read lookback, indel, unreachable; max incoming transitions of one destination;
// exact-DP source updates applied (after pruning); 0.
int ge_solve(void* gp, int n, const int64_t* anchors, const int32_t* item, const Params* pp, int32_t* hs,
             int32_t* hsrc, uint64_t* hmask, int32_t* hrec, int32_t* C, int32_t* Pj, int16_t* Pg, int32_t* B,
             double* stats) {
    try {
        auto* g = static_cast<Graph*>(gp);
        const Params& p = *pp;
        const int V = g->V, W = g->W, H = g->H;
        std::vector<uint64_t> amask((size_t)n * W);
        for (int i = 0; i < n; ++i) {
            const int64_t* a = anchors + (size_t)i * 6;
            if (a[0] < 0 || a[0] >= V || a[1] < 1 || a[2] < a[1] || a[2] > g->len[a[0]] || a[3] < 1 || a[4] < a[3]) {
                g_error = "bad anchor"; return 1;
            }
            if (i > 0) {
                const int64_t* b = anchors + (size_t)(i - 1) * 6;
                if (b[3] > a[3] || (b[3] == a[3] && b[4] < a[4])) { g_error = "anchors not sorted by (c asc, d desc)"; return 1; }
            }
            for (int w = 0; w < W; ++w) amask[(size_t)i * W + w] = g->vmask[(size_t)a[0] * W + w];
        }
        DP dp(n, H, W, item, amask.data(), p, hs, hsrc, hmask, hrec, C, Pj, Pg, B);
        std::vector<int32_t> tsrc, tj0;
        double total = 0, scanned = 0, rej_look = 0, rej_gap = 0, rej_read = 0, rej_indel = 0, rej_unreach = 0, mx = 0;
        const uint64_t L = (uint64_t)p.base_seed_length;
        for (int i = 0; i < n; ++i) {
            const int64_t* a = anchors + (size_t)i * 6;
            const int64_t vi = a[0], xi = a[1], ci = a[3];
            tsrc.clear(); tj0.clear();
            for (int j = 0; j < i; ++j) {
                const int64_t* s = anchors + (size_t)j * 6;
                const int64_t vj = s[0], xj = s[1], yj = s[2], dj = s[4];
                if (dj >= ci) continue;  // read: source must end before destination starts (d_j < c_i)
                scanned += 1;
                const int32_t D = g->dist[(size_t)vj * V + vi];
                if (D < 0) { rej_unreach += 1; continue; }
                const int64_t Z = (int64_t)D + xi - xj;          // zip-tree distance between seed starts
                if (Z <= 0 || Z > p.graph_lookback) { rej_look += 1; continue; }
                const int64_t graph_gap = Z - (yj - xj + 1);      // remove source length (end_hint_offset)
                if (graph_gap < 0) { rej_gap += 1; continue; }
                const int64_t read_gap = ci - 1 - dj;             // to.read_start - from.read_end (0-based)
                if (read_gap > p.read_lookback) { rej_read += 1; continue; }
                const int64_t indel = read_gap > graph_gap ? read_gap - graph_gap : graph_gap - read_gap;
                if (indel > p.indel_limit) { rej_indel += 1; continue; }
                tsrc.push_back(j);
                tj0.push_back(jump0_of((uint64_t)indel, L, p.gap_scale));
            }
            total += (double)tsrc.size();
            mx = std::max(mx, (double)tsrc.size());
            dp.process(i, tsrc.size(), tsrc.data(), tj0.data());
        }
        dp.finish();
        if (stats) {
            stats[0] = total; stats[1] = scanned; stats[2] = rej_look; stats[3] = rej_gap;
            stats[4] = rej_read; stats[5] = rej_indel; stats[6] = rej_unreach; stats[7] = mx;
            stats[8] = dp.vector_updates; stats[9] = 0;
        }
        return 0;
    } catch (const std::exception& e) {
        g_error = e.what();
        return 1;
    }
}

// Same DPs on an explicit transition list (hand-built examples, brute-force checks).
// trans_off[n+1]; trans_src / trans_jump0 list the incoming transitions of each destination.
int ge_solve_explicit(int n, int H, int W, const int32_t* item, const uint64_t* amask, const int64_t* trans_off,
                      const int32_t* trans_src, const int32_t* trans_jump0, const Params* pp, int32_t* hs,
                      int32_t* hsrc, uint64_t* hmask, int32_t* hrec, int32_t* C, int32_t* Pj, int16_t* Pg,
                      int32_t* B) {
    try {
        DP dp(n, H, W, item, amask, *pp, hs, hsrc, hmask, hrec, C, Pj, Pg, B);
        for (int i = 0; i < n; ++i) {
            for (int64_t t = trans_off[i]; t < trans_off[i + 1]; ++t)
                if (trans_src[t] < 0 || trans_src[t] >= i) { g_error = "explicit transition must go forward"; return 1; }
            dp.process(i, (size_t)(trans_off[i + 1] - trans_off[i]), trans_src + trans_off[i], trans_jump0 + trans_off[i]);
        }
        dp.finish();
        return 0;
    } catch (const std::exception& e) {
        g_error = e.what();
        return 1;
    }
}

// Transition list of an instance (for independent checks): returns count; fills up to cap entries.
int64_t ge_transitions(void* gp, int n, const int64_t* anchors, const Params* pp, int64_t cap, int32_t* out_src,
                       int32_t* out_dst, int32_t* out_jump0, int32_t* out_indel) {
    auto* g = static_cast<Graph*>(gp);
    const Params& p = *pp;
    const int V = g->V;
    const uint64_t L = (uint64_t)p.base_seed_length;
    int64_t k = 0;
    for (int i = 0; i < n; ++i) {
        const int64_t* a = anchors + (size_t)i * 6;
        for (int j = 0; j < i; ++j) {
            const int64_t* s = anchors + (size_t)j * 6;
            if (s[4] >= a[3]) continue;
            const int32_t D = g->dist[(size_t)s[0] * V + a[0]];
            if (D < 0) continue;
            const int64_t Z = (int64_t)D + a[1] - s[1];
            if (Z <= 0 || Z > p.graph_lookback) continue;
            const int64_t gg = Z - (s[2] - s[1] + 1);
            if (gg < 0) continue;
            const int64_t rg = a[3] - 1 - s[4];
            if (rg > p.read_lookback) continue;
            const int64_t indel = rg > gg ? rg - gg : gg - rg;
            if (indel > p.indel_limit) continue;
            if (k < cap) {
                out_src[k] = j; out_dst[k] = i;
                out_jump0[k] = jump0_of((uint64_t)indel, L, p.gap_scale);
                out_indel[k] = (int32_t)indel;
            }
            ++k;
        }
    }
    return k;
}

}  // extern "C"
