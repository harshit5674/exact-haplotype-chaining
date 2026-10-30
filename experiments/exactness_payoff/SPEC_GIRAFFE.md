# vg Giraffe recombination-aware chaining: specification for emulation

Status: **specification from source reading** (2026-10-06). Nothing here was obtained by running vg; vg does not run on this macOS arm64 host. Every rule below cites the source line that implements it. The emulation built from this specification is in [`giraffe_emul.cpp`](giraffe_emul.cpp); its faithfulness limits are listed in section 8.

## 1. Sources read

| Component | Version read | Pin |
| --- | --- | --- |
| vg | release **v1.77.0** (2026-09-09), the latest release on 2026-10-06 | commit `02244e42aac24a8375d71e7c3c7f07817921eab1` |
| gbwtgraph | the submodule commit pinned by vg v1.77.0 (`deps/gbwtgraph`) | commit `7fac0af212b0502d40cca3d2b2c3a5fd7d85c540` |
| vg, for history | v1.69.0, v1.70.0, v1.71.0 and v1.75.0 copies of `chain_items.cpp` | tags |

`src/algorithms/chain_items.cpp` on vg master has no commit after 2026-08-24 (GitHub commits API, read 2026-10-06), so master equals v1.77.0 for this file. The v1.77.0 release note says its chaining changes ("memory & speed optimizations") have "no effect on output"; the v1.75.0 DP code (lines 559–690) is the same rule, with the indel size computed inside the DP instead of during transition generation.

Short names used below (all at tag v1.77.0 unless stated):
- **CI** = [src/algorithms/chain_items.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/algorithms/chain_items.cpp) (969 lines)
- **CH** = [src/algorithms/chain_items.hpp](https://github.com/vgteam/vg/blob/v1.77.0/src/algorithms/chain_items.hpp) (613 lines)
- **MFC** = [src/minimizer_mapper_from_chains.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/minimizer_mapper_from_chains.cpp)
- **MM** = [src/minimizer_mapper.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/minimizer_mapper.cpp), **MMH** = [src/minimizer_mapper.hpp](https://github.com/vgteam/vg/blob/v1.77.0/src/minimizer_mapper.hpp)
- **GM** = [src/subcommand/giraffe_main.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/subcommand/giraffe_main.cpp)
- **ZT** = [src/zip_code_tree.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/zip_code_tree.cpp), **ZTH** = [src/zip_code_tree.hpp](https://github.com/vgteam/vg/blob/v1.77.0/src/zip_code_tree.hpp)
- **AS** = [src/alignment_scorer.hpp](https://github.com/vgteam/vg/blob/v1.77.0/src/alignment_scorer.hpp), **IR** = [src/index_registry.cpp](https://github.com/vgteam/vg/blob/v1.77.0/src/index_registry.cpp)
- **GI** = gbwtgraph [src/internal.cpp](https://github.com/jltsiren/gbwtgraph/blob/7fac0af212b0502d40cca3d2b2c3a5fd7d85c540/src/internal.cpp), **GX** = gbwtgraph [src/index.cpp](https://github.com/jltsiren/gbwtgraph/blob/7fac0af212b0502d40cca3d2b2c3a5fd7d85c540/src/index.cpp), **GXH** = gbwtgraph [include/gbwtgraph/index.h](https://github.com/jltsiren/gbwtgraph/blob/7fac0af212b0502d40cca3d2b2c3a5fd7d85c540/include/gbwtgraph/index.h)

A citation "CI 48–57" means lines 48 to 57 of CI, e.g. `https://github.com/vgteam/vg/blob/v1.77.0/src/algorithms/chain_items.cpp#L48-L57`.

## 2. Haplotype masks and path IDs

- **Mask type.** `path_flags_t = uint64_t`, `MAX_PATHS = 64` (CH 41–42). An anchor stores `start_paths` and `end_paths`; `add_path` crashes unless the ID is below 64 (CH 161–212, crash check CH 186).
- **Path IDs.** `PathIdMap` sorts all GBWT paths by (sample, phase, contig, count) and assigns consecutive IDs to distinct keys, trying the key (sample, contig, haplotype), then (sample, haplotype), then (sample). If a key type needs more than `MAX_HAPLOTYPES` = 64 IDs it tries the next; if all fail the map is empty (GI 263–322). With an empty map "all paths map to 0" (GXH 51–61). Generic and reference paths are counted (they are GBWT paths).
- **Mask of a minimizer hit.** For each minimizer occurrence found on a haplotype window, the k-mer's node path is extracted (GI 344–371) and its GBWT search state is located with an r-index; the mask is the OR of `1 << id(path)` over all located sequences, cached per search state (GX 192–211). Hits colliding on the same (key, position) from different k-mer paths are OR-merged (GX 214–226). For a k-mer inside one node the mask is therefore exactly the set of haplotype IDs whose paths visit that node.
- **Seeds and anchors.** With path payloads, each seed's `paths` is the payload's last word (MM 4475–4477); otherwise it stays 0 (`Seed::paths = 0`, [snarl_seed_clusterer.hpp](https://github.com/vgteam/vg/blob/v1.77.0/src/snarl_seed_clusterer.hpp) line 79). `to_anchor` copies `seed.paths` into both `start_paths` and `end_paths` (MFC 3990, 4037; Anchor constructor CH 230–243).
- **When the mode is on.** `default_rec_mode` = index has `PAYLOAD_ZIPCODES_WITH_PATHS` and reads are long (GM 1831); `--rec-mode`/`--no-rec-mode` override (GM 1398–1403); short reads are refused (GM 1580–1583); `use_payload_paths &= rec_mode` (GM 2084).

For the project's panels: `hprc60` and the pilot's MHC graph each have 60 haplotypes with distinct sample names (`HG02148#1`, ...), so every key type gives 60 ≤ 64 IDs: one bit per haplotype. `hprc463` has 463 distinct sample names, so all three key types exceed 64 and the map is empty; the recombination-aware mode cannot represent it.

## 3. Anchors, items and the transition set

**Anchor from a seed** (MFC 3978–4038). The anchor is the part of the minimizer match on the node where its anchoring end falls; the rest of the minimizer becomes `margin_before`/`margin_after`. Its score is `score_exact_match` over the whole minimizer (MFC 4036), i.e. `match × k` with the default `match = 1` (AS 18). `base_seed_length() = margin_before + length + margin_after` (CH 240). Long-read indexes use minimizers with k = 31, w = 50 (IR 121–122). Gapless extension is off for hifi (`gapless-extension-limit 0`, GM 972; used at MFC 1219, 1590), so the chained items are the seed anchors.

**Item points** = `score() + item_bonus` (CI 434, 453).

**Order.** Anchors of one zip-code tree are sorted by read start ascending and read end descending with `std::sort` (CI 96–106, called at MFC 1596). Ties (identical read interval) have unspecified order. The **index in this sorted view** is the anchor's `source` number used in tie-breaking (section 5).

**Transition enumeration** (CI 108–268). For every destination seed, walking the zip-code tree left to right, a distance iterator walks right-to-left and yields every source seed whose running graph distance `Z` satisfies `0 < Z ≤ graph_lookback_limit`:
- a seed is yielded only if the running distance is positive and its orientation matches (ZT 2027–2049);
- exceeding the limit halts the scan on the top-level chain or skips the current chain inside a snarl (ZT 2090–2109);
- inside a DAG snarl the distances to the chains on the left are taken from the snarl's stored minimum-distance matrix (ZT 1690–1707, 1866–1949), and leaving a chain inside a snarl does not continue into other chains (ZT 2079–2089), so each source is reported once, at its minimum distance;
- the tree's edge values are "the number of nucleotides on the minimum-length path in the variation graph" and seeds stand for the first nucleotide of the match; two seeds at the same position are at distance 0 (ZTH 72–86).

Distances are therefore **haplotype-agnostic minimum distances** between seed positions. Each yielded (source, destination) pair is passed to `add_transition_if_legal` (CI 270–355), which keeps it iff
1. the source ends no later than the destination starts in the read, `read_distance = to.read_start − from.read_end ≥ 0` (CI 282–289, 960–965);
2. `read_distance ≤ max_read_lookback` (CI 291–298);
3. the exclusion zones do not overlap, `from.read_exclusion_end ≤ to.read_exclusion_start` (CI 300–307);
4. `graph_distance = Z − (to.start_hint_offset + from.end_hint_offset) ≥ 0` (CI 319–336); for a forward seed anchor `start_hint_offset = 0` and `end_hint_offset = length` (CH 239), so `graph_distance` is the number of graph bases strictly between the two anchors;
5. `indel_size = |read_distance − graph_distance| ≤ max_indel_bases` (CI 343–352).

The surviving transitions are sorted by destination read start (CI 143–148) and fed to the DP (CI 150–152). Same-orientation pairs only: mixed orientations crash (CI 258–262); reverse-orientation pairs are flipped so the DP always runs along the read (CI 236–257).

**Limits** (MFC 1624–1628):
`graph_lookback_limit = max(max_graph_lookback_bases, (size_t)(max_graph_lookback_bases_per_base × read_length))`, and the same form for the read lookback and the indel limit.

| Parameter | Library default (MMH) | hifi preset (GM 957–1008) |
| --- | --- | --- |
| `max-graph-lookback-bases` / per base | 3000 / 0.3 (MMH 271–275) | **20000 / 0.10501002120802233** (GM 980–981) |
| `max-read-lookback-bases` / per base | `SIZE_MAX` / 1.0 (MMH 277–281) | not set: no read limit |
| `max-indel-bases` / per base | 2000 / 0.2 (MMH 297–301) | **5000 / 2.45** (GM 982–983) |
| `item-bonus` | 0 (MMH 284) | **2** (GM 984) |
| `gap-scale` | 1.0 (MMH 288) | **0.27579** (GM 985) |
| `rec-penalty` ρ | 0 (MMH 291) | **2** (GM 986) |
| `rec-consistency-bonus` | 0 (MMH 294; −1 means "use ρ", MFC 1644) | **12** (GM 987) |
| `max-chains-per-tree` | 1 (MMH 314) | **3** (GM 991) |
| `rec-penalty-aln` | 0 (MMH 378) | 28, alignment stage only (GM 1007; MFC 901–916) |
| `zipcode-tree-scale` | 2.0 (MMH 227) | 2.0 (GM 977) |

For the project's reads (17–32 kb), the hifi graph lookback is 20000 bp, the read lookback is unlimited and the indel limit is 2.45 × read length (≥ 41 kb).

## 4. Transition score

For a transition `j → i` with `indel = indel_size` (CI 509–519):

```
score_chain_gap(d, L) = 0                                   if d = 0
                      = (int)(0.01 * L * d + 0.5 * log2(d)) otherwise      (CI 365–373; double arithmetic, truncation)
jump0(j,i)  = (int)( -score_chain_gap(indel, L) * gap_scale )             (CI 509; int × double, truncation toward 0)
jump(j,i)   = jump0(j,i) - rho * [ record_j.paths & M_i == 0 ]            (CI 512, 377–383)
candidate   = record_j.score + jump(j,i) + item_i                         (CI 515–519)
```

`L = base_seed_length` is the integer mean of `base_seed_length()` over the tree's anchors (CI 408–415). `record_j.paths` is the **running mask of j's DP record**, not j's own anchor mask. Scores are `int` (CH 348).

## 5. The single DP record and how it is updated

**Record** (`TracedScore`, CH 296–354): `score` (int), `source` (index of the predecessor, or `nowhere() = SIZE_MAX`), `paths` (one 64-bit mask), `rec_num`. Comparison operators order by `(score, source)` only (CH 333–340).

**Initialisation** (CI 431–435): every anchor starts as "from nowhere": `{item_i, nowhere, end_paths_i}`, with evaluation bonus `consistency_bonus`.

**Mask update** (`set_shared_paths`, CI 71–94). For a seed anchor `start_paths = end_paths = M_i`:
- if `record_j.paths & M_i == 0`: the new mask is `M_i` and `rec_num + 1` (a reset);
- otherwise the new mask is `record_j.paths & M_i`.

(An anchor with `start_paths ≠ end_paths`, possible only for merged extension anchors, resets to `end_paths` without counting; not used by hifi.)

**Consistency bonus** (CI 419–431, 521–533). It is used only to choose among candidates and is never added to the stored score (comment at CI 422–425; MMH 293):
- from nowhere: `bonus = consistency_bonus` (CI 461–465);
- from `j`, when `consistency_bonus > 0`, `pre = popcount(record_j.paths) > 0` and `record_j.paths & M_i ≠ 0`: `bonus = (consistency_bonus × popcount(new mask)) / pre` (integer division);
- from `j` after a reset: `bonus = 0`.

**Selection and ties** (CI 462–474, 535–548). Let `eval = score + bonus`. Each transition into `i` first re-offers the from-nowhere candidate, then the candidate from `j`. A candidate replaces the current record iff its `eval` is larger, or `eval` ties and its `(score, source)` is lexicographically larger. `nowhere = SIZE_MAX` is the largest source, so on a full tie a fresh start wins. Because all transitions into `i` are processed after every transition into their sources (CH 469–470; sources end before `i` starts in the read), and the replacement rule is a maximum under a fixed total order, **the final record is independent of the order in which transitions into `i` arrive**: it is the candidate with the lexicographically largest `(eval, score, source)`. Stored scores never fall below `item_i`, since every candidate bonus is at most the fresh-start bonus.

**Traceback** (CI 646–727, 729–857). Items are sorted by `(score, source)` descending; tracebacks start from unused items in that order, follow `source`, and stop early at a used item, with a score correction (CI 689–699). At most `max_chains_per_tree` (hifi: 3) disjoint tracebacks are returned in score order. The **best chain score is the largest stored score**, and its traceback is the chain the DP built. `find_best_chains` re-derives the recombination positions along each chain with the same intersect-or-reset scan (forward pass CI 773–808), adds a symmetric backward pass, and pairs them into intervals `[left, right]` for each recombination (CI 810–848). The count `rec_positions.size()` is the chain's recombination count (MFC 1767), later used only to penalise alignments by `rec-penalty-aln` (MFC 901–916).

## 6. What the single-record DP approximates

Notation: transitions `T` (section 3), items `item_i`, gap terms `jump0(j,i) ≤ 0`, masks `M_i` (the haplotype IDs whose paths contain the seed), penalty ρ.

**Exact objective.** Over chains `i_1 → … → i_k` with every consecutive pair in `T`, and labels `h_t ∈ M_{i_t}`, maximise

```
Σ_t item(i_t) + Σ_{t≥2} jump0(i_{t-1}, i_t) − ρ · #{ t ≥ 2 : h_t ≠ h_{t−1} }.
```

This is the project's basic objective ([MATHEMATICAL_MODEL.md](../../docs/MATHEMATICAL_MODEL.md)) with anchor weights `item_i`, the uniform penalty ρ, a label-independent transition cost `−jump0(j,i)`, and the transition set `T` (a subset of the strict precedence relation: conditions 1 and 4 of section 3 give `d_j < c_i`, and `y_j < x_i` on a shared vertex). The cost is not separable (`|read gap − graph gap|` inside a concave function), so the project's range-maximum engines do not apply; the exact DP enumerates `T`.

**Proposition (per chain, proved here).** Along a fixed chain with nonempty masks `M_1, …, M_k`, the number of resets `r` of the intersect-or-reset scan equals the minimum number of label changes over labelings `h_t ∈ M_t`.
*Proof.* The scan splits the chain into `r + 1` blocks whose masks have a common element (the running intersection at the block's end), so labeling each block with one such element gives `r` changes. Conversely, take any labeling with `c` changes and blocks ending at `e_1 < … < e_{c+1} = k`, and let the scan's blocks end at `E_1 < … < E_{r+1} = k`. By induction `E_m ≥ e_m`: the labeling's block `m` has a common label, so if the scan's block `m` starts at or before `e_m` (it starts at `E_{m−1} + 1 ≤ e_m` unless `E_{m−1} ≥ e_m` already), its running intersection stays nonempty through `e_m`. Hence the scan needs no more blocks than the labeling: `r ≤ c`. ∎

**Consequences.**
- The score stored at anchor `i` equals the exact objective value of the chain obtained by following `source` pointers from `i`, with an optimal labeling of that chain. So `stored_i ≤ OPT_i` (best exact value of chains ending at `i`) and Giraffe's best chain score is `≤ OPT = max_i OPT_i`. The emulation asserts `stored_i ≤ B_i` on every anchor.
- The rule is a heuristic for the **joint** choice of chain and labels: a record keeps one predecessor and one running mask, so a predecessor that scores more now but carries a mask that forces a later reset displaces one that scores less but keeps the right haplotypes. The novelty audit's four-anchor example (A{h1} 10, B{h2} 11, C{h1,h2} 10, D{h1} 10; A→C, B→C, C→D; ρ = 5) gives 26 against an optimum of 30; it is executed by `run_phase2.py --counterexample` and the result is in [phase2_counterexample.json](phase2_counterexample.json).
- With the consistency bonus on, stored scores are still exact values of their traced chains; only predecessor choice changes (bonus-preferred candidates may score less).

**The penalty is per reset of the running intersection, not per pairwise disjointness.** In v1.69.0 and v1.70.0 the rule was different: `add_points` returned `{score + adjustment, source}`, dropping the mask (v1.69.0 [CI 60–62](https://github.com/vgteam/vg/blob/v1.69.0/src/algorithms/chain_items.cpp#L60-L62); same in v1.70.0), so `set_shared_paths` always saw an empty running mask and stored the destination's own mask. The penalty then fired iff two **consecutive anchors' own masks** were disjoint. That is an exact DP for a pairwise objective that undercounts label changes (A{h1} → B{h1,h2} → C{h2} pays nothing, yet needs one change). v1.71.0 keeps the mask in `add_points` ([CI 60–61](https://github.com/vgteam/vg/blob/v1.71.0/src/algorithms/chain_items.cpp#L60-L61); PR #4753, "fixes a bug that reset the paths every time an anchor was added"). v1.69.0 also had no consistency bonus, an `item_scale` and `points_per_possible_match` term ([CI 569–652](https://github.com/vgteam/vg/blob/v1.69.0/src/algorithms/chain_items.cpp#L569-L652)), and the mode was opt-in; v1.75.0 made it the hifi/r10 default with ρ = 2 and the separate bonus (release notes). The rule specified here is v1.71.0 and later, with the bonus from v1.75.0.

## 7. The exact DP used as reference

Over the same `T`, items, `jump0` and masks, for `h ∈ M_i`:

```
C[i,h] = item_i + max( 0,  max_{(j,i) ∈ T} ( max( C[j,h] if h ∈ M_j,  B_j − ρ ) + jump0(j,i) ) ),
B_i = max_{h ∈ M_i} C[i,h],   OPT = max_i B_i.
```

The `0` is the fresh start, matching Giraffe's from-nowhere candidate. Including `g = h` in `B_j − ρ` is harmless for ρ ≥ 0 ([MATHEMATICAL_MODEL.md](../../docs/MATHEMATICAL_MODEL.md)). Cost `O(|T| · H)` with explicit transitions. Pointers store `(j, g)` per state.

**Exact pruning (implementation, 2026-10-06).** For a destination `i`, `u_j = B_j + jump0(j,i)` bounds every label's candidate from source `j`. After the source with the largest `u` is applied, every label holds at least `u_max − ρ` (its switch term), so a source with `u_j ≤ u_max − ρ` cannot strictly improve any label and is skipped. All `C`, `B` values are unchanged and only the choice among tied pointers can differ. On large true-strand instances this leaves about 2–4 full `O(H)` source updates per anchor (1.6–4.3 measured on six instances) instead of 700–1,200 transitions per anchor. `--no-prune` (`Params.prune = 0`) restores the original all-sources loop. The first four Phase 2 outputs were produced with it; 48 re-solves of stored rows, with and without pruning, reproduced their optima, best Giraffe scores and exact-chain lengths and switches exactly ([RESULTS.md](RESULTS.md)).

## 8. Emulation on the project's anchors: faithful parts and approximations

Project anchors are `(v, x, y, c, d, w)`, 1-based closed, inside one vertex, from pinned Minichain (k = 17, w = 11; every anchor has span 17).

| Item | Emulation | Status |
| --- | --- | --- |
| Mask | `M_i = A(v_i)`, the haplotypes visiting `v_i` | Faithful for in-node k-mers (GX 192–211); one bit per haplotype at H = 60 (GI 263–322) |
| Item points | `17 × match(1) + item_bonus(2) = 19` | Faithful formula; **Minichain seeds (k = 17), not Giraffe's k = 31, w = 50 minimizers** |
| `L` | 17 (all margins 0) | Faithful formula on these anchors |
| Graph distance | Exact minimum distance on the panel graph between seed start positions: `Z = D(v_j, v_i) + x_i − x_j`, `D` = shortest start-to-start distance; `graph_distance = Z − len_j` | Faithful to the zip-tree semantics (section 3); computed directly instead of through zip codes. Every panel edge is used by some haplotype, so the GBWTGraph has the same edges |
| Limits | graph lookback 20000 (or 0.105 × read length if larger), no read limit, indel ≤ max(5000, 2.45 × read length) | Faithful (hifi preset) |
| Gap score, ρ, bonus, ties | Sections 4–5, with the same integer truncations | Faithful |
| Tie order of identical read intervals | `std::sort` leaves it unspecified; emulation breaks ties by graph position `(v, x)`, then input order; a sensitivity run uses the reverse | Approximation; sensitivity reported |
| Zip-code forest | One chaining problem per (read, strand) instance | Approximation; transitions never cross trees, because the 2 × read length tree split exceeds the 20 kb lookback, so per-anchor scores match; tree filtering by coverage/score is not modelled |
| Seed filtering | None (all Minichain hits) | Approximation; Giraffe downsamples minimizers and applies hit caps (GM 962–970) |
| Alignment stage | Not modelled (`rec-penalty-aln`, base-level alignment) | Out of scope |

**Observation, not verified by execution.** When the recombination-aware mode is off (`use_payload_paths = false`), seed masks stay 0 (`Seed::paths = 0`). `check_recombination` (CI 377–383) then returns 1 on every transition and the forward recount (CI 790–795) counts every anchor after the first as a recombination. I found no place that zeroes `rec_penalty` when the mode is off (searched GM, MM, MFC for `rec_penalty`). If that is right, the hifi preset with `--no-rec-mode`, or on an index without path payloads, charges ρ on every transition and `rec-penalty-aln × (chain length − 1)` at alignment. This may be handled somewhere not read here; it is irrelevant to the emulation, which always has masks.
