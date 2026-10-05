// map_merge_filter.cpp
// Combined mapping, merging, and density filtering in one optimized pass.
// No intermediate files - all processing done in memory with parallel execution.
// Optionally emits a second output for 2D epistasis (with bitstrings, no filter).
//
// Usage: map_merge_filter <index.umap> <reference.fa> <k> <output.tsv>
//        [--threads N] [--min-logp X] [--chunk-size M]
//        [--window-size W] [--min-density D]
//        [--emit-2d <path> --min-logp-2d X]

#include <bits/stdc++.h>
#include <filesystem>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#include "map_utils.hpp"

struct Hit {
    size_t contig_idx;   // which contig
    size_t pos;          // position within contig
    double effect;
    double pval;
    std::string bitstring;
};

struct Prec {
    static constexpr int EFFECT = 9;
    static constexpr int PVAL   = 9;
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 5) {
        cerr <<
        "Map k-mers to reference, merge chromosomes, and filter by density.\n"
        "All in one pass with no intermediate files.\n\n"
        "Usage:\n"
        "  " << argv[0] << " <index.umap> <reference.fa> <k> <output.tsv>\n"
        "    [--threads N] [--min-logp X] [--chunk-size M]\n"
        "    [--window-size W] [--min-density D]\n"
        "    [--emit-2d <path> --min-logp-2d X]\n\n"
        "Options:\n"
        "  --threads N       OpenMP threads (default: OMP max)\n"
        "  --min-logp X      Keep hits with -log10(p) >= X (default: 0 = no filter)\n"
        "  --chunk-size M    Split long contigs for intra-contig parallelism\n"
        "  --window-size W   Density filter window size in bp (default: 100)\n"
        "  --min-density D   Min fraction of positions with hits (default: 0.95)\n"
        "  --emit-2d PATH    Also write 2D input (pos, bitstring, effect, pval) to PATH\n"
        "  --min-logp-2d X   Min -log10(p) for 2D output (default: 15)\n";
        return 1;
    }

    const string umap_path = argv[1];
    const string ref_path  = argv[2];
    const size_t K         = (size_t)stoull(argv[3]);
    const string out_path  = argv[4];

    int    threads      = 0;
    double min_logp     = 0.0;
    size_t chunk_size   = 0;
    size_t window_size  = 100;
    double min_density  = 0.95;
    string emit_2d_path;
    double min_logp_2d  = 15.0;

    for (int i = 5; i < argc; ++i) {
        string a = argv[i];
        if (a == "--threads" && i+1 < argc)         threads = stoi(argv[++i]);
        else if (a == "--min-logp" && i+1 < argc)   min_logp = stod(argv[++i]);
        else if (a == "--chunk-size" && i+1 < argc)  chunk_size = (size_t)stoull(argv[++i]);
        else if (a == "--window-size" && i+1 < argc) window_size = (size_t)stoull(argv[++i]);
        else if (a == "--min-density" && i+1 < argc) min_density = stod(argv[++i]);
        else if (a == "--emit-2d" && i+1 < argc)    emit_2d_path = argv[++i];
        else if (a == "--min-logp-2d" && i+1 < argc) min_logp_2d = stod(argv[++i]);
        else cerr << "[WARN] Unknown arg: " << a << "\n";
    }

    const bool do_2d = !emit_2d_path.empty();

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    const int used_threads = (threads > 0 ? threads : omp_get_max_threads());
#else
    const int used_threads = 1;
#endif

    const double p_thresh = (min_logp <= 0.0) ? 1.0 : pow(10.0, -min_logp);
    const double p_thresh_2d = pow(10.0, -min_logp_2d);
    const size_t min_count = (size_t)(window_size * min_density);

    cerr << "[INFO] Threads: " << used_threads
         << ", min-logp: " << min_logp
         << ", window: " << window_size
         << ", min-density: " << min_density
         << " (" << min_count << " hits)";
    if (do_2d) cerr << ", 2D: min-logp-2d=" << min_logp_2d
                     << " -> " << emit_2d_path;
    cerr << "\n";

    // =========================================================================
    // Step 1: Load sharded d4 index + bloom filters
    // Compact (no bitstrings) when .bits absent → fast 1D path
    // Full (with bitstrings) when .bits present → needed for 2D output
    // =========================================================================
    vector<ShardMapCompact> shards;       // 1D path (compact)
    vector<ShardMap>        shards_full;  // 2D path (with bitstrings)
    vector<BloomFilter> blooms;
    uint32_t M_shards = 0;
    bool use_full = false;
    {
        if (!std::filesystem::is_directory(umap_path)) {
            cerr << "[ERROR] Index path is not a directory: " << umap_path << "\n";
            return 2;
        }
        auto t0 = chrono::high_resolution_clock::now();
        uint32_t idx_K; uint64_t idx_N;
        if (!read_index_meta(umap_path, idx_K, M_shards, idx_N)) return 2;

        use_full = std::filesystem::exists(umap_path + "/shard_000.bits");

        if (use_full) shards_full.resize(M_shards);
        else          shards.resize(M_shards);
        blooms.resize(M_shards);

        bool load_ok = true;
        bool have_blooms = true;
        #pragma omp parallel for schedule(dynamic, 1)
        for (uint32_t i = 0; i < M_shards; ++i) {
            bool ok = use_full
                ? load_shard(umap_path, i, idx_K, shards_full[i], true)
                : load_shard_compact(umap_path, i, idx_K, shards[i]);
            if (!ok) {
                #pragma omp critical
                { load_ok = false; }
            }
            if (!load_bloom(umap_path, i, blooms[i])) {
                #pragma omp critical
                { have_blooms = false; }
            }
        }
        if (!load_ok) return 2;

        // Build bloom filters in-memory if not on disk
        if (!have_blooms) {
            cerr << "[INFO] No bloom filters on disk, building in-memory...\n";
            for (uint32_t i = 0; i < M_shards; ++i) {
                if (use_full) {
                    blooms[i].init(shards_full[i].size(), 10);
                    for (const auto& [k, v] : shards_full[i]) blooms[i].insert(k);
                } else {
                    blooms[i].init(shards[i].size(), 10);
                    for (const auto& [k, v] : shards[i]) blooms[i].insert(k);
                }
            }
        }

        uint64_t total = 0;
        size_t bloom_mb = 0;
        for (uint32_t i = 0; i < M_shards; ++i) {
            total += use_full ? shards_full[i].size() : shards[i].size();
            bloom_mb += blooms[i].size_bytes();
        }
        auto t1 = chrono::high_resolution_clock::now();
        cerr << "[INFO] Loaded index: " << total << " kmers ("
             << M_shards << " shards, bloom " << (bloom_mb >> 20) << " MB"
             << (use_full ? ", with bitstrings" : ", compact") << ") in "
             << chrono::duration<double>(t1-t0).count() << " s\n";
    }

    // =========================================================================
    // Step 2: Load reference
    // =========================================================================
    vector<Contig> contigs;
    {
        auto t0 = chrono::high_resolution_clock::now();
        if (!read_reference_fasta(ref_path, contigs)) return 3;
        auto t1 = chrono::high_resolution_clock::now();
        size_t total_bp = 0;
        for (auto& c : contigs) total_bp += c.seq.size();
        cerr << "[INFO] Loaded reference: " << contigs.size() << " contigs, "
             << total_bp << " bp in " << chrono::duration<double>(t1-t0).count() << " s\n";
    }
    if (contigs.empty()) { cerr << "[WARN] No contigs.\n"; return 0; }

    // =========================================================================
    // Precompute contig offsets for genome-wide coordinates
    // =========================================================================
    vector<size_t> contig_offsets(contigs.size());
    {
        size_t offset = 0;
        for (size_t i = 0; i < contigs.size(); ++i) {
            contig_offsets[i] = offset;
            offset += contigs[i].seq.size();
        }
        cerr << "[INFO] Genome size: " << offset << " bp\n";
    }

    // =========================================================================
    // Step 3: Map + density filter per contig (fused, parallel)
    // =========================================================================
    // Each thread: map contig → density filter → collect 1D/2D hits
    // Hits within a contig are already position-sorted from the L→R scan.
    // Density filter uses two-pointer on sorted hits: O(hits) not O(genome).
    //
    // Template lambda handles both compact (1D) and full (2D with bitstrings)
    // shard types without duplicating the mapping loop.

    vector<vector<Hit>> thread_1d(used_threads);
    vector<vector<Hit>> thread_2d(used_threads);

    auto T0 = chrono::high_resolution_clock::now();

    const EncodedKmer mask = KMER_MASK(K);

    // Lookup helper: works for both ShardMapCompact and ShardMap
    // Returns {effect, pval, bitstring} or nullopt
    auto map_and_filter = [&](auto& shard_vec, vector<Hit>& local_1d, vector<Hit>& local_2d,
                               vector<Hit>& contig_hits, size_t cidx) {
        const string& seq = contigs[cidx].seq;
        if (seq.size() < K) return;
        const size_t gw_offset = contig_offsets[cidx];

        // --- Map: sliding-window encoded lookup with bloom gate ---
        contig_hits.clear();
        EncodedKmer enc = 0;
        size_t valid = 0;

        for (size_t i = 0; i < seq.size(); ++i) {
            uint8_t code = BASE_TO_2BIT.v[(unsigned char)seq[i]];
            if (code == ENCODE_INVALID) { valid = 0; enc = 0; continue; }
            enc = ((enc << 2) | code) & mask;
            if (++valid < K) continue;

            EncodedKmer canonical = canonical_encoded(enc, K);
            uint32_t si = static_cast<uint64_t>(canonical) & (M_shards - 1);
            if (!blooms[si].test(canonical)) continue;
            auto it = shard_vec[si].find(canonical);
            if (it == shard_vec[si].end()) continue;
            if (it->second.pval > p_thresh) continue;

            size_t pos = i - K + 1;
            if constexpr (std::is_same_v<std::decay_t<decltype(it->second)>, EntryD4>) {
                contig_hits.push_back({cidx, pos, it->second.effect, it->second.pval,
                    (do_2d && it->second.pval <= p_thresh_2d) ? it->second.bitstring : string{}});
            } else {
                contig_hits.push_back({cidx, pos, it->second.effect, it->second.pval, {}});
            }
        }

        if (contig_hits.empty()) return;

        // --- 2D output: all hits with bitstrings, no density filter ---
        if (do_2d) {
            for (auto& h : contig_hits) {
                if (!h.bitstring.empty()) {
                    Hit h2d = h;
                    h2d.pos += gw_offset;
                    local_2d.push_back(std::move(h2d));
                }
            }
        }

        // --- Density filter: two-pointer on position-sorted hits ---
        size_t n = contig_hits.size();
        size_t good_until = 0;
        size_t right = 0;

        for (size_t left = 0; left < n; ++left) {
            while (right < n &&
                   contig_hits[right].pos < contig_hits[left].pos + window_size)
                ++right;
            if (right - left >= min_count)
                good_until = max(good_until,
                                 contig_hits[left].pos + window_size - 1);

            if (contig_hits[left].pos <= good_until) {
                Hit h = contig_hits[left];
                h.pos += gw_offset;
                h.bitstring.clear();
                local_1d.push_back(std::move(h));
            }
        }
    };

#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        vector<Hit>& local_1d = thread_1d[tid];
        vector<Hit>& local_2d = thread_2d[tid];
        local_1d.reserve(100000);
        vector<Hit> contig_hits;

        #pragma omp for schedule(dynamic, 1)
        for (size_t cidx = 0; cidx < contigs.size(); ++cidx) {
            if (use_full)
                map_and_filter(shards_full, local_1d, local_2d, contig_hits, cidx);
            else
                map_and_filter(shards, local_1d, local_2d, contig_hits, cidx);
        }
    }
#else
    {
        vector<Hit>& local_1d = thread_1d[0];
        vector<Hit>& local_2d = thread_2d[0];
        vector<Hit> contig_hits;

        for (size_t cidx = 0; cidx < contigs.size(); ++cidx) {
            if (use_full)
                map_and_filter(shards_full, local_1d, local_2d, contig_hits, cidx);
            else
                map_and_filter(shards, local_1d, local_2d, contig_hits, cidx);
        }
    }
#endif

    // Merge thread buffers
    size_t total_raw = 0, total_kept = 0, total_2d = 0;
    vector<Hit> kept_hits, hits_2d;
    {
        for (auto& v : thread_1d) total_kept += v.size();
        for (auto& v : thread_2d) total_2d += v.size();
        kept_hits.reserve(total_kept);
        for (auto& v : thread_1d) {
            kept_hits.insert(kept_hits.end(),
                make_move_iterator(v.begin()), make_move_iterator(v.end()));
            total_raw += v.size(); // approximate; exact raw counted below
            v.clear(); v.shrink_to_fit();
        }
        if (do_2d) {
            hits_2d.reserve(total_2d);
            for (auto& v : thread_2d) {
                hits_2d.insert(hits_2d.end(),
                    make_move_iterator(v.begin()), make_move_iterator(v.end()));
                v.clear(); v.shrink_to_fit();
            }
        }
    }

    auto T1 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Map+filter: " << total_kept << " kept in "
         << chrono::duration<double>(T1-T0).count() << " s\n";

    // =========================================================================
    // Step 4: Write 1D output (density-filtered, no bitstrings)
    // =========================================================================
    {
        auto T6 = chrono::high_resolution_clock::now();
        ofstream out(out_path);
        if (!out) {
            cerr << "[ERROR] Cannot open output: " << out_path << "\n";
            return 4;
        }

        out << std::setprecision(Prec::EFFECT) << std::defaultfloat;
        for (const auto& h : kept_hits) {
            out << h.pos << '\t'
                << h.effect << '\t'
                << std::setprecision(Prec::PVAL) << std::scientific << h.pval
                << std::setprecision(Prec::EFFECT) << std::defaultfloat << '\n';
        }

        auto T7 = chrono::high_resolution_clock::now();
        cerr << "[INFO] Wrote " << kept_hits.size() << " hits to " << out_path
             << " in " << chrono::duration<double>(T7-T6).count() << " s\n";
    }

    // =========================================================================
    // Step 5: Write 2D output (no density filter, with bitstrings)
    // =========================================================================
    if (do_2d) {
        auto T8 = chrono::high_resolution_clock::now();
        ofstream out2d(emit_2d_path);
        if (!out2d) {
            cerr << "[ERROR] Cannot open 2D output: " << emit_2d_path << "\n";
            return 4;
        }

        out2d << std::setprecision(Prec::EFFECT) << std::defaultfloat;
        for (const auto& h : hits_2d) {
            out2d << h.pos << '\t'
                  << h.bitstring << '\t'
                  << h.effect << '\t'
                  << std::setprecision(Prec::PVAL) << std::scientific << h.pval
                  << std::setprecision(Prec::EFFECT) << std::defaultfloat << '\n';
        }

        auto T9 = chrono::high_resolution_clock::now();
        cerr << "[INFO] Wrote " << hits_2d.size() << " 2D hits to " << emit_2d_path
             << " in " << chrono::duration<double>(T9-T8).count() << " s\n";
    }

    auto Tend = chrono::high_resolution_clock::now();
    cerr << "[INFO] Total time: " << chrono::duration<double>(Tend-T0).count() << " s\n";

    return 0;
}
