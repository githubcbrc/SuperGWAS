#include <bits/stdc++.h>
#include <filesystem>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
namespace fs = std::filesystem;
#include "map_utils.hpp"

// --- list files by extension in a directory ---
static std::vector<std::string> list_files_by_ext(const std::string& dir_path, const std::string& ext) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            std::string name = entry.path().filename().string();
            if (name.size() > ext.size() && name.substr(name.size() - ext.size()) == ext) {
                files.push_back(entry.path().string());
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

// --- parse TSV, encode k-mers, append to per-shard vectors ---
// Supports 3-column (kmer \t effect \t p_value) and
// 4-column (kmer \t bitstring \t effect \t p_value) formats.
static size_t load_tsv_encoded(const std::string& path, size_t K,
                                uint32_t M,
                                vector<vector<pair<EncodedKmer, EntryD4>>>& shard_entries) {
    ifstream in(path);
    if (!in.is_open()) {
        cerr << "[ERROR] Could not open: " << path << "\n";
        return 0;
    }

    string line;
    bool first = true;
    size_t inserted = 0, seen = 0, skipped = 0;

    while (getline(in, line)) {
        if (line.empty()) continue;

        if (first) {
            first = false;
            if (looks_like_header(line)) continue;
        }

        size_t t1 = line.find('\t');
        if (t1 == string::npos) { ++skipped; ++seen; continue; }
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == string::npos) { ++skipped; ++seen; continue; }
        size_t t3 = line.find('\t', t2 + 1);

        const char* kmer_start = line.data();
        size_t kmer_len = t1;
        string bitstring;
        const char* eff_start;
        const char* pv_start;

        if (t3 != string::npos && line.find('\t', t3 + 1) == string::npos) {
            // 4-column: kmer \t bitstring \t effect \t pval
            bitstring.assign(line.data() + t1 + 1, t2 - t1 - 1);
            eff_start = line.data() + t2 + 1;
            pv_start = line.data() + t3 + 1;
        } else if (t3 == string::npos) {
            // 3-column: kmer \t effect \t pval
            eff_start = line.data() + t1 + 1;
            pv_start = line.data() + t2 + 1;
        } else {
            ++skipped; ++seen; continue;
        }

        if (kmer_len == 0 || kmer_len != K) { ++skipped; ++seen; continue; }

        // Check for invalid bases
        bool valid = true;
        for (size_t i = 0; i < kmer_len; ++i) {
            if (BASE_TO_2BIT.v[(unsigned char)kmer_start[i]] == ENCODE_INVALID) {
                valid = false; break;
            }
        }
        if (!valid) { ++skipped; ++seen; continue; }

        EncodedKmer enc = encode_kmer(kmer_start, K);
        EncodedKmer canonical = canonical_encoded(enc, K);
        uint32_t shard_id = static_cast<uint32_t>(static_cast<uint64_t>(canonical) & (M - 1));

        double eff = strtod(eff_start, nullptr);
        double pv = strtod(pv_start, nullptr);

        shard_entries[shard_id].emplace_back(canonical, EntryD4{eff, pv, move(bitstring)});
        ++inserted; ++seen;
    }

    cerr << "  [LOAD] rows=" << seen
         << " inserted=" << inserted
         << " skipped=" << skipped << "\n";
    return inserted;
}


static void print_usage(const char* prog) {
    cerr << "Usage:\n"
         << "  " << prog << " <input.tsv|input_dir> <out_dir> [--shards M] [--reserve N]\n\n"
         << "Input can be a single TSV file or a directory containing .tsv files.\n"
         << "Output is a directory of sharded binary index files (d4 format).\n"
         << "Default: 16 shards. M must be a power of 2.\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) { print_usage(argv[0]); return 1; }

    const string input_path = argv[1];
    const string out_dir    = argv[2];

    uint32_t M = 16;
    size_t reserve_hint = 10'000'000ULL;
    for (int i = 3; i < argc; ++i) {
        string a = argv[i];
        if (a == "--shards" && i + 1 < argc)       M = static_cast<uint32_t>(stoul(argv[++i]));
        else if (a == "--reserve" && i + 1 < argc)  reserve_hint = stoull(argv[++i]);
        else cerr << "[WARN] Unknown arg: " << a << "\n";
    }

    // Validate M is power of 2
    if (M == 0 || (M & (M - 1)) != 0) {
        cerr << "[ERROR] --shards must be a power of 2, got " << M << "\n";
        return 1;
    }

    // Determine input files
    vector<string> input_files;
    if (!fs::exists(input_path)) {
        cerr << "[ERROR] Path does not exist: " << input_path << "\n";
        return 1;
    }
    if (fs::is_directory(input_path)) {
        input_files = list_files_by_ext(input_path, ".tsv");
        if (input_files.empty()) {
            cerr << "[ERROR] No .tsv files in: " << input_path << "\n";
            return 1;
        }
        cerr << "[INFO] Found " << input_files.size() << " .tsv files in " << input_path << "\n";
    } else {
        input_files.push_back(input_path);
    }

    // Infer K from first file
    size_t K = infer_K(input_files[0]);
    if (!K) { cerr << "[ERROR] Could not infer k-mer length.\n"; return 1; }
    if (2 * K > 128) {
        cerr << "[ERROR] K=" << K << " too large for 128-bit encoding (max 64).\n";
        return 1;
    }
    cerr << "[INFO] Inferred k=" << K << ", shards=" << M << "\n";

    // Create output directory
    fs::create_directories(out_dir);

    // =========================================================================
    // Parallel TSV parsing
    // =========================================================================
    auto t0 = chrono::high_resolution_clock::now();

    // Global shard vectors
    vector<vector<pair<EncodedKmer, EntryD4>>> global_shards(M);
    for (auto& s : global_shards) s.reserve(reserve_hint / M);

    size_t total_inserted = 0;

#ifdef _OPENMP
    #pragma omp parallel
    {
        // Per-thread local shard vectors
        vector<vector<pair<EncodedKmer, EntryD4>>> local_shards(M);

        #pragma omp for schedule(dynamic, 1) reduction(+:total_inserted)
        for (size_t fi = 0; fi < input_files.size(); ++fi) {
            cerr << "[INFO] Loading: " << input_files[fi] << "\n";
            total_inserted += load_tsv_encoded(input_files[fi], K, M, local_shards);
        }

        // Merge local into global
        #pragma omp critical
        {
            for (uint32_t s = 0; s < M; ++s) {
                global_shards[s].insert(global_shards[s].end(),
                    make_move_iterator(local_shards[s].begin()),
                    make_move_iterator(local_shards[s].end()));
                local_shards[s].clear();
            }
        }
    }
#else
    for (const auto& file : input_files) {
        cerr << "[INFO] Loading: " << file << "\n";
        total_inserted += load_tsv_encoded(file, K, M, global_shards);
    }
#endif

    auto t1 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Parsed " << total_inserted << " k-mers in "
         << chrono::duration<double>(t1-t0).count() << " s\n";

    // =========================================================================
    // Dedup each shard: sort by encoded k-mer, keep first occurrence
    // =========================================================================
    auto t2 = chrono::high_resolution_clock::now();
    uint64_t total_unique = 0;

    for (uint32_t s = 0; s < M; ++s) {
        auto& vec = global_shards[s];
        sort(vec.begin(), vec.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
        auto last = unique(vec.begin(), vec.end(),
            [](const auto& a, const auto& b) { return a.first == b.first; });
        vec.erase(last, vec.end());
        total_unique += vec.size();
    }

    auto t3 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Dedup: " << total_unique << " unique k-mers in "
         << chrono::duration<double>(t3-t2).count() << " s\n";

    // =========================================================================
    // Write shards
    // =========================================================================
    auto t4 = chrono::high_resolution_clock::now();

    for (uint32_t s = 0; s < M; ++s) {
        if (!dump_shard(out_dir, s, static_cast<uint32_t>(K), global_shards[s])) {
            cerr << "[ERROR] Failed to write shard " << s << "\n";
            return 2;
        }
    }

    if (!write_index_meta(out_dir, static_cast<uint32_t>(K), M, total_unique)) {
        cerr << "[ERROR] Failed to write index.meta\n";
        return 2;
    }

    auto t5 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Wrote " << M << " shards to " << out_dir
         << " in " << chrono::duration<double>(t5-t4).count() << " s\n";

    // =========================================================================
    // Build and write bloom filters
    // =========================================================================
    auto t6 = chrono::high_resolution_clock::now();

    for (uint32_t s = 0; s < M; ++s) {
        BloomFilter bf;
        bf.init(global_shards[s].size(), 10);
        for (const auto& [enc, entry] : global_shards[s])
            bf.insert(enc);
        if (!dump_bloom(out_dir, s, bf)) {
            cerr << "[ERROR] Failed to write bloom filter for shard " << s << "\n";
            return 2;
        }
    }

    auto t7 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Wrote " << M << " bloom filters in "
         << chrono::duration<double>(t7-t6).count() << " s\n";
    cerr << "[INFO] Total time: " << chrono::duration<double>(t7-t0).count() << " s\n";

    return 0;
}
