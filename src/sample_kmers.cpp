// sample_kmers_distributed.cpp
// Reservoir-sample k-mers from multiple shards in a single job.
// Designed for SLURM array jobs where each job handles a batch of shards.
//
// Usage: sample_kmers_distributed <shards_dir> <accessions.txt> <output_dir>
//
// Environment variables:
//   RESERVOIR_SIZE   Total sample size across ALL jobs (default: 10000000)
//   SHARDS_PER_JOB   Number of shards this job processes (default: 50)
//   JOB_INDEX        This job's index, 0-based (default: SLURM_ARRAY_TASK_ID or 0)
//   TOTAL_JOBS       Total number of jobs (auto-calculated if not set)
//   SEED             Random seed (default: 42)
//
// Example with 1500 shards, 50 per job = 30 jobs:
//   SHARDS_PER_JOB=50 sbatch --array=0-29 sample_kmers.sbatch
//
// Output: <output_dir>/part_<JOB_INDEX>.tsv.gz
// Final step: cat all parts with concat_samples script

#include <bits/stdc++.h>
#include <zlib.h>
#include <filesystem>

namespace fs = std::filesystem;
using std::string; using std::vector; using std::size_t;

[[noreturn]] static void die(const string& m){ std::cerr << "[ERROR] " << m << "\n"; std::exit(1); }

static size_t env_size(const char* name, size_t def){
    const char* v = std::getenv(name);
    return v ? (size_t)std::stoull(v) : def;
}

static string env_str(const char* name, const string& def){
    const char* v = std::getenv(name);
    return v ? string(v) : def;
}

// ---------- Gz helpers ----------
struct GzReader {
    gzFile f = nullptr;
    explicit GzReader(const string& path){
        f = gzopen(path.c_str(), "rb");
        if(!f) die("gzopen(read): " + path);
        gzbuffer(f, 1<<20);
    }
    ~GzReader(){ if(f) gzclose(f); }
    GzReader(const GzReader&) = delete;
    bool getline(string& out){
        out.clear();
        int c;
        while(true){
            c = gzgetc(f);
            if(c == -1) return !out.empty();
            if(c == '\n') return true;
            out.push_back((char)c);
        }
    }
};

struct GzWriter {
    gzFile f = nullptr;
    explicit GzWriter(const string& path){
        f = gzopen(path.c_str(), "wb");
        if(!f) die("gzopen(write): " + path);
        gzbuffer(f, 1<<20);
    }
    ~GzWriter(){ if(f) gzclose(f); }
    GzWriter(const GzWriter&) = delete;
    void write(const char* data, size_t len){
        if(gzwrite(f, data, (unsigned)len) <= 0)
            die("gzwrite failed");
    }
};

// K-mer length: auto-detected from first k-mer, or from env K if set
static size_t KMER_LEN = 0;

int main(int argc, char** argv){
    if(argc < 4){
        std::cerr <<
            "Usage: sample_kmers_distributed <shards_dir> <accessions.txt> <output_dir>\n\n"
            "Environment variables:\n"
            "  RESERVOIR_SIZE   Total sample size (default: 10000000)\n"
            "  SHARDS_PER_JOB   Shards per job (default: 50)\n"
            "  JOB_INDEX        This job's index (default: SLURM_ARRAY_TASK_ID or 0)\n"
            "  TOTAL_JOBS       Total jobs (auto-calculated if not set)\n"
            "  SEED             Random seed (default: 42)\n\n"
            "Output: <output_dir>/part_<JOB_INDEX>.tsv.gz\n";
        return 1;
    }

    const string shards_dir      = argv[1];
    const string accessions_file = argv[2];
    const string output_dir      = argv[3];

    const size_t RESERVOIR_SIZE = env_size("RESERVOIR_SIZE", 10000000);
    const size_t SHARDS_PER_JOB = env_size("SHARDS_PER_JOB", 50);
    const size_t SEED           = env_size("SEED", 42);

    // Get job index from env or SLURM
    size_t JOB_INDEX = 0;
    {
        const char* v = std::getenv("JOB_INDEX");
        if (v) {
            JOB_INDEX = (size_t)std::stoull(v);
        } else {
            v = std::getenv("SLURM_ARRAY_TASK_ID");
            if (v) JOB_INDEX = (size_t)std::stoull(v);
        }
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    // Load accessions to know expected bitstring length
    size_t expected_bits = 0;
    {
        std::ifstream in(accessions_file);
        if(!in) die("cannot open " + accessions_file);
        string line;
        while(std::getline(in, line)){
            while(!line.empty() && isspace((unsigned char)line.back())) line.pop_back();
            if(!line.empty()) expected_bits++;
        }
    }

    // Enumerate all shards
    vector<fs::path> all_shards;
    for(auto& e : fs::directory_iterator(shards_dir)){
        auto fn = e.path().filename().string();
        if(fn.size() >= 7 && fn.substr(fn.size()-7) == ".tsv.gz")
            all_shards.push_back(e.path());
    }
    std::sort(all_shards.begin(), all_shards.end());

    if(all_shards.empty()) die("no .tsv.gz files in " + shards_dir);

    const size_t total_shards = all_shards.size();
    const size_t total_jobs = (total_shards + SHARDS_PER_JOB - 1) / SHARDS_PER_JOB;

    // Calculate this job's shard range
    const size_t shard_start = JOB_INDEX * SHARDS_PER_JOB;
    const size_t shard_end = std::min(shard_start + SHARDS_PER_JOB, total_shards);

    if (shard_start >= total_shards) {
        std::cerr << "[INFO] Job " << JOB_INDEX << " has no shards to process (start="
                  << shard_start << " >= total=" << total_shards << ")\n";
        return 0;
    }

    const size_t my_shard_count = shard_end - shard_start;

    // Per-job quota: divide total reservoir by number of jobs
    const size_t per_job_quota = (RESERVOIR_SIZE + total_jobs - 1) / total_jobs;

    std::cerr << "[INFO] Job " << JOB_INDEX << "/" << total_jobs
              << ": shards " << shard_start << "-" << (shard_end-1)
              << " (" << my_shard_count << " shards)\n"
              << "[INFO] Accessions: " << expected_bits
              << ", Total reservoir: " << RESERVOIR_SIZE
              << ", Per-job quota: " << per_job_quota << "\n";

    // Create output directory
    fs::create_directories(output_dir);

    // Reservoir sampling across all assigned shards
    std::mt19937_64 rng(SEED + JOB_INDEX);
    vector<string> reservoir;
    reservoir.reserve(std::min(per_job_quota, (size_t)2000000));

    size_t total_seen = 0;

    for (size_t si = shard_start; si < shard_end; ++si) {
        const auto& shard_path = all_shards[si];
        std::cerr << "[INFO] Processing shard " << si << ": " << shard_path.filename() << "\n";

        GzReader gz(shard_path.string());
        string line;

        while(gz.getline(line)){
            if(line.empty()) continue;
            size_t t = line.find('\t');
            if(t == string::npos) t = line.find(' ');
            if(t == string::npos) continue;

            // Auto-detect K from first valid k-mer
            if(KMER_LEN == 0){
                KMER_LEN = t;
                // Check against env K if set
                const char* env_k = std::getenv("K");
                if(env_k){
                    size_t expected_k = (size_t)std::stoull(env_k);
                    if(expected_k != KMER_LEN){
                        std::cerr << "[WARN] Detected k-mer length " << KMER_LEN
                                  << " differs from K=" << expected_k << " in environment\n";
                    }
                }
                std::cerr << "[INFO] Auto-detected k-mer length: " << KMER_LEN << "\n";
            }
            if(t != KMER_LEN) continue;

            size_t bits_len = line.size() - t - 1;
            while(bits_len > 0 && isspace((unsigned char)line[t + 1 + bits_len - 1])) --bits_len;
            if(bits_len != expected_bits) continue;

            if(total_seen < per_job_quota){
                reservoir.push_back(line.substr(0, t + 1 + bits_len));
            } else {
                std::uniform_int_distribution<size_t> dist(0, total_seen);
                size_t j = dist(rng);
                if(j < per_job_quota)
                    reservoir[j] = line.substr(0, t + 1 + bits_len);
            }
            total_seen++;
        }
    }

    // Write output partition
    string part_name = "part_" + std::to_string(JOB_INDEX) + ".tsv.gz";
    string output_path = (fs::path(output_dir) / part_name).string();

    size_t kept = reservoir.size();
    {
        GzWriter out(output_path);
        for(auto& rec : reservoir){
            rec.push_back('\n');
            out.write(rec.data(), rec.size());
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "[INFO] Total seen: " << total_seen << ", kept: " << kept << "\n"
              << "[INFO] Wrote " << output_path << " in " << secs << "s\n";

    // Write stats file
    string stats_path = (fs::path(output_dir) / ("part_" + std::to_string(JOB_INDEX) + ".stats")).string();
    {
        std::ofstream sf(stats_path);
        sf << "job_index\t" << JOB_INDEX << "\n"
           << "shard_start\t" << shard_start << "\n"
           << "shard_end\t" << shard_end << "\n"
           << "total_seen\t" << total_seen << "\n"
           << "kept\t" << kept << "\n"
           << "time_seconds\t" << secs << "\n";
    }

    return 0;
}
