// kgwas_2d.cpp — 2D epistasis: all-pairs (i<j) interaction testing on index k-mers
// Loads k-mer bitstrings directly from the sharded d4 index (no reference needed).
// Outputs pairs with k-mer indices + a companion kmers file for downstream mapping.
//
// Usage: kgwas_2d samples.txt phenotype.txt U2D.bin INDEX_DIR out_prefix
//        [log10_thresh] [--emit-threads N] [--out-chunk N]

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <omp.h>
#include <fcntl.h>
#include <unistd.h>

extern "C" void openblas_set_num_threads(int);
#include "kgwas_utils.hpp"
#include "map_utils.hpp"
#include "chisq.hpp"

using std::string;
using std::vector;
using std::size_t;
using namespace kgwas;

// ---------------- Types ----------------
struct Config {
  string samples_file, pheno_file;
  string U_file;
  string index_dir;
  string out_prefix;
  double log10_thresh = 6.0;
  int BLAS_THREADS = 112;
  int EMIT_THREADS = 112;
  int OUT_CHUNK    = 1000;
  int task_id      = 0;
  int num_tasks    = 1;
};

struct Pair {
  int64_t idx_a, idx_b;
  double effect, pval;
};

// ---------------- Arg parsing ----------------
static Config parse_args(int argc, char** argv) {
  if (argc < 6) {
    std::cerr
      << "Usage: " << argv[0]
      << " samples.txt phenotype.txt U2D.bin INDEX_DIR out_prefix"
      << " [log10_thresh] [--emit-threads N] [--out-chunk N]\n"
      << "\nOutputs:\n"
      << "  out_prefix.tsv   — idx_a, idx_b, effect, pvalue\n"
      << "  out_prefix.kmers — one k-mer sequence per line (line# = index)\n";
    std::exit(1);
  }
  Config cfg;
  cfg.samples_file = argv[1];
  cfg.pheno_file   = argv[2];
  cfg.U_file       = argv[3];
  cfg.index_dir    = argv[4];
  cfg.out_prefix   = argv[5];
  if (argc > 6 && argv[6][0] != '-') cfg.log10_thresh = std::stod(argv[6]);

  int start = (argc > 6 && argv[6][0] != '-') ? 7 : 6;
  for (int i = start; i < argc; ++i) {
    string a = argv[i];
    if (a == "--emit-threads" && i+1 < argc)    cfg.EMIT_THREADS = std::stoi(argv[++i]);
    else if (a == "--out-chunk" && i+1 < argc)  cfg.OUT_CHUNK    = std::stoi(argv[++i]);
    else if (a == "--task-id" && i+1 < argc)    cfg.task_id      = std::stoi(argv[++i]);
    else if (a == "--num-tasks" && i+1 < argc)  cfg.num_tasks    = std::stoi(argv[++i]);
  }

  cfg.BLAS_THREADS = env_int("BLAS_THREADS", cfg.BLAS_THREADS);
  cfg.EMIT_THREADS = env_int("EMIT_THREADS", cfg.EMIT_THREADS);
  return cfg;
}

// -------- Load X matrix directly from index shards --------
static void load_X_from_index(const Config& cfg,
                               const std::vector<int>& overlap_idxs,
                               int n,
                               std::vector<std::string>& kmer_seqs_out,
                               MatRow& Xrm_out) {
  Timer t;
  kmer_seqs_out.clear();

  // --- Load sharded d4 index WITH bitstrings ---
  vector<ShardMap> shards_full;
  uint32_t M_shards = 0;
  uint32_t idx_K = 0;
  {
    if (!std::filesystem::is_directory(cfg.index_dir))
      die("Index path is not a directory: " + cfg.index_dir);

    uint64_t idx_N;
    if (!read_index_meta(cfg.index_dir, idx_K, M_shards, idx_N))
      die("Failed to read index meta from " + cfg.index_dir);

    shards_full.resize(M_shards);

    bool load_ok = true;
    #pragma omp parallel for schedule(dynamic, 1)
    for (uint32_t i = 0; i < M_shards; ++i) {
      bool ok = load_shard(cfg.index_dir, i, idx_K, shards_full[i], true);
      if (!ok) {
        #pragma omp critical
        { load_ok = false; }
      }
    }
    if (!load_ok) die("Failed to load index shards");

    uint64_t total = 0;
    for (uint32_t i = 0; i < M_shards; ++i) total += shards_full[i].size();
    std::cerr << "[INFO] Loaded index: " << total << " kmers ("
              << M_shards << " shards, K=" << idx_K << ") in "
              << t.sec() << " s\n";
  }

  // --- Collect all entries with bitstrings, sorted by encoded k-mer for determinism ---
  struct IndexEntry {
    EncodedKmer enc;
    std::string bitstring;
  };

  vector<IndexEntry> entries;
  {
    size_t total_with_bits = 0;
    for (uint32_t i = 0; i < M_shards; ++i) {
      for (const auto& [enc, entry] : shards_full[i]) {
        if (!entry.bitstring.empty()) ++total_with_bits;
      }
    }
    entries.reserve(total_with_bits);

    for (uint32_t i = 0; i < M_shards; ++i) {
      for (auto& [enc, entry] : shards_full[i]) {
        if (!entry.bitstring.empty()) {
          entries.push_back({enc, std::move(entry.bitstring)});
        }
      }
    }

    // Sort by encoded k-mer for deterministic ordering
    std::sort(entries.begin(), entries.end(),
              [](const IndexEntry& a, const IndexEntry& b) { return a.enc < b.enc; });

    std::cerr << "[INFO] " << entries.size() << " k-mers with bitstrings\n";
  }

  // Free shards (bitstrings moved out)
  shards_full.clear();
  shards_full.shrink_to_fit();

  if (entries.empty()) die("No k-mers with bitstrings in index");

  // --- Decode k-mer sequences for the kmers file ---
  kmer_seqs_out.resize(entries.size());
  for (size_t i = 0; i < entries.size(); ++i) {
    kmer_seqs_out[i].resize(idx_K);
    decode_kmer(entries[i].enc, idx_K, kmer_seqs_out[i].data());
  }

  // --- Build X matrix from bitstrings ---
  int maxidx = -1;
  if (!overlap_idxs.empty())
    maxidx = *std::max_element(overlap_idxs.begin(), overlap_idxs.end());

  std::vector<float> data;
  data.reserve(static_cast<std::size_t>(n) * entries.size());

  // Track which entries are usable (valid bitstring length)
  vector<bool> valid(entries.size(), false);

  for (size_t j = 0; j < entries.size(); ++j) {
    const int B = static_cast<int>(entries[j].bitstring.size());
    if (B == n) {
      for (int i = 0; i < n; ++i)
        data.push_back((entries[j].bitstring[i] == '1') ? 1.0f : 0.0f);
      valid[j] = true;
    } else if (maxidx >= 0 && B > maxidx) {
      for (int i = 0; i < n; ++i)
        data.push_back((entries[j].bitstring[overlap_idxs[i]] == '1') ? 1.0f : 0.0f);
      valid[j] = true;
    }
  }

  // Compact: remove invalid entries from kmer_seqs_out
  if (std::count(valid.begin(), valid.end(), true) < (int64_t)entries.size()) {
    vector<std::string> compact_seqs;
    compact_seqs.reserve(entries.size());
    for (size_t j = 0; j < entries.size(); ++j) {
      if (valid[j]) compact_seqs.push_back(std::move(kmer_seqs_out[j]));
    }
    kmer_seqs_out = std::move(compact_seqs);
  }

  const int64_t N = static_cast<int64_t>(kmer_seqs_out.size());
  if (N <= 0) die("No valid k-mers after bitstring parsing");

  Xrm_out.resize(n, static_cast<int>(N));
  for (int64_t j = 0; j < N; ++j) {
    const float* src = &data[static_cast<std::size_t>(j * n)];
    for (int i = 0; i < n; ++i) Xrm_out(i, static_cast<int>(j)) = src[i];
  }

  std::cerr << "[INFO] Built X from index (N=" << N << ", n=" << n
            << ", " << t.sec() << " s)\n";
}

// --------------- Stats caches for self-2D ---------------
struct Stats2D {
  double yyM = 0.0;
  int n = 0, r = 0, N = 0;
  Eigen::MatrixXf Uf;
  Eigen::VectorXf yf;
  Eigen::VectorXf uyf;
  Eigen::MatrixXf UX;
  Eigen::VectorXf Xy;
  Eigen::ArrayXf  xx;
};

static Stats2D build_stats2d_self(const Eigen::MatrixXd& U,
                                  const Eigen::VectorXd& y,
                                  const MatRow& Xrm,
                                  int BLAS_THREADS) {
  Timer t;
  omp_set_num_threads(1);
  openblas_set_num_threads(BLAS_THREADS);

  const int n  = (int)U.rows();
  const int r  = (int)U.cols();
  const int N  = (int)Xrm.cols();

  Eigen::MatrixXf Uf = U.cast<float>();
  Eigen::VectorXf yf = y.cast<float>();
  Eigen::VectorXf uyf = Uf.transpose() * yf;

  Eigen::MatrixXf UX = Uf.transpose() * Xrm;
  Eigen::VectorXf Xy = Xrm.transpose() * yf;
  Eigen::ArrayXf  xx = Xrm.colwise().squaredNorm().array();

  const Eigen::VectorXd Uy = U.transpose()*y;
  const double yyM = y.squaredNorm() - Uy.squaredNorm();
  if(!(yyM>1e-12)) warn("yyM ~ 0");

  Stats2D st;
  st.n=n; st.r=r; st.N=N;
  st.yyM = yyM;
  st.Uf  = std::move(Uf);
  st.yf  = std::move(yf);
  st.uyf = std::move(uyf);
  st.UX  = std::move(UX);
  st.Xy  = std::move(Xy);
  st.xx  = std::move(xx);
  std::cerr << "[INFO] self 2D caches ("<<t.sec()<<" s, BLAS_THREADS="<<BLAS_THREADS<<")\n";
  return st;
}

static void emit_pairs_self(const vector<std::string>& kmer_seqs,
                            const MatRow& Xrm,
                            const Stats2D& st,
                            const Config& cfg) {
  Timer t;
  openblas_set_num_threads(1);

  #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(cfg.EMIT_THREADS);
    const int used_threads = omp_get_max_threads();
  #else
    const int used_threads = 1;
  #endif

  const double p_thresh = std::pow(10.0, -cfg.log10_thresh);
  const double eps = 1e-12;

  const int64_t N  = (int64_t)kmer_seqs.size();
  const int     n  = st.n;
  const int64_t CH = std::max<int>(1, cfg.OUT_CHUNK);

  if (N < 2) die("emit_pairs_self: need at least 2 columns");

  // Row range for this task
  const int64_t rows_per_task = (N + cfg.num_tasks - 1) / cfg.num_tasks;
  const int64_t row_start = (int64_t)cfg.task_id * rows_per_task;
  const int64_t row_end   = std::min(row_start + rows_per_task, N);

  if (row_start >= N) {
    std::cerr << "[INFO] task " << cfg.task_id << ": no rows to process\n";
    return;
  }

  // Count pairs in this task's range
  std::uint64_t task_pairs = 0;
  for (int64_t i = row_start; i < row_end; ++i) task_pairs += (N - i - 1);

  const auto t0 = std::chrono::steady_clock::now();
  std::atomic<std::uint64_t> pairs_done{0};
  std::atomic<int64_t>       chunks_done{0};
  const int64_t total_chunks = (row_end - row_start + CH - 1) / CH;

  auto pairs_in_chunk = [&](int64_t s, int64_t e) -> std::uint64_t {
    const long long m = (long long)(e - s);
    if (m <= 0) return 0ULL;
    const __int128 tri = (__int128)m * (N - 1) - ((__int128)(s + (e - 1)) * m) / 2;
    return (tri > 0) ? (std::uint64_t)tri : 0ULL;
  };

  std::cerr << "[INFO] emit: N=" << N << " rows=[" << row_start << "," << row_end << ")"
            << " n=" << n << " chunks=" << total_chunks << " CH=" << CH
            << " threads=" << used_threads
            << " task=" << cfg.task_id << "/" << cfg.num_tasks
            << "  (p <= " << std::scientific << p_thresh << std::defaultfloat << ")\n";

  auto progress_print = [&](bool force=false){
    static std::atomic<int>         last_pct{-1};
    static std::atomic<long long>   last_ms{0};

    const auto  now = std::chrono::steady_clock::now();
    const long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - t0).count();
    const std::uint64_t done = pairs_done.load(std::memory_order_relaxed);
    const int pct = (int)std::floor(100.0 * (double)done / (double)task_pairs);

    bool do_print = force;
    int expected = last_pct.load(std::memory_order_relaxed);
    if (!do_print && pct >= expected + 1)
      do_print = last_pct.compare_exchange_strong(expected, pct);
    long long lm = last_ms.load(std::memory_order_relaxed);
    if (!do_print && ms - lm >= 2000)
      do_print = last_ms.compare_exchange_strong(lm, ms);

    if (do_print) {
      const double secs = ms / 1000.0;
      const double rate = secs > 0 ? (double)done / secs : 0.0;
      const double rem  = rate > 0 ? ((double)task_pairs - (double)done) / rate : 0.0;
      const int64_t cd  = chunks_done.load(std::memory_order_relaxed);
      std::cerr << "[PROG] " << pct << "%  chunks " << cd << "/" << total_chunks
                << "  pairs " << done << "/" << task_pairs
                << "  rate " << std::fixed << std::setprecision(2) << rate << "/s"
                << "  ETA "  << std::setprecision(1) << rem << " s"
                << "        \r" << std::flush;
    }
  };

  // --- Collect all passing pairs into per-thread vectors ---
  vector<vector<Pair>> thread_pairs(used_threads);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided,1)
  #endif
  for (int64_t s = row_start; s < row_end; s += CH) {
    const int64_t e = std::min<int64_t>(s + CH, row_end);
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif
    auto& local = thread_pairs[tid];

    Eigen::VectorXf x12, t12;
    x12.resize(n); t12.resize(st.r);

    for (int64_t i = s; i < e; ++i) {
      const float a11 = st.xx[(int)i];
      if (!(a11 > eps)) continue;
      const Eigen::VectorXf x1col = Xrm.col((int)i);

      for (int64_t j = i + 1; j < N; ++j) {
        const float a22 = st.xx[(int)j];
        if (!(a22 > eps)) continue;
        const Eigen::VectorXf x2col = Xrm.col((int)j);

        x12.array() = x1col.array() * x2col.array();
        const float a33 = x12.squaredNorm();
        if (!(a33 > eps)) continue;

        const float b3 = x12.dot(st.yf);
        t12.noalias() = st.Uf.transpose() * x12;

        const double xxM = (double)a33 - (double)t12.squaredNorm();
        if (xxM <= eps || st.yyM <= eps) continue;

        const double xyM = (double)b3 - (double)t12.dot(st.uyf);
        double r2 = (xyM * xyM) / (xxM * st.yyM);
        if (r2 < 0) r2 = 0; if (r2 >= 1) r2 = 1 - 1e-15;

        const double T = (double)n * std::log(1.0 / (1.0 - r2));
        const double p = chi2_sf_k1(T);
        if (p > p_thresh) continue;

        const double effect = xyM / (xxM * st.yyM);

        local.push_back({i, j, effect, p});
      }
    }

    pairs_done.fetch_add(pairs_in_chunk(s, e), std::memory_order_relaxed);
    chunks_done.fetch_add(1, std::memory_order_relaxed);
    progress_print(false);
  }

  progress_print(true);
  std::cerr << "\n";

  // --- Merge thread buffers and sort ---
  vector<Pair> all_pairs;
  {
    size_t total = 0;
    for (auto& v : thread_pairs) total += v.size();
    all_pairs.reserve(total);
    for (auto& v : thread_pairs) {
      all_pairs.insert(all_pairs.end(),
        std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
      v.clear(); v.shrink_to_fit();
    }
    std::sort(all_pairs.begin(), all_pairs.end(), [](const Pair& a, const Pair& b) {
      if (a.idx_a != b.idx_a) return a.idx_a < b.idx_a;
      return a.idx_b < b.idx_b;
    });
  }

  std::cerr << "[INFO] found " << all_pairs.size() << " passing pairs in "
            << t.sec() << " s, threads=" << used_threads << "\n";

  // --- Write pairs output (per-task file if parallel) ---
  {
    Timer tw;
    string pairs_path = cfg.out_prefix;
    if (cfg.num_tasks > 1)
      pairs_path += ".part_" + std::to_string(cfg.task_id);
    pairs_path += ".tsv";
    std::ofstream out(pairs_path, std::ios::binary);
    if (!out) die("open out: " + pairs_path);

    std::ostringstream oss;
    oss.setf(std::ios::fmtflags(0), std::ios::floatfield);
    for (const auto& p : all_pairs) {
      oss << p.idx_a << '\t' << p.idx_b << '\t'
          << std::setprecision(9) << std::defaultfloat << p.effect << '\t'
          << std::setprecision(9) << std::scientific   << p.pval   << '\n';
    }
    const std::string buf = std::move(oss).str();
    out.write(buf.data(), (std::streamsize)buf.size());
    out.close();
    if (!out) die("write out failed");

    std::cerr << "[INFO] wrote " << all_pairs.size() << " pairs to " << pairs_path
              << " in " << tw.sec() << " s\n";
  }

  // --- Write kmers file (only from task 0) ---
  if (cfg.task_id == 0) {
    Timer tw;
    string kmers_path = cfg.out_prefix + ".kmers";
    std::ofstream out(kmers_path);
    if (!out) die("open out: " + kmers_path);
    for (const auto& seq : kmer_seqs) out << seq << '\n';
    out.close();
    std::cerr << "[INFO] wrote " << kmer_seqs.size() << " k-mers to " << kmers_path
              << " in " << tw.sec() << " s\n";
  }
}

// ---------------- main ----------------
int main(int argc, char** argv) {
  Config cfg = parse_args(argc, argv);
  Timer tall;

  Overlap ov = build_overlap(cfg.samples_file, cfg.pheno_file);
  UData   ud = load_Ubin(cfg.U_file, (int)ov.y.size(), ov.y);

  vector<std::string> kmer_seqs;
  MatRow Xrm;
  load_X_from_index(cfg, ov.idxs, (int)ov.y.size(), kmer_seqs, Xrm);

  Stats2D st = build_stats2d_self(ud.U, ov.y, Xrm, cfg.BLAS_THREADS);

  emit_pairs_self(kmer_seqs, Xrm, st, cfg);

  std::cerr << "[TIME] total " << tall.sec() << " s\n";
  return 0;
}
