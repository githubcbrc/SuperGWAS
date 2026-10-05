// kgwas.cpp — Unified 1D k-mer GWAS (replaces kgwas_cols + kgwas_rows)
// Auto-detects mode: directory input → COLS, file input → ROWS.
// Override with --mode cols|rows.
//
// Usage: kgwas samples.txt pheno.txt INPUT U.bin output.tsv [log10_thresh] [--mode cols|rows]
//
// COLS mode: INPUT is a batch directory containing bits/ and kmers.txt.gz
// ROWS mode: INPUT is a .tsv or .tsv.gz file with kmer\tbitstring rows

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <omp.h>
#include <fcntl.h>
#include <unistd.h>

extern "C" void openblas_set_num_threads(int);

#include "kgwas_utils.hpp"

using std::string;
using std::vector;
using std::size_t;
using namespace kgwas;

// ----- Types -----
struct Config {
  string samples_file, pheno_file, input_path, U_file, out_file;
  double log10_thresh = 6.0;
  double log10_thresh_2d = 15.0;
  double corr_threshold = 0.2;
  int LOAD_THREADS = 24;
  int BLAS_THREADS = 112;
  int EMIT_THREADS = 112;
  string mode; // "cols" or "rows", empty = auto-detect
};

struct HitStats {
  double effect;
  double T;
  double pval;
};

// ----- Shared stat computation -----
static inline HitStats compute_hit(double xxM, double xyM, double yyM, int n) {
  const double eps = 1e-12;
  HitStats h{0.0, 0.0, 1.0};
  if (xxM > eps && yyM > eps) {
    double r2 = (xyM * xyM) / (xxM * yyM);
    if (r2 < 0) r2 = 0;
    if (r2 >= 1) r2 = 1 - 1e-15;
    h.T = (double)n * std::log(1.0 / (1.0 - r2));
    h.pval = std::erfc(std::sqrt(0.5 * h.T));
    h.effect = xyM / (xxM * yyM);
  }
  return h;
}

// ----- Reconstruct bitstring from column matrix -----
static string reconstruct_bitstring(const MatRow& Xrm, int col_idx, int n) {
  string bits(n, '0');
  for (int i = 0; i < n; ++i)
    bits[i] = (Xrm(i, col_idx) > 0.5f) ? '1' : '0';
  return bits;
}

// ----- Arg parsing -----
static Config parse_args(int argc, char** argv) {
  if (argc < 6) {
    std::cerr
      << "Usage: " << argv[0]
      << " samples.txt phenotype.txt INPUT U.bin output.tsv [log10_thresh] [--mode cols|rows]\n"
      << "\n"
      << "  INPUT: directory (batch with bits/ + kmers.txt.gz) for COLS mode,\n"
      << "         or a .tsv[.gz] file for ROWS mode.\n"
      << "  Auto-detects mode from INPUT type unless --mode is specified.\n"
      << "\n"
      << "  Env: LOAD_THREADS, BLAS_THREADS, EMIT_THREADS, CORR_THRESHOLD, BATCH_SIZE\n";
    std::exit(1);
  }
  Config cfg;
  cfg.samples_file = argv[1];
  cfg.pheno_file   = argv[2];
  cfg.input_path   = argv[3];
  cfg.U_file       = argv[4];
  cfg.out_file     = argv[5];

  // Parse remaining positional/named args
  for (int i = 6; i < argc; ++i) {
    string a = argv[i];
    if (a == "--mode" && i + 1 < argc) {
      cfg.mode = argv[++i];
    } else if (a == "--log10-thresh-2d" && i + 1 < argc) {
      cfg.log10_thresh_2d = std::stod(argv[++i]);
    } else {
      try { cfg.log10_thresh = std::stod(a); }
      catch (...) { die("invalid argument: " + a); }
      if (cfg.log10_thresh < 0.0) die("threshold must be >= 0");
    }
  }

  // env overrides
  cfg.LOAD_THREADS = env_int("LOAD_THREADS", cfg.LOAD_THREADS);
  cfg.BLAS_THREADS = env_int("BLAS_THREADS", cfg.BLAS_THREADS);
  cfg.EMIT_THREADS = env_int("EMIT_THREADS", cfg.EMIT_THREADS);
  {
    const char* v = std::getenv("CORR_THRESHOLD");
    if (v) cfg.corr_threshold = std::stod(v);
  }
  {
    const char* v = std::getenv("LOG10_THRESH_2D_INPUT");
    if (v) cfg.log10_thresh_2d = std::stod(v);
  }

  // Auto-detect mode if not specified
  if (cfg.mode.empty()) {
    struct stat st{};
    if (::stat(cfg.input_path.c_str(), &st) != 0)
      die("cannot stat INPUT: " + cfg.input_path);
    cfg.mode = S_ISDIR(st.st_mode) ? "cols" : "rows";
    std::cerr << "[INFO] auto-detected mode: " << cfg.mode << "\n";
  }
  if (cfg.mode != "cols" && cfg.mode != "rows")
    die("--mode must be 'cols' or 'rows', got: " + cfg.mode);

  return cfg;
}

// =====================================================================
// COLS mode
// =====================================================================

struct ColsBulkStats {
  Eigen::ArrayXf xxM;
  Eigen::ArrayXf xyM;
  double yyM = 0.0;
};

static ColsBulkStats compute_cols_stats(const Eigen::MatrixXd& U,
                                        const Eigen::VectorXd& y,
                                        const MatRow& Xrm,
                                        int BLAS_THREADS) {
  Timer t;
  omp_set_num_threads(1);
  openblas_set_num_threads(BLAS_THREADS);

  Eigen::MatrixXf Uf = U.cast<float>();
  Eigen::VectorXf yf = y.cast<float>();

  Eigen::VectorXf Uyf = Uf.transpose() * yf;
  Eigen::MatrixXf UX  = Uf.transpose() * Xrm;
  Eigen::VectorXf Xy  = Xrm.transpose() * yf;
  Eigen::ArrayXf xx   = Xrm.colwise().squaredNorm().array();
  Eigen::ArrayXf sUX  = UX.colwise().squaredNorm().array();
  Eigen::VectorXf UyUX = UX.transpose() * Uyf;

  ColsBulkStats st;
  st.xyM = Xy.array() - UyUX.array();
  st.xxM = xx - sUX;

  const Eigen::VectorXd Uy = U.transpose() * y;
  st.yyM = y.squaredNorm() - Uy.squaredNorm();
  if (!(st.yyM > 1e-12)) warn("yyM ~ 0");

  std::cerr << "[INFO] LA (" << t.sec() << " s, BLAS_THREADS=" << BLAS_THREADS << ")\n";
  return st;
}

static void run_cols(const Config& cfg) {
  Timer tall;
  Overlap ov  = build_overlap(cfg.samples_file, cfg.pheno_file);
  UData ud    = load_Ubin(cfg.U_file, (int)ov.y.size(), ov.y);

  // input_path is a batch directory
  string bits_dir = cfg.input_path + "/bits";
  string kmers_gz = cfg.input_path + "/kmers.txt.gz";

  size_t N    = infer_N_from_bits(bits_dir, ov.idxs.front());
  auto kmers  = load_kmers(kmers_gz, N);

  MatRow Xrm  = load_matrix_rows(bits_dir, ov.idxs, (int)ov.y.size(), N, cfg.LOAD_THREADS);

  // --- Pearson correlation pre-filter ---
  {
    Timer tcor;
    const size_t total_kmers = kmers.size();
    const int n = (int)ov.y.size();
    const Eigen::VectorXf yf = ov.y.cast<float>();
    const float ymean = yf.mean();
    const Eigen::VectorXf yc = yf.array() - ymean;
    const float ynorm = yc.norm();
    const float thresh = (float)cfg.corr_threshold;

    std::vector<int> keep; keep.reserve(total_kmers);
    for (size_t j = 0; j < total_kmers; ++j) {
      Eigen::VectorXf xc = Xrm.col((int)j).array() - Xrm.col((int)j).mean();
      float xnorm = xc.norm();
      if (xnorm < 1e-12f || ynorm < 1e-12f) continue;
      float cor = xc.dot(yc) / (xnorm * ynorm);
      if (std::fabs(cor) >= thresh) keep.push_back((int)j);
    }

    const size_t passed = keep.size();
    std::cerr << "[INFO] CORR_THRESHOLD=" << cfg.corr_threshold
              << "  total kmers: " << total_kmers
              << ", passed correlation filter: " << passed << "\n";

    if (passed < total_kmers) {
      MatRow Xfilt(n, (int)passed);
      vector<string> kfilt(passed);
      for (size_t i = 0; i < passed; ++i) {
        Xfilt.col((int)i) = Xrm.col(keep[i]);
        kfilt[i] = std::move(kmers[(size_t)keep[i]]);
      }
      Xrm = std::move(Xfilt);
      kmers = std::move(kfilt);
    }
    std::cerr << "[INFO] Correlation filter (" << tcor.sec() << " s)\n";
  }

  ColsBulkStats st = compute_cols_stats(ud.U, ov.y, Xrm, cfg.BLAS_THREADS);

  // --- Emit hits ---
  // 1D output (3-column: kmer \t effect \t pval) — all significant hits, no bitstrings
  // 2D output (4-column: kmer \t bitstring \t effect \t pval) — only if --log10-thresh-2d set
  {
    Timer t;
    openblas_set_num_threads(1);
    omp_set_num_threads(cfg.EMIT_THREADS);

    const double p_thresh = std::pow(10.0, -cfg.log10_thresh);
    const bool do_2d = (cfg.log10_thresh_2d > 0.0);
    const double p_thresh_2d = do_2d ? std::pow(10.0, -cfg.log10_thresh_2d) : 0.0;
    const int64_t N_k = (int64_t)kmers.size();
    const int n = (int)Xrm.rows();
    const int64_t CH = 1'000'000;
    std::vector<std::string> buffers_1d((N_k + CH - 1) / CH);
    std::vector<std::string> buffers_2d(do_2d ? (N_k + CH - 1) / CH : 0);

    #pragma omp parallel for schedule(static)
    for (int64_t s = 0; s < N_k; s += CH) {
      int64_t e = std::min<int64_t>(s + CH, N_k);
      std::ostringstream oss_1d;
      oss_1d.setf(std::ios::fmtflags(0), std::ios::floatfield);
      std::ostringstream oss_2d;
      if (do_2d) oss_2d.setf(std::ios::fmtflags(0), std::ios::floatfield);

      for (int64_t b = s; b < e; ++b) {
        HitStats h = compute_hit((double)st.xxM[(int)b], (double)st.xyM[(int)b], st.yyM, n);
        if (h.pval <= p_thresh) {
          oss_1d << kmers[(size_t)b] << '\t' << h.effect << '\t' << h.pval << '\n';
          if (do_2d && h.pval <= p_thresh_2d) {
            string bits = reconstruct_bitstring(Xrm, (int)b, n);
            oss_2d << kmers[(size_t)b] << '\t' << bits << '\t' << h.effect << '\t' << h.pval << '\n';
          }
        }
      }
      buffers_1d[(size_t)(s / CH)] = std::move(oss_1d).str();
      if (do_2d) buffers_2d[(size_t)(s / CH)] = std::move(oss_2d).str();
    }

    // Write 1D output
    std::ofstream out(cfg.out_file, std::ios::binary);
    if (!out) die("open out: " + cfg.out_file);
    for (auto& s : buffers_1d) out.write(s.data(), (std::streamsize)s.size());
    out.close();
    if (!out) die("write out failed");

    std::cerr << "[INFO] wrote 1D hits at -log10(p) >= " << cfg.log10_thresh
              << " (" << t.sec() << " s, EMIT_THREADS=" << cfg.EMIT_THREADS << ")\n";

    // Write 2D output
    if (do_2d) {
      string out_2d_path = cfg.out_file;
      size_t dot = out_2d_path.rfind('.');
      if (dot != string::npos)
        out_2d_path = out_2d_path.substr(0, dot) + ".2d" + out_2d_path.substr(dot);
      else
        out_2d_path += ".2d";

      std::ofstream out2d(out_2d_path, std::ios::binary);
      if (!out2d) die("open 2D out: " + out_2d_path);
      for (auto& s : buffers_2d) out2d.write(s.data(), (std::streamsize)s.size());
      out2d.close();
      if (!out2d) die("write 2D out failed");

      std::cerr << "[INFO] wrote 2D hits at -log10(p) >= " << cfg.log10_thresh_2d
                << " to " << out_2d_path << "\n";
    }
  }

  // Write count file
  {
    std::string countfile = cfg.out_file + ".count";
    std::string shardname = cfg.out_file;
    size_t lastslash = shardname.find_last_of("/\\");
    if (lastslash != std::string::npos) shardname = shardname.substr(lastslash + 1);
    size_t lastdot = shardname.find_last_of('.');
    if (lastdot != std::string::npos) shardname = shardname.substr(0, lastdot);
    std::ofstream countout(countfile);
    if (countout.is_open()) {
      countout << shardname << "\t" << kmers.size() << "\n";
      std::cerr << "[INFO] Count written to " << countfile << "\n";
    }
  }

  std::cerr << "[TIME] total " << tall.sec() << " s\n";
}

// =====================================================================
// ROWS mode (streaming with batched BLAS)
// =====================================================================

static inline bool has_suffix(const std::string& s, const std::string& suf) {
  if (s.size() < suf.size()) return false;
  return std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

static bool gz_getline(gzFile f, std::string& out, const size_t chunk = (1u << 20)) {
  out.clear();
  if (f == nullptr) return false;
  thread_local std::unique_ptr<char[]> buf(new char[chunk]);
  bool got_any = false;
  for (;;) {
    char* res = gzgets(f, buf.get(), static_cast<int>(chunk));
    if (!res) {
      if (!got_any) return false;
      return true;
    }
    got_any = true;
    size_t len = std::strlen(res);
    if (len == 0) continue;
    if (res[len - 1] == '\n') {
      size_t end = len - 1;
      if (end > 0 && res[end - 1] == '\r') --end;
      out.append(res, res + end);
      return true;
    } else {
      out.append(res, res + len);
    }
  }
}

static void run_rows(const Config& cfg) {
  auto t0 = std::chrono::high_resolution_clock::now();
  auto elapsed = [](auto a, auto b) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
  };

  // Load genome IDs
  auto t_ids0 = std::chrono::high_resolution_clock::now();
  const auto all_genome_ids = load_lines(cfg.samples_file);
  std::cerr << "[TIME] Loaded all_genome_ids in " << elapsed(t_ids0, std::chrono::high_resolution_clock::now())
            << "s (" << all_genome_ids.size() << " IDs)\n";

  // Load phenotype and align
  auto t_ph0 = std::chrono::high_resolution_clock::now();
  const auto pheno_map = load_pheno(cfg.pheno_file);
  std::vector<int> idxs; idxs.reserve(all_genome_ids.size());
  Eigen::VectorXd y(pheno_map.size());
  int m = 0;
  for (int i = 0; i < (int)all_genome_ids.size(); ++i) {
    auto it = pheno_map.find(all_genome_ids[i]);
    if (it != pheno_map.end()) { idxs.push_back(i); y(m++) = it->second; }
  }
  y.conservativeResize(m);
  const int n = (int)idxs.size();
  if (n == 0) die("No overlapping genomes in phenotype.");
  std::cerr << "[TIME] Loaded & aligned phenotype in " << elapsed(t_ph0, std::chrono::high_resolution_clock::now())
            << "s (" << n << " genomes)\n";

  // Load U
  auto t_u0 = std::chrono::high_resolution_clock::now();
  std::ifstream uin(cfg.U_file, std::ios::binary);
  if (!uin) die("Cannot open " + cfg.U_file);
  int n_u = 0, r = 0;
  uin.read((char*)&n_u, sizeof(int));
  uin.read((char*)&r, sizeof(int));
  if (n_u != n) die("U.bin row count (" + std::to_string(n_u) + ") does not match phenotype size (" + std::to_string(n) + ")");
  Eigen::MatrixXd U(n, r);
  uin.read((char*)U.data(), sizeof(double) * n * r);
  std::cerr << "[TIME] Loaded U (" << n << " x " << r << ") in " << elapsed(t_u0, std::chrono::high_resolution_clock::now()) << "s\n";

  // Correlation pre-filter threshold
  const double CORR_THRESHOLD = cfg.corr_threshold;
  std::cerr << "[INFO] CORR_THRESHOLD=" << CORR_THRESHOLD << "\n";

  const Eigen::VectorXd y_centered = y.array() - y.mean();
  const double y_centered_norm = y_centered.norm();

  // Precompute Uy, yyM
  const Eigen::VectorXd Uy = U.transpose() * y;
  const double yyM = y.squaredNorm() - Uy.squaredNorm();
  if (!(yyM > 1e-12)) std::cerr << "[WARN] yyM ~ 0. Tests unstable.\n";

  // Open input
  const bool gz_mode = has_suffix(cfg.input_path, ".gz");
  std::ifstream kin_plain;
  gzFile kin_gz = nullptr;
  if (gz_mode) {
    kin_gz = gzopen(cfg.input_path.c_str(), "rb");
    if (kin_gz == nullptr) die("Cannot open gz file: " + cfg.input_path);
    gzbuffer(kin_gz, 1u << 20);
  } else {
    kin_plain.open(cfg.input_path);
    if (!kin_plain) die("Cannot open kmer file: " + cfg.input_path);
  }

  // Open outputs
  std::ofstream out(cfg.out_file);
  if (!out) die("Cannot open output: " + cfg.out_file);

  const bool do_2d = (cfg.log10_thresh_2d > 0.0);
  const double p_thresh_2d = do_2d ? std::pow(10.0, -cfg.log10_thresh_2d) : 0.0;
  std::ofstream out2d;
  string out_2d_path;
  if (do_2d) {
    out_2d_path = cfg.out_file;
    size_t dot = out_2d_path.rfind('.');
    if (dot != string::npos)
      out_2d_path = out_2d_path.substr(0, dot) + ".2d" + out_2d_path.substr(dot);
    else
      out_2d_path += ".2d";
    out2d.open(out_2d_path);
    if (!out2d) die("Cannot open 2D output: " + out_2d_path);
  }

  const double p_thresh = std::pow(10.0, -cfg.log10_thresh);
  const double eps = 1e-12;
  size_t line_count = 0;
  size_t passed_filter = 0;
  auto t_start_scan = std::chrono::high_resolution_clock::now();

  // Stream k-mers
  std::string line;
  while (true) {
    bool ok = false;
    if (gz_mode) ok = gz_getline(kin_gz, line);
    else ok = static_cast<bool>(std::getline(kin_plain, line));
    if (!ok) break;
    if (line.empty()) continue;
    ++line_count;

    std::string kmer, bits;
    {
      size_t tab = line.find('\t');
      if (tab == std::string::npos) {
        tab = line.find(' ');
        if (tab == std::string::npos) continue;
      }
      kmer = line.substr(0, tab);
      bits = line.substr(tab + 1);
    }

    // Build x on the phenotype subset
    Eigen::VectorXd x(n);
    for (int j = 0; j < n; ++j) {
      const char c = bits[idxs[j]];
      x(j) = (c == '1') ? 1.0 : 0.0;
    }

    // Raw Pearson correlation pre-filter
    {
      const Eigen::VectorXd x_centered = x.array() - x.mean();
      const double x_centered_norm = x_centered.norm();
      if (x_centered_norm < eps || y_centered_norm < eps) continue;
      const double cor = x_centered.dot(y_centered) / (x_centered_norm * y_centered_norm);
      if (std::fabs(cor) < CORR_THRESHOLD) continue;
    }
    ++passed_filter;

    // Low-rank trick
    const Eigen::VectorXd Ux = U.transpose() * x;
    const double xxM = x.squaredNorm() - Ux.squaredNorm();
    const double xyM = x.dot(y) - Ux.dot(Uy);

    HitStats h = compute_hit(xxM, xyM, yyM, n);
    if (h.pval <= p_thresh) {
      // 1D output: 3-column, no bitstring
      out << kmer << '\t' << h.effect << '\t' << h.pval << '\n';

      // 2D output: 4-column with bitstring, stricter threshold
      if (do_2d && h.pval <= p_thresh_2d) {
        string obits(n, '0');
        for (int j = 0; j < n; ++j)
          obits[j] = bits[idxs[j]];
        out2d << kmer << '\t' << obits << '\t' << h.effect << '\t' << h.pval << '\n';
      }
    }

    if ((line_count % 1000000) == 0) {
      std::cerr << "[INFO] Processed " << line_count << " kmers in "
                << elapsed(t_start_scan, std::chrono::high_resolution_clock::now()) << "s\n";
    }
  }

  if (gz_mode) gzclose(kin_gz);

  std::cerr << "[INFO] Total kmers: " << line_count
            << ", passed correlation filter: " << passed_filter << "\n";
  std::cerr << "[TIME] Total scan " << elapsed(t_start_scan, std::chrono::high_resolution_clock::now()) << "s\n";
  if (do_2d) {
    out2d.close();
    std::cerr << "[INFO] 2D hits written to " << out_2d_path << "\n";
  }

  // Write count file
  {
    std::string countfile = cfg.out_file + ".count";
    std::string shardname = cfg.out_file;
    size_t lastslash = shardname.find_last_of("/\\");
    if (lastslash != std::string::npos) shardname = shardname.substr(lastslash + 1);
    size_t lastdot = shardname.find_last_of('.');
    if (lastdot != std::string::npos) shardname = shardname.substr(0, lastdot);
    std::ofstream countout(countfile);
    if (countout.is_open()) {
      countout << shardname << "\t" << passed_filter << "\n";
      std::cerr << "[INFO] Count written to " << countfile << "\n";
    }
  }

  std::cerr << "[TIME] End-to-end " << elapsed(t0, std::chrono::high_resolution_clock::now()) << "s\n";
}

// =====================================================================
// main — dispatch to cols or rows
// =====================================================================
int main(int argc, char** argv) {
  Config cfg = parse_args(argc, argv);
  if (cfg.mode == "cols") run_cols(cfg);
  else                    run_rows(cfg);
  return 0;
}
