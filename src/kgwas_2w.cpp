// kgwas_2w.cpp — Two-genome all-pairs sweep with grid-structured LRT
//
// Usage: kgwas_2w [--interaction_only]
//        U.txt Th_cols.txt Tp_cols.txt phenotype.txt host_matrix.txt pathogen_matrix.txt out.txt

#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>

#include "chisq.hpp"
#include "io_2w.hpp"
#include "la_2w.hpp"
#include "stats_2w.hpp"

using la::Mat;
using la::Vec;
using la::Row;

namespace prog {
using clock_t = std::chrono::steady_clock;

inline void print(std::size_t done, std::size_t total, double elapsed_s) {
    const double rate  = done  / std::max(elapsed_s, 1e-12);
    const double eta_s = (total > done) ? (total - done) / std::max(rate, 1e-12) : 0.0;
    const int pct = total ? int((100.0 * done) / total + 0.5) : 100;
    std::cerr << "\r[progress] " << std::setw(3) << pct << "% ("
              << done << "/" << total << ")  "
              << "elapsed " << std::fixed << std::setprecision(1) << elapsed_s << "s  "
              << "rate "    << std::setprecision(1) << rate   << "/s  "
              << "ETA "     << std::setprecision(1) << eta_s  << "s" << std::flush;
}
} // namespace prog

int main(int argc, char** argv) {
    bool interaction_only = false;
    std::vector<std::string> pos; pos.reserve(argc);

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--interaction_only") {
            interaction_only = true;
        } else {
            pos.push_back(std::move(a));
        }
    }

    if (pos.size() != 7) {
        std::cerr <<
        "Usage:\n"
        "  " << argv[0] << " [--interaction_only]\n"
        "        U.txt Th_cols.txt Tp_cols.txt phenotype.txt host_matrix.txt pathogen_matrix.txt out.txt\n"
        "Notes:\n"
        "  - Phenotype file lines: <host_id> <path_id> <y>\n"
        "  - host_matrix / pathogen_matrix are BITSTRING files with lines:\n"
        "        <label><space or tab><bitstring-of-0/1>\n"
        "  - Th_cols.txt: r floats per line, one line per host k-mer (mh lines)\n"
        "  - Tp_cols.txt: r floats per line, one line per patho k-mer (mp lines)\n";
        return 1;
    }

    const std::string U_path     = pos[0];
    const std::string Th_path    = pos[1];
    const std::string Tp_path    = pos[2];
    const std::string phen_path  = pos[3];
    const std::string host_path  = pos[4];
    const std::string patho_path = pos[5];
    const std::string out_path   = pos[6];

    try {
        Mat U = io::read_matrix_txt(U_path);
        const int n = static_cast<int>(U.rows());
        const int r = static_cast<int>(U.cols());

        std::vector<int> host_id, path_id; host_id.reserve(n); path_id.reserve(n);
        Vec y(n);
        {
            std::ifstream in(phen_path);
            if (!in) throw std::runtime_error("Cannot open: " + phen_path);
            std::string line; long long ih, ip; double v; int t = 0;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream ss(line);
                if (ss >> ih >> ip >> v) {
                    host_id.push_back(static_cast<int>(ih));
                    path_id.push_back(static_cast<int>(ip));
                    if (t < n) y[t] = v;
                    ++t;
                }
            }
            if (static_cast<int>(host_id.size()) != n)
                throw std::runtime_error("phenotype rows != U.rows()");
        }

        int nh = 0, np = 0;
        for (int h : host_id) nh = std::max(nh, h + 1);
        for (int p : path_id) np = std::max(np, p + 1);
        if (n != nh * np) throw std::runtime_error("Full grid assumption failed: n != nh*np");

        auto Hfull = io::load_kmer_matrix_full(host_path);
        auto Pfull = io::load_kmer_matrix_full(patho_path);
        if (Hfull.G.rows() != nh) throw std::runtime_error("host_matrix rows != nh");
        if (Pfull.G.rows() != np) throw std::runtime_error("pathogen_matrix rows != np");
        const std::size_t mh = static_cast<std::size_t>(Hfull.G.cols());
        const std::size_t mp = static_cast<std::size_t>(Pfull.G.cols());

        const Mat& Gch = Hfull.Gc;
        const Mat& Gcp = Pfull.Gc;

        Mat Th = io::load_T_cols(Th_path, r, mh);
        Mat Tp = io::load_T_cols(Tp_path, r, mp);

        FullGridCaches C = build_fullgrid_caches(U, host_id, path_id, y, nh, np);
        C.Th = Th;
        C.Tp = Tp;

        std::ofstream out(out_path);
        if (!out) throw std::runtime_error("Cannot write: " + out_path);
        out.setf(std::ios::fmtflags(0), std::ios::floatfield);
        out.precision(10);

        const std::size_t total_pairs = mh * mp;
        std::atomic<std::size_t> done{0};
        const auto t0 = prog::clock_t::now();

        Eigen::setNbThreads(1);
        constexpr std::size_t FLUSH_EVERY = 4096;
        constexpr int OMP_CHUNK = 64;

        #pragma omp parallel
        {
            std::ostringstream local; local.setf(std::ios::fmtflags(0), std::ios::floatfield);
            local.precision(10);
            std::size_t local_count = 0;

            Vec gh(nh), gp(np), th_cache, tp_cache;

            #pragma omp for schedule(dynamic, OMP_CHUNK) collapse(2)
            for (std::size_t a = 0; a < mh; ++a) {
                for (std::size_t b = 0; b < mp; ++b) {
                    gh = Gch.col(static_cast<Eigen::Index>(a));
                    gp = Gcp.col(static_cast<Eigen::Index>(b));

                    th_cache = C.Th.col(static_cast<Eigen::Index>(a));
                    tp_cache = C.Tp.col(static_cast<Eigen::Index>(b));

                    const std::string& hlabel = Hfull.labels[a];
                    const std::string& plabel = Pfull.labels[b];

                    StatsOut stats = interaction_only
                        ? interaction_only_stat_full(C, gh, gp)
                        : joint_3df_stat_full(C, gh, gp, th_cache, tp_cache);

                    local << a << ' ' << b << ' ' << hlabel << ' ' << plabel << ' '
                          << stats.ratio << ' ' << stats.delta << ' '
                          << stats.LRT << ' ' << stats.pchi << '\n';

                    if (++local_count % FLUSH_EVERY == 0) {
                        #pragma omp critical(out_flush)
                        { out << local.str(); }
                        local.str(std::string()); local.clear();
                    }

                    const std::size_t d = done.fetch_add(1, std::memory_order_relaxed) + 1;
                    if ((d % 1000ull == 0ull) || d == total_pairs) {
                        const auto now = prog::clock_t::now();
                        const double elapsed_s =
                            std::chrono::duration<double>(now - t0).count();
                        #pragma omp critical(progress)
                        { prog::print(d, total_pairs, elapsed_s); }
                    }
                }
            }
            #pragma omp critical(out_flush)
            { out << local.str(); }
        }

        const auto t1 = prog::clock_t::now();
        const double total_s = std::chrono::duration<double>(t1 - t0).count();
        const double rate = total_pairs / std::max(total_s, 1e-12);
        std::cerr << "\n[done] wrote " << out_path
                  << " | pairs=" << total_pairs
                  << " | total=" << std::fixed << std::setprecision(3) << total_s << "s"
                  << " | rate="  << std::setprecision(1)  << rate << "/s"
                  << " | threads=" << omp_get_max_threads() << "\n";

    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
