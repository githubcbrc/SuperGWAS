#include <bits/stdc++.h>
#include <filesystem>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#include "map_utils.hpp" 

struct Prec {
    static constexpr int EFFECT = 9;
    static constexpr int PVAL   = 9;
};

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 5) {
        cerr <<
        "Map a multi-FASTA reference against a k-mer significance .umap index.\n"
        "Writes one file per contig in <out_dir>, each line: pos\\teffect\\tp_value (no header). Skips empty.\n\n"
        "Usage:\n"
        "  " << argv[0] << " <index.umap> <reference.fa> <k> <out_dir>\n"
        "    [--threads N] [--dual] [--min-logp X] [--chunk-size M]\n\n"
        "Options:\n"
        "  --threads N     OpenMP threads (default: OMP max setting)\n"
        "  --dual          Lookup both k-mer and reverse complement (non-canonical)\n"
        "  --min-logp X    Keep hits with -log10(p) >= X (i.e., p <= 10^-X). Default 0 = no filter\n"
        "  --chunk-size M  Split long contigs into ~M-sized k-mer start blocks for intra-contig parallelism\n";
        return 1;
    }

    const string umap_path  = argv[1];
    const string ref_path   = argv[2];
    const size_t K          = (size_t)stoull(argv[3]);
    const string out_dir    = argv[4];

    int    threads     = 0;
    bool   dual_lookup = false;
    double min_logp    = 0.0;
    size_t chunk_size  = 0; // 0 => no intra-contig splitting

    for (int i=5; i<argc; ++i) {
        string a = argv[i];
        if (a == "--threads" && i+1 < argc)      threads = stoi(argv[++i]);
        else if (a == "--dual")                  dual_lookup = true;
        else if (a == "--min-logp" && i+1 < argc) min_logp = stod(argv[++i]);
        else if (a == "--chunk-size" && i+1 < argc) chunk_size = (size_t)stoull(argv[++i]);
        else cerr << "[WARN] Unknown/ignored arg: " << a << "\n";
    }

    std::filesystem::create_directories(out_dir);

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    const int used_threads = (threads>0 ? threads : omp_get_max_threads());
#else
    const int used_threads = 1;
    if (threads>0) cerr << "[WARN] OpenMP not enabled; running single-threaded.\n";
#endif

    // Threshold p <= 10^-min_logp
    const double p_thresh = (min_logp <= 0.0) ? 1.0 : pow(10.0, -min_logp);

    // Load index (double precision)
    KmerMap db;
    {
        auto t0 = chrono::high_resolution_clock::now();
        if (!load_umap(umap_path, db)) return 2;
        auto t1 = chrono::high_resolution_clock::now();
        cerr << "[INFO] Loaded index: " << db.size() << " kmers in "
             << chrono::duration<double>(t1-t0).count() << " s; threads=" << used_threads
             << " | min-logp (for mapping only) =" << min_logp << " (p<=" << std::scientific << p_thresh << std::defaultfloat << ")\n";
    }

    // Read reference contigs
    vector<Contig> contigs;
    {
        auto t0 = chrono::high_resolution_clock::now();
        if (!read_reference_fasta(ref_path, contigs)) return 3;
        auto t1 = chrono::high_resolution_clock::now();
        size_t total_bp = 0; for (auto &c: contigs) total_bp += c.seq.size();
        cerr << "[INFO] Loaded reference: " << contigs.size() << " contigs, "
             << total_bp << " bp in " << chrono::duration<double>(t1-t0).count() << " s\n";
    }
    if (contigs.empty()) { cerr << "[WARN] No contigs found.\n"; return 0; }

    auto T0 = chrono::high_resolution_clock::now();

#ifdef _OPENMP
    // Use OpenMP tasks so big contigs can be split and processed cooperatively
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for (size_t cidx = 0; cidx < contigs.size(); ++cidx) {
                #pragma omp task firstprivate(cidx) shared(contigs, db)
                {
                    const string& name = contigs[cidx].name;
                    const string& seq  = contigs[cidx].seq;
                    const string out_path = (std::filesystem::path(out_dir) / (name + ".tsv")).string();

                    if (seq.size() < K) {
                    } else {
                        const size_t max_start = seq.size() - K + 1;

                        auto emit_block = [&](size_t beg, size_t end, std::string& buf, size_t& hits){
                            std::ostringstream local;
                            local.setf(std::ios::fmtflags(0), std::ios::floatfield);
                            for (size_t i = beg; i < end; ++i) {
                                if (has_N(seq, i, K)) continue;

                                string kmer(seq.data() + i, K);
                                string rc = kmer;
                                revcomp_inplace(rc);

                                const Entry* val = nullptr;
                                if (!dual_lookup) {
                                    string key = kmer;
                                    if (rc < key) key.swap(rc);
                                    auto it = db.find(key);
                                    if (it != db.end()) val = &it->second;
                                } else {
                                    auto it = db.find(kmer);
                                    if (it == db.end()) it = db.find(rc);
                                    if (it != db.end()) val = &it->second;
                                }
                                if (!val) continue;
                                if (val->pval > p_thresh) continue;

                                local << i << '\t'
                                      << std::setprecision(Prec::EFFECT) << std::defaultfloat << val->effect << '\t'
                                      << std::setprecision(Prec::PVAL)   << std::scientific   << val->pval   << '\n';
                                ++hits;
                            }
                            buf = std::move(local).str();
                        };

                        if (chunk_size == 0 || max_start <= chunk_size) {
                            std::string buf; size_t hits=0;
                            emit_block(0, max_start, buf, hits);
                            if (hits > 0) {
                                ofstream out(out_path);
                                if (out.is_open()) out << buf;
                                else {
                                    #pragma omp critical
                                    { cerr << "[WARN] cannot open " << out_path << " for write\n"; }
                                }
                            }
                        } else {
                            const size_t nb = (max_start + chunk_size - 1) / chunk_size;
                            vector<string> block_buf(nb);
                            vector<size_t> block_hits(nb, 0);

                            for (size_t b = 0; b < nb; ++b) {
                                #pragma omp task firstprivate(b) shared(block_buf, block_hits, seq, db)
                                {
                                    const size_t beg = b * chunk_size;
                                    const size_t end = std::min(max_start, beg + chunk_size);
                                    emit_block(beg, end, block_buf[b], block_hits[b]);
                                }
                            }
                            #pragma omp taskwait

                            size_t total_hits = 0;
                            for (auto h : block_hits) total_hits += h;
                            if (total_hits > 0) {
                                ofstream out(out_path);
                                if (out.is_open()) {
                                    for (auto &bb : block_buf) out << bb;
                                } else {
                                    #pragma omp critical
                                    { cerr << "[WARN] cannot open " << out_path << " for write\n"; }
                                }
                            }
                        }
                    }
                } // contig task
            } // for contigs
        } // single
    } // parallel
#else
    // No OpenMP: sequential
    for (const auto& c : contigs) {
        const string out_path = (std::filesystem::path(out_dir) / (c.name + ".tsv")).string();
        if (c.seq.size() < K) continue;
        const size_t max_start = c.seq.size() - K + 1;

        std::ostringstream local;
        local.setf(std::ios::fmtflags(0), std::ios::floatfield);
        size_t hits = 0;
        for (size_t i = 0; i < max_start; ++i) {
            if (has_N(c.seq, i, K)) continue;

            string kmer(c.seq.data() + i, K);
            string rc = kmer; revcomp_inplace(rc);

            const Entry* val = nullptr;
            if (!dual_lookup) {
                string key = kmer;
                if (rc < key) key.swap(rc);
                auto it = db.find(key);
                if (it != db.end()) val = &it->second;
            } else {
                auto it = db.find(kmer);
                if (it == db.end()) it = db.find(rc);
                if (it != db.end()) val = &it->second;
            }
            if (!val) continue;
            if (val->pval > p_thresh) continue;

            local << i << '\t'
                  << std::setprecision(Prec::EFFECT) << std::defaultfloat << val->effect << '\t'
                  << std::setprecision(Prec::PVAL)   << std::scientific   << val->pval   << '\n';
            ++hits;
        }
        if (hits > 0) {
            ofstream out(out_path);
            if (out.is_open()) out << local.str();
            else cerr << "[WARN] cannot open " << out_path << " for write\n";
        }
    }
#endif

    auto T1 = chrono::high_resolution_clock::now();
    cerr << "[INFO] Completed mapping in "
         << chrono::duration<double>(T1 - T0).count() << " s\n";
    return 0;
}

