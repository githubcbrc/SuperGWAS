// sort_matrix.cpp
// - Parallel load *.tsv.gz
// - Parallel LSD radix sort by (bitstring bytes || k-mer bytes)
// - Emit: kmers.txt (one per line) + bits/col_XXXXXX.txt (each is ONE line of '0'/'1')

#include <bits/stdc++.h>
#include <zlib.h>
#include <filesystem>
#include <thread>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace fs = std::filesystem;
using std::string; using std::vector; using std::size_t; using std::uint32_t;

[[noreturn]] static void die(const string& m){ std::cerr << "[ERROR] " << m << "\n"; std::exit(1); }
static void warn(const string& m){ std::cerr << "[WARN] " << m << "\n"; }

// K-mer length: auto-detected from first k-mer
static size_t KMER_LEN = 0;

struct GzReader {
    gzFile f = nullptr;
    explicit GzReader(const string& path){
        f = gzopen(path.c_str(), "rb");
        if(!f) die("gzopen(read): " + path);
        gzbuffer(f, 1<<20);
    }
    ~GzReader(){ if(f) gzclose(f); }
    bool getline(string& out){
        out.clear();
        int c;
        while(true){
            c = gzgetc(f);
            if(c == -1){
                int err = 0; gzerror(f, &err);
                return !out.empty();
            }
            if(c == '\n') return true;
            out.push_back((char)c);
        }
    }
};

struct Row { string bit; string kmer; };

struct LoadResult {
    vector<Row> rows;
    size_t bit_len = SIZE_MAX;
};

static LoadResult read_all_rows_parallel(const vector<fs::path>& files, int threads){
    threads = std::max(1, std::min<int>(threads, (int)files.size()));
    std::atomic<size_t> next{0};
    vector<vector<Row>> buckets(threads);
    vector<size_t> seen_len(threads, SIZE_MAX);

    auto worker = [&](int tid){
        auto& local = buckets[tid];
        local.reserve(200000);
        size_t local_len = SIZE_MAX;

        for(;;){
            size_t idx = next.fetch_add(1, std::memory_order_relaxed);
            if(idx >= files.size()) break;
            const auto& p = files[idx];

            GzReader in(p.string());
            string line; size_t ln=0;
            while(in.getline(line)){
                ++ln; if(line.empty()) continue;
                size_t t = line.find('\t');
                if(t == string::npos){
                    t = line.find(' ');
                    if(t == string::npos){ warn(p.filename().string()+": no delimiter on line "+std::to_string(ln)); continue; }
                }
                string kmer = line.substr(0, t);
                string bits = line.substr(t+1);
                while(!bits.empty() && isspace((unsigned char)bits.back())) bits.pop_back();
                // Auto-detect K from first k-mer
                if(KMER_LEN == 0){
                    KMER_LEN = kmer.size();
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
                if(kmer.size() != KMER_LEN){ warn("bad kmer length at "+p.filename().string()+":"+std::to_string(ln)); continue; }
                if(local_len == SIZE_MAX) local_len = bits.size();
                if(bits.size() != local_len) die("inconsistent bitstring length inside "+p.filename().string());
                for(char c: bits) if(c!='0' && c!='1') die("non-binary bit at "+p.filename().string()+":"+std::to_string(ln));
                local.push_back(Row{std::move(bits), std::move(kmer)});
            }
            seen_len[tid] = (local_len==SIZE_MAX ? seen_len[tid] : local_len);
        }
    };

    vector<std::thread> pool; pool.reserve(threads);
    for(int t=0;t<threads;t++) pool.emplace_back(worker, t);
    for(auto& th: pool) th.join();

    size_t L = SIZE_MAX;
    for(size_t s: seen_len) if(s != SIZE_MAX){
        if(L == SIZE_MAX) L = s; else if(L != s) die("bitstring length mismatch across files");
    }
    if(L == SIZE_MAX) die("failed to infer bitstring length");

    LoadResult out;
    out.bit_len = L;
    size_t total = 0; for(auto& v: buckets) total += v.size();
    out.rows.reserve(total);
    for(auto& v: buckets){
        for(const auto& r : v) if(r.bit.size()!=L) die("bitstring length mismatch");
        out.rows.insert(out.rows.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
        vector<Row>().swap(v);
    }
    return out;
}

// ---------------- key blobs for radix (0/1 bytes for bits, raw kmer bytes) ----------------
struct KeyBlobs {
    std::unique_ptr<unsigned char[]> bits_blob; // N*L
    std::unique_ptr<unsigned char[]> kmer_blob; // N*KMER_LEN
    size_t N=0, L=0;
    static KeyBlobs make_from_rows(const vector<Row>& rows, size_t L){
        KeyBlobs kb; kb.N = rows.size(); kb.L = L;
        if(kb.N==0) return kb;
        kb.bits_blob.reset(new unsigned char[kb.N*L]);
        kb.kmer_blob.reset(new unsigned char[kb.N*KMER_LEN]);
        #pragma omp parallel for schedule(static)
        for (ptrdiff_t i=0;i<(ptrdiff_t)kb.N;++i){
            const auto& b = rows[(size_t)i].bit;
            const auto& k = rows[(size_t)i].kmer;
            auto* bp = kb.bits_blob.get() + (size_t)i*L;
            for (size_t j=0;j<L;++j) bp[j] = (unsigned char)(b[j] & 1u);
            memcpy(kb.kmer_blob.get() + (size_t)i*KMER_LEN, k.data(), KMER_LEN);
        }
        return kb;
    }
    inline unsigned char key_byte(size_t row, size_t pos, size_t K) const {
        if (pos < L) return bits_blob[row*L + pos];
        size_t off = pos - L; return kmer_blob[row*KMER_LEN + off];
    }
};

// ---------------- parallel LSD radix-256 on (bits||kmer) ----------------
static void radix_sort_rows_by_combined_key(vector<size_t>& idx, const KeyBlobs& kb){
    const size_t N = idx.size(); if(N<=1) return;
    const size_t K = kb.L + KMER_LEN;
    vector<size_t> out(N);
    const int T = std::max(1, (int)std::thread::hardware_concurrency());

    for (size_t pass=0; pass<K; ++pass){
        size_t byte_pos = K - 1 - pass;

        vector<std::array<uint32_t,256>> H(T);
        #pragma omp parallel num_threads(T)
        {
            int t=0; 
	    #ifdef _OPENMP
            t = omp_get_thread_num();
            #endif
            auto& h = H[t]; h.fill(0);
            size_t chunk = (N + T - 1)/T;
            size_t a = (size_t)t*chunk, b = std::min(N, a+chunk);
            for(size_t i=a;i<b;++i){
                unsigned char c = kb.key_byte(idx[i], byte_pos, K);
                h[c]++;
            }
        }

        std::array<uint32_t,256> G{}; G.fill(0);
        for(int t=0;t<T;++t) for(int c=0;c<256;++c) G[c]+=H[t][c];

        std::array<uint32_t,256> pref{}; pref[0]=0;
        for(int c=1;c<256;++c) pref[c]=pref[c-1]+G[c-1];

        vector<std::array<uint32_t,256>> S(T);
        std::array<uint32_t,256> acc{}; acc.fill(0);
        for(int t=0;t<T;++t){
            for(int c=0;c<256;++c){ S[t][c]=pref[c]+acc[c]; acc[c]+=H[t][c]; }
        }

        #pragma omp parallel num_threads(T)
        {
            int t=0; 
	    #ifdef _OPENMP
            t = omp_get_thread_num();
            #endif
            auto offs = S[t];
            size_t chunk = (N + T - 1)/T;
            size_t a = (size_t)t*chunk, b = std::min(N, a+chunk);
            for(size_t i=a;i<b;++i){
                unsigned char c = kb.key_byte(idx[i], byte_pos, K);
                out[offs[c]++] = idx[i];
            }
        }
        idx.swap(out);
    }
}

// ---------------- emit kmers ----------------
static void write_kmers_plain(const fs::path& out_dir, const vector<string>& kmers){
    fs::create_directories(out_dir);
    std::ofstream out(out_dir / "kmers.txt", std::ios::binary);
    if(!out) die("open kmers.txt");
    for(const auto& k: kmers){ out.write(k.data(), (std::streamsize)k.size()); out.put('\n'); }
    out.close(); if(!out) die("write kmers.txt");
}

// ---------------- write columns as ONE-LINE "0/1" strings ----------------

// Full-buffer mode: allocate N bytes per column; single write per column.
static void write_columns_text_fullbuffer(const fs::path& bits_dir,
                                          const vector<Row>& rows,
                                          size_t N, size_t L,
                                          int jobs)
{
    fs::create_directories(bits_dir);
    vector<const char*> rowp(N);
    for(size_t i=0;i<N;++i) rowp[i] = rows[i].bit.data();

    std::atomic<size_t> next{0};
    auto worker = [&](){
        static thread_local vector<char> tbuf(1<<20);
        for(;;){
            size_t j = next.fetch_add(1, std::memory_order_relaxed);
            if(j >= L) break;
            // Allocate and fill the column line
            std::unique_ptr<char[]> col(new char[N+1]);
            char* dst = col.get();
            for(size_t i=0;i<N;++i){
                dst[i] = char( (rowp[i][j] & 1u) + '0' ); // branchless
            }
            dst[N] = '\n';

            char name[64]; std::snprintf(name, sizeof(name), "col_%06zu.txt", j);
            fs::path path = bits_dir / name;

            std::ofstream out(path, std::ios::binary);
            if(!out){ std::cerr << "[ERROR] open " << path << "\n"; std::exit(1); }
            out.rdbuf()->pubsetbuf(tbuf.data(), (std::streamsize)tbuf.size());
            out.write(col.get(), (std::streamsize)(N+1));
            out.close(); if(!out){ std::cerr << "[ERROR] write " << path << "\n"; std::exit(1); }
        }
    };

    jobs = std::max(1, jobs);
    vector<std::thread> pool; pool.reserve(jobs);
    for(int t=0;t<jobs;++t) pool.emplace_back(worker);
    for(auto& th: pool) th.join();
}


// ---------------- main ----------------
int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if(argc != 2 && argc != 3){
        std::cerr << "Usage: " << argv[0] << " <chunks_dir> [out_dir]\n";
        return 2;
    }
    fs::path chunks_dir = argv[1];
    fs::path out_dir = (argc==3) ? fs::path(argv[2]) : fs::path("sorted_cols_text");
    if(!fs::is_directory(chunks_dir)) die("not a directory: " + chunks_dir.string());
    fs::create_directories(out_dir);
    fs::create_directories(out_dir / "bits");

    // collect input files
    vector<fs::path> gz_files;
    for(const auto& e : fs::directory_iterator(chunks_dir)){
        if(e.is_regular_file()){
            auto p = e.path();
            if(p.extension()==".gz" && p.filename().string().find(".tsv.gz") != string::npos)
                gz_files.push_back(p);
        }
    }
    if(gz_files.empty()) die("no *.tsv.gz in " + chunks_dir.string());

    // 1) load
    int thr = std::max(1u, std::thread::hardware_concurrency());
    auto loaded = read_all_rows_parallel(gz_files, thr);
    auto& rows = loaded.rows;
    const size_t L = loaded.bit_len;
    const size_t N = rows.size();
    std::cerr << "[INFO] Loaded rows=" << N << "  bit_length="<< L << "  threads="<< thr << "\n";
    if(N==0){ std::cerr << "[INFO] nothing to do\n"; return 0; }

    // 2) radix key blobs
    std::cerr << "[INFO] Building key blobs...\n";
    struct KeyBlobs {
        std::unique_ptr<unsigned char[]> bits_blob, kmer_blob; size_t N=0,L=0;
        static KeyBlobs make(const vector<Row>& rows, size_t L){
            KeyBlobs kb; kb.N=rows.size(); kb.L=L;
            kb.bits_blob.reset(new unsigned char[kb.N*L]);
            kb.kmer_blob.reset(new unsigned char[kb.N*KMER_LEN]);
            #pragma omp parallel for schedule(static)
            for (ptrdiff_t i=0;i<(ptrdiff_t)kb.N;++i){
                const auto& b = rows[(size_t)i].bit;
                const auto& k = rows[(size_t)i].kmer;
                auto* bp = kb.bits_blob.get() + (size_t)i*L;
                for(size_t j=0;j<L;++j) bp[j] = (unsigned char)(b[j] & 1u);
                memcpy(kb.kmer_blob.get() + (size_t)i*KMER_LEN, k.data(), KMER_LEN);
            }
            return kb;
        }
        inline unsigned char key_byte(size_t row, size_t pos, size_t K) const {
            if (pos < L) return bits_blob[row*L + pos];
            size_t off = pos - L; return kmer_blob[row*KMER_LEN + off];
        }
    } kb = KeyBlobs::make(rows, L);

    auto radix_sort_rows_by_combined_key = [&](vector<size_t>& idx){
        const size_t K = L + KMER_LEN, Nloc = idx.size();
        vector<size_t> out(Nloc);
        const int T = std::max(1, (int)std::thread::hardware_concurrency());
        for (size_t pass=0; pass<K; ++pass){
            size_t byte_pos = K - 1 - pass;
            vector<std::array<uint32_t,256>> H(T);
            #pragma omp parallel num_threads(T)
            {
                int t=0; 
		#ifdef _OPENMP
                t = omp_get_thread_num();
                #endif
                auto& h = H[t]; h.fill(0);
                size_t chunk = (Nloc + T - 1)/T;
                size_t a = (size_t)t*chunk, b = std::min(Nloc, a+chunk);
                for(size_t i=a;i<b;++i){
                    unsigned char c = kb.key_byte(idx[i], byte_pos, K);
                    h[c]++;
                }
            }
            std::array<uint32_t,256> G{}; G.fill(0);
            for(int t=0;t<T;++t) for(int c=0;c<256;++c) G[c]+=H[t][c];
            std::array<uint32_t,256> pref{}; pref[0]=0;
            for(int c=1;c<256;++c) pref[c]=pref[c-1]+G[c-1];
            vector<std::array<uint32_t,256>> S(T);
            std::array<uint32_t,256> acc{}; acc.fill(0);
            for(int t=0;t<T;++t){ for(int c=0;c<256;++c){ S[t][c]=pref[c]+acc[c]; acc[c]+=H[t][c]; } }
            #pragma omp parallel num_threads(T)
            {
                int t=0; 
		#ifdef _OPENMP
                t = omp_get_thread_num();
                #endif
                auto offs = S[t];
                size_t chunk = (Nloc + T - 1)/T;
                size_t a = (size_t)t*chunk, b = std::min(Nloc, a+chunk);
                for(size_t i=a;i<b;++i){
                    unsigned char c = kb.key_byte(idx[i], byte_pos, K);
                    out[offs[c]++] = idx[i];
                }
            }
            idx.swap(out);
        }
    };

    std::cerr << "[INFO] Radix sorting by (bits||kmer)...\n";
    vector<size_t> idx(N); for(size_t i=0;i<N;++i) idx[i]=i;
    radix_sort_rows_by_combined_key(idx);

    std::cerr << "[INFO] Reordering rows...\n";
    vector<Row> sorted(N);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t r=0;r<(ptrdiff_t)N;++r) sorted[(size_t)r] = std::move(rows[idx[(size_t)r]]);
    vector<Row>().swap(rows);

    std::cerr << "[INFO] Writing kmers.txt...\n";
    vector<string> kmers(N);
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i=0;i<(ptrdiff_t)N;++i) kmers[(size_t)i] = std::move(sorted[(size_t)i].kmer);
    write_kmers_plain(out_dir, kmers);

    std::cerr << "[INFO] Writing " << L << " column bitstrings (one line each) in parallel...\n";
    int jobs = std::max(1u, std::thread::hardware_concurrency());
    write_columns_text_fullbuffer(out_dir / "bits", sorted, N, L, jobs);

    // meta
    std::ofstream meta(out_dir / "meta.txt");
    if(meta){
        meta << "format=column_text_bitstrings\n";
        meta << "num_rows=" << N << "\n";
        meta << "num_columns=" << L << "\n";
        meta << "column_file_pattern=bits/col_%06d.txt\n";
        meta << "each_column=single_line_of_0_or_1_followed_by_newline\n";
        meta << "sorted_by=bits_then_kmer_LSD_radix256\n";
    }

    std::cerr << "[INFO] Done. Output at: " << out_dir << "\n";
    return 0;
}

