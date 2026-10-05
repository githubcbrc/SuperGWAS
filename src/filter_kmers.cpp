// sample_matrix.cpp
// Read sampled k-mers, filter by MAF on the phenotyped panel, subsample,
// and output a matrix suitable for precompute.
//
// Usage: sample_matrix <sampled_kmers> <samples_file> <pheno_file> <output_tsv>
// Env:   SUBSAMPLE_SIZE (default 200000), MAF_THRESHOLD (default 5, percent), SEED (default 42)

#include <bits/stdc++.h>
#include <zlib.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using std::string; using std::vector; using std::size_t;

[[noreturn]] static void die(const string& m){ std::cerr << "[ERROR] " << m << "\n"; std::exit(1); }
static void warn(const string& m){ std::cerr << "[WARN] " << m << "\n"; }

static size_t env_size(const char* name, size_t def){
    const char* v = std::getenv(name);
    return v ? (size_t)std::stoull(v) : def;
}

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

struct Record {
    string kmer;
    vector<char> bits;  // only phenotyped positions, '0'/'1'
};

int main(int argc, char** argv){
    if(argc < 5){
        std::cerr << "Usage: sample_matrix <sampled_kmers> <samples_file> <pheno_file> <output_tsv>\n"
                  << "Env: SUBSAMPLE_SIZE (default 200000), MAF_THRESHOLD (default 5), SEED (default 42)\n";
        return 1;
    }
    const string sampled_file = argv[1];
    const string samples_file = argv[2];
    const string pheno_file   = argv[3];
    const string output_tsv   = argv[4];

    const size_t SUBSAMPLE_SIZE = env_size("SUBSAMPLE_SIZE", 200000);
    const size_t MAF_THRESHOLD  = env_size("MAF_THRESHOLD", 5);
    const size_t SEED           = env_size("SEED", 42);
    const double maf_min = MAF_THRESHOLD / 100.0;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Load ordered sample IDs
    vector<string> samples;
    {
        std::ifstream in(samples_file);
        if(!in) die("cannot open " + samples_file);
        string line;
        while(std::getline(in, line)){
            while(!line.empty() && isspace((unsigned char)line.back())) line.pop_back();
            if(!line.empty()) samples.push_back(line);
        }
    }
    std::cerr << "[INFO] Total samples in panel: " << samples.size() << "\n";

    // Load phenotype IDs
    std::unordered_set<string> pheno_ids;
    {
        std::ifstream in(pheno_file);
        if(!in) die("cannot open " + pheno_file);
        string id; double v;
        while(in >> id >> v) pheno_ids.insert(id);
    }
    std::cerr << "[INFO] Phenotyped accessions: " << pheno_ids.size() << "\n";

    // Compute phenotyped bit positions and ordered IDs
    vector<size_t> pheno_positions;
    vector<string> phenotyped_ids;
    for(size_t i = 0; i < samples.size(); i++){
        if(pheno_ids.count(samples[i])){
            pheno_positions.push_back(i);
            phenotyped_ids.push_back(samples[i]);
        }
    }
    const size_t n_pheno = pheno_positions.size();
    if(n_pheno == 0) die("no overlap between samples and phenotype file");
    std::cerr << "[INFO] Phenotyped positions in bitstring: " << n_pheno << "\n";

    // Stream through sampled k-mers, filter by MAF
    vector<Record> passing;
    {
        GzReader gz(sampled_file);
        string line;
        size_t total = 0, kept = 0;
        while(gz.getline(line)){
            if(line.empty()) continue;
            total++;

            size_t t = line.find('\t');
            if(t == string::npos) t = line.find(' ');
            if(t == string::npos) continue;

            size_t bits_start = t + 1;
            size_t bits_len = line.size() - bits_start;
            while(bits_len > 0 && isspace((unsigned char)line[bits_start + bits_len - 1])) --bits_len;

            // Extract phenotyped bits and compute MAF
            size_t ones = 0;
            vector<char> pbits(n_pheno);
            bool valid = true;
            for(size_t j = 0; j < n_pheno; j++){
                size_t pos = pheno_positions[j];
                if(pos >= bits_len){ valid = false; break; }
                char c = line[bits_start + pos];
                pbits[j] = c;
                if(c == '1') ones++;
            }
            if(!valid) continue;

            double freq = (double)ones / n_pheno;
            double maf = std::min(freq, 1.0 - freq);
            if(maf < maf_min) continue;

            Record rec;
            rec.kmer.assign(line.data(), t);
            rec.bits = std::move(pbits);
            passing.push_back(std::move(rec));
            kept++;
        }
        std::cerr << "[INFO] Scanned " << total << " k-mers, " << kept
                  << " passed MAF >= " << MAF_THRESHOLD << "% filter\n";
    }

    // Subsample if needed
    if(passing.size() > SUBSAMPLE_SIZE){
        std::mt19937_64 rng(SEED);
        std::shuffle(passing.begin(), passing.end(), rng);
        passing.resize(SUBSAMPLE_SIZE);
        std::cerr << "[INFO] Subsampled to " << SUBSAMPLE_SIZE << " records\n";
    } else {
        std::cerr << "[WARN] Only " << passing.size()
                  << " records passed filter (requested " << SUBSAMPLE_SIZE << ")\n";
    }

    // Write output: rows = samples, columns = kmers
    // precompute.cpp:load_X_matrix expects:
    //   header line (skipped)
    //   each data line: <sample_id>\t<val1>\t<val2>\t...
    {
        std::ofstream out(output_tsv);
        if(!out) die("cannot open " + output_tsv + " for writing");

        // Header
        out << "sample_id";
        for(size_t k = 0; k < passing.size(); k++)
            out << '\t' << passing[k].kmer;
        out << '\n';

        // One row per phenotyped sample
        for(size_t s = 0; s < n_pheno; s++){
            out << phenotyped_ids[s];
            for(size_t k = 0; k < passing.size(); k++)
                out << '\t' << (passing[k].bits[s] == '1' ? '1' : '0');
            out << '\n';
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cerr << "[INFO] Wrote " << n_pheno << " x " << passing.size()
              << " matrix to " << output_tsv << " in " << secs << "s\n";
    return 0;
}
