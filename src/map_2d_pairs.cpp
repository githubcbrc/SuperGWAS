// map_2d_pairs.cpp — Map k-mer pair indices to genome positions per reference
// Reads kgwas_2d output (idx_a, idx_b, effect, pval) + kmers file,
// scans reference FASTA to find genome positions, translates pairs.
//
// Usage: map_2d_pairs kmers.txt pairs.tsv REF_FASTA K out.tsv [--threads N]

#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "map_utils.hpp"

using std::string;
using std::vector;

struct Pair {
  int64_t idx_a, idx_b;
  double effect, pval;
};

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);

  if (argc < 6) {
    std::cerr << "Usage: " << argv[0]
              << " kmers.txt pairs.tsv REF_FASTA K out.tsv [--threads N]\n";
    return 1;
  }

  const string kmers_path = argv[1];
  const string pairs_path = argv[2];
  const string ref_path   = argv[3];
  const size_t K          = (size_t)std::stoull(argv[4]);
  const string out_path   = argv[5];

  int threads = 0;
  for (int i = 6; i < argc; ++i) {
    if (string(argv[i]) == "--threads" && i+1 < argc)
      threads = std::stoi(argv[++i]);
  }

#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif

  // --- 1. Load k-mer sequences ---
  vector<string> kmer_seqs;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::ifstream in(kmers_path);
    if (!in) { std::cerr << "[ERROR] cannot open " << kmers_path << "\n"; return 1; }
    string line;
    while (std::getline(in, line)) {
      if (!line.empty()) kmer_seqs.push_back(line);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "[INFO] Loaded " << kmer_seqs.size() << " k-mers in "
              << std::chrono::duration<double>(t1-t0).count() << " s\n";
  }

  // --- 2. Build encoded k-mer → index lookup ---
  std::unordered_map<EncodedKmer, int64_t, EncodedKmerHash> kmer_to_idx;
  kmer_to_idx.reserve(kmer_seqs.size() * 2);
  for (size_t i = 0; i < kmer_seqs.size(); ++i) {
    EncodedKmer enc = encode_kmer(kmer_seqs[i].c_str(), K);
    EncodedKmer can = canonical_encoded(enc, K);
    kmer_to_idx[can] = (int64_t)i;
  }

  // --- 3. Load reference FASTA ---
  vector<Contig> contigs;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    if (!read_reference_fasta(ref_path, contigs)) return 2;
    auto t1 = std::chrono::high_resolution_clock::now();
    size_t total_bp = 0;
    for (auto& c : contigs) total_bp += c.seq.size();
    std::cerr << "[INFO] Loaded reference: " << contigs.size() << " contigs, "
              << total_bp << " bp in " << std::chrono::duration<double>(t1-t0).count() << " s\n";
  }

  // Contig offsets for genome-wide coordinates
  vector<size_t> contig_offsets(contigs.size());
  {
    size_t offset = 0;
    for (size_t i = 0; i < contigs.size(); ++i) {
      contig_offsets[i] = offset;
      offset += contigs[i].seq.size();
    }
  }

  // --- 4. Scan reference, map k-mers to genome positions ---
  // For each k-mer index, store the first genome position found
  vector<int64_t> idx_to_pos(kmer_seqs.size(), -1);
  const EncodedKmer mask = KMER_MASK(K);

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    size_t mapped = 0;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:mapped)
    #endif
    for (size_t cidx = 0; cidx < contigs.size(); ++cidx) {
      const string& seq = contigs[cidx].seq;
      if (seq.size() < K) continue;
      const size_t gw_offset = contig_offsets[cidx];

      EncodedKmer enc = 0;
      size_t valid = 0;

      for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t code = BASE_TO_2BIT.v[(unsigned char)seq[i]];
        if (code == ENCODE_INVALID) { valid = 0; enc = 0; continue; }
        enc = ((enc << 2) | code) & mask;
        if (++valid < K) continue;

        EncodedKmer canonical = canonical_encoded(enc, K);
        auto it = kmer_to_idx.find(canonical);
        if (it == kmer_to_idx.end()) continue;

        int64_t idx = it->second;
        // First hit wins
        if (idx_to_pos[idx] == -1) {
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          {
            if (idx_to_pos[idx] == -1) {
              idx_to_pos[idx] = (int64_t)(gw_offset + i - K + 1);
              ++mapped;
            }
          }
        }
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "[INFO] Mapped " << mapped << "/" << kmer_seqs.size()
              << " k-mers to genome in "
              << std::chrono::duration<double>(t1-t0).count() << " s\n";
  }

  // --- 5. Load pairs and translate ---
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::ifstream in(pairs_path);
    if (!in) { std::cerr << "[ERROR] cannot open " << pairs_path << "\n"; return 1; }

    std::ofstream out(out_path);
    if (!out) { std::cerr << "[ERROR] cannot open " << out_path << "\n"; return 1; }

    size_t total = 0, written = 0;
    string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      ++total;

      // Parse: idx_a \t idx_b \t effect \t pval
      size_t t1p = line.find('\t');
      size_t t2p = line.find('\t', t1p+1);
      if (t1p == string::npos || t2p == string::npos) continue;

      int64_t ia = std::stoll(line.substr(0, t1p));
      int64_t ib = std::stoll(line.substr(t1p+1, t2p-(t1p+1)));

      if (ia < 0 || ia >= (int64_t)kmer_seqs.size() ||
          ib < 0 || ib >= (int64_t)kmer_seqs.size()) continue;

      int64_t pa = idx_to_pos[ia];
      int64_t pb = idx_to_pos[ib];
      if (pa < 0 || pb < 0) continue;  // k-mer(s) don't map to this ref

      // Ensure pos_a < pos_b
      if (pa > pb) std::swap(pa, pb);

      // Output: pos_a \t pos_b \t effect \t pval (rest of line from t2p)
      out << pa << '\t' << pb << line.substr(t2p) << '\n';
      ++written;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cerr << "[INFO] Translated " << written << "/" << total << " pairs in "
              << std::chrono::duration<double>(t1-t0).count() << " s\n";
  }

  return 0;
}
