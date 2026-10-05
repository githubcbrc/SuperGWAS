#pragma once
#include <bits/stdc++.h>
#include <filesystem>


#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

// ===============================
// Types & constants
// ===============================

struct Entry { double effect; double pval; std::string bitstring; };
struct Contig {std::string name; std::string seq; };

using KmerMap = std::unordered_map<std::string, Entry>;

// ===============================
// 2-bit k-mer encoding
// ===============================

using EncodedKmer = __uint128_t;

static constexpr uint8_t ENCODE_INVALID = 0xFF;

struct Base2BitTable {
    uint8_t v[256];
    constexpr Base2BitTable() : v{} {
        for (int i = 0; i < 256; ++i) v[i] = ENCODE_INVALID;
        v['A'] = 0; v['a'] = 0;
        v['C'] = 1; v['c'] = 1;
        v['G'] = 2; v['g'] = 2;
        v['T'] = 3; v['t'] = 3;
    }
};
static inline constexpr Base2BitTable BASE_TO_2BIT{};

static inline EncodedKmer KMER_MASK(size_t K) {
    if (2 * K >= 128) return ~EncodedKmer(0);
    return (EncodedKmer(1) << (2 * K)) - 1;
}

static inline EncodedKmer encode_kmer(const char* seq, size_t K) {
    EncodedKmer enc = 0;
    for (size_t i = 0; i < K; ++i)
        enc = (enc << 2) | BASE_TO_2BIT.v[(unsigned char)seq[i]];
    return enc;
}

static inline void decode_kmer(EncodedKmer enc, size_t K, char* out) {
    static constexpr char DECODE[4] = {'A', 'C', 'G', 'T'};
    for (size_t i = K; i > 0; --i) {
        out[i - 1] = DECODE[enc & 3];
        enc >>= 2;
    }
}

static inline EncodedKmer revcomp_encoded(EncodedKmer enc, size_t K) {
    uint64_t lo = static_cast<uint64_t>(enc);
    uint64_t hi = static_cast<uint64_t>(enc >> 64);

    // Byte-swap each half
    lo = __builtin_bswap64(lo);
    hi = __builtin_bswap64(hi);

    // Reverse 2-bit pairs within each byte
    lo = ((lo >> 2) & 0x3333333333333333ULL) | ((lo & 0x3333333333333333ULL) << 2);
    lo = ((lo >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((lo & 0x0F0F0F0F0F0F0F0FULL) << 4);
    hi = ((hi >> 2) & 0x3333333333333333ULL) | ((hi & 0x3333333333333333ULL) << 2);
    hi = ((hi >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((hi & 0x0F0F0F0F0F0F0F0FULL) << 4);

    // Swap halves and reconstruct
    EncodedKmer reversed = (static_cast<EncodedKmer>(lo) << 64) | hi;

    // Align: shift right to put K bases in lowest 2K bits
    reversed >>= (128 - 2 * K);

    // Complement: XOR all 2K bits (A↔T, C↔G)
    reversed ^= KMER_MASK(K);

    return reversed;
}

static inline EncodedKmer canonical_encoded(EncodedKmer enc, size_t K) {
    EncodedKmer rc = revcomp_encoded(enc, K);
    return (enc < rc) ? enc : rc;
}

struct EncodedKmerHash {
    std::size_t operator()(EncodedKmer key) const noexcept {
        uint64_t x = static_cast<uint64_t>(key) ^ static_cast<uint64_t>(key >> 64);
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return static_cast<std::size_t>(x);
    }
};

// ===============================
// Sharded index types (d4 format)
// ===============================

struct EntryD4 { double effect; double pval; std::string bitstring; };
using ShardMap = std::unordered_map<EncodedKmer, EntryD4, EncodedKmerHash>;

struct EntryCompact { double effect; double pval; };
using ShardMapCompact = std::unordered_map<EncodedKmer, EntryCompact, EncodedKmerHash>;

static inline constexpr char UMAP_MAGIC_D4[8] = {'U','M','A','P','d','4','\0','\0'};

// ===============================
// Bloom filter (2-hash, power-of-2 size)
// ===============================

struct BloomFilter {
    std::vector<uint64_t> data;
    uint64_t mask;  // m - 1, where m = number of bits

    BloomFilter() : mask(0) {}

    void init(uint64_t num_elements, uint64_t bits_per_element = 10) {
        uint64_t m = 1;
        uint64_t target = num_elements * bits_per_element;
        if (target == 0) target = 64;
        while (m < target) m <<= 1;
        mask = m - 1;
        data.assign((m + 63) / 64, 0);
    }

    static std::pair<uint64_t, uint64_t> hash2(EncodedKmer key) {
        uint64_t x = static_cast<uint64_t>(key) ^ static_cast<uint64_t>(key >> 64);
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        uint64_t h1 = x;
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27; x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return {h1, x};
    }

    void insert(EncodedKmer key) {
        auto [h1, h2] = hash2(key);
        uint64_t i1 = h1 & mask, i2 = h2 & mask;
        data[i1 >> 6] |= (1ULL << (i1 & 63));
        data[i2 >> 6] |= (1ULL << (i2 & 63));
    }

    bool test(EncodedKmer key) const {
        auto [h1, h2] = hash2(key);
        uint64_t i1 = h1 & mask, i2 = h2 & mask;
        if (!(data[i1 >> 6] & (1ULL << (i1 & 63)))) return false;
        return (data[i2 >> 6] & (1ULL << (i2 & 63))) != 0;
    }

    std::size_t size_bytes() const { return data.size() * 8; }
};

// ===============================
// Helpers
// ===============================

static inline bool looks_like_header(const std::string& line) {
    std::string t = line;
    for (char& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return t.find("kmer") != std::string::npos && t.find("effect") != std::string::npos;
}

static std::size_t infer_K(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[ERROR] Could not open: " << path << "\n";
        return 0;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (looks_like_header(line)) continue;

        std::size_t i1 = line.find_first_of("\t ");
        if (i1 == std::string::npos) continue; // need at least kmer + one value

        std::string kmer = line.substr(0, i1);
        if (!kmer.empty()) return kmer.size();
    }
    return 0;
}

// ===============================
// Loader (.umap, double)
// Format (little-endian):
//   magic[8] = "UMAPd1\0"
//   uint32_t K
//   uint64_t N
//   N × { char[K], double effect, double pval }
// ===============================

static inline constexpr char UMAP_MAGIC_D1[8] = {'U','M','A','P','d','1','\0','\0'};
static inline constexpr char UMAP_MAGIC_D2[8] = {'U','M','A','P','d','2','\0','\0'};
static inline constexpr char UMAP_MAGIC_D3[8] = {'U','M','A','P','d','3','\0','\0'};
static inline constexpr char BITS_MAGIC[8]    = {'U','M','P','b','\0','\0','\0','\0'};
namespace detail {
    inline bool read_exact(std::istream& in, char* buf, std::streamsize n) {
        in.read(buf, n);
        return in.good();
    }
    inline bool write_exact(std::ostream& out, const char* buf, std::streamsize n) {
        out.write(buf, n);
        return static_cast<bool>(out);
    }
    inline bool is_little_endian() {
        const uint16_t x = 0x0102;
        return *(reinterpret_cast<const unsigned char*>(&x)) == 0x02;
    }
}

enum class UmapVersion { d1, d2, d3 };

inline bool load_umap(const std::string& path,
                             KmerMap& db,
                             std::size_t expected_K = 0,
                             std::size_t reserve_hint = 0,
                             UmapVersion* out_version = nullptr)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::cerr << "[ERROR] cannot open " << path << "\n";
        return false;
    }
    if (!detail::is_little_endian()) {
        std::cerr << "[ERROR] big-endian host unsupported for UMAPd1 (add byte-swap).\n";
        return false;
    }

    char magic[8];
    if (!detail::read_exact(in, magic, sizeof(magic))) {
        std::cerr << "[ERROR] short read (magic): " << path << "\n";
        return false;
    }
    bool is_d1 = (std::memcmp(magic, UMAP_MAGIC_D1, sizeof(magic)) == 0);
    bool is_d2 = (std::memcmp(magic, UMAP_MAGIC_D2, sizeof(magic)) == 0);
    bool is_d3 = (std::memcmp(magic, UMAP_MAGIC_D3, sizeof(magic)) == 0);
    if (!is_d1 && !is_d2 && !is_d3) {
        std::cerr << "[ERROR] incompatible .umap (need UMAPd1, UMAPd2 or UMAPd3): " << path << "\n";
        return false;
    }
    if (out_version) {
        *out_version = is_d3 ? UmapVersion::d3 : (is_d2 ? UmapVersion::d2 : UmapVersion::d1);
    }

    std::uint32_t K = 0;
    if (!detail::read_exact(in, reinterpret_cast<char*>(&K), sizeof(K)) || K == 0) {
        std::cerr << "[ERROR] invalid or short read (K): " << path << "\n";
        return false;
    }
    if (expected_K && expected_K != static_cast<std::size_t>(K)) {
        std::cerr << "[ERROR] K mismatch: file K=" << K << " expected " << expected_K << "\n";
        return false;
    }

    std::uint64_t N = 0;
    if (!detail::read_exact(in, reinterpret_cast<char*>(&N), sizeof(N))) {
        std::cerr << "[ERROR] short read (N): " << path << "\n";
        return false;
    }

    db.clear();
    if (reserve_hint == 0) reserve_hint = static_cast<std::size_t>(N * 1.3);
    db.reserve(reserve_hint);
    db.max_load_factor(0.70f);

    std::string k; k.resize(K);
    for (std::uint64_t i = 0; i < N; ++i) {
        if (!detail::read_exact(in, k.data(), static_cast<std::streamsize>(K))) {
            std::cerr << "[ERROR] short read (kmer@" << i << "): " << path << "\n";
            return false;
        }
        Entry e{};
        if (!detail::read_exact(in, reinterpret_cast<char*>(&e.effect), sizeof(double)) ||
            !detail::read_exact(in, reinterpret_cast<char*>(&e.pval),   sizeof(double))) {
            std::cerr << "[ERROR] short read (payload@" << i << "): " << path << "\n";
            return false;
        }
        if (is_d2) {
            std::uint32_t blen = 0;
            if (!detail::read_exact(in, reinterpret_cast<char*>(&blen), sizeof(blen))) {
                std::cerr << "[ERROR] short read (bitstring len@" << i << "): " << path << "\n";
                return false;
            }
            e.bitstring.resize(blen);
            if (blen > 0 && !detail::read_exact(in, e.bitstring.data(), blen)) {
                std::cerr << "[ERROR] short read (bitstring@" << i << "): " << path << "\n";
                return false;
            }
        }
        db.emplace(k, std::move(e));
    }
    return true;
}

// ===============================
// Sharded d4 index I/O
// ===============================

inline bool write_index_meta(const std::string& dir,
                              std::uint32_t K, std::uint32_t M, std::uint64_t N) {
    std::ofstream out(dir + "/index.meta");
    if (!out) { std::cerr << "[ERROR] cannot write index.meta in " << dir << "\n"; return false; }
    out << "K=" << K << "\nM=" << M << "\nN=" << N << "\n";
    return static_cast<bool>(out);
}

inline bool read_index_meta(const std::string& dir,
                             std::uint32_t& K, std::uint32_t& M, std::uint64_t& N) {
    std::ifstream in(dir + "/index.meta");
    if (!in) { std::cerr << "[ERROR] cannot read index.meta in " << dir << "\n"; return false; }
    std::string line;
    K = 0; M = 0; N = 0;
    while (std::getline(in, line)) {
        if (line.compare(0, 2, "K=") == 0) K = static_cast<std::uint32_t>(std::stoul(line.substr(2)));
        else if (line.compare(0, 2, "M=") == 0) M = static_cast<std::uint32_t>(std::stoul(line.substr(2)));
        else if (line.compare(0, 2, "N=") == 0) N = std::stoull(line.substr(2));
    }
    if (K == 0 || M == 0) {
        std::cerr << "[ERROR] invalid index.meta (K=" << K << " M=" << M << ")\n";
        return false;
    }
    return true;
}

inline bool dump_shard(const std::string& dir, std::uint32_t shard_id,
                        std::uint32_t K,
                        const std::vector<std::pair<EncodedKmer, EntryD4>>& entries) {
    char fname[32];
    std::snprintf(fname, sizeof(fname), "shard_%03u", shard_id);
    std::string umap_path = dir + "/" + fname + ".umap";
    std::string bits_path = dir + "/" + fname + ".bits";

    // .umap: header + fixed 32-byte records
    std::ofstream uout(umap_path, std::ios::binary);
    if (!uout) { std::cerr << "[ERROR] cannot write " << umap_path << "\n"; return false; }

    detail::write_exact(uout, UMAP_MAGIC_D4, 8);
    detail::write_exact(uout, reinterpret_cast<const char*>(&K), sizeof(K));
    std::uint64_t N = entries.size();
    detail::write_exact(uout, reinterpret_cast<const char*>(&N), sizeof(N));

    for (const auto& [enc, entry] : entries) {
        detail::write_exact(uout, reinterpret_cast<const char*>(&enc), sizeof(EncodedKmer));
        detail::write_exact(uout, reinterpret_cast<const char*>(&entry.effect), sizeof(double));
        detail::write_exact(uout, reinterpret_cast<const char*>(&entry.pval), sizeof(double));
    }

    // .bits: only written if any entry has a bitstring
    bool has_bits = false;
    for (const auto& [enc, entry] : entries) {
        if (!entry.bitstring.empty()) { has_bits = true; break; }
    }

    if (has_bits) {
        std::ofstream bout(bits_path, std::ios::binary);
        if (!bout) { std::cerr << "[ERROR] cannot write " << bits_path << "\n"; return false; }

        detail::write_exact(bout, reinterpret_cast<const char*>(&N), sizeof(N));
        for (const auto& [enc, entry] : entries) {
            std::uint32_t blen = static_cast<std::uint32_t>(entry.bitstring.size());
            detail::write_exact(bout, reinterpret_cast<const char*>(&blen), sizeof(blen));
            if (blen > 0) detail::write_exact(bout, entry.bitstring.data(), blen);
        }
        if (!bout) return false;
    }

    return static_cast<bool>(uout);
}

inline bool load_shard(const std::string& dir, std::uint32_t shard_id,
                        std::uint32_t expected_K, ShardMap& shard,
                        bool load_bits) {
    char fname[32];
    std::snprintf(fname, sizeof(fname), "shard_%03u", shard_id);
    std::string umap_path = dir + "/" + fname + ".umap";

    // Bulk-read entire .umap file
    std::ifstream in(umap_path, std::ios::binary);
    if (!in) { std::cerr << "[ERROR] cannot open " << umap_path << "\n"; return false; }

    in.seekg(0, std::ios::end);
    std::size_t file_size = static_cast<std::size_t>(in.tellg());
    in.seekg(0, std::ios::beg);

    std::vector<char> buf(file_size);
    in.read(buf.data(), static_cast<std::streamsize>(file_size));
    if (!in) { std::cerr << "[ERROR] short read: " << umap_path << "\n"; return false; }
    in.close();

    // Parse header (20 bytes: magic[8] + K[4] + N[8])
    if (file_size < 20) { std::cerr << "[ERROR] file too small: " << umap_path << "\n"; return false; }
    const char* ptr = buf.data();

    if (std::memcmp(ptr, UMAP_MAGIC_D4, 8) != 0) {
        std::cerr << "[ERROR] bad magic in " << umap_path << "\n";
        return false;
    }
    ptr += 8;

    std::uint32_t K;
    std::memcpy(&K, ptr, sizeof(K)); ptr += sizeof(K);
    if (expected_K && K != expected_K) {
        std::cerr << "[ERROR] K mismatch in " << umap_path << ": got " << K
                  << " expected " << expected_K << "\n";
        return false;
    }

    std::uint64_t N;
    std::memcpy(&N, ptr, sizeof(N)); ptr += sizeof(N);

    shard.reserve(static_cast<std::size_t>(N * 1.3));
    shard.max_load_factor(0.70f);

    // Parse records: 32-byte stride (16B key + 8B effect + 8B pval)
    for (std::uint64_t i = 0; i < N; ++i) {
        EncodedKmer enc;
        std::memcpy(&enc, ptr, sizeof(EncodedKmer)); ptr += sizeof(EncodedKmer);
        double effect, pval;
        std::memcpy(&effect, ptr, sizeof(double)); ptr += sizeof(double);
        std::memcpy(&pval, ptr, sizeof(double)); ptr += sizeof(double);
        shard.emplace(enc, EntryD4{effect, pval, {}});
    }

    // Load bitstrings if requested
    if (load_bits) {
        std::string bits_path = dir + "/" + fname + ".bits";
        std::ifstream bin(bits_path, std::ios::binary);
        if (!bin) { std::cerr << "[ERROR] cannot open " << bits_path << "\n"; return false; }

        bin.seekg(0, std::ios::end);
        std::size_t bfile_size = static_cast<std::size_t>(bin.tellg());
        bin.seekg(0, std::ios::beg);

        std::vector<char> bbuf(bfile_size);
        bin.read(bbuf.data(), static_cast<std::streamsize>(bfile_size));
        if (!bin) { std::cerr << "[ERROR] short read: " << bits_path << "\n"; return false; }
        bin.close();

        const char* bptr = bbuf.data();
        std::uint64_t bN;
        std::memcpy(&bN, bptr, sizeof(bN)); bptr += sizeof(bN);
        if (bN != N) {
            std::cerr << "[ERROR] .bits N=" << bN << " != .umap N=" << N << "\n";
            return false;
        }

        // Re-read keys from .umap buffer to match bitstrings by file order
        const char* kptr = buf.data() + 20; // skip header
        for (std::uint64_t i = 0; i < N; ++i) {
            EncodedKmer enc;
            std::memcpy(&enc, kptr, sizeof(EncodedKmer));
            kptr += 32; // stride to next record

            std::uint32_t blen;
            std::memcpy(&blen, bptr, sizeof(blen)); bptr += sizeof(blen);

            if (blen > 0) {
                auto it = shard.find(enc);
                if (it != shard.end()) {
                    it->second.bitstring.assign(bptr, blen);
                }
            }
            bptr += blen;
        }
    }

    return true;
}

inline bool load_sharded_index(const std::string& dir,
                                std::vector<ShardMap>& shards,
                                bool load_bits) {
    std::uint32_t K, M;
    std::uint64_t N;
    if (!read_index_meta(dir, K, M, N)) return false;

    shards.resize(M);

    bool ok = true;
    #pragma omp parallel for schedule(dynamic, 1)
    for (std::uint32_t i = 0; i < M; ++i) {
        if (!load_shard(dir, i, K, shards[i], load_bits)) {
            #pragma omp critical
            { ok = false; std::cerr << "[ERROR] Failed to load shard " << i << "\n"; }
        }
    }

    if (ok) {
        std::uint64_t total = 0;
        for (const auto& s : shards) total += s.size();
        std::cerr << "[INFO] Loaded sharded index: " << total << " kmers ("
                  << M << " shards) from " << dir << "\n";
    }
    return ok;
}

// ===============================
// Compact shard loading (no bitstrings)
// ===============================

inline bool load_shard_compact(const std::string& dir, std::uint32_t shard_id,
                                std::uint32_t expected_K, ShardMapCompact& shard) {
    char fname[32];
    std::snprintf(fname, sizeof(fname), "shard_%03u", shard_id);
    std::string umap_path = dir + "/" + fname + ".umap";

    std::ifstream in(umap_path, std::ios::binary);
    if (!in) { std::cerr << "[ERROR] cannot open " << umap_path << "\n"; return false; }

    in.seekg(0, std::ios::end);
    std::size_t file_size = static_cast<std::size_t>(in.tellg());
    in.seekg(0, std::ios::beg);

    std::vector<char> buf(file_size);
    in.read(buf.data(), static_cast<std::streamsize>(file_size));
    if (!in) { std::cerr << "[ERROR] short read: " << umap_path << "\n"; return false; }
    in.close();

    if (file_size < 20) return false;
    const char* ptr = buf.data();

    if (std::memcmp(ptr, UMAP_MAGIC_D4, 8) != 0) return false;
    ptr += 8;

    std::uint32_t K;
    std::memcpy(&K, ptr, sizeof(K)); ptr += sizeof(K);
    if (expected_K && K != expected_K) return false;

    std::uint64_t N;
    std::memcpy(&N, ptr, sizeof(N)); ptr += sizeof(N);

    shard.reserve(static_cast<std::size_t>(N * 1.3));
    shard.max_load_factor(0.70f);

    for (std::uint64_t i = 0; i < N; ++i) {
        EncodedKmer enc;
        std::memcpy(&enc, ptr, sizeof(EncodedKmer)); ptr += sizeof(EncodedKmer);
        double effect, pval;
        std::memcpy(&effect, ptr, sizeof(double)); ptr += sizeof(double);
        std::memcpy(&pval, ptr, sizeof(double)); ptr += sizeof(double);
        shard.emplace(enc, EntryCompact{effect, pval});
    }
    return true;
}

// ===============================
// Bloom filter I/O
// ===============================

inline bool dump_bloom(const std::string& dir, std::uint32_t shard_id,
                        const BloomFilter& bf) {
    char fname[32];
    std::snprintf(fname, sizeof(fname), "shard_%03u.bloom", shard_id);
    std::ofstream out(dir + "/" + fname, std::ios::binary);
    if (!out) return false;
    std::uint64_t m = bf.mask + 1;
    std::uint64_t nwords = bf.data.size();
    out.write(reinterpret_cast<const char*>(&m), sizeof(m));
    out.write(reinterpret_cast<const char*>(&nwords), sizeof(nwords));
    out.write(reinterpret_cast<const char*>(bf.data.data()),
              static_cast<std::streamsize>(nwords * 8));
    return static_cast<bool>(out);
}

inline bool load_bloom(const std::string& dir, std::uint32_t shard_id,
                        BloomFilter& bf) {
    char fname[32];
    std::snprintf(fname, sizeof(fname), "shard_%03u.bloom", shard_id);
    std::ifstream in(dir + "/" + fname, std::ios::binary);
    if (!in) return false;
    std::uint64_t m, nwords;
    in.read(reinterpret_cast<char*>(&m), sizeof(m));
    in.read(reinterpret_cast<char*>(&nwords), sizeof(nwords));
    bf.mask = m - 1;
    bf.data.resize(nwords);
    in.read(reinterpret_cast<char*>(bf.data.data()),
            static_cast<std::streamsize>(nwords * 8));
    return static_cast<bool>(in);
}

// ===============================
// FASTA / sequence helpers
// ===============================

// Fast A/C/G/T -> complement (uppercase). Others -> 'N'
inline char comp_base_fast(char c) {
    // Uppercase ASCII expected. Map directly.
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

// In-place reverse complement;
inline void revcomp_inplace(std::string &s) {
    if (s.empty()) return;
    size_t i = 0, j = s.size() - 1;
    while (i < j) {
        char a = comp_base_fast(s[i]);
        char b = comp_base_fast(s[j]);
        s[i] = b; s[j] = a;
        ++i; --j;
    }
    if (i == j) s[i] = comp_base_fast(s[i]); // middle char for odd length
}

// Return true if any character in s[i..i+k) is 'N'
inline bool has_N(const std::string &s, size_t i, size_t k) {
    const char* p = s.data() + i;
    for (size_t j = 0; j < k; ++j) if (p[j] == 'N') return true;
    return false;
}

static bool read_reference_fasta(const std::string& path, std::vector<Contig>& out){
    std::ifstream in(path);
    if (!in) { std::cerr << "[ERROR] cannot open reference: " << path << "\n"; return false; }

    auto to_base = [](unsigned char c)->char {
        c = (unsigned char)std::toupper(c);
        switch (c) { case 'A': case 'C': case 'G': case 'T': return (char)c; default: return 'N'; }
    };

    out.clear();
    std::string line, cur_name, cur_seq;
    cur_seq.reserve(1<<20);

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!cur_name.empty()) {
                out.push_back({cur_name, cur_seq});
                cur_seq.clear();
            }
            std::string tok = line.substr(1);
            size_t sp = tok.find_first_of(" \t");
            if (sp != std::string::npos) tok = tok.substr(0, sp);
            cur_name = tok;
        } else {
            for (unsigned char c : line) {
                if (!std::isspace(c)) cur_seq.push_back(to_base(c));
            }
        }
    }
    if (!cur_name.empty()) out.push_back({cur_name, cur_seq});
    return true;
}








