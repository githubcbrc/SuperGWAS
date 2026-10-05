#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <zlib.h>

namespace kgwas {

using std::string;
using std::vector;

struct Timer {
  std::chrono::high_resolution_clock::time_point t0;
  Timer(){ reset(); }
  void reset(){ t0 = std::chrono::high_resolution_clock::now(); }
  double sec() const {
    using clk = std::chrono::high_resolution_clock;
    return std::chrono::duration_cast<std::chrono::duration<double>>(clk::now() - t0).count();
  }
};

// ---------- small utils ----------
[[noreturn]] inline void die(const string& m){
  std::cerr << "[ERROR] " << m << "\n";
  std::exit(1);
}

inline void warn(const string& m){
  std::cerr << "[WARN] " << m << "\n";
}

inline int env_int(const char* k, int dflt){
  const char* v = std::getenv(k);
  if(!v || !*v) return dflt;
  try { return std::max(1, std::stoi(v)); } catch(...) { return dflt; }
}


// ---------- I/O helpers ----------
inline vector<string> load_lines(const string& file){
  std::ifstream in(file);
  if(!in) die("open: " + file);
  vector<string> v; v.reserve(1024);
  string s;
  while(std::getline(in,s)) if(!s.empty()){
    if(!s.empty() && s.back()=='\r') s.pop_back();
    v.push_back(s);
  }
  return v;
}

inline std::unordered_map<string,double> load_pheno(const string& file){
  std::ifstream in(file);
  if(!in) die("open: " + file);
  std::unordered_map<string,double> m;
  string id; double val;
  while(in>>id>>val) m[id]=val;
  return m;
}

inline size_t probe_N_txt(const string& path){
  int fd = ::open(path.c_str(), O_RDONLY);
  if(fd < 0) die("open: " + path);
  struct stat st{}; if(fstat(fd,&st)!=0){ ::close(fd); die("fstat: " + path); }
  size_t sz = (size_t)st.st_size;
  if(sz == 0){ ::close(fd); die("empty: " + path); }
  void* p = mmap(nullptr, sz, PROT_READ, MAP_PRIVATE, fd, 0);
  if(p==MAP_FAILED){ ::close(fd); die("mmap: " + path); }
  madvise(p, sz, MADV_SEQUENTIAL);
  const char* b = (const char*)p;
  const void* nlv = memchr(b, '\n', sz);
  size_t N = nlv ? (size_t)((const char*)nlv - b) : sz; // tolerate missing NL
  munmap(p, sz); ::close(fd);
  if(N==0) die("no bytes before newline: " + path);
  return N;
}

inline bool ends_with(const string& s, const string& suf){
  return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

inline vector<string> load_kmers(const string& path, size_t N){
  vector<string> kmers; kmers.reserve(N);
  if(ends_with(path, ".gz")){
    gzFile f = gzopen(path.c_str(), "rb");
    if(!f) die("gzopen: " + path);
    const int BUFSZ = 1<<15;
    std::unique_ptr<char[]> buf(new char[BUFSZ]);
    string line; line.reserve(256);
    while(kmers.size() < N){
      char* r = gzgets(f, buf.get(), BUFSZ);
      if(!r){
        int err=0; const char* msg = gzerror(f,&err);
        if(err == Z_STREAM_END) break;  // normal EOF
        if(err != Z_OK) { gzclose(f); die(string("gzgets error: ") + (msg?msg:"")); }
        break;
      }
      char* nl = strchr(r, '\n');
      if(nl){ *nl = '\0'; line.append(r); kmers.push_back(std::move(line)); line.clear(); }
      else { line.append(r); continue; }
    }
    gzclose(f);
  } else {
    std::ifstream in(path);
    if(!in) die("open: " + path);
    string s;
    while(kmers.size() < N && std::getline(in, s)){
      if(!s.empty() && s.back()=='\r') s.pop_back();
      kmers.push_back(std::move(s));
    }
  }
  if(kmers.size() != N){
    die("kmers lines != N (got " + std::to_string(kmers.size()) + ", expected " + std::to_string(N) + ")");
  }
  return kmers;
}

struct UData {
  Eigen::MatrixXd U;               // (n x r), aligned to Overlap order
  int r = 0;                       // rank
  double yyM = 0.0;                // y'y - (U'y)^2
};

static UData load_Ubin(const string& U_file, int n, const Eigen::VectorXd& y) {
  Timer t;
  std::ifstream uin(U_file, std::ios::binary);
  if(!uin) die("open: " + U_file);
  int n_u=0, r=0;
  uin.read((char*)&n_u, sizeof(int));
  uin.read((char*)&r, sizeof(int));
  if(!uin) die("short U.bin header");
  if(n_u!=n) die("U.bin n mismatch: file n="+std::to_string(n_u)+" vs overlap n="+std::to_string(n));

  UData ud;
  ud.U.resize(n, r);
  uin.read((char*)ud.U.data(), sizeof(double)*n*r);
  if(!uin) die("short U.bin body");
  ud.r = r;

  const Eigen::VectorXd Uy = ud.U.transpose()*y;
  ud.yyM = y.squaredNorm() - Uy.squaredNorm();
  if(!(ud.yyM>1e-12)) warn("yyM ~ 0");

  std::cerr << "[INFO] loaded U (n="<<n<<", r="<<r<<")  ("<<t.sec()<<" s)\n";
  return ud;
}


using MatRow = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

static inline bool file_exists(const std::string& p){
  struct stat st{}; return ::stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

static size_t probe_N_gz(const std::string& path){
  gzFile f = gzopen(path.c_str(), "rb");
  if(!f) die("gzopen: " + path);
  gzbuffer(f, 1<<20); // 1 MiB zlib buffer helps throughput

  const size_t CHUNK = 1<<20;
  std::unique_ptr<unsigned char[]> buf(new unsigned char[CHUNK]);
  size_t total = 0;

  for(;;){
    int n = gzread(f, buf.get(), (unsigned)CHUNK);
    if(n < 0){
      int err=0; const char* msg = gzerror(f,&err);
      gzclose(f);
      die(std::string("gzread error: ") + (msg?msg:""));
    }
    if(n == 0) break; // EOF

    void* nl = memchr(buf.get(), '\n', (size_t)n);
    if(nl){
      const unsigned char* p = (const unsigned char*)buf.get();
      size_t off = (size_t)((const unsigned char*)nl - p);
      total += off;
      gzclose(f);
      return total;
    }
    total += (size_t)n;
  }
  gzclose(f);
  if(total == 0) die("empty gz column: " + path);
  return total; 
}


static void read_gz_exact_N(const std::string& path, unsigned char* dst, size_t N){
  gzFile f = gzopen(path.c_str(), "rb");
  if(!f) die("gzopen: " + path);
  gzbuffer(f, 1<<20); // bigger zlib buffer helps throughput

  size_t got = 0;
  while(got < N){
    int n = gzread(f, dst + got, (unsigned int)std::min<size_t>(N - got, 1<<20));
    if(n <= 0){
      int err=0; const char* msg = gzerror(f, &err);
      gzclose(f);
      if(err == Z_STREAM_END) die("short gz column: " + path + " (got " + std::to_string(got) + " < N " + std::to_string(N) + ")");
      die(std::string("gzread error: ") + (msg?msg:""));
    }
    got += (size_t)n;
  }
  unsigned char extra;
  int n = gzread(f, &extra, 1);
  if(n == 1 && extra != '\n'){
    gzclose(f);
    die("gz column longer than expected (N+ extra non-newline): " + path);
  }
  gzclose(f);
}

static MatRow load_matrix_rows(const std::string& bits_dir,
                               const std::vector<int>& idxs,
                               int n, size_t N,
                               int LOAD_THREADS)
{
  Timer t;
  MatRow Xrm(n, (int)N);

  std::atomic<int> done{0};
  omp_set_num_threads(LOAD_THREADS);

  #pragma omp parallel for schedule(dynamic,1)
  for(int k=0;k<n;++k){
    char base[64]; std::snprintf(base,sizeof(base),"col_%06d.txt", idxs[k]);
    const std::string p_txt = bits_dir + "/" + base;
    const std::string p_gz  = p_txt + ".gz";

    float* row = Xrm.row(k).data();

    if(file_exists(p_txt)){
      int fd = ::open(p_txt.c_str(), O_RDONLY);
      if(fd<0) die("open: " + p_txt);
      struct stat st{}; if(fstat(fd,&st)!=0){ ::close(fd); die("fstat: " + p_txt); }
      size_t sz = (size_t)st.st_size;
      if(sz < N){ ::close(fd); die("short column: " + p_txt + " (size " + std::to_string(sz) + " < N " + std::to_string(N) + ")"); }

      void* p = mmap(nullptr, sz, PROT_READ, MAP_PRIVATE, fd, 0);
      if(p==MAP_FAILED){ ::close(fd); die("mmap: " + p_txt); }
      madvise(p, sz, MADV_SEQUENTIAL);

      const unsigned char* bytes = (const unsigned char*)p;
      #pragma omp simd
      for(size_t i=0;i<N;++i) row[i] = float(bytes[i] & 1u);

      munmap(p, sz); ::close(fd);
    }
    else if(file_exists(p_gz)){
      std::unique_ptr<unsigned char[]> buf(new unsigned char[N]);
      read_gz_exact_N(p_gz, buf.get(), N);

      #pragma omp simd
      for(size_t i=0;i<N;++i) row[i] = float(buf[i] & 1u);
    }
    else{
      die("missing column: " + p_txt + " (and " + p_gz + ")");
    }

    int d = ++done;
    if((d % 8) == 0){
      #pragma omp critical
      std::cerr << "\r[LOAD] " << d << "/" << n << std::flush;
    }
  }

  std::cerr << "\r[LOAD] " << n << "/" << n
            << "  ("<<t.sec()<<" s, LOAD_THREADS="<<LOAD_THREADS<<")\n";
  return Xrm;
}

// ---------- Overlap: unified sample/phenotype loader ----------
struct Overlap {
  vector<int> idxs;                // positions in samples that exist in phenotype
  vector<string> samples;          // all samples (original order)
  Eigen::VectorXd y;               // overlapped phenotype vector (n)
};

static Overlap build_overlap(const string& samples_file, const string& pheno_file) {
  Timer t;
  Overlap ov;
  ov.samples = load_lines(samples_file);
  const int L = (int)ov.samples.size();
  if(L<=0) die("samples empty");

  auto ph = load_pheno(pheno_file);
  ov.idxs.reserve(L);
  std::vector<double> yv; yv.reserve(L);

  for(int j=0;j<L;++j){
    auto it = ph.find(ov.samples[j]);
    if(it!=ph.end()){
      ov.idxs.push_back(j);
      yv.push_back(it->second);
    }
  }
  if(ov.idxs.empty()) die("no phenotype overlap");
  ov.y = Eigen::Map<Eigen::VectorXd>(yv.data(), (int)yv.size());

  std::cerr << "[INFO] n="<<ov.y.size()<<" (overlap)  ("<<t.sec()<<" s)\n";
  return ov;
}

// ---------- infer N from bits directory ----------
static size_t infer_N_from_bits(const std::string& bits_dir, int first_idx) {
  Timer t;
  char base[64]; std::snprintf(base, sizeof(base), "col_%06d.txt", first_idx);
  const std::string p_txt = bits_dir + "/" + base;
  const std::string p_gz  = p_txt + ".gz";

  size_t N = 0;
  if(file_exists(p_txt))      N = probe_N_txt(p_txt);
  else if(file_exists(p_gz))  N = probe_N_gz(p_gz);
  else die("missing column: " + p_txt + " (and " + p_gz + ")");

  if(N == 0) die("N=0");
  std::cerr << "[INFO] inferred N=" << N << " rows  (" << t.sec() << " s)"
            << (file_exists(p_txt) ? " [from .txt]" : " [from .txt.gz]") << "\n";
  return N;
}

} // namespace kgwas

