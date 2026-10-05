// compute_structure_all.cpp
// PCA -> Ph, Khatri-Rao -> Phh, SVD of [1, Ph, Phh] -> U2D.bin
// Emits both U.bin (1D structure) and U2D.bin (2D structure)
//
// Usage: compute_structure_all X_sample.tsv phenotype.txt r outfolder zscore(0|1)

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <unordered_map>
#include <cmath>
#include <limits>

using Clock = std::chrono::high_resolution_clock;
inline double elapsed(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
}

struct MatrixWithIDs { Eigen::MatrixXd data; std::vector<std::string> ids; };

static inline std::string trim(const std::string& s) {
    size_t i = s.find_first_not_of(" \t\r\n"); if (i == std::string::npos) return "";
    size_t j = s.find_last_not_of(" \t\r\n");  return s.substr(i, j - i + 1);
}

static MatrixWithIDs load_X_matrix(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open " + filename);

    std::string line;
    std::vector<std::vector<double>> rows;
    std::vector<std::string> ids;

    if (!std::getline(file, line)) throw std::runtime_error("Empty X file");

    size_t expected_cols = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok; std::vector<double> row;

        std::getline(ss, tok, '\t');
        ids.push_back(trim(tok));

        while (std::getline(ss, tok, '\t')) {
            tok = trim(tok);
            row.push_back(tok.empty() ? 0.0 : std::stod(tok));
        }
        if (rows.empty()) expected_cols = row.size();
        if (row.size() != expected_cols)
            throw std::runtime_error("Inconsistent column count in X at id " + ids.back());
        rows.push_back(std::move(row));
    }
    if (rows.empty()) throw std::runtime_error("No data rows in X");

    Eigen::MatrixXd mat(rows.size(), rows[0].size());
    for (size_t i=0;i<rows.size();++i)
        for (size_t j=0;j<rows[0].size();++j)
            mat((Eigen::Index)i,(Eigen::Index)j) = rows[i][j];

    return {std::move(mat), std::move(ids)};
}

static std::vector<std::string> load_pheno_ids_ordered(const std::string& fn) {
    std::ifstream in(fn);
    if (!in.is_open()) throw std::runtime_error("Cannot open " + fn);
    std::vector<std::string> ids;
    std::string id; double v;
    while (in >> id >> v) ids.push_back(id);
    if (ids.empty()) throw std::runtime_error("Phenotype file has no rows");
    return ids;
}

static void save_Ubin(const std::string& path, const Eigen::MatrixXd& U) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write " + path);
    int n = (int)U.rows();
    int k = (int)U.cols();
    out.write((char*)&n, sizeof(int));
    out.write((char*)&k, sizeof(int));
    out.write((char*)U.data(), sizeof(double)*n*k);
}

static void zscore_columns(Eigen::MatrixXd& M) {
    const Eigen::RowVectorXd mu = M.colwise().mean();
    M.rowwise() -= mu;
    Eigen::RowVectorXd sd = ((M.array().square().colwise().sum())/(M.rows()-1)).sqrt().matrix();
    for (Eigen::Index j=0;j<sd.size();++j) if (sd[j]==0.0) sd[j]=1.0;
    M.array().rowwise() /= sd.array();
}

static Eigen::MatrixXd khatri_rao_upper(const Eigen::MatrixXd& A) {
    const Eigen::Index n = A.rows(), r = A.cols();
    const Eigen::Index R = r*(r+1)/2;
    Eigen::MatrixXd B(n, R);
    Eigen::Index k = 0;
    for (Eigen::Index a=0; a<r; ++a)
        for (Eigen::Index b=a; b<r; ++b)
            B.col(k++) = A.col(a).cwiseProduct(A.col(b));
    return B;
}

static Eigen::MatrixXd orthonormal_basis(const Eigen::MatrixXd& X) {
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU);
    return svd.matrixU();
}

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " X_sample.tsv phenotype.txt r outfolder zscore(0|1)\n";
        return 1;
    }

    auto t0 = Clock::now();
    const std::string Xfile     = argv[1];
    const std::string phenofile = argv[2];
    int r = std::stoi(argv[3]);
    const std::string outFolder = argv[4];
    const bool do_zscore = std::stoi(argv[5]) != 0;

    auto X_with_ids = load_X_matrix(Xfile);
    Eigen::MatrixXd X_all = std::move(X_with_ids.data);
    std::vector<std::string> ids_all = std::move(X_with_ids.ids);
    const int n_all = (int)ids_all.size();
    const int p = (int)X_all.cols();

    auto pheno_ids = load_pheno_ids_ordered(phenofile);
    std::unordered_map<std::string,int> pos; pos.reserve(ids_all.size()*2);
    for (int i=0;i<n_all;++i) pos[ids_all[i]] = i;

    std::vector<int> keep; keep.reserve(pheno_ids.size());
    for (auto& id : pheno_ids) {
        auto it = pos.find(id);
        if (it != pos.end()) keep.push_back(it->second);
    }
    if (keep.empty()) throw std::runtime_error("No overlapping genomes between X and phenotype.");

    const int n = (int)keep.size();
    if (r <= 0 || r > n) throw std::runtime_error("r must be in [1, n]");

    Eigen::MatrixXd X(n, p);
    for (int i=0;i<n;++i) X.row(i) = X_all.row(keep[i]);

    X.rowwise() -= X.colwise().mean();

    auto t1 = Clock::now();
    std::cerr << "[TIME] Loaded & subset X (" << n << "/" << n_all
              << " genomes, " << p << " features) in " << elapsed(t0,t1) << "s\n";

    // PCA in sample space
    Eigen::MatrixXd G = X * X.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(G);
    if (eig.info() != Eigen::Success) throw std::runtime_error("Eigen decomposition failed");

    Eigen::VectorXd evals = eig.eigenvalues();
    Eigen::MatrixXd evecs = eig.eigenvectors();
    Eigen::MatrixXd U_r = evecs.rightCols(r);
    Eigen::VectorXd lam_r(r);
    for (int i=0;i<r;++i) lam_r(i) = evals(evals.size()-r+i);

    Eigen::VectorXd sig_r = lam_r.cwiseMax(0.0).cwiseSqrt();
    Eigen::MatrixXd Ph = U_r * sig_r.asDiagonal();

    if (do_zscore) zscore_columns(Ph);

    Eigen::MatrixXd Phh = khatri_rao_upper(Ph);
    if (do_zscore) zscore_columns(Phh);

    // 1D null: [1, Ph]
    Eigen::MatrixXd Z1D(n, 1 + Ph.cols());
    Z1D.col(0).setOnes();
    Z1D.block(0,1,n,Ph.cols()) = Ph;

    // 2D null: [1, Ph, Phh]
    Eigen::MatrixXd Z2D(n, 1 + Ph.cols() + Phh.cols());
    Z2D.col(0).setOnes();
    Z2D.block(0,1,n,Ph.cols()) = Ph;
    Z2D.block(0,1+Ph.cols(), n, Phh.cols()) = Phh;

    Eigen::MatrixXd U1D  = orthonormal_basis(Z1D);
    Eigen::MatrixXd U2D  = orthonormal_basis(Z2D);

    auto t2 = Clock::now();
    std::cerr << "[TIME] Built U1D (rank=" << U1D.cols()
              << "), U2D (rank=" << U2D.cols()
              << ") in " << elapsed(t1,t2) << "s\n";

    save_Ubin(outFolder + "/U.bin",   U1D);
    save_Ubin(outFolder + "/U2D.bin", U2D);

    auto t3 = Clock::now();
    std::cerr << "[TIME] Saved U.bin and U2D.bin in " << elapsed(t2,t3) << "s\n";
    return 0;
}
