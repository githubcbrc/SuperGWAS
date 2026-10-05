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

MatrixWithIDs load_X_matrix(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open " + filename);

    std::string line;
    std::vector<std::vector<double>> rows;
    std::vector<std::string> ids;

    if (!std::getline(file, line)) throw std::runtime_error("Empty X file"); // skip header

    size_t expected_cols = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok; std::vector<double> row;

        std::getline(ss, tok, '\t');  // genome ID
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

// SVD -> thin U (orthonormal basis of span(X))
static Eigen::MatrixXd orthonormal_basis(const Eigen::MatrixXd& X) {
    Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU);
    return svd.matrixU(); // n x rank
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: ./precompute X_sample.tsv phenotype.txt r outfolder\n";
        return 1;
    }

    auto t0 = Clock::now();
    const std::string Xfile = argv[1];
    const std::string phenofile = argv[2];
    int r = std::stoi(argv[3]);
    const std::string outFolder = argv[4];

    // 1) Load X (all genomes)
    auto X_with_ids = load_X_matrix(Xfile);
    Eigen::MatrixXd X_all = std::move(X_with_ids.data);
    std::vector<std::string> ids_all = std::move(X_with_ids.ids);
    const int n_all = (int)ids_all.size();
    const int p = (int)X_all.cols();

    // 2) Subset X to phenotype genomes, in phenotype order
    auto pheno_ids = load_pheno_ids_ordered(phenofile);
    std::unordered_map<std::string,int> pos; pos.reserve(ids_all.size()*2);
    for (int i=0;i<n_all;++i) pos[ids_all[i]] = i;

    std::vector<int> keep; keep.reserve(pheno_ids.size());
    for (auto& id : pheno_ids) {
        auto it = pos.find(id);
        if (it != pos.end()) { keep.push_back(it->second); }
    }
    if (keep.empty()) throw std::runtime_error("No overlapping genomes between X and phenotype.");

    const int n = (int)keep.size();
    if (r <= 0 || r > n) throw std::runtime_error("r must be in [1, n]");

    Eigen::MatrixXd X(n, p);
    for (int i=0;i<n;++i) X.row(i) = X_all.row(keep[i]);

    // 3) Center columns on the analysis subset
    X.rowwise() -= X.colwise().mean();

    auto t1 = Clock::now();
    std::cerr << "[TIME] Loaded & subset X (" << n << "/" << n_all
              << " genomes, " << p << " features) in " << elapsed(t0,t1) << "s\n";

    // 4) PCA in sample space: eig of X X^T
    Eigen::MatrixXd G = X * X.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(G);
    if (eig.info() != Eigen::Success) throw std::runtime_error("Eigen decomposition failed");

    Eigen::VectorXd evals = eig.eigenvalues();      
    Eigen::MatrixXd evecs = eig.eigenvectors();     
    Eigen::MatrixXd U_r = evecs.rightCols(r);      
    Eigen::VectorXd lam_r(r);
    for (int i=0;i<r;++i) lam_r(i) = evals(evals.size()-r+i);

    // 5) Scores S = U_r * diag(sqrt(lam_r))  (n x r)
    Eigen::VectorXd sig_r = lam_r.cwiseMax(0.0).cwiseSqrt();
    Eigen::MatrixXd S = U_r * sig_r.asDiagonal();

    // 6) Build Z = [1, S] and compute its SVD to get an orthonormal basis of span(Z) --> Uz
    Eigen::MatrixXd Z(n, r + 1);
    Z.col(0).setOnes();
    Z.block(0, 1, n, r) = S;

    Eigen::MatrixXd U_final = orthonormal_basis(Z);  
    int rankZ = U_final.cols();

    auto t2 = Clock::now();
    std::cerr << "[TIME] Built U via SVD (rank=" << rankZ << ") in "
              << elapsed(t1,t2) << "s\n";

    // 7) Save outputs
    {
        // Save U_final as U.bin (includes intercept + PCs)
        std::ofstream uout(outFolder + "/U.bin", std::ios::binary);
        if (!uout) throw std::runtime_error("Cannot write "+outFolder+"/U.bin");
        int k = U_final.cols();
        uout.write((char*)&n, sizeof(int));
        uout.write((char*)&k, sizeof(int));
        uout.write((char*)U_final.data(), sizeof(double)*n*k);
    }

    auto t3 = Clock::now();
    std::cerr << "[TIME] Saved outputs in " << elapsed(t2,t3) << "s\n";
    return 0;
}

