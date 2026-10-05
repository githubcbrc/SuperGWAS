#pragma once
// Contents:
//   - Types: Mat / Vec / Row
//   - Utils: trim()
//   - Matrix I/O: read_matrix_txt(), write_matrix_txt(), write_vector_row()
//   - Phenotypes: read_phenotype_ids()
//   - K-mers: load_kmer_matrix_full()  -> {labels, G, Gc, mu}
//   - Th/Tp columns: load_T_cols()
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace io {

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using Row = Eigen::RowVectorXd;

inline std::string trim(std::string s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    s.erase(b);
    s.erase(0, a);
    return s;
}

inline Mat read_matrix_txt(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    std::vector<double> data;
    data.reserve(1 << 20);

    std::string line;
    std::size_t ncols = 0, nrows = 0;

    while (std::getline(in, line)) {
        line = trim(line);
        std::istringstream ss(line);
        std::size_t cols_this = 0;
        for (double x; ss >> x; ) {
            data.push_back(x);
            ++cols_this;
        }
        if (cols_this == 0) continue;

        if (ncols == 0) ncols = cols_this;
        else if (cols_this != ncols)
            throw std::runtime_error("Inconsistent columns in " + path +
                                     " at row " + std::to_string(nrows));
        ++nrows;
    }

    if (ncols == 0) return Mat();

    using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    RowMajorMat tmp(static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(ncols));
    std::copy(data.begin(), data.end(), tmp.data());
    return Mat(tmp);
}

inline void write_matrix_txt(const std::string& path, const Mat& X) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write: " + path);
    out.setf(std::ios::fmtflags(0), std::ios::floatfield);
    out.precision(10);

    const Eigen::Index R = X.rows(), C = X.cols();
    for (Eigen::Index i = 0; i < R; ++i) {
        for (Eigen::Index j = 0; j < C; ++j) {
            if (j) out.put(' ');
            out << X(i, j);
        }
        out.put('\n');
    }
}

inline void write_vector_row(std::ostream& out, const Vec& v) {
    for (Eigen::Index j = 0; j < v.size(); ++j) {
        if (j) out.put(' ');
        out << v[j];
    }
    out.put('\n');
}

inline void read_phenotype_ids(const std::string& path,
                               std::vector<int>& host_ids,
                               std::vector<int>& path_ids) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    host_ids.clear();
    path_ids.clear();

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        std::istringstream iss(line);
        int h, p;
        if (!(iss >> h >> p))
            throw std::runtime_error("Bad phenotype line (need host_id pathogen_id): " + line);

        host_ids.push_back(h);
        path_ids.push_back(p);
    }
}

struct KmerMatrixFull {
    std::vector<std::string> labels;
    Mat G;
    Mat Gc;
    Row mu;
};

inline KmerMatrixFull load_kmer_matrix_full(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    std::vector<std::string> labels; labels.reserve(1 << 16);
    std::vector<std::string> bitcols; bitcols.reserve(1 << 16);

    std::string line;
    int bitlen = -1;

    while (std::getline(in, line)) {
        line = trim(line);
        std::istringstream ss(line);
        std::string lab, bits;
        if (!(ss >> lab >> bits)) continue;

        if (bitlen < 0) bitlen = static_cast<int>(bits.size());
        else if (static_cast<int>(bits.size()) != bitlen)
            throw std::runtime_error("Bitstring length mismatch in " + path);

        labels.push_back(std::move(lab));
        bitcols.push_back(std::move(bits));
    }

    if (bitlen < 0)
        throw std::runtime_error("No data lines in: " + path);

    const int n = bitlen;
    const int m = static_cast<int>(bitcols.size());

    Mat G(n, m);
    for (int j = 0; j < m; ++j) {
        const std::string& bits = bitcols[j];
        for (int i = 0; i < n; ++i)
            G(i, j) = (bits[static_cast<std::size_t>(i)] == '1') ? 1.0 : 0.0;
    }

    Row mu = G.colwise().mean();
    Mat Gc = G; Gc.rowwise() -= mu;

    return KmerMatrixFull{std::move(labels), std::move(G), std::move(Gc), std::move(mu)};
}

inline Mat load_T_cols(const std::string& path, int r_expected, std::size_t m_expected) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open: " + path);

    std::vector<double> data;
    data.reserve((r_expected > 0 ? r_expected : 128) * (m_expected ? m_expected : 1));

    std::string line;
    std::size_t rows_seen = 0;
    std::size_t r = 0;

    while (std::getline(in, line)) {
        line = trim(line);
        std::istringstream ss(line);
        std::size_t cnt = 0;
        for (double x; ss >> x; ) { data.push_back(x); ++cnt; }
        if (cnt == 0) continue;

        if (r == 0) r = cnt;
        else if (cnt != r)
            throw std::runtime_error("Inconsistent r across lines in " + path);

        ++rows_seen;
    }

    if (m_expected && rows_seen != m_expected)
        throw std::runtime_error("Th/Tp line count != number of k-mers in " + path +
                                 " (" + std::to_string(rows_seen) + " vs " + std::to_string(m_expected) + ")");

    if (r_expected > 0 && r != static_cast<std::size_t>(r_expected))
        throw std::runtime_error("Th/Tp column length (r) mismatch in " + path);

    if (rows_seen == 0) return Mat();

    using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    RowMajorMat tmp(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(rows_seen));
    std::copy(data.begin(), data.end(), tmp.data());
    return Mat(tmp);
}

} // namespace io
