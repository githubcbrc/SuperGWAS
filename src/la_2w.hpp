#pragma once
#include <Eigen/Dense>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <sstream>

namespace la {

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using Row = Eigen::RowVectorXd;

inline void print_shape(const char* name, Eigen::Index rows, Eigen::Index cols) {
    std::ostringstream oss;
    oss << "  [shape] " << name << ": " << rows << " x " << cols << '\n';
    std::cerr << oss.str();
}

template <typename Derived>
inline void log_matrix_info(const char* name, const Eigen::MatrixBase<Derived>& M) {
    std::ostringstream oss;
    oss << "[dims] " << name << ": " << M.rows() << " x " << M.cols();
    if (M.size() > 0) {
        const auto mn = M.minCoeff();
        const auto mx = M.maxCoeff();
        oss << "  (min=" << mn << ", max=" << mx << ")";
    }
    oss << (M.allFinite() ? "  [finite]\n" : "  [NON-FINITE]\n");
    std::cerr << oss.str();
}

inline void center_columns(Mat& M) {
    const Row mu = M.colwise().mean();
    M.rowwise() -= mu;
    std::cerr << "[center] columns centered\n";
}

inline Mat pc_scores(const Mat& M_in, int r, bool centered=true) {
    if (r <= 0) throw std::runtime_error("pc_scores: r must be positive");
    Mat M = M_in;
    if (!centered) center_columns(M);

    const Eigen::Index n = M.rows();
    const Eigen::Index p = M.cols();
    const int r_eff = std::min<int>(r, static_cast<int>(std::min(n, p)));
    if (r_eff == 0) throw std::runtime_error("pc_scores: effective rank is zero");

    std::cerr << "[PCA] computing thin SVD of centered M\n";
    print_shape("M", n, p);
    Eigen::BDCSVD<Mat> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svd.info() != Eigen::Success) throw std::runtime_error("pc_scores: SVD failed");

    const Mat U_r = svd.matrixU().leftCols(r_eff);
    const Vec S_r = svd.singularValues().head(r_eff);

    Mat Scores = U_r * S_r.asDiagonal();
    return Scores;
}

inline Mat build_kron_matrix(const Mat& Ph, const Mat& Pp) {
    const Eigen::Index n_h = Ph.rows(), r_h = Ph.cols();
    const Eigen::Index n_p = Pp.rows(), r_p = Pp.cols();
    Mat Php(n_h * n_p, r_h * r_p);

    for (Eigen::Index a = 0; a < r_h; ++a) {
        const auto u = Ph.col(a);
        for (Eigen::Index b = 0; b < r_p; ++b) {
            const auto v = Pp.col(b);
            const Eigen::Index col = a * r_p + b;
            for (Eigen::Index ih = 0; ih < n_h; ++ih) {
                Php.block(ih * n_p, col, n_p, 1).noalias() = u(ih) * v;
            }
        }
    }
    return Php;
}

inline Mat build_X0(const Mat& Ph, const Mat& Pp, const Mat& Php) {
    const Eigen::Index n_h = Ph.rows(),  r_h = Ph.cols();
    const Eigen::Index n_p = Pp.rows(),  r_p = Pp.cols();
    if (Php.rows() != n_h * n_p) throw std::runtime_error("build_X0: Php.rows() mismatch");
    if (Php.cols() != r_h * r_p) throw std::runtime_error("build_X0: Php.cols() mismatch");

    const Eigen::Index n_rows = Php.rows();
    const Eigen::Index n_cols = 1 + r_h + r_p + r_h * r_p;
    Mat X0(n_rows, n_cols);
    X0.setZero();

    const Eigen::Index off_1  = 0;
    const Eigen::Index off_Ph = off_1  + 1;
    const Eigen::Index off_Pp = off_Ph + r_h;
    const Eigen::Index off_K  = off_Pp + r_p;

    X0.col(off_1).setOnes();

    Eigen::Index row = 0;
    for (Eigen::Index ih = 0; ih < n_h; ++ih) {
        for (Eigen::Index ip = 0; ip < n_p; ++ip, ++row) {
            X0.block(row, off_Ph, 1, r_h) = Ph.row(ih);
            X0.block(row, off_Pp, 1, r_p) = Pp.row(ip);
            X0.block(row, off_K,  1, r_h * r_p) = Php.row(row);
        }
    }
    return X0;
}

inline Mat left_U_from_X0(const Mat& X0, int r) {
    if (r <= 0) throw std::runtime_error("left_U_from_X0: r must be positive");
    Mat Xw = X0;
    center_columns(Xw);

    const Eigen::Index n = Xw.rows();
    const Eigen::Index p = Xw.cols();
    const int r_eff = std::min<int>(r, static_cast<int>(std::min(n, p)));
    if (r_eff == 0) throw std::runtime_error("left_U_from_X0: effective rank is zero");

    Eigen::BDCSVD<Mat> svd(Xw, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svd.info() != Eigen::Success) throw std::runtime_error("left_U_from_X0: SVD failed");

    return svd.matrixU().leftCols(r_eff);
}

inline void build_Sh_Sp(const Mat& U,
                        const std::vector<int>& host_ids,
                        const std::vector<int>& path_ids,
                        int nh, int np,
                        Mat& Sh, Mat& Sp) {
    const int n = static_cast<int>(U.rows());
    const int r = static_cast<int>(U.cols());
    if (static_cast<int>(host_ids.size()) != n || static_cast<int>(path_ids.size()) != n)
        throw std::runtime_error("build_Sh_Sp: host_ids/path_ids size mismatch");

    Sh = Mat::Zero(r, nh);
    Sp = Mat::Zero(r, np);
    for (int t = 0; t < n; ++t) {
        const int ih = host_ids[t];
        const int ip = path_ids[t];
        if (ih < 0 || ih >= nh) throw std::runtime_error("host_id out of range at obs " + std::to_string(t));
        if (ip < 0 || ip >= np) throw std::runtime_error("pathogen_id out of range at obs " + std::to_string(t));
        Sh.col(ih).noalias() += U.row(t).transpose();
        Sp.col(ip).noalias() += U.row(t).transpose();
    }
}

inline void zscore_columns(Mat& M) {
    if (M.cols() == 0) return;
    const Row mu = M.colwise().mean();
    Mat C = M.rowwise() - mu;
    Row sd(M.cols());
    for (Eigen::Index j = 0; j < M.cols(); ++j) {
        double s2 = C.col(j).array().square().mean();
        double s  = std::sqrt(std::max(0.0, s2));
        sd(j) = (s > 0.0 ? s : 1.0);
    }
    M = C.array().rowwise() / sd.array();
}

inline Mat khatri_rao_upper(const Mat& Ph, bool cross_only=false) {
    const Eigen::Index n = Ph.rows();
    const Eigen::Index r = Ph.cols();
    if (r == 0) return Mat(n, 0);

    const Eigen::Index R = cross_only ? (r * (r - 1)) / 2 : (r * (r + 1)) / 2;
    Mat Phh(n, R);

    Eigen::Index c = 0;
    for (Eigen::Index a = 0; a < r; ++a) {
        for (Eigen::Index b = a; b < r; ++b) {
            if (cross_only && a == b) continue;
            Phh.col(c++) = Ph.col(a).array() * Ph.col(b).array();
        }
    }
    return Phh;
}

inline Mat build_X0_from_Ph_Phh(const Mat& Ph, const Mat& Phh) {
    const Eigen::Index n  = Ph.rows();
    if (Phh.rows() != 0 && Phh.rows() != n)
        throw std::runtime_error("build_X0_from_Ph_PhH: row mismatch");

    const Eigen::Index r  = Ph.cols();
    const Eigen::Index R  = Phh.cols();
    Mat X0(n, 1 + r + R);
    X0.col(0).setOnes();
    if (r > 0)  X0.block(0, 1,           n, r) = Ph;
    if (R > 0)  X0.block(0, 1 + r,       n, R) = Phh;
    return X0;
}

} // namespace la
