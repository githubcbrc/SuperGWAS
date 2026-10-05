#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <iostream>

#include "chisq.hpp"

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

struct FullGridCaches {
    int nh=0, np=0, n=0, r=0;
    Mat Y;
    std::vector<Mat> Uplanes;
    Vec uy; double RSS0=0.0;
    Vec yh;
    Vec yp;
    Mat Th;
    Mat Tp;
};

inline FullGridCaches build_fullgrid_caches(const Mat& U,
                                            const std::vector<int>& host_id,
                                            const std::vector<int>& path_id,
                                            const Vec& y,
                                            int nh, int np) {
    const int n = static_cast<int>(U.rows());
    const int r = static_cast<int>(U.cols());
    if (y.size() != n) throw std::runtime_error("build_fullgrid_caches: y.size()!=U.rows()");
    if (static_cast<int>(host_id.size()) != n || static_cast<int>(path_id.size()) != n)
        throw std::runtime_error("build_fullgrid_caches: id vectors must have length n");
    if (n != nh * np) throw std::runtime_error("Full grid assumption failed: n != nh*np");

    Mat Y = Mat::Zero(nh, np);
    std::vector<Mat> Uplanes;
    Uplanes.reserve(r);
    for (int j = 0; j < r; ++j) Uplanes.emplace_back(Mat::Zero(nh, np));

    for (int t = 0; t < n; ++t) {
        const int h = host_id[t];
        const int p = path_id[t];
        if (h < 0 || h >= nh || p < 0 || p >= np)
            throw std::runtime_error("build_fullgrid_caches: phenotype id out of range");
        Y(h, p) = y[t];
        for (int j = 0; j < r; ++j) Uplanes[j](h, p) = U(t, j);
    }

    const Vec uy = U.transpose() * y;
    double RSS0 = y.squaredNorm() - uy.squaredNorm();
    if (RSS0 < 0 && RSS0 > -1e-12) RSS0 = 0.0;

    const Vec yh = Y.rowwise().sum();
    const Vec yp = Y.colwise().sum().transpose();

    FullGridCaches C;
    C.nh = nh; C.np = np; C.n = n; C.r = r;
    C.Y = std::move(Y);
    C.Uplanes = std::move(Uplanes);
    C.uy = uy;
    C.RSS0 = RSS0;
    C.yh = yh;
    C.yp = yp;

    return C;
}

struct StatsOut {
    double ratio  = std::numeric_limits<double>::quiet_NaN();
    double delta  = std::numeric_limits<double>::quiet_NaN();
    double LRT    = std::numeric_limits<double>::quiet_NaN();
    double pchi   = std::numeric_limits<double>::quiet_NaN();
};

inline StatsOut interaction_only_stat_full(const FullGridCaches& C,
                                           const Vec& gh, const Vec& gp) {
    const Vec Ygp = C.Y * gp;
    const double xhp_y = gh.dot(Ygp);

    Vec thp(C.r);
    for (int j = 0; j < C.r; ++j) {
        const Vec Ujgp = C.Uplanes[j] * gp;
        thp[j] = gh.dot(Ujgp);
    }

    double qx = gh.squaredNorm() * gp.squaredNorm() - thp.squaredNorm();
    if (qx < 0 && qx > -1e-12) qx = 0.0;

    const double sx = xhp_y - thp.dot(C.uy);

    StatsOut out;
    const double denom = qx * C.RSS0;
    double R2 = 0.0;

    if (denom > 0.0) {
        R2 = (sx * sx) / denom;
        if (R2 < 0 && R2 > -1e-15) R2 = 0.0;
        if (R2 > 1 && R2 < 1 + 1e-12) R2 = 1.0;
    } else {
        out.ratio = std::numeric_limits<double>::infinity();
        out.delta = C.RSS0;
        out.LRT   = std::numeric_limits<double>::infinity();
        out.pchi  = 0.0;
        return out;
    }

    if (R2 >= 1.0) {
        out.ratio = std::numeric_limits<double>::infinity();
        out.delta = C.RSS0;
        out.LRT   = std::numeric_limits<double>::infinity();
        out.pchi  = 0.0;
        return out;
    }

    out.delta = C.RSS0 * R2;
    out.ratio = 1.0 / (1.0 - R2);

    if (std::isfinite(out.ratio) && out.ratio > 0.0) {
        out.LRT = double(C.n) * std::log(out.ratio);
        if (out.LRT < 0 && out.LRT > -1e-12) out.LRT = 0.0;
        out.pchi = chisq_sf(out.LRT, 1.0);
        if (std::isfinite(out.pchi)) {
            if (out.pchi < 0.0) out.pchi = 0.0;
            else if (out.pchi > 1.0) out.pchi = 1.0;
        }
    } else {
        out.LRT  = std::numeric_limits<double>::infinity();
        out.pchi = 0.0;
    }
    return out;
}

static StatsOut joint_3df_stat_full(const FullGridCaches& C,
                                    const Vec& gh, const Vec& gp,
                                    const Vec& th, const Vec& tp) {
    const double xh_y  = gh.dot(C.yh);
    const double xp_y  = gp.dot(C.yp);
    const double xhp_y = gh.dot(C.Y * gp);

    Vec thp(C.r);
    for (int j=0;j<C.r;++j) {
        Vec Ujgp = C.Uplanes[j] * gp;
        thp[j] = gh.dot(Ujgp);
    }

    Eigen::Vector3d svec;
    svec << (xh_y  - th.dot(C.uy)),
            (xp_y  - tp.dot(C.uy)),
            (xhp_y - thp.dot(C.uy));

    const double a = double(C.np) * gh.squaredNorm();
    const double d = double(C.nh) * gp.squaredNorm();
    const double f = gh.squaredNorm() * gp.squaredNorm();

    Eigen::Matrix3d XTX;
    XTX << a, 0.0, 0.0,
           0.0, d, 0.0,
           0.0, 0.0, f;

    const double th_th  = th.squaredNorm();
    const double tp_tp  = tp.squaredNorm();
    const double th_tp  = th.dot(tp);
    const double th_thp = th.dot(thp);
    const double tp_thp = tp.dot(thp);
    const double thp_thp= thp.squaredNorm();

    Eigen::Matrix3d G;
    G << th_th, th_tp, th_thp,
         th_tp, tp_tp, tp_thp,
         th_thp, tp_thp, thp_thp;

    Eigen::Matrix3d S = XTX - G;

    Eigen::FullPivLU<Eigen::Matrix3d> lu(S);
    StatsOut out;
    if (!lu.isInvertible()) return out;

    Eigen::Vector3d x = lu.solve(svec);
    out.delta = svec.dot(x);
    const double RSS1 = C.RSS0 - out.delta;

    if (RSS1 > 0.0) {
        out.ratio = C.RSS0 / RSS1;
        if (out.ratio > 0.0 && std::isfinite(out.ratio)) {
            out.LRT  = double(C.n) * std::log(out.ratio);
            if (out.LRT < 0.0 && out.LRT > -1e-12) out.LRT = 0.0;
            out.pchi = chisq_sf(out.LRT, 3.0);
            if (std::isfinite(out.pchi)) {
                if (out.pchi < 0.0) out.pchi = 0.0;
                else if (out.pchi > 1.0) out.pchi = 1.0;
            }
        } else { out.LRT = std::numeric_limits<double>::infinity(); out.pchi = 0.0; }
    } else {
        out.ratio = std::numeric_limits<double>::infinity();
        out.LRT   = std::numeric_limits<double>::infinity();
        out.pchi  = 0.0;
    }
    return out;
}
