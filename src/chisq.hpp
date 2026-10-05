#pragma once
#include <cmath>
#include <limits>

// Fast 1-df chi-square survival function (used by 1D and 2D)
inline double chi2_sf_k1(double T) {
    if (T <= 0.0) return 1.0;
    return std::erfc(std::sqrt(0.5 * T));
}

// Fast 3-df chi-square survival function (used by 2D full interaction)
inline double chi2_sf_k3(double T) {
    if (T <= 0.0) return 1.0;
    double h = std::sqrt(0.5 * T);
    double e = std::exp(-0.5 * T);
    return std::erfc(h) + e * h * (2.0 / std::sqrt(M_PI));
}

// General incomplete gamma series
inline void gser(double& gamser, double a, double x, double& gln) {
    const int ITMAX = 200;
    const double EPS = 3e-14;
    gln = std::lgamma(a);
    if (x <= 0.0) { gamser = 0.0; return; }

    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= ITMAX; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::fabs(del) < std::fabs(sum) * EPS) break;
    }
    gamser = sum * std::exp(-x + a * std::log(x) - gln);
}

// General incomplete gamma continued fraction
inline void gcf(double& gammcf, double a, double x, double& gln) {
    const int ITMAX = 200;
    const double EPS = 3e-14;
    const double FPMIN =
        std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();

    gln = std::lgamma(a);

    double b = x + 1.0 - a;
    if (std::fabs(b) < FPMIN) b = (b >= 0 ? FPMIN : -FPMIN);
    double c = b;
    double d = 1.0 / b;
    double h = d;

    for (int i = 1; i <= ITMAX; ++i) {
        double an = -static_cast<double>(i) * (i - a);
        b += 2.0;
        double denom = an * d + b;
        if (std::fabs(denom) < FPMIN) denom = (denom >= 0 ? FPMIN : -FPMIN);
        d = 1.0 / denom;
        double cden = c;
        if (std::fabs(cden) < FPMIN) cden = (cden >= 0 ? FPMIN : -FPMIN);
        c = b + an / cden;
        if (std::fabs(c) < FPMIN) c = (c >= 0 ? FPMIN : -FPMIN);
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }
    gammcf = std::exp(-x + a * std::log(x) - gln) * h;
}

// Regularized upper incomplete gamma Q(a,x)
inline double gammq(double a, double x) {
    if (!(x >= 0.0) || a <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    double gln, gamser, gammcf;
    if (x < a + 1.0) { gser(gamser, a, x, gln); return 1.0 - gamser; }
    else             { gcf(gammcf, a, x, gln); return        gammcf; }
}

// General chi-square survival function: P(X > x) where X ~ chi2(k)
inline double chisq_sf(double x, double k) {
    if (!(x >= 0.0) || k <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return gammq(0.5 * k, 0.5 * x);
}
