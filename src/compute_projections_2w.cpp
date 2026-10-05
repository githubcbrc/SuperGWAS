// compute_projections_2w.cpp — Factored T projection via build_Sh_Sp
//
// Usage: compute_projections_2w <U0.txt> <phenotype.txt> <host_bits_full.txt> <patho_bits_full.txt> <out_prefix>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "io_2w.hpp"
#include "la_2w.hpp"

using la::Mat;
using la::Vec;

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr
            << "Usage: " << argv[0]
            << " <U0.txt> <phenotype.txt> <host_bits_full.txt> <patho_bits_full.txt> <out_prefix>\n"
            << "Notes:\n"
            << "  - U0 is from compute_structure_2w (rows must match phenotype rows).\n"
            << "  - host_bits_full / patho_bits_full are FULL bitstring files:\n"
            << "        <label><ws><bitstring-of-0/1>\n";
        return 1;
    }

    const std::string U0_path   = argv[1];
    const std::string phen_path = argv[2];
    const std::string host_path = argv[3];
    const std::string patho_path= argv[4];
    const std::string out_prefix= argv[5];

    try {
        Mat U0 = io::read_matrix_txt(U0_path);
        std::vector<int> host_ids, path_ids;
        io::read_phenotype_ids(phen_path, host_ids, path_ids);

        const int n  = static_cast<int>(U0.rows());
        const int r0 = static_cast<int>(U0.cols());
        if (static_cast<int>(host_ids.size()) != n)
            throw std::runtime_error("phenotype rows != U0.rows()");
        la::log_matrix_info("U0", U0);

        auto H = io::load_kmer_matrix_full(host_path);
        auto P = io::load_kmer_matrix_full(patho_path);
        const int nh = static_cast<int>(H.G.rows());
        const int np = static_cast<int>(P.G.rows());
        const int mh = static_cast<int>(H.G.cols());
        const int mp = static_cast<int>(P.G.cols());

        Mat Sh, Sp;
        la::build_Sh_Sp(U0, host_ids, path_ids, nh, np, Sh, Sp);

        Vec u_sum = U0.colwise().sum().transpose();
        Mat Th = Sh * H.G;
        Mat Tp = Sp * P.G;
        Th.noalias() -= (u_sum * H.mu);
        Tp.noalias() -= (u_sum * P.mu);

        la::log_matrix_info("Th (r×mh)", Th);
        la::log_matrix_info("Tp (r×mp)", Tp);

        const std::string th_path = out_prefix + "Th.txt";
        const std::string tp_path = out_prefix + "Tp.txt";
        std::ofstream fh(th_path), fp(tp_path);
        if (!fh) throw std::runtime_error("Cannot write: " + th_path);
        if (!fp) throw std::runtime_error("Cannot write: " + tp_path);
        fh.setf(std::ios::fmtflags(0), std::ios::floatfield); fh.precision(10);
        fp.setf(std::ios::fmtflags(0), std::ios::floatfield); fp.precision(10);

        for (int k = 0; k < mh; ++k) io::write_vector_row(fh, Th.col(k));
        for (int k = 0; k < mp; ++k) io::write_vector_row(fp, Tp.col(k));

        std::cerr << "[done]\n"
                  << "  Wrote:\n"
                  << "   " << th_path << "\n"
                  << "   " << tp_path << "\n";

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 2;
    }
    return 0;
}
