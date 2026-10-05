// compute_structure_2w.cpp — Dual PCA + Kronecker Ph⊗Pp → U0.txt
//
// Usage: compute_structure_2w <host_bits_sub.txt> <r_host> <patho_bits_sub.txt> <r_patho> <phenotype.txt> <out_prefix>

#include <algorithm>
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
    if (argc != 7) {
        std::cerr
            << "Usage: " << argv[0]
            << " <host_bits_sub.txt> <r_host> <patho_bits_sub.txt> <r_patho> <phenotype.txt> <out_prefix>\n"
            << "Notes:\n"
            << "  - host_bits_sub / patho_bits_sub are *subsampled* bitstring files:\n"
            << "        <label><ws><bitstring-of-0/1>\n"
            << "  - phenotype.txt has lines: <host_id> <pathogen_id> (0-based)\n";
        return 1;
    }

    const std::string host_sub_path = argv[1];
    const int         r_host        = std::stoi(argv[2]);
    const std::string patho_sub_path= argv[3];
    const int         r_patho       = std::stoi(argv[4]);
    const std::string phen_path     = argv[5];
    const std::string out_prefix    = argv[6];

    try {
        auto Hs = io::load_kmer_matrix_full(host_sub_path);
        auto Ps = io::load_kmer_matrix_full(patho_sub_path);

        const int nh = static_cast<int>(Hs.G.rows());
        const int np = static_cast<int>(Ps.G.rows());

        la::log_matrix_info("G_h_sub (nh×mh_sub)", Hs.G);
        la::log_matrix_info("G_p_sub (np×mp_sub)", Ps.G);

        std::cerr << "[Host] PCA scores (top " << r_host << ")...\n";
        Mat Ph = la::pc_scores(Hs.Gc, r_host);
        std::cerr << "[Pathogen] PCA scores (top " << r_patho << ")...\n";
        Mat Pp = la::pc_scores(Ps.Gc, r_patho);

        std::cerr << "[Build] Php = Ph ⊗ Pp ...\n";
        Mat Php = la::build_kron_matrix(Ph, Pp);
        std::cerr << "[Build] X0 = [1 Ph Pp Php] ...\n";
        Mat X0  = la::build_X0(Ph, Pp, Php);

        const int r0 = static_cast<int>(std::min<Eigen::Index>(X0.rows(), X0.cols()));
        std::cerr << "[U] Computing left singular vectors of X0 (r0=" << r0 << ")...\n";
        Mat U0 = la::left_U_from_X0(X0, r0);

        std::vector<int> host_ids, path_ids;
        io::read_phenotype_ids(phen_path, host_ids, path_ids);
        if (static_cast<int>(host_ids.size()) != U0.rows())
            throw std::runtime_error("phenotype rows != U0.rows()");

        const std::string u0_path = out_prefix + "U0.txt";
        std::cerr << "[Output] Writing U0 -> " << u0_path << "\n";
        io::write_matrix_txt(u0_path, U0);

        std::cerr << "[done] wrote " << u0_path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 2;
    }
    return 0;
}
