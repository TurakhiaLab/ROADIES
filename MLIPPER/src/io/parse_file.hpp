#pragma once

#include <string>
#include <vector>

namespace parse {

struct Alignment {
    // names[i] and sequences[i] describe one row; every sequence has sites
    // characters after input validation and preprocessing.
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    size_t sites = 0;
};

struct ModelConfig {
    // freqs uses state order A,C,G,T for DNA. rates contains the six reversible
    // GTR exchangeabilities AC, AG, AT, CG, CT, GT. ncat/alpha describe the
    // mean-discretized Gamma categories consumed by MLIPPER.
    int states = 4;
    std::string subst_model;
    int ncat = 1;
    double alpha = 1.0;
    double pinv = 0.0;
    std::vector<double> freqs;
    std::vector<double> rates;
    bool per_rate_scaling = true;
};

// Read FASTA or sequential/interleaved PHYLIP based on file content. This is a
// syntactic reader; workflow-specific name, symbol, and model checks are done
// by input_validation.
Alignment read_alignment_file(const std::string& path);

} // namespace parse
