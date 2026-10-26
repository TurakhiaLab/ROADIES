#include "io/parse_file.hpp"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string/trim.hpp>
#include <libpll/pll.h>

namespace parse {

namespace {

using PllMsaPtr = std::unique_ptr<pll_msa_t, decltype(&pll_msa_destroy)>;

Alignment copy_pll_alignment(const pll_msa_t& msa)
{
    if (msa.count <= 0 || msa.length <= 0 || !msa.label || !msa.sequence) {
        throw std::runtime_error("Alignment contains no sequences or sites.");
    }

    Alignment alignment;
    alignment.names.reserve(static_cast<size_t>(msa.count));
    alignment.sequences.reserve(static_cast<size_t>(msa.count));
    for (int index = 0; index < msa.count; ++index) {
        if (!msa.label[index] || !msa.sequence[index]) {
            throw std::runtime_error("Alignment contains a missing name or sequence.");
        }
        alignment.names.emplace_back(msa.label[index]);
        alignment.sequences.emplace_back(msa.sequence[index]);
    }
    alignment.sites = static_cast<size_t>(msa.length);
    return alignment;
}

PllMsaPtr load_phylip_alignment(const std::string& path)
{
    // libpll requires the caller to select sequential or interleaved parsing.
    // Try both so the public reader supports both standard PHYLIP layouts.
    PllMsaPtr msa(
        pll_phylip_load(path.c_str(), PLL_FALSE),
        &pll_msa_destroy);
    if (!msa) {
        msa.reset(pll_phylip_load(path.c_str(), PLL_TRUE));
    }
    return msa;
}

} // namespace

Alignment read_alignment_file(const std::string& path) {
    std::ifstream probe(path);
    if (!probe) {
        throw std::runtime_error("Cannot open alignment file: " + path);
    }

    std::string first_line;
    std::string trimmed;
    while (std::getline(probe, first_line)) {
        trimmed = boost::algorithm::trim_copy(first_line);
        if (!trimmed.empty()) break;
    }
    if (trimmed.empty()) {
        throw std::runtime_error("Alignment file is empty.");
    }

    PllMsaPtr msa(
        trimmed.front() == '>'
            ? pll_fasta_load(path.c_str())
            : nullptr,
        &pll_msa_destroy);
    if (trimmed.front() != '>') {
        msa = load_phylip_alignment(path);
    }
    if (!msa) {
        const std::string parser_message = pll_errmsg[0]
            ? std::string(pll_errmsg)
            : std::string("unknown parser error");
        throw std::runtime_error(
            "Cannot parse alignment file '" + path + "': " + parser_message);
    }
    return copy_pll_alignment(*msa);
}

} // namespace parse
