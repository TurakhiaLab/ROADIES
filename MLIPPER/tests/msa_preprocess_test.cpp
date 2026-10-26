#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "util/msa_preprocess.hpp"

std::vector<mlipper::SequenceRecord> build_placement_query(
    const std::vector<std::string>& names,
    const std::vector<std::string>& sequences)
{
    if (names.size() != sequences.size()) {
        throw std::runtime_error("Alignment names/sequences size mismatch.");
    }
    std::vector<mlipper::SequenceRecord> queries;
    queries.reserve(names.size());
    for (size_t index = 0; index < names.size(); ++index) {
        mlipper::SequenceRecord query;
        query.name = names[index];
        query.sequence = sequences[index];
        queries.push_back(std::move(query));
    }
    return queries;
}

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void testCombinedReferenceAndQueryPatterns()
{
    const parse::Alignment reference{
        {"r1", "r2"},
        {"AACGT", "CCGTT"},
        5};
    const parse::Alignment queries{{"q1"}, {"GGAAA"}, 5};

    const PreprocessedAlignments result =
        preprocess_alignments(reference, queries);

    require(result.tree_alignment.sites == 4, "duplicate site was not removed");
    require(result.query_alignment.sites == 4, "query site count was not updated");
    require(result.tree_alignment.sequences ==
            std::vector<std::string>({"ACGT", "CGTT"}),
            "reference columns were not restored to original order");
    require(result.query_alignment.sequences ==
            std::vector<std::string>({"GAAA"}),
            "query columns do not match compressed reference columns");
    require(result.pattern_weights ==
            std::vector<unsigned>({2u, 1u, 1u, 1u}),
            "duplicate-site weight was not accumulated");
}

void testQueryPatternPreventsReferenceOnlyMerge()
{
    const parse::Alignment reference{
        {"r1", "r2"},
        {"AACGT", "CCGTT"},
        5};
    const parse::Alignment queries{
        {"q1", "q2"},
        {"GGAAA", "GTCCC"},
        5};

    const PreprocessedAlignments result =
        preprocess_alignments(reference, queries);
    require(result.tree_alignment.sites == 5,
            "columns with different query states were incorrectly merged");
    require(result.pattern_weights ==
            std::vector<unsigned>({1u, 1u, 1u, 1u, 1u}),
            "unexpected weights for unique combined patterns");
}

void testInvalidInputsFail()
{
    bool threw = false;
    try {
        std::vector<std::string> rows{"AA"};
        std::vector<mlipper::SequenceRecord> queries;
        std::vector<unsigned> weights{1u};
        size_t sites = 2;
        remove_repetitive_columns(rows, queries, weights, sites);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "pattern-weight size mismatch was accepted");

    threw = false;
    try {
        std::vector<std::string> rows{"A"};
        std::vector<mlipper::SequenceRecord> queries;
        std::vector<unsigned> weights{1u, 1u};
        size_t sites = 2;
        remove_repetitive_columns(rows, queries, weights, sites);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "short reference sequence was accepted");

    threw = false;
    try {
        std::vector<std::string> rows{"AA"};
        std::vector<mlipper::SequenceRecord> queries{{"q1", "A"}};
        std::vector<unsigned> weights{1u, 1u};
        size_t sites = 2;
        remove_repetitive_columns(rows, queries, weights, sites);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "short query sequence was accepted");

    threw = false;
    try {
        const parse::Alignment empty{{"r1"}, {""}, 0};
        (void)preprocess_alignments(empty);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "zero-site alignment was accepted");
}

} // namespace

int main()
{
    try {
        testCombinedReferenceAndQueryPatterns();
        testQueryPatternPreventsReferenceOnlyMerge();
        testInvalidInputsFail();
        std::cout << "msa_preprocess_test: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "msa_preprocess_test: FAIL: " << error.what() << "\n";
        return 1;
    }
}
