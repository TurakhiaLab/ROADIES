#include <cassert>
#include <string>

#include "io/input_validation.hpp"
#include "io/parse_file.hpp"

namespace {

template <typename Function>
void expect_validation_error(Function&& function)
{
    try {
        function();
    } catch (const mlipper::input::ValidationError&) {
        return;
    }
    assert(false && "Expected input validation to fail");
}

} // namespace

int main()
{
    parse::ModelConfig model;
    model.states = 4;
    model.ncat = 4;
    model.subst_model = "GTR";
    model.alpha = 0.5;
    model.pinv = 0.0;
    model.freqs = {0.25, 0.25, 0.25, 0.25};
    model.rates = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    mlipper::input::validate_model_inputs(model);

    model.alpha = 0.019;
    expect_validation_error([&model] {
        mlipper::input::validate_model_inputs(model);
    });
    model.alpha = 0.5;

    model.states = 5;
    expect_validation_error([&model] {
        mlipper::input::validate_model_inputs(model);
    });
    model.states = 4;

    model.subst_model = "HKY";
    expect_validation_error([&model] {
        mlipper::input::validate_model_inputs(model);
    });
    model.subst_model = "GTR";

    model.pinv = 0.1;
    expect_validation_error([&model] {
        mlipper::input::validate_model_inputs(model);
    });

    parse::Alignment dna_alignment;
    dna_alignment.names = {"sample"};
    dna_alignment.sequences = {"ACGTURYSWKMBDHVN-.?X"};
    dna_alignment.sites = dna_alignment.sequences.front().size();
    mlipper::input::validate_alignment_symbols(
        dna_alignment, 4, "--tree-alignment");
    dna_alignment.sequences.front().back() = '0';
    expect_validation_error([&dna_alignment] {
        mlipper::input::validate_alignment_symbols(
            dna_alignment, 4, "--tree-alignment");
    });

    parse::Alignment mismatched_alignment;
    mismatched_alignment.names = {"sample"};
    mismatched_alignment.sequences = {};
    mismatched_alignment.sites = 1;
    expect_validation_error([&mismatched_alignment] {
        mlipper::input::validate_alignment_symbols(
            mismatched_alignment, 4, "--tree-alignment");
    });

    parse::Alignment short_alignment;
    short_alignment.names = {"sample"};
    short_alignment.sequences = {"A"};
    short_alignment.sites = 2;
    expect_validation_error([&short_alignment] {
        mlipper::input::validate_alignment_symbols(
            short_alignment, 4, "--tree-alignment");
    });

    return 0;
}
