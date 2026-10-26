#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>

namespace parse {
struct Alignment;
struct ModelConfig;
}

namespace mlipper {
namespace input {

// User-facing configuration/input failures. CLI entry points catch this family
// and report the message without treating it as an internal program defect.
class CliError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ValidationError : public CliError {
public:
    ValidationError(const std::string& option_name, const std::string& message)
        : CliError(option_name + ": " + message) {}
};

// Resolve raw_path against base, returning a normalized absolute path. An empty
// raw path remains empty so the option-specific validator can report it.
std::filesystem::path normalize_cli_path(
    const std::filesystem::path& base,
    const std::string& raw_path);

// Check that an output names a file and that its nearest existing ancestor is
// a directory. The function validates only; it does not create directories.
void validate_output_path(
    const std::filesystem::path& base,
    const std::string& option_name,
    const std::string& raw_path);

// Require a one-to-one names/sequences mapping with non-empty unique names.
void validate_alignment_names(
    const parse::Alignment& alignment,
    const std::string& option_name);

// Require equal row lengths and symbols supported by the selected state model.
void validate_alignment_symbols(
    const parse::Alignment& alignment,
    int states,
    const std::string& option_name);

// Enforce the public MLIPPER model scope before any CPU/GPU allocation.
void validate_model_inputs(const parse::ModelConfig& model);

// Reject taxa present in both alignments because commit mode requires each
// query to create a new tree tip.
void validate_query_reference_name_overlap(
    const parse::Alignment& tree_alignment,
    const parse::Alignment& query_alignment,
    const std::string& option_name);

} // namespace input
} // namespace mlipper
