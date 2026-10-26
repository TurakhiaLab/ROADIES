#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include "util/model_utils.hpp"

namespace {

std::filesystem::path write_model_file(
    const std::filesystem::path& directory,
    const std::string& filename,
    const std::string& model)
{
    const std::filesystem::path path = directory / filename;
    std::ofstream output(path);
    output << model << '\n';
    output.close();
    assert(output.good());
    return path;
}

void expect_parse_error(
    const std::filesystem::path& path,
    const std::string& expected_message)
{
    try {
        (void)mlipper::model::parse_best_model_file(path);
    } catch (const std::runtime_error& error) {
        assert(std::string(error.what()).find(expected_message) != std::string::npos);
        return;
    }
    assert(false && "Expected bestModel parsing to fail");
}

} // namespace

int main()
{
    const std::filesystem::path test_directory =
        std::filesystem::temp_directory_path() /
        ("mlipper_model_utils_test_" + std::to_string(::getpid()));
    std::filesystem::create_directories(test_directory);

    const auto gamma_path = write_model_file(
        test_directory,
        "gamma.bestModel",
        "GTR{1/2/3/4/5/6}+FO{0.1/0.2/0.3/0.4}+G4m{0.5}");
    const auto gamma = mlipper::model::parse_best_model_file(gamma_path);
    assert(gamma.model.subst_model == "GTR");
    assert(gamma.model.ncat == 4);
    assert(gamma.model.alpha == 0.5);

    const auto free_rate_path = write_model_file(
        test_directory,
        "free_rate.bestModel",
        "GTR{1/2/3/4/5/6}+R4{0.1/0.5/1/2}{0.25/0.25/0.25/0.25}");
    expect_parse_error(free_rate_path, "FreeRate (+R) models are not supported");

    const auto mixture_path = write_model_file(
        test_directory,
        "mixture.bestModel",
        "GTR{1/2/3/4/5/6}+MIXTURE{model-a/model-b}");
    expect_parse_error(mixture_path, "Mixture models are not supported");

    const auto unknown_modifier_path = write_model_file(
        test_directory,
        "unknown_modifier.bestModel",
        "GTR{1/2/3/4/5/6}+FO{0.1/0.2/0.3/0.4}+ASC+G4m{0.5}");
    expect_parse_error(
        unknown_modifier_path, "Unsupported bestModel syntax or modifier");

    const auto invalid_rate_path = write_model_file(
        test_directory,
        "invalid_rate.bestModel",
        "GTR{1/2/3/4/5/nan}+G4m{0.5}");
    expect_parse_error(invalid_rate_path, "GTR rates");

    const auto invalid_alpha_path = write_model_file(
        test_directory,
        "invalid_alpha.bestModel",
        "GTR{1/2/3/4/5/6}+G4m{0}");
    expect_parse_error(invalid_alpha_path, "Gamma alpha");

    for (double alpha : {0.02, 0.1, 0.5, 1.0, 10.0}) {
        for (int categories : {1, 4, 8}) {
            const auto rates =
                mlipper::model::build_gamma_rate_categories(alpha, categories);
            assert(rates.size() == static_cast<size_t>(categories));
            double sum = 0.0;
            for (double rate : rates) {
                assert(std::isfinite(rate));
                assert(rate > 0.0);
                sum += rate;
            }
            assert(std::fabs(sum / categories - 1.0) < 1.0e-12);
        }
    }

    const auto weights = mlipper::model::build_discrete_gamma_weights(4);
    assert(weights.size() == 4);
    for (double weight : weights) assert(weight == 0.25);

    std::filesystem::remove(gamma_path);
    std::filesystem::remove(free_rate_path);
    std::filesystem::remove(mixture_path);
    std::filesystem::remove(unknown_modifier_path);
    std::filesystem::remove(invalid_rate_path);
    std::filesystem::remove(invalid_alpha_path);
    std::filesystem::remove(test_directory);
    return 0;
}
