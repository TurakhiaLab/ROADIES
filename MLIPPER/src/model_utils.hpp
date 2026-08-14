#pragma once
	
#include <filesystem>
#include <string>
#include <vector>
	
#include "io/parse_file.hpp"
	
namespace mlipper {
namespace model {

struct BestModelConfig {
    parse::ModelConfig model;
    bool empirical_freqs = false;
    std::string raw_text;
};

BestModelConfig parse_best_model_file(const std::filesystem::path& path);

std::vector<double> estimate_empirical_pi(const parse::Alignment& alignment, int states);
std::vector<double> ensure_normalized_pi(std::vector<double> pi, int states);

} // namespace model
} // namespace mlipper
