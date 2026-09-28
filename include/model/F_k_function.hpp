#pragma once
#include <string>
#include <vector>
#include "utils/logger.hpp"
#include "utils/cli.hpp"
#include "model/KmerCount.hpp"

#include <unordered_map>

class FKFunction {
public:
    FKFunction(const std::string& seq1, const std::string& seq2, int pattern_length, int min_k_min, const Options& opts, Logger& logger);
    double calculate_p_hat(const std::vector<std::vector<KmerCount>> &kmer_counts1,const std::vector<std::vector<KmerCount>> &kmer_counts2);
    void compute_fk(const std::vector<std::vector<KmerCount>> &kmer_counts1, const std::vector<std::vector<KmerCount>> &kmer_counts2);
    const std::vector<double>& get_FkLog() const { return Fk_log; }
    std::vector<int> k_vals;
    std::vector<double> Fk;
    std::vector<double> Fk_log;

private:
    std::string seq1, seq2;
    int L1, L2;
    int k_min, k_max;
    int min_k_min;
    int pattern_length;
    double L_avg;
    bool use_one_to_one_matching;
    bool use_background_matches;
    Logger& logger;
};
