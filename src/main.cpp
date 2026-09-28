#include "utils/logger.hpp"
#include "utils/cli.hpp"
#include "model/F_k_function.hpp"
#include "model/match_tool.hpp"
#include "utils/evolution_model.hpp"
#include "utils/file_system.hpp"
#include "utils/progress_bar.hpp"
#include "model/ry_sampling_word_sets.hpp"
#include "utils/sequence_tool.hpp"
#include "model/KmerCount.hpp"

#include <unordered_map>
#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>
#include <string>
#include <set>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <numeric>

#ifdef USE_OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    // === parse command line arguments ===
    Options opts = parse_args(argc, argv);

#ifdef USE_OPENMP
    omp_set_num_threads(opts.threads);
#endif
    const bool use_parallel = opts.threads > 1;

    Logger logger(opts.output_directory);
    logger.info("Program started");
    std::string output_dir = logger.get_log_directory();
    logger.info("Output directory: " + output_dir);
    logger.info("Output file: " + logger.get_log_filename());

    // === start logging performance info ===
    std::string thread_log_path = output_dir + "/performance_info.log";
    std::ofstream thread_log(thread_log_path);
    if (!thread_log) {
        logger.error("Failed to open performance info log file: " + thread_log_path);
        return 1;
    }

    // === record run parameters ===
    std::string params_path = output_dir + "/run_parameters.txt";
    std::ofstream params_out(params_path);
    if (params_out) {
        params_out << format_options(opts);
        params_out.close();
    } else {
        logger.warn("Failed to write run parameters to: " + params_path);
    }
    logger.info("Sampling method: " + opts.sampling_method);
    logger.info("Threads: " + std::to_string(opts.threads));

    // === load sequences ===
    logger.info("Loading sequences from directory: " + opts.input_directory);
    auto [names, sequences] = load_sequences_from_multiple_files(opts.input_directory);
    size_t N = sequences.size();

    if (N < 2) {
        logger.error("At least two sequences are required, found " + std::to_string(N) + " in: " + opts.input_directory);
        return 1;
    }

    logger.info("Loaded " + std::to_string(N) + " sequences");

    // === compute distance matrix ===
    logger.info("Computing pairwise distances...");
    
    std::vector<std::vector<float>> distance_matrix(N, std::vector<float>(N, 0.0f));

    ProgressBar progress_bar(static_cast<size_t>(N * (N - 1) / 2), "Computing pairwise distances", opts.use_progress_bar);
    progress_bar.start();
    auto start_time = std::chrono::high_resolution_clock::now();  // start timing

    const auto& patterns_vec = get_patterns(opts.sampling_method);
    int pattern_length = get_pattern_length(opts.sampling_method);
    std::set<std::string> patterns(patterns_vec.begin(), patterns_vec.end());

    // map the DNA k-mer to a number
    char dna_kmer_to_num[128] = {0};
    dna_kmer_to_num['A'] = 0;
    dna_kmer_to_num['C'] = 1;
    dna_kmer_to_num['G'] = 2;
    dna_kmer_to_num['T'] = 3;

    // pattern lookup table
    size_t table_size = 1ULL << pattern_length; // 2^pattern_length
    std::vector<char> table(table_size, 0);
    for (const auto& pattern : patterns) {
        size_t pattern_num = 0;
        for (char j : pattern) {
            pattern_num = (pattern_num << 1) + (j == 'y');
        }
        table[pattern_num] = 1;
    }

    // calculate min and max k-mer length
    size_t min_seq_len = -1;  // maximum possible value of size_t
    size_t max_seq_len = 0;
    for (size_t i = 0; i < N; ++i) {
        min_seq_len = std::min(min_seq_len, sequences[i].size());
        max_seq_len = std::max(max_seq_len, sequences[i].size());
    }
    size_t min_L_avg = min_seq_len;
    size_t max_L_avg = max_seq_len;
    int min_k_min = std::max(static_cast<int>(std::ceil((std::log(min_L_avg) + 0.69) / 0.875)), pattern_length);
    int max_k_max = std::max(static_cast<int>(std::floor(std::log(max_L_avg) / 0.634)), pattern_length);

    // store k-mer counts for each sequence
    std::vector< std::vector< std::vector<KmerCount> > > kmerCounts(N);
    #pragma omp parallel for schedule(dynamic) collapse(1) if(use_parallel)
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < sequences[i].size(); ++j) {
            sequences[i][j] = dna_kmer_to_num[sequences[i][j]];
        }
    
        // get reverse complement of sequence
        std::string rc_seq = reverse_complement(sequences[i]);
        std::vector<std::string> seq_list = {sequences[i], rc_seq};

        std::vector< std::vector<size_t> > kmer_lists = extract_kmers_with_pattern(seq_list, min_k_min, max_k_max, dna_kmer_to_num, pattern_length, table);

        kmerCounts[i].resize(max_k_max - min_k_min + 1);

        for (int k = min_k_min; k <= max_k_max; ++k) {
            std::vector<size_t> &kmer_list = kmer_lists[k - min_k_min];
            std::vector<KmerCount> &counts = kmerCounts[i][k - min_k_min];

            sort(kmer_list.begin(), kmer_list.end());

            size_t oldPos = 0;
            size_t newPos = 0;

            while (newPos < kmer_list.size()) {
                ++newPos;
                if (newPos == kmer_list.size() || kmer_list[newPos] > kmer_list[oldPos]) {
                    int count_of_old_kmer = newPos - oldPos;
                    KmerCount kc = {kmer_list[oldPos], count_of_old_kmer};
                    counts.push_back(kc);
                    oldPos = newPos;
                }
            }
        }
    }



    #pragma omp parallel for schedule(dynamic) collapse(1) if(use_parallel)
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {  // only compute upper triangle part
            FKFunction fk(sequences[i], sequences[j], pattern_length, min_k_min, opts, logger);
            double p_hat = fk.calculate_p_hat(kmerCounts[i], kmerCounts[j]);

            // === compute evolutionary distance ===
            float distance = estimate_jukes_cantor_distance(static_cast<float>(p_hat), logger);
            
            distance_matrix[i][j] = distance;
            distance_matrix[j][i] = distance;  // symmetric assignment
            
            progress_bar.increment();
        }
    }

    progress_bar.finish();

    auto end_time = std::chrono::high_resolution_clock::now();  // end timing
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    thread_log << "Distance matrix computation took "
           << duration.count() << " milliseconds\n";
    thread_log << "\n";
    thread_log.close();  // close log

    // === save PHYLIP file ===
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::stringstream time_ss;
    time_ss << std::put_time(std::localtime(&now_c), "%Y%m%d_%H%M%S");
    std::string timestamp = time_ss.str();

    std::string phylip_path = output_dir + "/distance_matrix_" + timestamp + ".phy";
    std::ofstream phy_out(phylip_path);
    if (!phy_out) {
        logger.error("Failed to write PHYLIP file to: " + phylip_path);
        return 1;
    }

    phy_out << N << "\n";
    for (size_t i = 0; i < N; ++i) {
        std::string name = names[i].substr(0, 10);
        name.resize(10, ' ');  // PHYLIP name field fixed width
        phy_out << name;

        for (size_t j = 0; j < N; ++j) {
            phy_out << std::fixed << std::setprecision(6) << " " << distance_matrix[i][j];
        }
        phy_out << "\n";
    }

    logger.info("PHYLIP file saved to: " + phylip_path);
    logger.info("Program finished successfully.");
    return 0;
}

