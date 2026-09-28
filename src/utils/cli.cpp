#include "utils/cli.hpp"
#include "model/ry_sampling_word_sets.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifndef RY_SLOPESEARCH_VERSION
#define RY_SLOPESEARCH_VERSION "unknown"
#endif

namespace fs = std::filesystem;

namespace {

void print_usage(std::ostream& os, const char* prog) {
    os << "RY-SlopeSearch " << RY_SLOPESEARCH_VERSION
       << " - slope-based alignment-free DNA sequence comparison\n\n"
       << "Usage: " << prog << " -i <input_dir> -o <output_dir> [options]\n\n"
       << "Required:\n"
       << "  -i, --input DIR       directory containing .fasta files (one genome per file)\n"
       << "  -o, --output DIR      output root directory (created if it does not exist)\n\n"
       << "Options:\n"
       << "  -m, --method NAME     sampling word set (default: start_ry_128_matches)\n"
       << "      --many-to-many    use many-to-many matching (default: one-to-one)\n"
       << "      --background      subtract expected background matches (default: off)\n"
       << "  -t, --threads N       number of threads (default: all available)\n"
       << "      --no-progress     disable the progress bar\n"
       << "  -h, --help            print this help and exit\n"
       << "  -v, --version         print version and exit\n\n"
       << "Available methods:\n";
    for (const auto& [name, patterns] : kmer_sampling_methods) {
        os << "  " << name << "\n";
    }
}

[[noreturn]] void fail(const char* prog, const std::string& msg) {
    std::cerr << "Error: " << msg << "\n\n";
    print_usage(std::cerr, prog);
    std::exit(1);
}

int default_threads() {
#ifdef USE_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace

Options parse_args(int argc, char** argv) {
    const char* prog = argv[0];
    Options opts;
    opts.threads = default_threads();

    std::ostringstream cmd;
    for (int i = 0; i < argc; ++i) {
        cmd << (i ? " " : "") << argv[i];
    }
    opts.command_line = cmd.str();

    auto next_value = [&](int& i, const std::string& flag) -> std::string {
        if (i + 1 >= argc) fail(prog, "missing value for " + flag);
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(std::cout, prog);
            std::exit(0);
        } else if (arg == "-v" || arg == "--version") {
            std::cout << "RY-SlopeSearch " << RY_SLOPESEARCH_VERSION << "\n";
            std::exit(0);
        } else if (arg == "-i" || arg == "--input") {
            opts.input_directory = next_value(i, arg);
        } else if (arg == "-o" || arg == "--output") {
            opts.output_directory = next_value(i, arg);
        } else if (arg == "-m" || arg == "--method") {
            opts.sampling_method = next_value(i, arg);
        } else if (arg == "--many-to-many") {
            opts.use_one_to_one_matching = false;
        } else if (arg == "--background") {
            opts.use_background_matches = true;
        } else if (arg == "-t" || arg == "--threads") {
            std::string value = next_value(i, arg);
            size_t pos = 0;
            int n = 0;
            try {
                n = std::stoi(value, &pos);
            } catch (...) {
                pos = 0;
            }
            if (pos != value.size() || n < 1) fail(prog, "--threads must be a positive integer, got '" + value + "'");
            opts.threads = n;
        } else if (arg == "--no-progress") {
            opts.use_progress_bar = false;
        } else {
            fail(prog, "unknown argument '" + arg + "'");
        }
    }

    if (opts.input_directory.empty()) fail(prog, "--input is required");
    if (opts.output_directory.empty()) fail(prog, "--output is required");
    if (!fs::is_directory(opts.input_directory)) fail(prog, "input directory does not exist: " + opts.input_directory);
    if (kmer_sampling_methods.find(opts.sampling_method) == kmer_sampling_methods.end()) {
        fail(prog, "unknown method '" + opts.sampling_method + "'");
    }
#ifndef USE_OPENMP
    if (opts.threads > 1) {
        std::cerr << "Warning: built without OpenMP, running single-threaded\n";
        opts.threads = 1;
    }
#endif

    opts.input_directory = fs::absolute(opts.input_directory).lexically_normal().string();
    opts.output_directory = fs::absolute(opts.output_directory).lexically_normal().string();
    return opts;
}

std::string format_options(const Options& opts) {
    std::ostringstream oss;
    oss << "version: " << RY_SLOPESEARCH_VERSION << "\n"
        << "command: " << opts.command_line << "\n"
        << "input_directory: " << opts.input_directory << "\n"
        << "output_directory: " << opts.output_directory << "\n"
        << "sampling_method: " << opts.sampling_method << "\n"
        << "matching: " << (opts.use_one_to_one_matching ? "one-to-one" : "many-to-many") << "\n"
        << "background_matches: " << (opts.use_background_matches ? "subtracted" : "not subtracted") << "\n"
        << "threads: " << opts.threads << "\n";
    return oss.str();
}
