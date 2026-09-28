#pragma once
#include <string>

struct Options {
    std::string input_directory;
    std::string output_directory;
    std::string sampling_method = "start_ry_128_matches";
    bool use_one_to_one_matching = true;
    bool use_background_matches = false;
    int threads = 1;
    bool use_progress_bar = true;
    std::string command_line;
};

// parse command line arguments; exits the program on --help, --version or invalid input
Options parse_args(int argc, char** argv);

// human-readable summary of the parsed options
std::string format_options(const Options& opts);
