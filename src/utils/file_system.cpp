#include "utils/file_system.hpp"
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/views/to_char.hpp>

namespace fs = std::filesystem;

std::vector<std::string> load_sequences(const std::string& file_path) {
    std::vector<std::string> sequences;
    seqan3::sequence_file_input fasta_input(file_path);

    for (auto & record : fasta_input) {
        auto seq_view = record.sequence() | seqan3::views::to_char;
        sequences.emplace_back(std::string(seq_view.begin(), seq_view.end()));
    }

    if (sequences.size() < 2) {
        throw std::runtime_error("At least two sequences are required");
    }

    return sequences;
}

std::tuple<std::vector<std::string>, std::vector<std::string>>
load_sequences_for_evaluation(const std::string& file_path) {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    seqan3::sequence_file_input fasta_input(file_path);

    for (auto & record : fasta_input) {
        auto id_view = record.id() | seqan3::views::to_char;
        auto seq_view = record.sequence() | seqan3::views::to_char;

        names.emplace_back(std::string(id_view.begin(), id_view.end()));
        sequences.emplace_back(std::string(seq_view.begin(), seq_view.end()));
    }

    return {names, sequences};
}

std::tuple<std::vector<std::string>, std::vector<std::string>>
load_sequences_from_multiple_files(const std::string& directory) {
    std::vector<std::string> names;
    std::vector<std::string> sequences;

    for (auto const& entry : fs::recursive_directory_iterator(directory)) {
        if (entry.path().extension() == ".fasta") {
            seqan3::sequence_file_input fasta_input(entry.path());

            for (auto & record : fasta_input) {
                auto id_view = record.id() | seqan3::views::to_char;
                auto seq_view = record.sequence() | seqan3::views::to_char;

                names.emplace_back(std::string(id_view.begin(), id_view.end()));
                sequences.emplace_back(std::string(seq_view.begin(), seq_view.end()));
            }
        }
    }

    return {names, sequences};
}
