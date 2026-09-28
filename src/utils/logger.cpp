#include "utils/logger.hpp"

Logger::Logger(const std::string& log_root_dir) {
    namespace fs = std::filesystem;

    // generate subdirectory based on timestamp
    std::string time_str = timestamp(true);
    log_directory = log_root_dir + "/" + time_str;

    // create directory
    fs::create_directories(log_directory);

    // build log file path
    log_filename = time_str + ".log";
    std::string log_path = log_directory + "/" + log_filename;

    // open log file
    logfile.open(log_path, std::ios::app);
}

Logger::~Logger() {
    if (logfile.is_open()) {
        logfile.close();
    }
}

std::string Logger::timestamp(bool for_path) {
    std::time_t now = std::time(nullptr);
    std::tm* tm = std::localtime(&now);
    char buf[64];

    if (for_path) {
        // for path: 2025-10-25_20-45-33
        std::strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", tm);
    } else {
        // for log line: 2025-10-25 20:45:33
        std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", tm);
    }

    return std::string(buf);
}