#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <sstream>
#include <filesystem>

class Logger {
public:
    explicit Logger(const std::string& log_root_dir);
    ~Logger();

    // log output function
    template<typename T>
    void info(const T& msg) { write("INFO", msg); }

    template<typename T>
    void warn(const T& msg) { write("WARN", msg); }

    template<typename T>
    void error(const T& msg) { write("ERROR", msg); }

    // get log directory and filename
    std::string get_log_directory() const { return log_directory; }
    std::string get_log_filename() const { return log_filename; }

private:
    std::ofstream logfile;
    std::string log_directory;
    std::string log_filename;

    std::string timestamp(bool for_path = false);

    template<typename T>
    void write(const std::string& level, const T& msg) {
        std::ostringstream oss;
        oss << "[" << timestamp(false) << "] [" << level << "] " << msg;
        std::string logline = oss.str();

        std::cout << logline << std::endl;
        if (logfile.is_open()) {
            logfile << logline << std::endl;
        }
    }
};