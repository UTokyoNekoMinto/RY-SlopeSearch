cmake_minimum_required(VERSION 3.20)
project(RY-SlopeSearch VERSION 1.0.0)

# === set C++ standard ===
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Ensure that the build type is Release ===
set(CMAKE_BUILD_TYPE Release)

# === set optimization level to -O3 for Release builds ===
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# # === Build type ===
# set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)

# # === Optimization + Debug symbols ===
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -fno-omit-frame-pointer")

# === OpenMP option switch ===
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# === add SeqAn3 path ===
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third_party/seqan3/cmake")

# === include range-v3 (SeqAn3's dependency) ===
include_directories(${CMAKE_SOURCE_DIR}/third_party/range-v3/include)

# === add SeqAn3 ===
find_package(seqan3 3.0 REQUIRED)

# === optional enable OpenMP ===
if(USE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)  # only find C++ part, avoid OpenMP_C error
    if(OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP for C++")
        add_compile_definitions(USE_OPENMP)  # make code can use #ifdef USE_OPENMP
    else()
        message(WARNING "OpenMP not found. Compiling without parallelization.")
        set(USE_OPENMP OFF)
    endif()
endif()

# === include header file path ===
include_directories(${CMAKE_SOURCE_DIR}/include)

# === find source files ===
file(GLOB_RECURSE SOURCES "src/*.cpp")

# === add executable file ===
add_executable(RY-SlopeSearch ${SOURCES})

# === link dependencies ===
target_link_libraries(RY-SlopeSearch PRIVATE seqan3::seqan3)

# === expose project version to the code ===
target_compile_definitions(RY-SlopeSearch PRIVATE RY_SLOPESEARCH_VERSION="${PROJECT_VERSION}")

# === link OpenMP (if available) ===
if(USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(RY-SlopeSearch PRIVATE OpenMP::OpenMP_CXX)
endif()

