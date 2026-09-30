# Wrapper build of the upstream C++ RSHash (branch multishapes) used as a test
# oracle. Builds from the clone without modifying it (the clone's
# test/CMakeLists.txt references a missing freq_mini.cpp).
#
#   module load gcc/14.2.0 cmake/3.29.2
#   cmake -S parity/cpp -B parity/work/cpp-build -DCMAKE_BUILD_TYPE=Release \
#         -DRSHASH_SRC=/path/to/rshash
#   cmake --build parity/work/cpp-build -j
cmake_minimum_required(VERSION 3.20...3.31)
project(rshash_parity CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(RSHASH_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../../rshash" CACHE PATH "upstream rshash clone")
list(APPEND CMAKE_PREFIX_PATH "${RSHASH_SRC}/seqan3/cmake" "${RSHASH_SRC}/sharg-parser/cmake")
find_package(seqan3 3.0 REQUIRED)
find_package(sharg 1.0 REQUIRED)

include(FetchContent)
FetchContent_Declare(gtl GIT_REPOSITORY https://github.com/greg7mdp/gtl.git GIT_TAG v1.2.0)
FetchContent_MakeAvailable(gtl)

set(S "${RSHASH_SRC}/source")
set(FLAGS -O3 -march=native -Wno-interference-size)

# the upstream binary, unchanged
add_executable(rshash ${S}/main.cpp ${S}/locate.cpp ${S}/lookup.cpp ${S}/build.cpp ${S}/io.cpp)
target_compile_options(rshash PRIVATE ${FLAGS})
target_include_directories(rshash PRIVATE ${RSHASH_SRC}/sux ${S})
target_compile_definitions(rshash PRIVATE PARSEARCH)
target_link_libraries(rshash PRIVATE seqan3::seqan3 sharg::sharg gtl)

# per-record dump driver (same library sources, our main)
add_executable(rshash_dump cpp_dump.cpp ${S}/locate.cpp ${S}/lookup.cpp ${S}/build.cpp ${S}/io.cpp)
target_compile_options(rshash_dump PRIVATE ${FLAGS})
target_include_directories(rshash_dump PRIVATE ${RSHASH_SRC}/sux ${S})
target_compile_definitions(rshash_dump PRIVATE PARSEARCH)
target_link_libraries(rshash_dump PRIVATE seqan3::seqan3 sharg::sharg gtl)

# the upstream hash-table oracle
add_executable(hashtable ${RSHASH_SRC}/test/hashtable.cpp)
target_compile_options(hashtable PRIVATE ${FLAGS})
target_include_directories(hashtable PRIVATE ${S})
target_link_libraries(hashtable PRIVATE seqan3::seqan3 sharg::sharg gtl)

# component micro-benchmark (mirrors crates/rshash/examples/components.rs)
add_executable(cpp_components cpp_components.cpp ${S}/locate.cpp ${S}/lookup.cpp ${S}/build.cpp ${S}/io.cpp)
target_compile_options(cpp_components PRIVATE ${FLAGS})
target_include_directories(cpp_components PRIVATE ${RSHASH_SRC}/sux ${S})
target_compile_definitions(cpp_components PRIVATE PARSEARCH)
target_link_libraries(cpp_components PRIVATE seqan3::seqan3 sharg::sharg gtl)
