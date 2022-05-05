cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(graph VERSION 1.0)
file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(graph ${SRC_FILES})
target_compile_options(graph PUBLIC -Werror -Wall -Wextra -g -O3 -Ofast)
add_compile_options(-Wall -Wextra -pedantic -Werror -g -O3 -Ofast)
target_link_libraries(graph PUBLIC -pthread)
target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)



