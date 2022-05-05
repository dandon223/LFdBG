cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_STANDARD 20)
project(tests)



message(STATUS "start running cmake...")

find_package(Boost 1.71.0 COMPONENTS unit_test_framework REQUIRED)
if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})

endif()
add_executable(tests unitTests.cpp )

if(Boost_FOUND)

    target_link_libraries(tests ${Boost_LIBRARIES} pthread)

endif()
