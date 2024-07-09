cmake_minimum_required(VERSION 3.10)

project(sampling)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -std=c++17 -Wall -g")

find_package(Boost COMPONENTS program_options REQUIRED)

set(ZLIB_USE_STATIC_LIBS "ON")
find_package(ZLIB REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

add_executable (sampling
    sampling.cpp)

target_link_libraries(sampling PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
#target_link_libraries(sampling PRIVATE ${Boost_LIBRARIES})

