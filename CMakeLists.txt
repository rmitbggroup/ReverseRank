cmake_minimum_required(VERSION 3.9)
project(whynot)

set(CMAKE_CXX_STANDARD 11)
include_directories(C:/Users/liang/Google\ Drive/Code/glpk-4.65/src)
add_executable(whynot main.cpp)
TARGET_LINK_LIBRARIES(whynot C:/Users/liang/Google\ Drive/Code/glpk-4.65/w64/glpk_4_65.lib)
